"""Tests for mundey.archive (chunked TSCAL/TZERO, S3 mock, --offline) and mundey CLI."""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from mundey.archive import (
    S3_TESS_PUBLIC_ROOT,
    iter_collateral_chunks,
    prefetch_collateral,
    resolve_collateral_file,
    resolve_static_models,
)
from mundey.cli import main as cli_main


class FakeS3FileSystem:
    """In-memory mock of anonymous ``s3fs.S3FileSystem`` for deterministic offline CI."""

    def __init__(self, files: dict[str, bytes]) -> None:
        self._files = dict(files)
        self.opened_keys: list[str] = []

    def ls(self, prefix: str) -> list[str]:
        clean = prefix.rstrip("/") + "/"
        return sorted(k for k in self._files if k.startswith(clean))

    def glob(self, pattern: str) -> list[str]:
        import fnmatch

        return sorted(k for k in self._files if fnmatch.fnmatch(k, pattern))

    def open(self, key: str, mode: str = "rb") -> io.BytesIO:
        assert mode == "rb"
        if key not in self._files:
            raise FileNotFoundError(key)
        self.opened_keys.append(key)
        return io.BytesIO(self._files[key])


def test_offline_guard_blocks_before_s3_instantiation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """offline=True must raise FileNotFoundError before touching fsspec / network."""

    def _forbid_s3(**_kwargs):
        raise AssertionError("S3 filesystem must NEVER be opened when offline=True")

    monkeypatch.setattr("mundey.archive._open_anon_s3_filesystem", _forbid_s3)

    with pytest.raises(FileNotFoundError, match="offline=True"):
        resolve_static_models(2, 3, cache_dir=tmp_path, offline=True)

    with pytest.raises(FileNotFoundError, match="offline=True"):
        resolve_collateral_file(
            sector=6,
            camera=2,
            ccd=3,
            output="B",
            kind="smrow",
            cache_dir=tmp_path,
            offline=True,
        )


def test_chunked_tscal_tzero_collateral_reader_and_sentinel_masking(tmp_path: Path) -> None:
    """iter_collateral_chunks must stream TSCALn/TZEROn unsigned integer tables in slices."""
    n_rows = 25
    # Unsigned uint16 stored via BZERO/TZERO=32768, TSCAL=1 with sentinel TNULL=-32768
    vals_u16 = (np.arange(n_rows * 4, dtype=np.uint16).reshape(n_rows, 2, 2) + np.uint16(100))
    vals_u16[7, 0, 1] = np.uint16(0)  # raw stored int16 = 0 - 32768 = -32768 == TNULL

    col_time = fits.Column(name="TIME", format="D", array=np.arange(n_rows, dtype=np.float64))
    col_cad = fits.Column(name="CADENCENO", format="J", array=1000 + np.arange(n_rows))
    col_raw = fits.Column(
        name="SMROW_RAW",
        format="4I",
        bscale=1.0,
        bzero=32768,
        null=-32768,
        dim="(2,2)",
        array=vals_u16,
    )
    col_path = tmp_path / "scaled_col.fits"
    tbl = fits.BinTableHDU.from_columns([col_time, col_cad, col_raw])
    fits.HDUList([fits.PrimaryHDU(), tbl]).writeto(col_path, overwrite=True)

    chunks = list(iter_collateral_chunks(col_path, "SMROW_RAW", chunk_rows=10))
    assert [(c[0], c[1]) for c in chunks] == [(0, 10), (10, 20), (20, 25)]

    stacked = np.concatenate([c[4] for c in chunks], axis=0)
    assert stacked.dtype == np.float64
    assert np.isnan(stacked[7, 0, 1])

    mask = np.ones_like(vals_u16, dtype=bool)
    mask[7, 0, 1] = False
    np.testing.assert_allclose(stacked[mask], vals_u16[mask].astype(np.float64))

    # Selective row_indices streaming must skip chunks [0, 10) and [10, 20) completely!
    sub_chunks = list(
        iter_collateral_chunks(
            col_path,
            "SMROW_RAW",
            chunk_rows=10,
            row_indices=np.array([21, 23], dtype=np.int64),
        )
    )
    assert [(c[0], c[1]) for c in sub_chunks] == [(20, 25)]


def test_prefetch_from_mock_s3_stages_jobfs_and_cleans_raw_col(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """prefetch_collateral must fetch from S3, stream chunks, publish -blf, and unlink raw stage."""
    src_dir = tmp_path / "src_gen"
    tpf_builder(
        src_dir,
        col0_1based=200,
        row0_1based=1500,
        nx=8,
        ny=10,
        n_cadences=12,
        sector=6,
        camera=2,
        ccd=3,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        seed=303,
    )

    # Populate FakeS3FileSystem using real MAST S3 path layout (including 1-s timestamp offset!)
    s3_objects: dict[str, bytes] = {
        f"{S3_TESS_PUBLIC_ROOT}/models/tess2018323-2-3-2dblack.fits": (
            src_dir / "tess2018323-2-3-2dblack.fits"
        ).read_bytes(),
        f"{S3_TESS_PUBLIC_ROOT}/models/tess2018323-2-3-flat.fits": (
            src_dir / "tess2018323-2-3-flat.fits"
        ).read_bytes(),
        f"{S3_TESS_PUBLIC_ROOT}/models/tess2018143203310-41006_100-linearity.xml": (
            src_dir / "tess2018143203310-41006_100-linearity.xml"
        ).read_bytes(),
        (
            f"{S3_TESS_PUBLIC_ROOT}/ffi/s0006/2018/349/2-3/"
            "tess2020001000001-s0006-smrow-2-3-a-0165-s_col.fits"
        ): next(src_dir.glob("*-smrow-2-3-a-*-s_col.fits")).read_bytes(),
        (
            f"{S3_TESS_PUBLIC_ROOT}/ffi/s0006/2018/349/2-3/"
            "tess2020001000001-s0006-tvcol-2-3-a-0165-s_col.fits"
        ): next(src_dir.glob("*-tvcol-2-3-a-*-s_col.fits")).read_bytes(),
    }
    fake_s3 = FakeS3FileSystem(s3_objects)

    cache_dir = tmp_path / "scratch_cache"
    jobfs_dir = tmp_path / "jobfs"

    published = prefetch_collateral(
        sector=6,
        camera=2,
        ccd=3,
        outputs=("A",),
        cadences=("2min",),
        cache_dir=cache_dir,
        stage_dir=jobfs_dir,
        chunk_rows=5,  # force 3 streaming chunks over 12 cadences
        cleanup_staged_raw=True,
        s3_fs=fake_s3,
    )
    assert len(published) == 5
    for p in published:
        assert p.is_file()

    # Raw staged _col.fits inside $JOBFS were unlinked after -blf.fits publication, and zero .tmp!
    assert list(jobfs_dir.glob("*-s_col.fits")) == []
    assert list(jobfs_dir.glob("*.tmp")) == []
    assert list((cache_dir / "s0006").glob("*.tmp")) == []

    # Subsequent offline recalibrate_tpf succeeds using ONLY the prefetched cache_dir
    out_cal = cli_main(
        [
            "calibrate",
            str(src_dir / "tess2020001000000-s0006-0000000123456789-0165-s_tp.fits"),
            "-o",
            str(tmp_path / "cli_cal_tp.fits"),
            "--cache-dir",
            str(cache_dir),
            "--smear",
            "standard",
            "--offline",
            "--overwrite",
        ]
    )
    assert out_cal == 0
    assert (tmp_path / "cli_cal_tp.fits").is_file()


def test_cli_prefetch_task_id_and_tpf_inference(tpf_builder, tmp_path: Path) -> None:
    """mundey prefetch CLI must support --task-id 0..15 and --tpf inference offline."""
    bundle = tpf_builder(
        tmp_path / "cli_pre",
        col0_1based=200,
        row0_1based=1500,
        nx=6,
        ny=8,
        n_cadences=4,
        sector=9,
        camera=2,  # task_id = (2-1)*4 + (3-1) = 6
        ccd=3,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        seed=12,
    )
    cache1 = tmp_path / "cache1"
    rc1 = cli_main(
        [
            "prefetch",
            "--sector",
            "9",
            "--task-id",
            "6",
            "--outputs",
            "A",
            "--ddir",
            str(bundle.ddir),
            "--cache-dir",
            str(cache1),
            "--chunk-rows",
            "2",
            "--offline",
        ]
    )
    assert rc1 == 0
    assert len(list((cache1 / "s0009").glob("*-blf.fits"))) == 2

    cache2 = tmp_path / "cache2"
    rc2 = cli_main(
        [
            "prefetch",
            "--tpf",
            str(bundle.tpf_path),
            "--ddir",
            str(bundle.ddir),
            "--cache-dir",
            str(cache2),
            "--offline",
        ]
    )
    assert rc2 == 0
    assert len(list((cache2 / "s0009").glob("*-blf.fits"))) == 2

    # --strict-hpc rejects /fred cache directory via CLI
    with pytest.raises(ValueError, match="under /fred"):
        cli_main(
            [
                "prefetch",
                "--sector",
                "9",
                "--task-id",
                "6",
                "--outputs",
                "A",
                "--cache-dir",
                "/fred/oz440/bad_cache",
                "--strict-hpc",
                "--offline",
            ]
        )
