"""Tests for HALOCVER = 1 FITS contract and Swinburne NT / OzSTAR HPC helpers."""

from __future__ import annotations

import ast
import os
from pathlib import Path

import astropy.units as u
import lightkurve as lk
import numpy as np
import pytest

from mundey import (
    HALOCVER,
    PixelQuality,
    _contract,
    append_smear_err_column,
    atomic_publish,
    atomic_temp_path,
    make_atomic_tmp_path,
    read_mundey_provenance,
    read_pixel_quality,
    read_smear_err_column,
    resolve_cache_dir,
    write_mundey_provenance,
    write_pixel_quality,
)


def test_contract_module_has_zero_cross_package_or_fcntl_imports() -> None:
    """_contract.py must not import tessbkgd, halophot, jax, or fcntl."""
    src = Path(_contract.__file__).read_text(encoding="utf-8")
    tree = ast.parse(src)
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imported.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])

    forbidden = {"tessbkgd", "halophot", "fcntl", "IPython", "pkg_resources"}
    assert imported.isdisjoint(forbidden)


def test_pixel_quality_bit_positions() -> None:
    """Verify exact bit positions specified by HALOCVER = 1."""
    assert HALOCVER == 1
    assert int(PixelQuality.STRAP_COL) == (1 << 0)
    assert int(PixelQuality.SAT_BLEED) == (1 << 1)
    assert int(PixelQuality.SMEAR_UNRESOLVED) == (1 << 2)
    assert int(PixelQuality.SOURCE_MASK) == (1 << 3)
    assert int(PixelQuality.BKG_USED) == (1 << 4)

    combo = PixelQuality.STRAP_COL | PixelQuality.SAT_BLEED | PixelQuality.BKG_USED
    assert int(combo) == 0b10011


def test_halocver_roundtrip_and_slicing_lockstep(
    synthetic_2output_tpf: lk.TessTargetPixelFile,
    tmp_path: Path,
) -> None:
    """Verify BKGSUB=False, PIXEL_QUALITY ImageHDU, and tpf[10:30] SMEAR_ERR slicing."""
    tpf = synthetic_2output_tpf
    n_raw = len(tpf.hdu[1].data)
    _, ny, nx = tpf.flux.shape
    assert n_raw > len(tpf.flux)  # fixture includes 1 quality-flagged cadence

    write_mundey_provenance(tpf, bkgsub=False)
    assert tpf.hdu[0].header["HALOCVER"] == 1
    assert tpf.hdu[1].header["HALOCVER"] == 1
    assert tpf.hdu[0].header["BKGSUB"] is False
    assert tpf.hdu[1].header["BKGSUB"] is False
    assert tpf.meta["HALOCVER"] == 1
    assert tpf.meta["BKGSUB"] is False
    assert read_mundey_provenance(tpf) == {"halocver": 1, "bkgsub": False}

    pq = np.zeros((ny, nx), dtype=np.uint16)
    pq[:, 0] = int(PixelQuality.STRAP_COL)
    pq[4:, nx // 2] = int(PixelQuality.SAT_BLEED | PixelQuality.SMEAR_UNRESOLVED)
    pq[0, -1] = int(PixelQuality.BKG_USED)
    write_pixel_quality(tpf, pq)
    # Idempotent overwrite must not duplicate PIXEL_QUALITY extension
    write_pixel_quality(tpf, pq)
    assert [h.name for h in tpf.hdu].count("PIXEL_QUALITY") == 1
    assert read_pixel_quality(tpf).dtype == np.uint16
    np.testing.assert_array_equal(read_pixel_quality(tpf), pq)

    # Append united Quantity SMEAR_ERR over all N_raw rows of HDU 1
    smear_err_raw = (
        np.arange(n_raw * nx, dtype=np.float32).reshape(n_raw, nx) * 0.25 + 1.0
    ) * (u.electron / u.s)
    append_smear_err_column(tpf, smear_err_raw)

    # Preserve critical HDU 1 headers and WCS
    assert tpf.hdu[1].name == "PIXELS"
    assert tpf.hdu[1].header["1CRV4P"] == 551
    assert tpf.hdu[1].header["BKGSUB"] is False
    assert np.isclose(tpf.hdu[1].header["INT_TIME"], 1.98)

    masked_smear_err = read_smear_err_column(tpf, quality_masked=True)
    assert masked_smear_err.shape == (len(tpf.flux), nx)
    np.testing.assert_allclose(masked_smear_err, smear_err_raw.value[tpf.quality_mask])

    # Cadence slicing tpf[10:30] must slice SMEAR_ERR in lockstep with tpf.flux
    sliced = tpf[10:30]
    assert sliced.flux.shape == (20, ny, nx)
    sliced_smear_err = read_smear_err_column(sliced)
    assert sliced_smear_err.shape == (20, nx)
    assert len(sliced.hdu[1].data["SMEAR_ERR"]) == len(sliced.flux) == 20
    np.testing.assert_allclose(sliced_smear_err, masked_smear_err[10:30])
    # Spatial PIXEL_QUALITY ImageHDU must remain unsliced (ny, nx)
    np.testing.assert_array_equal(read_pixel_quality(sliced), pq)

    # Full FITS disk round-trip via tpf.to_fits -> lk.read
    out_fits = tmp_path / "roundtrip_tpf.fits"
    tpf.to_fits(str(out_fits), overwrite=True)
    reloaded = lk.read(out_fits)

    assert reloaded.meta["HALOCVER"] == 1
    assert reloaded.meta["BKGSUB"] is False
    np.testing.assert_array_equal(read_pixel_quality(reloaded), pq)
    np.testing.assert_allclose(read_smear_err_column(reloaded), masked_smear_err)

    reloaded_slice = reloaded[10:30]
    assert reloaded_slice.flux.shape[0] == 20
    assert read_smear_err_column(reloaded_slice).shape == (20, nx)
    np.testing.assert_allclose(read_smear_err_column(reloaded_slice), masked_smear_err[10:30])


def test_contract_validation_errors(synthetic_1output_tpf: lk.TessTargetPixelFile) -> None:
    """Invalid PIXEL_QUALITY or SMEAR_ERR shapes/dtypes must raise clear errors."""
    tpf = synthetic_1output_tpf
    _, ny, nx = tpf.flux.shape

    with pytest.raises(TypeError, match="integer or boolean dtype"):
        write_pixel_quality(tpf, np.zeros((ny, nx), dtype=np.float32))

    with pytest.raises(ValueError, match="does not match TPF spatial shape"):
        write_pixel_quality(tpf, np.zeros((ny + 1, nx), dtype=np.uint16))

    with pytest.raises(ValueError, match="does not match HDU\\[1\\]"):
        append_smear_err_column(tpf, np.ones((len(tpf.flux) - 3, nx), dtype=np.float32))


def test_resolve_cache_dir_scratch_default_and_fred_hpc_guard(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """resolve_cache_dir must use $SCRATCH/.cache/tess_bright and block /fred in strict_hpc."""
    monkeypatch.delenv("MUNDEY_CACHE_DIR", raising=False)
    scratch_dir = tmp_path / "aphid_scratch"
    monkeypatch.setenv("SCRATCH", str(scratch_dir))

    resolved = resolve_cache_dir(strict_hpc=True, create=True)
    assert resolved == scratch_dir / ".cache" / "tess_bright"
    assert resolved.is_dir()

    # Direct /fred paths must fail in strict_hpc=True
    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir("/fred/oz440/user/.cache", strict_hpc=True)

    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir("/fred", strict_hpc=True)

    # Traversal into /fred must also fail
    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir("/tmp/../fred/oz440/cache", strict_hpc=True)

    # Symlink into /fred must also fail
    fred_link = tmp_path / "symlink_to_fred"
    fred_link.symlink_to("/fred/oz440/cache")
    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir(fred_link, strict_hpc=True)

    # Non-/fred prefix like /frederick must NOT be falsely rejected
    assert resolve_cache_dir("/frederick/cache", strict_hpc=True) == Path("/frederick/cache")


def test_atomic_publish_same_dir_cross_dir_and_fsync(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """atomic_publish must fsync and os.replace cleanly without leaving .tmp files."""
    dst = tmp_path / "cache" / "sector06-blf.fits"
    tmp_candidate = make_atomic_tmp_path(dst)
    assert tmp_candidate.parent == dst.parent
    assert tmp_candidate.name.startswith("sector06-blf.fits.")
    assert f".{os.getpid()}." in tmp_candidate.name
    assert tmp_candidate.name.endswith(".tmp")

    fsync_fds: list[int] = []
    real_fsync = os.fsync

    def tracked_fsync(fd: int) -> None:
        fsync_fds.append(fd)
        real_fsync(fd)

    monkeypatch.setattr(os, "fsync", tracked_fsync)

    # 1. Same-directory atomic_temp_path publish
    with atomic_temp_path(dst) as staged:
        staged.write_bytes(b"initial-payload")
    assert dst.read_bytes() == b"initial-payload"
    assert len(fsync_fds) >= 1
    assert list(dst.parent.glob("*.tmp")) == []

    # 2. Cross-directory ($JOBFS -> $SCRATCH) atomic_publish overwrite
    jobfs_tmp = tmp_path / "jobfs" / "worker_out.bin"
    jobfs_tmp.parent.mkdir(parents=True)
    jobfs_tmp.write_bytes(b"updated-from-jobfs")

    published = atomic_publish(jobfs_tmp, dst)
    assert published == dst
    assert dst.read_bytes() == b"updated-from-jobfs"
    assert not jobfs_tmp.exists()
    assert list(dst.parent.glob("*.tmp")) == []

    # 3. Failed write inside atomic_temp_path cleans up temp file and preserves dst
    with pytest.raises(RuntimeError, match="abort write"):
        with atomic_temp_path(dst) as staged_fail:
            staged_fail.write_bytes(b"corrupt-partial")
            raise RuntimeError("abort write")
    assert dst.read_bytes() == b"updated-from-jobfs"
    assert list(dst.parent.glob("*.tmp")) == []
