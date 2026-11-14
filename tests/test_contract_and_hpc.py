"""Tests for HALOCVER = 1 FITS contract and Swinburne NT / OzSTAR HPC helpers."""

from __future__ import annotations

import ast
import os
import shutil
import subprocess
from pathlib import Path

import astropy.units as u
import lightkurve as lk
import numpy as np
import pytest

import mundey
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
    read_smear_err,
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

    assert tpf.hdu[1].header["BACKAPP"] is True
    write_mundey_provenance(tpf, bkgsub=False)
    assert tpf.hdu[0].header["HALOCVER"] == 1
    assert tpf.hdu[1].header["HALOCVER"] == 1
    assert tpf.hdu[0].header["BKGSUB"] is False
    assert tpf.hdu[1].header["BKGSUB"] is False
    assert tpf.hdu[1].header["BACKAPP"] is False
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
    """Invalid PIXEL_QUALITY, SMEAR_ERR, or provenance headers must raise clear errors."""
    tpf = synthetic_1output_tpf
    _, ny, nx = tpf.flux.shape

    with pytest.raises(TypeError, match="integer or boolean dtype"):
        write_pixel_quality(tpf, np.zeros((ny, nx), dtype=np.float32))

    with pytest.raises(TypeError, match="integer or boolean dtype"):
        write_pixel_quality(tpf, np.full((ny, nx), "0", dtype="<U2"))

    with pytest.raises(ValueError, match="negative"):
        write_pixel_quality(tpf, np.full((ny, nx), -1, dtype=np.int32))

    overflow = np.zeros((ny, nx), dtype=np.uint32)
    overflow[0, 0] = 1 << 16
    with pytest.raises(ValueError, match="exceeding uint16 max"):
        write_pixel_quality(tpf, overflow)

    with pytest.raises(ValueError, match="does not match TPF spatial shape"):
        write_pixel_quality(tpf, np.zeros((ny + 1, nx), dtype=np.uint16))

    with pytest.raises(ValueError, match="does not match HDU\\[1\\]"):
        append_smear_err_column(tpf, np.ones((len(tpf.flux) - 3, nx), dtype=np.float32))

    write_mundey_provenance(tpf, bkgsub=False)
    tpf.hdu[1].header["BKGSUB"] = True
    with pytest.raises(ValueError, match="Conflicting BKGSUB"):
        read_mundey_provenance(tpf)

    tpf.hdu[1].header["BKGSUB"] = "F"
    with pytest.raises(ValueError, match="Invalid BKGSUB"):
        read_mundey_provenance(tpf)

    tpf.hdu[1].header["BKGSUB"] = False
    tpf.hdu[1].header["HALOCVER"] = 2
    with pytest.raises(ValueError, match="Conflicting HALOCVER"):
        read_mundey_provenance(tpf)


def test_resolve_cache_dir_scratch_default_and_fred_hpc_guard(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """resolve_cache_dir must enforce precedence, strict_hpc=True default, and /fred guards."""
    monkeypatch.delenv("MUNDEY_CACHE_DIR", raising=False)
    monkeypatch.delenv("TESS_BRIGHT_CACHE_DIR", raising=False)
    monkeypatch.delenv("SCRATCH", raising=False)
    monkeypatch.delenv("XDG_CACHE_HOME", raising=False)

    # 1. Fallback to ~/.cache/tess_bright (create=False by default does not mkdir)
    fake_home = tmp_path / "fake_home"
    fake_home.mkdir()
    monkeypatch.setenv("HOME", str(fake_home))
    home_resolved = resolve_cache_dir()
    assert home_resolved == (fake_home / ".cache" / "tess_bright").resolve(strict=False)
    assert not home_resolved.exists()

    # 2. $XDG_CACHE_HOME/tess_bright overrides ~/.cache/tess_bright
    xdg_dir = tmp_path / "xdg_cache"
    monkeypatch.setenv("XDG_CACHE_HOME", str(xdg_dir))
    assert resolve_cache_dir() == (xdg_dir / "tess_bright").resolve(strict=False)

    # 3. $SCRATCH/.cache/tess_bright overrides $XDG_CACHE_HOME; create=True creates directory
    scratch_dir = tmp_path / "aphid_scratch"
    monkeypatch.setenv("SCRATCH", str(scratch_dir))
    resolved = resolve_cache_dir(create=True)
    assert resolved == (scratch_dir / ".cache" / "tess_bright").resolve(strict=False)
    assert resolved.is_dir()

    # 4. $TESS_BRIGHT_CACHE_DIR overrides $SCRATCH, and $MUNDEY_CACHE_DIR overrides all env vars
    tb_dir = tmp_path / "tess_bright_env"
    monkeypatch.setenv("TESS_BRIGHT_CACHE_DIR", str(tb_dir))
    assert resolve_cache_dir() == tb_dir.resolve(strict=False)

    mundey_dir = tmp_path / "mundey_env"
    monkeypatch.setenv("MUNDEY_CACHE_DIR", str(mundey_dir))
    assert resolve_cache_dir() == mundey_dir.resolve(strict=False)

    # Explicit cache_dir argument overrides all env vars
    explicit_dir = tmp_path / "explicit_cache"
    assert resolve_cache_dir(explicit_dir) == explicit_dir.resolve(strict=False)

    # 5. Direct /fred paths must fail by default (strict_hpc=True is the default!)
    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir("/fred/oz440/user/.cache")

    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir("/fred")

    # Traversal into /fred must also fail
    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir("/tmp/../fred/oz440/cache")

    # Symlink into /fred must fail by default (including via $TESS_BRIGHT_CACHE_DIR / $HOME)
    fred_link = tmp_path / "symlink_to_fred"
    fred_link.symlink_to("/fred/oz440/cache")
    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir(fred_link)

    monkeypatch.delenv("MUNDEY_CACHE_DIR", raising=False)
    monkeypatch.setenv("TESS_BRIGHT_CACHE_DIR", str(fred_link))
    with pytest.raises(ValueError, match="Refusing cache directory .* under /fred"):
        resolve_cache_dir()

    # 6. strict_hpc=False emits UserWarning on /fred (direct and symlink) and proceeds
    with pytest.warns(UserWarning, match="under /fred"):
        warned_direct = resolve_cache_dir("/fred/oz440/user/.cache", strict_hpc=False)
    assert warned_direct == Path("/fred/oz440/user/.cache").resolve(strict=False)

    with pytest.warns(UserWarning, match="under /fred"):
        warned_link = resolve_cache_dir(fred_link, strict_hpc=False)
    assert warned_link == Path("/fred/oz440/cache").resolve(strict=False)

    # Non-/fred prefix like /frederick must NOT be falsely rejected
    assert resolve_cache_dir("/frederick/cache") == Path("/frederick/cache").resolve(strict=False)


def test_version_2_0_0_and_read_smear_err_export_symmetry(
    synthetic_1output_bundle,
) -> None:
    """mundey v2.0.0 must export read_smear_err alias and SMEARMOD provenance card."""
    assert mundey.__version__ == "2.0.0"
    assert "read_smear_err" in mundey.__all__
    assert read_smear_err is read_smear_err_column

    tpf = lk.read(synthetic_1output_bundle.tpf_path)
    n_raw = len(tpf.hdu[1].data)
    nx = tpf.shape[2]
    vals = np.arange(n_raw * nx, dtype=np.float32).reshape(n_raw, nx)
    append_smear_err_column(tpf, vals)
    write_mundey_provenance(tpf, bkgsub=False, smear_mode="standard")

    assert tpf.hdu[0].header["SMEARMOD"] == "standard"
    assert tpf.meta["SMEARMOD"] == "standard"
    np.testing.assert_array_equal(
        read_smear_err(tpf),
        read_smear_err_column(tpf),
    )
    np.testing.assert_array_equal(
        read_smear_err(tpf, apply_quality_mask=False),
        vals,
    )


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


def test_hpc_sbatch_redirects_apptainer_cache_and_tmp_off_home() -> None:
    """Slurm template must configure NT bash, Apptainer, scratch logs, and trevor guidance."""
    repo_root = Path(__file__).resolve().parent.parent
    sbatch_path = repo_root / "hpc" / "mundey_prefetch.sbatch"
    text = sbatch_path.read_text(encoding="utf-8")

    assert text.startswith("#!/bin/bash\n")
    assert "set -eo pipefail" in text
    assert "set -euo pipefail" not in text
    assert ': "${SECTOR:?Set SECTOR=<int> via sbatch --export=ALL,SECTOR=...}"' in text

    assert "module -q load apptainer" in text
    assert "$MYFRED/containers/halo.sif" in text
    assert "virgil.sif" not in text

    assert "/aphid/scratch-3month/%u/logs/%x-%A_%a.out" in text
    assert "/aphid/scratch-3month/%u/logs/%x-%A_%a.err" in text
    assert "trevor" in text

    assert (
        'export TESS_BRIGHT_CACHE_DIR="${TESS_BRIGHT_CACHE_DIR:-$SCRATCH/.cache/tess_bright}"'
        in text
    )
    assert 'export MUNDEY_CACHE_DIR="${MUNDEY_CACHE_DIR:-$TESS_BRIGHT_CACHE_DIR}"' in text
    assert 'export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$SCRATCH/.apptainer/cache}"' in text
    assert 'export APPTAINER_TMPDIR="${JOBFS:-$STAGE_DIR}"' in text
    assert 'export APPTAINER_BIND="/fred,/aphid/scratch-3month${JOBFS:+,$JOBFS}"' in text
    assert 'export APPTAINER_BINDPATH="/fred,/aphid/scratch-3month${JOBFS:+,$JOBFS}"' in text
    assert '"$APPTAINER_CACHEDIR"' in text
    assert 'mkdir -p "$STAGE_DIR" "$APPTAINER_TMPDIR" "$OUT"' in text

    bash_bin = shutil.which("bash")
    if bash_bin is not None:
        subprocess.run([bash_bin, "-n", str(sbatch_path)], check=True)

    readme_text = (repo_root / "README.md").read_text(encoding="utf-8")
    assert "trevor" in readme_text
    assert "/aphid/scratch-3month/$USER/logs" in readme_text

    pyproject_text = (repo_root / "pyproject.toml").read_text(encoding="utf-8")
    assert "zarr" not in pyproject_text
