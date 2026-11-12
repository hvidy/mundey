"""Versioned FITS contract (``HALOCVER = 1``) and HPC filesystem helpers.

This module has zero cross-package imports: ``mundey``, ``tessbkgd``, and
``halophot`` communicate strictly through standard ``lightkurve.TargetPixelFile``
and ``astropy.io.fits.HDUList`` structures conforming to ``HALOCVER = 1``.
"""

from __future__ import annotations

import contextlib
import errno
import os
import shutil
import socket
import uuid
from collections.abc import Iterator
from enum import IntFlag
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits

HALOCVER: int = 1
PIXEL_QUALITY_EXTNAME: str = "PIXEL_QUALITY"
SMEAR_ERR_COLNAME: str = "SMEAR_ERR"


class PixelQuality(IntFlag):
    """Spatial pixel-quality bitmask flags for ``PIXEL_QUALITY`` (``uint16``)."""

    STRAP_COL = 1 << 0
    SAT_BLEED = 1 << 1
    SMEAR_UNRESOLVED = 1 << 2
    SOURCE_MASK = 1 << 3
    BKG_USED = 1 << 4


def _unwrap_hdulist(target: fits.HDUList | Any) -> fits.HDUList:
    """Return the underlying ``fits.HDUList`` from an HDUList or TargetPixelFile."""
    if isinstance(target, fits.HDUList):
        return target
    if hasattr(target, "hdu") and isinstance(target.hdu, fits.HDUList):
        return target.hdu
    raise TypeError(
        f"Expected astropy.io.fits.HDUList or object with .hdu HDUList, got {type(target)!r}."
    )


def _infer_cutout_shape(hdul: fits.HDUList) -> tuple[int, int, int] | None:
    """Infer ``(n_raw, ny, nx)`` from HDU 1 ``PIXELS`` table columns if present."""
    if len(hdul) < 2 or getattr(hdul[1], "data", None) is None:
        return None
    data = hdul[1].data
    names = getattr(data, "names", ()) or ()
    for col_name in ("FLUX", "RAW_CNTS", "FLUX_ERR"):
        if col_name in names:
            shape = data[col_name].shape
            if len(shape) == 3:
                return int(shape[0]), int(shape[1]), int(shape[2])
    return None


def write_mundey_provenance(
    target: fits.HDUList | Any,
    *,
    bkgsub: bool = False,
    halocver: int = HALOCVER,
) -> fits.HDUList:
    """Write ``HALOCVER`` and ``BKGSUB`` provenance cards to ``hdu[0]`` and ``hdu[1]``.

    Writing to ``hdu[0].header`` exposes ``tpf.meta['HALOCVER']`` and
    ``tpf.meta['BKGSUB']`` automatically on ``lightkurve`` v2 ``TargetPixelFile``
    instances without mutating the read-only ``MappingProxyType``. When SPOC's
    standard ``BACKAPP`` card is present, it is synchronized to ``bool(bkgsub)``
    so FITS consumers never observe contradictory background-subtraction flags.
    """
    hdul = _unwrap_hdulist(target)
    if len(hdul) < 2:
        raise ValueError("Expected HDUList with at least primary HDU[0] and PIXELS HDU[1].")
    bkgsub_bool = bool(bkgsub)
    halocver_int = int(halocver)
    for idx in (0, 1):
        hdr = hdul[idx].header
        hdr["HALOCVER"] = (
            halocver_int,
            "Halo photometry FITS contract version",
        )
        hdr["BKGSUB"] = (
            bkgsub_bool,
            "Whether sky background has been subtracted",
        )
        if "BACKAPP" in hdr:
            hdr["BACKAPP"] = bkgsub_bool
    meta = getattr(target, "meta", None)
    meta_dict = getattr(meta, "_dict", None)
    if isinstance(meta_dict, dict):
        meta_dict["HALOCVER"] = halocver_int
        meta_dict["BKGSUB"] = bkgsub_bool
        if "BACKAPP" in meta_dict:
            meta_dict["BACKAPP"] = bkgsub_bool
    return hdul


def read_mundey_provenance(target: fits.HDUList | Any) -> dict[str, int | bool | None]:
    """Read and validate ``HALOCVER`` and ``BKGSUB`` across ``hdu[0]`` and ``hdu[1]``."""
    hdul = _unwrap_hdulist(target)
    halocver: int | None = None
    bkgsub: bool | None = None
    for idx in range(min(2, len(hdul))):
        hdr = hdul[idx].header
        if "HALOCVER" in hdr:
            raw_ver = hdr["HALOCVER"]
            if not isinstance(raw_ver, (int, np.integer)) or isinstance(
                raw_ver, (bool, np.bool_)
            ):
                raise ValueError(
                    f"Invalid HALOCVER header value in HDU[{idx}] (expected int): {raw_ver!r}."
                )
            ver_int = int(raw_ver)
            if halocver is None:
                halocver = ver_int
            elif halocver != ver_int:
                raise ValueError(
                    f"Conflicting HALOCVER values between HDU[0] ({halocver}) "
                    f"and HDU[{idx}] ({ver_int})."
                )
        if "BKGSUB" in hdr:
            raw_bkg = hdr["BKGSUB"]
            if not isinstance(raw_bkg, (bool, np.bool_)):
                raise ValueError(
                    f"Invalid BKGSUB header value in HDU[{idx}] (expected bool): {raw_bkg!r}."
                )
            bkg_bool = bool(raw_bkg)
            if bkgsub is None:
                bkgsub = bkg_bool
            elif bkgsub != bkg_bool:
                raise ValueError(
                    f"Conflicting BKGSUB values between HDU[0] ({bkgsub}) "
                    f"and HDU[{idx}] ({bkg_bool})."
                )
    return {"halocver": halocver, "bkgsub": bkgsub}


def write_pixel_quality(
    target: fits.HDUList | Any,
    mask: np.ndarray,
) -> fits.HDUList:
    """Write or update the spatial ``PIXEL_QUALITY`` ``ImageHDU`` (``uint16``, ``(ny, nx)``)."""
    hdul = _unwrap_hdulist(target)
    raw = np.asarray(mask)
    if raw.ndim != 2:
        raise ValueError(f"PIXEL_QUALITY mask must be 2D (ny, nx), got shape {raw.shape}.")
    if not (np.issubdtype(raw.dtype, np.integer) or np.issubdtype(raw.dtype, np.bool_)):
        raise TypeError(f"PIXEL_QUALITY mask must have integer or boolean dtype, got {raw.dtype}.")
    if np.issubdtype(raw.dtype, np.integer):
        if np.any(raw < 0):
            raise ValueError("PIXEL_QUALITY mask cannot contain negative bitmask values.")
        uint16_max = int(np.iinfo(np.uint16).max)
        if np.any(raw > uint16_max):
            raise ValueError(
                f"PIXEL_QUALITY mask contains values exceeding uint16 max ({uint16_max})."
            )

    inferred = _infer_cutout_shape(hdul)
    if inferred is not None:
        _, ny, nx = inferred
        if raw.shape != (ny, nx):
            raise ValueError(
                f"PIXEL_QUALITY shape {raw.shape} does not match TPF spatial shape {(ny, nx)}."
            )

    arr = np.asarray(raw, dtype=np.uint16)
    img_hdu = fits.ImageHDU(data=arr, name=PIXEL_QUALITY_EXTNAME)
    img_hdu.header["HALOCVER"] = (HALOCVER, "Halo photometry FITS contract version")

    if PIXEL_QUALITY_EXTNAME in hdul:
        idx = hdul.index_of(PIXEL_QUALITY_EXTNAME)
        hdul[idx] = img_hdu
    else:
        hdul.append(img_hdu)
    return hdul


def read_pixel_quality(target: fits.HDUList | Any) -> np.ndarray:
    """Read the ``PIXEL_QUALITY`` spatial bitmask as a ``uint16`` ``(ny, nx)`` array."""
    hdul = _unwrap_hdulist(target)
    if PIXEL_QUALITY_EXTNAME not in hdul:
        raise KeyError(f"Extension {PIXEL_QUALITY_EXTNAME!r} not found in HDUList.")
    return np.asarray(hdul[PIXEL_QUALITY_EXTNAME].data, dtype=np.uint16)


def _strip_smear_units(smear_err: Any) -> np.ndarray:
    """Convert an array or Astropy Quantity of ``e-/s`` to a plain ``float32`` ndarray."""
    if isinstance(smear_err, u.Quantity):
        if smear_err.unit == u.dimensionless_unscaled:
            raw = smear_err.value
        else:
            raw = smear_err.to_value(u.electron / u.s)
    else:
        raw = smear_err
    return np.asarray(raw, dtype=np.float32)


def append_smear_err_column(
    target: fits.HDUList | Any,
    smear_err: Any,
) -> fits.HDUList:
    """Append or update ``SMEAR_ERR`` (``float32``, ``(N_raw, nx)``) in HDU 1 ``PIXELS``.

    Storing ``SMEAR_ERR`` directly inside ``hdu[1]`` ensures ``tpf[mask]``,
    ``tpf.quality_mask``, and ``tpf.to_fits()`` slice and serialize column-coherent
    smear uncertainties in 1-to-1 lockstep with ``TIME`` and ``FLUX``.
    """
    hdul = _unwrap_hdulist(target)
    if len(hdul) < 2 or not isinstance(hdul[1], fits.BinTableHDU) or hdul[1].data is None:
        raise ValueError("HDU[1] must be a populated PIXELS BinTableHDU.")

    arr = _strip_smear_units(smear_err)
    if arr.ndim != 2:
        raise ValueError(f"SMEAR_ERR must be 2D (N_raw, nx), got shape {arr.shape}.")

    hdu1 = hdul[1]
    n_raw = len(hdu1.data)
    inferred = _infer_cutout_shape(hdul)
    if inferred is not None:
        _, _, nx = inferred
        expected = (n_raw, nx)
    else:
        nx = arr.shape[1]
        expected = (n_raw, nx)

    if arr.shape != expected:
        raise ValueError(
            f"SMEAR_ERR shape {arr.shape} does not match HDU[1] (N_raw, nx) = {expected}."
        )

    existing_names = list(hdu1.columns.names or ())
    if SMEAR_ERR_COLNAME in existing_names:
        col_shape = hdu1.data[SMEAR_ERR_COLNAME].shape
        if col_shape == arr.shape:
            hdu1.data[SMEAR_ERR_COLNAME] = arr
            return hdul
        kept_cols = [col for col in hdu1.columns if col.name != SMEAR_ERR_COLNAME]
    else:
        kept_cols = list(hdu1.columns)

    new_col = fits.Column(
        name=SMEAR_ERR_COLNAME,
        format=f"{nx}E",
        unit="e-/s",
        dim=f"({nx})",
        array=arr,
    )
    preserved_header = hdu1.header.copy()
    ext_name = hdu1.name or "PIXELS"
    new_hdu1 = fits.BinTableHDU.from_columns(
        fits.ColDefs([*kept_cols, new_col]),
        header=preserved_header,
        name=ext_name,
    )
    hdul[1] = new_hdu1
    return hdul


def read_smear_err_column(
    target: fits.HDUList | Any,
    *,
    quality_masked: bool = True,
) -> np.ndarray:
    """Read ``SMEAR_ERR`` as ``float32`` ``(N, nx)``, optionally applying ``quality_mask``."""
    hdul = _unwrap_hdulist(target)
    if len(hdul) < 2 or hdul[1].data is None:
        raise KeyError("HDU[1] PIXELS table is missing or empty.")
    names = hdul[1].data.names or ()
    if SMEAR_ERR_COLNAME not in names:
        raise KeyError(f"Column {SMEAR_ERR_COLNAME!r} not found in HDU[1].")

    raw = np.asarray(hdul[1].data[SMEAR_ERR_COLNAME], dtype=np.float32)
    if raw.ndim == 1:
        raw = raw[:, np.newaxis]
    if quality_masked and hasattr(target, "quality_mask"):
        return np.asarray(raw[target.quality_mask], dtype=np.float32)
    return raw


def _is_fred_path(path: Path) -> bool:
    """Return True if ``path`` or its resolved target resides under ``/fred``."""
    lexical = Path(os.path.abspath(os.fspath(path)))
    resolved = path.resolve(strict=False)
    return lexical.parts[:2] == ("/", "fred") or resolved.parts[:2] == ("/", "fred")


def resolve_cache_dir(
    cache_dir: str | os.PathLike[str] | None = None,
    *,
    strict_hpc: bool = False,
    create: bool = False,
) -> Path:
    """Resolve the ``mundey`` cache directory, guarding against ``/fred`` on HPC.

    Default resolution order when ``cache_dir`` is ``None``:
    1. ``$MUNDEY_CACHE_DIR`` if set and non-empty.
    2. ``$SCRATCH/.cache/tess_bright`` when ``$SCRATCH`` is set and non-empty.
    3. ``$XDG_CACHE_HOME/tess_bright`` when ``$XDG_CACHE_HOME`` is set and non-empty.
    4. ``~/.cache/tess_bright`` otherwise.

    When ``strict_hpc=True``, paths pointing under ``/fred`` (either directly or
    via symlinks) raise ``ValueError`` to protect Swinburne NT / OzSTAR project
    inode quotas.
    """
    if cache_dir is not None:
        candidate = Path(cache_dir).expanduser()
    elif os.environ.get("MUNDEY_CACHE_DIR", "").strip():
        candidate = Path(os.environ["MUNDEY_CACHE_DIR"].strip()).expanduser()
    elif os.environ.get("SCRATCH", "").strip():
        candidate = Path(os.environ["SCRATCH"].strip()).expanduser() / ".cache" / "tess_bright"
    elif os.environ.get("XDG_CACHE_HOME", "").strip():
        candidate = Path(os.environ["XDG_CACHE_HOME"].strip()).expanduser() / "tess_bright"
    else:
        candidate = Path.home() / ".cache" / "tess_bright"

    if strict_hpc and _is_fred_path(candidate):
        raise ValueError(
            f"Refusing cache directory {candidate!s} under /fred in strict_hpc=True mode "
            "(violates Swinburne NT/OzSTAR shared inode quota); use $SCRATCH/.cache/tess_bright."
        )

    if create:
        candidate.mkdir(parents=True, exist_ok=True)
    return candidate


def make_atomic_tmp_path(dst: str | os.PathLike[str]) -> Path:
    """Construct a collision-free sibling temp path ``{dst}.{hostname}.{pid}.{uuid}.tmp``."""
    dst_path = Path(dst)
    hostname = socket.gethostname().split(".")[0] or "host"
    pid = os.getpid()
    token = uuid.uuid4().hex
    return dst_path.parent / f"{dst_path.name}.{hostname}.{pid}.{token}.tmp"


def _fsync_file(path: Path) -> None:
    """Flush file data and metadata to stable storage."""
    try:
        fh = path.open("r+b")
    except OSError:
        fd = os.open(os.fspath(path), os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
        return
    with fh:
        os.fsync(fh.fileno())


def _fsync_parent_dir(parent: Path) -> None:
    """Best-effort fsync of the parent directory entry after atomic rename."""
    if not hasattr(os, "O_DIRECTORY"):
        return
    try:
        dir_fd = os.open(os.fspath(parent), os.O_RDONLY | os.O_DIRECTORY)
    except OSError:
        return
    try:
        os.fsync(dir_fd)
    except OSError:
        pass
    finally:
        os.close(dir_fd)


def atomic_publish(
    src_tmp: str | os.PathLike[str],
    dst: str | os.PathLike[str],
) -> Path:
    """Atomically publish ``src_tmp`` to ``dst`` without ``fcntl`` locks.

    Stages via a same-directory ``{dst}.{hostname}.{pid}.{uuid4().hex}.tmp``
    sibling when ``src_tmp`` resides in another directory or filesystem (such as
    Slurm ``$JOBFS`` staging into ``$SCRATCH``), calls ``os.fsync``, and executes
    atomic POSIX ``os.replace``.
    """
    src_path = Path(src_tmp)
    dst_path = Path(dst)
    if not src_path.exists():
        raise FileNotFoundError(f"Temporary source file does not exist: {src_path!s}")

    dst_path.parent.mkdir(parents=True, exist_ok=True)
    same_dir = os.path.abspath(os.fspath(src_path.parent)) == os.path.abspath(
        os.fspath(dst_path.parent)
    )

    if same_dir:
        _fsync_file(src_path)
        try:
            os.replace(src_path, dst_path)
            _fsync_parent_dir(dst_path.parent)
            return dst_path
        except OSError as exc:
            if exc.errno != errno.EXDEV:
                raise

    stage_tmp = make_atomic_tmp_path(dst_path)
    try:
        with src_path.open("rb") as src_fh, stage_tmp.open("wb") as dst_fh:
            shutil.copyfileobj(src_fh, dst_fh, length=1024 * 1024)
            dst_fh.flush()
            os.fsync(dst_fh.fileno())
        os.replace(stage_tmp, dst_path)
        _fsync_parent_dir(dst_path.parent)
    except BaseException:
        with contextlib.suppress(OSError):
            stage_tmp.unlink()
        raise

    with contextlib.suppress(OSError):
        src_path.unlink()
    return dst_path


@contextlib.contextmanager
def atomic_temp_path(dst: str | os.PathLike[str]) -> Iterator[Path]:
    """Yield a sibling ``{dst}.{hostname}.{pid}.{uuid}.tmp`` path and publish on exit."""
    dst_path = Path(dst)
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = make_atomic_tmp_path(dst_path)
    try:
        yield tmp_path
        if tmp_path.exists():
            atomic_publish(tmp_path, dst_path)
    except BaseException:
        with contextlib.suppress(OSError):
            tmp_path.unlink()
        raise
