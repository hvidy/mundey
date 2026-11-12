"""MAST S3 archive discovery, 4,096-row chunked collateral streaming, and ``-blf`` caching.

Supports dynamic resolution of static Camera/CCD calibration models (``2dblack``,
``flat``, ``linearity``) and both 2-minute (``-s_col.fits``) and 20-second
(``-a_fast-col.fits``) collateral pixel tables from local directories or
anonymous MAST S3 (``s3://stpubdata/tess/public/``), while enforcing strict
offline execution on HPC compute nodes when ``offline=True``.
"""

from __future__ import annotations

import contextlib
import os
import shutil
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from astropy.io import fits

from ._contract import atomic_publish, atomic_temp_path, make_atomic_tmp_path, resolve_cache_dir

S3_TESS_PUBLIC_ROOT: str = "stpubdata/tess/public"
DEFAULT_CHUNK_ROWS: int = 4096
CadenceKind = Literal["2min", "20s"]


@dataclass(frozen=True)
class StaticModelPaths:
    """Resolved local paths for static TESS Camera/CCD calibration models."""

    black2d_path: Path
    flat_path: Path
    linearity_path: Path


def infer_cadence_kind(
    *,
    tpf_path: str | os.PathLike[str] | None = None,
    hdu0_header: fits.Header | None = None,
    hdu1_header: fits.Header | None = None,
) -> CadenceKind:
    """Infer whether a TPF requires 2-minute (``2min``) or 20-second (``20s``) collateral data."""
    if tpf_path is not None:
        name = Path(tpf_path).name.lower()
        if "fast-tp" in name or "fast_tp" in name or "-a_fast" in name:
            return "20s"
        if "-s_tp" in name:
            return "2min"

    if hdu1_header is not None:
        if "TIMEDEL" in hdu1_header:
            try:
                dt_sec = float(hdu1_header["TIMEDEL"]) * 86400.0
                if 0.0 < dt_sec < 45.0:
                    return "20s"
            except (TypeError, ValueError):
                pass
        if "NUM_FRM" in hdu1_header:
            try:
                num_frm = int(hdu1_header["NUM_FRM"])
                crmiten = bool(
                    hdu0_header.get("CRMITEN", False)
                    if hdu0_header is not None and "CRMITEN" in hdu0_header
                    else hdu1_header.get("CRMITEN", True)
                )
                if num_frm <= 15 and not crmiten:
                    return "20s"
            except (TypeError, ValueError):
                pass
    return "2min"


def collateral_suffix(cadence: CadenceKind) -> str:
    """Return the MAST collateral filename tag for ``cadence``."""
    if cadence == "20s":
        return "a_fast-col"
    if cadence == "2min":
        return "s_col"
    raise ValueError(f"Unsupported cadence kind {cadence!r}; expected '2min' or '20s'.")


def _candidate_search_dirs(
    ddir: str | os.PathLike[str] | None,
    cache_dir: Path | None,
    *,
    subdir: str | None = None,
) -> list[Path]:
    """Return ordered local directories to search for reference and collateral files."""
    dirs: list[Path] = []
    if ddir is not None:
        dpath = Path(ddir).expanduser()
        if subdir:
            dirs.append(dpath / subdir)
        dirs.append(dpath)
    if cache_dir is not None:
        if subdir:
            dirs.append(cache_dir / subdir)
        dirs.append(cache_dir)
    deduped: list[Path] = []
    seen: set[str] = set()
    for d in dirs:
        key = os.path.abspath(os.fspath(d))
        if key not in seen:
            seen.add(key)
            deduped.append(d)
    return deduped


def _find_local_match(search_dirs: list[Path], patterns: list[str]) -> Path | None:
    """Find the first matching existing non-temporary file across ``search_dirs``."""
    for directory in search_dirs:
        if not directory.is_dir():
            continue
        for pat in patterns:
            if "*" not in pat and "?" not in pat:
                exact = directory / pat
                if exact.is_file():
                    return exact
            else:
                matches = sorted(
                    p
                    for p in directory.glob(pat)
                    if p.is_file() and not p.name.endswith(".tmp")
                )
                if matches:
                    return matches[0]
    return None


def _open_anon_s3_filesystem(*, s3_fs: Any = None) -> Any:
    """Open an anonymous ``s3fs`` filesystem handle lazily."""
    if s3_fs is not None:
        return s3_fs
    try:
        import fsspec
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "Anonymous MAST S3 downloads require 'fsspec' and 's3fs'. "
            "Install them via `pip install 'mundey[hpc]'` or pre-stage files locally."
        ) from exc
    return fsspec.filesystem("s3", anon=True)


def _download_s3_to_path(
    fs: Any,
    s3_key: str,
    dst: Path,
    *,
    stage_dir: Path | None = None,
) -> Path:
    """Stream an S3 object to local staging and publish atomically via ``atomic_publish``."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if stage_dir is not None:
        stage_dir.mkdir(parents=True, exist_ok=True)
        tmp_path = make_atomic_tmp_path(stage_dir / dst.name)
    else:
        tmp_path = make_atomic_tmp_path(dst)

    try:
        with fs.open(s3_key, "rb") as src_fh, tmp_path.open("wb") as dst_fh:
            shutil.copyfileobj(src_fh, dst_fh, length=4 * 1024 * 1024)
            dst_fh.flush()
            os.fsync(dst_fh.fileno())
        return atomic_publish(tmp_path, dst)
    except BaseException:
        with contextlib.suppress(OSError):
            tmp_path.unlink()
        raise


def _guard_offline(offline: bool, description: str, searched: list[Path]) -> None:
    """Raise ``FileNotFoundError`` immediately before any network access when ``offline=True``."""
    if offline:
        locs = ", ".join(str(d) for d in searched) if searched else "<none>"
        raise FileNotFoundError(
            f"offline=True: refusing network access for missing {description} "
            f"(searched local directories: {locs}). Run `mundey prefetch` first."
        )


def resolve_static_models(
    camera: int,
    ccd: int,
    *,
    ddir: str | os.PathLike[str] | None = None,
    cache_dir: Path | None = None,
    offline: bool = False,
    stage_dir: Path | None = None,
    s3_fs: Any = None,
    verbose: bool = False,
) -> StaticModelPaths:
    """Locate or fetch ``2dblack``, ``flat``, and ``linearity`` models for ``(camera, ccd)``."""
    if not (1 <= int(camera) <= 4 and 1 <= int(ccd) <= 4):
        raise ValueError(f"Expected camera, ccd in 1..4, got camera={camera}, ccd={ccd}.")

    search_dirs = _candidate_search_dirs(ddir, cache_dir, subdir="models")

    b2d_patterns = [
        f"tess2018323-{camera}-{ccd}-2dblack.fits",
        f"tess*-{camera}-{ccd}-2dblack.fits",
        f"*{camera}-{ccd}-2dblack.fits",
    ]
    flat_patterns = [
        f"tess2018323-{camera}-{ccd}-flat.fits",
        f"tess*-{camera}-{ccd}-flat.fits",
        f"*{camera}-{ccd}-flat.fits",
    ]
    lin_patterns = [
        "tess2018143203310-41006_100-linearity.xml",
        "tess*-linearity.xml",
        "*-linearity.xml",
    ]

    black2d_path = _find_local_match(search_dirs, b2d_patterns)
    flat_path = _find_local_match(search_dirs, flat_patterns)
    linearity_path = _find_local_match(search_dirs, lin_patterns)

    if black2d_path is not None and flat_path is not None and linearity_path is not None:
        return StaticModelPaths(
            black2d_path=black2d_path,
            flat_path=flat_path,
            linearity_path=linearity_path,
        )

    missing_labels: list[str] = []
    if black2d_path is None:
        missing_labels.append(f"2dblack ({camera}-{ccd})")
    if flat_path is None:
        missing_labels.append(f"flat ({camera}-{ccd})")
    if linearity_path is None:
        missing_labels.append("linearity.xml")
    _guard_offline(offline, ", ".join(missing_labels), search_dirs)

    if cache_dir is not None:
        target_models_dir = cache_dir / "models"
    elif ddir is not None:
        target_models_dir = Path(ddir).expanduser()
    else:
        target_models_dir = resolve_cache_dir(create=True) / "models"
    fs = _open_anon_s3_filesystem(s3_fs=s3_fs)
    s3_models = sorted(fs.ls(f"{S3_TESS_PUBLIC_ROOT}/models"))

    def _pick_s3_model(suffix_candidates: tuple[str, ...], preferred: str) -> str:
        for key in s3_models:
            if Path(key).name == preferred:
                return key
        matches = [k for k in s3_models if Path(k).name.endswith(suffix_candidates)]
        if not matches:
            raise FileNotFoundError(
                f"Could not find static model matching {suffix_candidates} under "
                f"s3://{S3_TESS_PUBLIC_ROOT}/models."
            )
        return matches[0]

    if black2d_path is None:
        key = _pick_s3_model(
            (f"-{camera}-{ccd}-2dblack.fits",),
            f"tess2018323-{camera}-{ccd}-2dblack.fits",
        )
        dst = target_models_dir / Path(key).name
        if verbose:
            print(f"Downloading S3 static model s3://{key} -> {dst}")
        black2d_path = _download_s3_to_path(fs, key, dst, stage_dir=stage_dir)

    if flat_path is None:
        key = _pick_s3_model(
            (f"-{camera}-{ccd}-flat.fits",),
            f"tess2018323-{camera}-{ccd}-flat.fits",
        )
        dst = target_models_dir / Path(key).name
        if verbose:
            print(f"Downloading S3 static model s3://{key} -> {dst}")
        flat_path = _download_s3_to_path(fs, key, dst, stage_dir=stage_dir)

    if linearity_path is None:
        key = _pick_s3_model(
            ("-linearity.xml",),
            "tess2018143203310-41006_100-linearity.xml",
        )
        dst = target_models_dir / Path(key).name
        if verbose:
            print(f"Downloading S3 static model s3://{key} -> {dst}")
        linearity_path = _download_s3_to_path(fs, key, dst, stage_dir=stage_dir)

    return StaticModelPaths(
        black2d_path=black2d_path,
        flat_path=flat_path,
        linearity_path=linearity_path,
    )


def collateral_blf_cache_path(
    base_dir: Path,
    *,
    sector: int,
    camera: int,
    ccd: int,
    output: str,
    kind: Literal["smrow", "tvcol"],
    cadence: CadenceKind = "2min",
    smear_mode: str = "standard",
    tpf_parts: list[str] | None = None,
    col_stem: str | None = None,
) -> Path:
    """Construct the canonical ``-blf.fits`` cache path for ``tvcol`` or ``smrow``."""
    out_lower = output.lower()
    col_tag = collateral_suffix(cadence)

    if kind == "smrow" and smear_mode == "alternate":
        if tpf_parts is not None and len(tpf_parts) >= 4:
            prefix, sec_str, tic_str, run_str = (
                tpf_parts[0],
                tpf_parts[1],
                tpf_parts[2],
                tpf_parts[3],
            )
            return base_dir / (
                f"{prefix}-{sec_str}-smrow-{camera}-{ccd}-"
                f"{out_lower}-{tic_str}-{run_str}-alt-{col_tag}-blf.fits"
            )
        stem = col_stem or f"tess-s{sector:04d}-smrow-{camera}-{ccd}-{out_lower}-{col_tag}"
        return base_dir / f"{stem}-alt-blf.fits"

    if col_stem is not None:
        return base_dir / f"{col_stem}-blf.fits"

    if tpf_parts is not None and len(tpf_parts) >= 4:
        return base_dir / (
            f"{tpf_parts[0]}-{tpf_parts[1]}-{kind}-{camera}-{ccd}-"
            f"{out_lower}-{tpf_parts[3]}-{col_tag}-blf.fits"
        )

    return base_dir / f"tess-s{sector:04d}-{kind}-{camera}-{ccd}-{out_lower}-{col_tag}-blf.fits"


def find_cached_blf(
    *,
    sector: int,
    camera: int,
    ccd: int,
    output: str,
    kind: Literal["smrow", "tvcol"],
    cadence: CadenceKind = "2min",
    smear_mode: str = "standard",
    ddir: str | os.PathLike[str] | None = None,
    cache_dir: Path | None = None,
    tpf_parts: list[str] | None = None,
    ticid: int | None = None,
) -> Path | None:
    """Locate an existing ``-blf.fits`` cache locally across ``ddir`` and ``cache_dir``."""
    out_lower = output.lower()
    col_tag = collateral_suffix(cadence)
    sec_sub = f"s{sector:04d}"
    search_dirs = _candidate_search_dirs(ddir, cache_dir, subdir=sec_sub)

    patterns: list[str] = []
    if kind == "smrow" and smear_mode == "alternate":
        if tpf_parts is not None and len(tpf_parts) >= 4:
            patterns.append(
                f"{tpf_parts[0]}-{tpf_parts[1]}-smrow-{camera}-{ccd}-"
                f"{out_lower}-{tpf_parts[2]}-{tpf_parts[3]}-alt-{col_tag}-blf.fits"
            )
        if ticid is not None:
            patterns.append(
                f"tess*-s{sector:04d}-smrow-{camera}-{ccd}-"
                f"{out_lower}-{int(ticid):016d}-*-alt-{col_tag}-blf.fits"
            )
    else:
        if tpf_parts is not None and len(tpf_parts) >= 4:
            patterns.append(
                f"{tpf_parts[0]}-{tpf_parts[1]}-{kind}-{camera}-{ccd}-"
                f"{out_lower}-{tpf_parts[3]}-{col_tag}-blf.fits"
            )
        patterns.append(
            f"tess*-s{sector:04d}-{kind}-{camera}-{ccd}-{out_lower}-*-{col_tag}-blf.fits"
        )
        patterns.append(
            f"tess-s{sector:04d}-{kind}-{camera}-{ccd}-{out_lower}-{col_tag}-blf.fits"
        )

    return _find_local_match(search_dirs, patterns)


def resolve_collateral_file(
    *,
    sector: int,
    camera: int,
    ccd: int,
    output: str,
    kind: Literal["smrow", "tvcol"],
    cadence: CadenceKind = "2min",
    ddir: str | os.PathLike[str] | None = None,
    cache_dir: Path | None = None,
    offline: bool = False,
    stage_dir: Path | None = None,
    keep_raw_in_stage_only: bool = False,
    tpf_parts: list[str] | None = None,
    s3_fs: Any = None,
    verbose: bool = False,
) -> Path:
    """Resolve a raw ``_col.fits`` / ``_fast-col.fits`` file locally or from MAST S3."""
    out_lower = output.lower()
    col_tag = collateral_suffix(cadence)
    sec_sub = f"s{sector:04d}"
    search_dirs = _candidate_search_dirs(ddir, cache_dir, subdir=sec_sub)
    if stage_dir is not None:
        search_dirs = [stage_dir, *search_dirs]

    patterns: list[str] = []
    if tpf_parts is not None and len(tpf_parts) >= 4:
        patterns.append(
            f"{tpf_parts[0]}-{tpf_parts[1]}-{kind}-{camera}-{ccd}-"
            f"{out_lower}-{tpf_parts[3]}-{col_tag}.fits"
        )
    patterns.append(f"tess*-s{sector:04d}-{kind}-{camera}-{ccd}-{out_lower}-*-{col_tag}.fits")

    found = _find_local_match(search_dirs, patterns)
    if found is not None:
        return found

    _guard_offline(
        offline,
        f"{kind} collateral ({sec_sub}, cam={camera}, ccd={ccd}, out={output}, {col_tag})",
        search_dirs,
    )

    fs = _open_anon_s3_filesystem(s3_fs=s3_fs)
    glob_pat = (
        f"{S3_TESS_PUBLIC_ROOT}/ffi/s{sector:04d}/*/*/{camera}-{ccd}/"
        f"tess*-s{sector:04d}-{kind}-{camera}-{ccd}-{out_lower}-*-{col_tag}.fits"
    )
    matches = sorted(fs.glob(glob_pat))
    if not matches:
        raise FileNotFoundError(
            f"Collateral file not found on MAST S3 matching s3://{glob_pat}."
        )
    s3_key = matches[0]
    filename = Path(s3_key).name

    if keep_raw_in_stage_only and stage_dir is not None:
        dst = stage_dir / filename
    elif cache_dir is not None:
        dst = cache_dir / sec_sub / filename
    elif ddir is not None:
        dst = Path(ddir).expanduser() / filename
    else:
        dst = resolve_cache_dir(create=True) / sec_sub / filename

    if verbose:
        print(f"Downloading S3 collateral file s3://{s3_key} -> {dst}")
    return _download_s3_to_path(fs, s3_key, dst, stage_dir=stage_dir)


def read_collateral_headers_and_index(
    col_path: str | os.PathLike[str],
) -> tuple[fits.Header, fits.Header, np.ndarray, np.ndarray]:
    """Read ``(hdr0, hdr1, time, cadenceno)`` without loading multi-GB pixel columns."""
    with fits.open(col_path, memmap=True, uint=False) as hdul:
        hdr0 = hdul[0].header.copy()
        hdr1 = hdul[1].header.copy()
        data = hdul[1].data
        time = np.asarray(data["TIME"], dtype=np.float64)
        cadenceno = np.asarray(data["CADENCENO"], dtype=np.int64)
    return hdr0, hdr1, time, cadenceno


def iter_collateral_chunks(
    col_path: str | os.PathLike[str],
    raw_colname: str,
    *,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
    row_indices: np.ndarray | None = None,
) -> Iterator[tuple[int, int, np.ndarray, np.ndarray, np.ndarray]]:
    """Stream ``TSCALn``/``TZEROn`` collateral tables in ``chunk_rows`` slices (< 2 GB RSS).

    Slices ``hdul[1].data[t0:t1]`` by row range *before* column extraction so
    Astropy never materializes the full-sector scaled ``float64`` table in RAM,
    and masks negative sentinel values (``-1`` / ``TNULLn``) to ``NaN`` in ``float64``.

    When ``row_indices`` is provided (sorted or arbitrary 1D ``int64`` indices into
    the collateral table), only chunks intersecting ``row_indices`` are read.
    """
    if chunk_rows <= 0:
        raise ValueError(f"chunk_rows must be positive, got {chunk_rows}.")

    with fits.open(col_path, memmap=True, uint=False) as hdul:
        hdu1 = hdul[1]
        data = hdu1.data
        n_total = len(data)

        # Inspect TNULL (and TSCAL/TZERO scaling) for raw_colname if present
        scaled_null: float | None = None
        for idx, col_def in enumerate(hdu1.columns, start=1):
            if col_def.name == raw_colname:
                tnull_val = col_def.null
                if tnull_val is None and f"TNULL{idx}" in hdu1.header:
                    with contextlib.suppress(TypeError, ValueError):
                        tnull_val = int(hdu1.header[f"TNULL{idx}"])
                if tnull_val is not None:
                    bscale = float(col_def.bscale) if col_def.bscale is not None else 1.0
                    bzero = float(col_def.bzero) if col_def.bzero is not None else 0.0
                    scaled_null = bscale * float(tnull_val) + bzero
                break

        if row_indices is not None:
            needed = np.unique(np.asarray(row_indices, dtype=np.int64))
        else:
            needed = None

        for t0 in range(0, n_total, chunk_rows):
            t1 = min(t0 + chunk_rows, n_total)
            if needed is not None:
                pos0 = int(np.searchsorted(needed, t0, side="left"))
                if pos0 >= len(needed) or int(needed[pos0]) >= t1:
                    continue

            chunk_rec = data[t0:t1]
            time_chunk = np.asarray(chunk_rec["TIME"], dtype=np.float64)
            cad_chunk = np.asarray(chunk_rec["CADENCENO"], dtype=np.int64)
            raw_arr = chunk_rec[raw_colname]
            raw_f64 = np.asarray(raw_arr, dtype=np.float64)

            bad = (raw_f64 < 0.0) | (~np.isfinite(raw_f64))
            if scaled_null is not None:
                bad = bad | (raw_f64 == scaled_null)
            if np.any(bad):
                raw_f64 = raw_f64.copy()
                raw_f64[bad] = np.nan

            yield t0, t1, time_chunk, cad_chunk, raw_f64


def write_onedblack_blf(
    dst_path: str | os.PathLike[str],
    *,
    time: np.ndarray,
    cadenceno: np.ndarray,
    black1d: np.ndarray,
    black1d_err: np.ndarray | None = None,
    n_eff: float | None = None,
    time_unit: str = "BJD - 2457000, days",
    stage_dir: Path | None = None,
) -> Path:
    """Atomically write a ``tvcol-...-blf.fits`` 1D-black cache with ``CADENCENO`` (``'J'``)."""
    dst = Path(dst_path)
    time_arr = np.asarray(time, dtype=np.float64)
    cad_arr = np.asarray(cadenceno, dtype=np.int32)
    b1d_arr = np.asarray(black1d, dtype=np.float64)
    nrows = b1d_arr.shape[1]

    cols = [
        fits.Column(
            name="TIME",
            format="D",
            unit=time_unit,
            disp="D14.7",
            array=time_arr,
        ),
        fits.Column(
            name="CADENCENO",
            format="J",
            disp="I10",
            array=cad_arr,
        ),
        fits.Column(
            name="1DBLACK",
            format=f"{nrows}D",
            unit="count",
            disp="F14.7",
            dim=f"({nrows})",
            array=b1d_arr,
        ),
    ]
    if black1d_err is not None:
        cols.append(
            fits.Column(
                name="1DBLACK_ERR",
                format=f"{nrows}E",
                unit="count",
                disp="E14.7",
                dim=f"({nrows})",
                array=np.asarray(black1d_err, dtype=np.float32),
            )
        )

    hdu1 = fits.BinTableHDU.from_columns(fits.ColDefs(cols), name="1DBLACK")
    hdu1.header["BLFVER"] = (2, "Mundey -blf cache schema version")
    if n_eff is not None:
        hdu1.header["NEFF"] = (float(n_eff), "Effective co-added frames")
    hdul = fits.HDUList([fits.PrimaryHDU(), hdu1])

    if stage_dir is not None:
        stage_dir.mkdir(parents=True, exist_ok=True)
        tmp = make_atomic_tmp_path(stage_dir / dst.name)
        try:
            hdul.writeto(tmp, overwrite=True)
            return atomic_publish(tmp, dst)
        except BaseException:
            with contextlib.suppress(OSError):
                tmp.unlink()
            raise

    with atomic_temp_path(dst) as tmp:
        hdul.writeto(tmp, overwrite=True)
    return dst


def write_smear_blf(
    dst_path: str | os.PathLike[str],
    *,
    time: np.ndarray,
    cadenceno: np.ndarray,
    smear: np.ndarray,
    smear_err: np.ndarray | None = None,
    smear_mode: str = "standard",
    cutout_bounds: tuple[int, int, int, int] | None = None,
    int_time: float | None = None,
    time_unit: str = "BJD - 2457000, days",
    stage_dir: Path | None = None,
) -> Path:
    """Atomically write a ``smrow-...-blf.fits`` cache with ``CADENCENO`` (``'J'``)."""
    dst = Path(dst_path)
    time_arr = np.asarray(time, dtype=np.float64)
    cad_arr = np.asarray(cadenceno, dtype=np.int32)
    sm_arr = np.asarray(smear, dtype=np.float64)
    ncols = sm_arr.shape[1]

    cols = [
        fits.Column(
            name="TIME",
            format="D",
            unit=time_unit,
            disp="D14.7",
            array=time_arr,
        ),
        fits.Column(
            name="CADENCENO",
            format="J",
            disp="I10",
            array=cad_arr,
        ),
        fits.Column(
            name="SMEAR",
            format=f"{ncols}D",
            unit="count",
            disp="F14.7",
            dim=f"({ncols})",
            array=sm_arr,
        ),
    ]
    if smear_err is not None:
        cols.append(
            fits.Column(
                name="SMEAR_ERR",
                format=f"{ncols}E",
                unit="count",
                disp="E14.7",
                dim=f"({ncols})",
                array=np.asarray(smear_err, dtype=np.float32),
            )
        )

    hdu1 = fits.BinTableHDU.from_columns(fits.ColDefs(cols), name="SMEAR")
    hdu1.header["BLFVER"] = (2, "Mundey -blf cache schema version")
    hdu1.header["SMEARMOD"] = (smear_mode, "Mundey smear estimation mode")
    if int_time is not None:
        hdu1.header["INT_TIME"] = (float(int_time), "Single-frame integration time (s)")
    if cutout_bounds is not None:
        crval1p, cutnx, crval2p, cutny = cutout_bounds
        hdu1.header["CRVAL1P"] = (int(crval1p), "Cutout start CCD column")
        hdu1.header["CUTNX"] = (int(cutnx), "Cutout column width")
        hdu1.header["CRVAL2P"] = (int(crval2p), "Cutout start CCD row")
        hdu1.header["CUTNY"] = (int(cutny), "Cutout row height")

    hdul = fits.HDUList([fits.PrimaryHDU(), hdu1])

    if stage_dir is not None:
        stage_dir.mkdir(parents=True, exist_ok=True)
        tmp = make_atomic_tmp_path(stage_dir / dst.name)
        try:
            hdul.writeto(tmp, overwrite=True)
            return atomic_publish(tmp, dst)
        except BaseException:
            with contextlib.suppress(OSError):
                tmp.unlink()
            raise

    with atomic_temp_path(dst) as tmp:
        hdul.writeto(tmp, overwrite=True)
    return dst


def prefetch_collateral(
    *,
    sector: int,
    camera: int,
    ccd: int,
    outputs: tuple[str, ...] = ("A", "B", "C", "D"),
    cadences: tuple[CadenceKind, ...] = ("2min",),
    cache_dir: str | os.PathLike[str] | None = None,
    ddir: str | os.PathLike[str] | None = None,
    stage_dir: str | os.PathLike[str] | None = None,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
    offline: bool = False,
    strict_hpc: bool = False,
    cleanup_staged_raw: bool = True,
    s3_fs: Any = None,
    verbose: bool = False,
) -> list[Path]:
    """Prefetch static models and build ``-blf.fits`` caches in ``chunk_rows`` slices.

    Streams raw ``_col.fits`` / ``_fast-col.fits`` tables in ``chunk_rows`` (default
    4,096) cadence chunks through the ``@in_float64`` calibration kernels, keeping
    memory below 2 GB RSS and publishing ``-blf.fits`` files locklessly via
    :func:`~mundey._contract.atomic_publish`.
    """
    from .calibration import build_sector_output_blf_caches

    resolved_cache = resolve_cache_dir(cache_dir, strict_hpc=strict_hpc, create=True)
    resolved_stage: Path | None
    if stage_dir is not None:
        resolved_stage = Path(stage_dir).expanduser()
    elif os.environ.get("JOBFS", "").strip():
        resolved_stage = Path(os.environ["JOBFS"].strip()).expanduser()
    else:
        resolved_stage = None

    models = resolve_static_models(
        camera,
        ccd,
        ddir=ddir,
        cache_dir=resolved_cache,
        offline=offline,
        stage_dir=resolved_stage,
        s3_fs=s3_fs,
        verbose=verbose,
    )

    published: list[Path] = [models.black2d_path, models.flat_path, models.linearity_path]
    normalized_outputs = tuple(o.upper() for o in outputs)
    for out in normalized_outputs:
        if out not in ("A", "B", "C", "D"):
            raise ValueError(f"Invalid CCD output {out!r}; expected 'A', 'B', 'C', or 'D'.")

    for cad in cadences:
        for out in normalized_outputs:
            tv_blf, sm_blf = build_sector_output_blf_caches(
                sector=sector,
                camera=camera,
                ccd=ccd,
                output=out,
                cadence=cad,
                models=models,
                cache_dir=resolved_cache,
                ddir=ddir,
                stage_dir=resolved_stage,
                chunk_rows=chunk_rows,
                offline=offline,
                cleanup_staged_raw=cleanup_staged_raw,
                s3_fs=s3_fs,
                verbose=verbose,
            )
            published.extend([tv_blf, sm_blf])

    return published
