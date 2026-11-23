"""MAST S3 archive discovery, 4,096-row chunked collateral streaming, and ``-blf`` caching.

Pure I/O and archive discovery layer (zero imports from :mod:`mundey.calibration`).
Supports dynamic resolution of static Camera/CCD calibration models (``2dblack``,
``flat``, ``linearity``) and both 2-minute (``-s_col.fits``) and 20-second
(``-a_fast-col.fits``) collateral pixel tables from local directories or anonymous
MAST S3 (``s3://stpubdata/tess/public/``), while enforcing strict offline execution
on HPC compute nodes when ``offline=True``.
"""

from __future__ import annotations

import contextlib
import os
import shutil
from collections.abc import Callable, Iterator
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


@dataclass(frozen=True)
class CollateralSpec:
    """Identifying metadata for a TESS sector/camera/CCD/output collateral stream."""

    sector: int
    camera: int
    ccd: int
    output: str
    cadence: CadenceKind = "2min"
    tpf_parts: tuple[str, ...] | None = None
    ticid: int | None = None

    @property
    def output_lower(self) -> str:
        return self.output.lower()

    @property
    def sector_subdir(self) -> str:
        return f"s{self.sector:04d}"

    @property
    def col_tag(self) -> str:
        return collateral_suffix(self.cadence)


def infer_cadence_kind(
    *,
    tpf_path: str | os.PathLike[str] | Any = None,
    hdu0_header: fits.Header | None = None,
    hdu1_header: fits.Header | None = None,
) -> CadenceKind:
    """Infer whether a TPF requires 2-minute (``2min``) or 20-second (``20s``) collateral data."""
    if isinstance(tpf_path, str | os.PathLike) and os.fspath(tpf_path):
        name = Path(tpf_path).name.lower()
        if "fast-tp" in name or "fast_tp" in name or "-a_fast" in name:
            return "20s"
        if "-s_tp" in name:
            return "2min"

    for hdr in (hdu1_header, hdu0_header):
        if hdr is not None and "TIMEDEL" in hdr:
            with contextlib.suppress(TypeError, ValueError):
                dt_sec = float(hdr["TIMEDEL"]) * 86400.0
                if 0.0 < dt_sec < 45.0:
                    return "20s"
    for hdr in (hdu1_header, hdu0_header):
        if hdr is not None and "NUM_FRM" in hdr:
            with contextlib.suppress(TypeError, ValueError):
                num_frm = int(hdr["NUM_FRM"])
                crmiten = bool(
                    hdu0_header.get("CRMITEN", False)
                    if hdu0_header is not None and "CRMITEN" in hdu0_header
                    else (
                        hdu1_header.get("CRMITEN", True)
                        if hdu1_header is not None
                        else True
                    )
                )
                if num_frm <= 15 and not crmiten:
                    return "20s"
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
    for base in (Path(ddir).expanduser() if ddir is not None else None, cache_dir):
        if base is None:
            continue
        if subdir:
            dirs.append(base / subdir)
        dirs.append(base)
    deduped: list[Path] = []
    seen: set[str] = set()
    for d in dirs:
        key = os.path.abspath(os.fspath(d))
        if key not in seen:
            seen.add(key)
            deduped.append(d)
    return deduped


def _find_all_local_matches(
    search_dirs: list[Path],
    patterns: list[str],
    *,
    predicate: Callable[[Path], bool] | None = None,
) -> list[Path]:
    """Return all matching existing non-temporary files across ``search_dirs`` in priority order."""
    results: list[Path] = []
    seen: set[str] = set()
    for directory in search_dirs:
        if not directory.is_dir():
            continue
        for pat in patterns:
            if "*" not in pat and "?" not in pat:
                exact = directory / pat
                candidates = [exact] if exact.is_file() else []
            else:
                candidates = sorted(
                    (
                        p
                        for p in directory.glob(pat)
                        if p.is_file() and not p.name.endswith(".tmp")
                    ),
                    reverse=True,
                )
            for p in candidates:
                if p.name.endswith(".tmp"):
                    continue
                if predicate is not None and not predicate(p):
                    continue
                key = os.path.abspath(os.fspath(p))
                if key not in seen:
                    seen.add(key)
                    results.append(p)
    return results


def _find_local_match(
    search_dirs: list[Path],
    patterns: list[str],
    *,
    predicate: Callable[[Path], bool] | None = None,
) -> Path | None:
    matches = _find_all_local_matches(search_dirs, patterns, predicate=predicate)
    return matches[0] if matches else None


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
    model_specs = (
        (
            f"2dblack ({camera}-{ccd})",
            [
                f"tess2018323-{camera}-{ccd}-2dblack.fits",
                f"tess*-{camera}-{ccd}-2dblack.fits",
                f"*{camera}-{ccd}-2dblack.fits",
            ],
            f"-{camera}-{ccd}-2dblack.fits",
            f"tess2018323-{camera}-{ccd}-2dblack.fits",
        ),
        (
            f"flat ({camera}-{ccd})",
            [
                f"tess2018323-{camera}-{ccd}-flat.fits",
                f"tess*-{camera}-{ccd}-flat.fits",
                f"*{camera}-{ccd}-flat.fits",
            ],
            f"-{camera}-{ccd}-flat.fits",
            f"tess2018323-{camera}-{ccd}-flat.fits",
        ),
        (
            "linearity.xml",
            [
                "tess2018143203310-41006_100-linearity.xml",
                "tess*-linearity.xml",
                "*-linearity.xml",
            ],
            "-linearity.xml",
            "tess2018143203310-41006_100-linearity.xml",
        ),
    )

    found = [_find_local_match(search_dirs, pats) for _, pats, _, _ in model_specs]
    if all(p is not None for p in found):
        return StaticModelPaths(found[0], found[1], found[2])  # type: ignore[arg-type]

    missing = [label for (label, _, _, _), p in zip(model_specs, found, strict=True) if p is None]
    _guard_offline(offline, ", ".join(missing), search_dirs)

    if cache_dir is not None:
        target_models_dir = cache_dir / "models"
    elif ddir is not None:
        target_models_dir = Path(ddir).expanduser()
    else:
        target_models_dir = (
            resolve_cache_dir(env_var="MUNDEY_CACHE_DIR", create=True) / "models"
        )

    fs = _open_anon_s3_filesystem(s3_fs=s3_fs)
    s3_models = sorted(fs.ls(f"{S3_TESS_PUBLIC_ROOT}/models"))

    resolved: list[Path] = []
    for (_, _, suffix, preferred), local_path in zip(model_specs, found, strict=True):
        if local_path is not None:
            resolved.append(local_path)
            continue
        key = next((k for k in s3_models if Path(k).name == preferred), None)
        if key is None:
            matches = [k for k in s3_models if Path(k).name.endswith(suffix)]
            if not matches:
                raise FileNotFoundError(
                    f"Could not find static model ending in {suffix!r} under "
                    f"s3://{S3_TESS_PUBLIC_ROOT}/models."
                )
            key = matches[-1]
        dst = target_models_dir / Path(key).name
        if verbose:
            print(f"Downloading S3 static model s3://{key} -> {dst}")
        resolved.append(_download_s3_to_path(fs, key, dst, stage_dir=stage_dir))

    return StaticModelPaths(resolved[0], resolved[1], resolved[2])


def collateral_blf_cache_path(
    base_dir: Path,
    spec: CollateralSpec,
    *,
    kind: Literal["smrow", "tvcol"],
    smear_mode: str = "standard",
    col_stem: str | None = None,
) -> Path:
    """Construct the canonical ``-blf.fits`` cache path for ``tvcol`` or ``smrow``."""
    parts = spec.tpf_parts
    out_l = spec.output_lower
    tag = spec.col_tag

    if kind == "smrow" and smear_mode == "alternate":
        if parts is not None and len(parts) >= 4:
            return base_dir / (
                f"{parts[0]}-{parts[1]}-smrow-{spec.camera}-{spec.ccd}-"
                f"{out_l}-{parts[2]}-{parts[3]}-alt-{tag}-blf.fits"
            )
        stem = col_stem or f"tess-s{spec.sector:04d}-smrow-{spec.camera}-{spec.ccd}-{out_l}-{tag}"
        return base_dir / f"{stem}-alt-blf.fits"

    if col_stem is not None:
        return base_dir / f"{col_stem}-blf.fits"

    if parts is not None and len(parts) >= 4:
        return base_dir / (
            f"{parts[0]}-{parts[1]}-{kind}-{spec.camera}-{spec.ccd}-"
            f"{out_l}-{parts[3]}-{tag}-blf.fits"
        )

    return (
        base_dir
        / f"tess-s{spec.sector:04d}-{kind}-{spec.camera}-{spec.ccd}-{out_l}-{tag}-blf.fits"
    )


def _is_alt_blf_name(path: Path) -> bool:
    name = path.name
    return "-alt-" in name or name.endswith("-alt-blf.fits")


def find_all_cached_blf(
    spec: CollateralSpec,
    *,
    kind: Literal["smrow", "tvcol"],
    smear_mode: str = "standard",
    ddir: str | os.PathLike[str] | None = None,
    cache_dir: Path | None = None,
) -> list[Path]:
    """Return all candidate ``-blf.fits`` cache paths across ``ddir`` and ``cache_dir``."""
    out_l = spec.output_lower
    tag = spec.col_tag
    sec_s = spec.sector_subdir
    parts = spec.tpf_parts
    search_dirs = _candidate_search_dirs(ddir, cache_dir, subdir=sec_s)

    patterns: list[str] = []
    if kind == "smrow" and smear_mode == "alternate":
        if parts is not None and len(parts) >= 4:
            patterns.extend(
                (
                    f"tess*-{sec_s}-smrow-{spec.camera}-{spec.ccd}-"
                    f"{out_l}-{parts[2]}-{parts[3]}-alt-{tag}-blf.fits",
                    f"{parts[0]}-{parts[1]}-smrow-{spec.camera}-{spec.ccd}-"
                    f"{out_l}-{parts[2]}-{parts[3]}-alt-{tag}-blf.fits",
                )
            )
        if spec.ticid is not None:
            tic = f"{int(spec.ticid):016d}"
            patterns.extend(
                (
                    f"tess*-{sec_s}-smrow-{spec.camera}-{spec.ccd}-{out_l}-{tic}-*-alt-{tag}-blf.fits",
                    f"tess-{sec_s}-smrow-{spec.camera}-{spec.ccd}-{out_l}-{tic}-{tag}-alt-blf.fits",
                )
            )
            if spec.cadence == "2min":
                patterns.append(
                    f"tess-{sec_s}-smrow-{spec.camera}-{spec.ccd}-{out_l}-{tic}-alt-blf.fits"
                )
        patterns.extend(
            (
                f"tess*-{sec_s}-smrow-{spec.camera}-{spec.ccd}-{out_l}-*-alt-{tag}-blf.fits",
                f"tess*-{sec_s}-smrow-{spec.camera}-{spec.ccd}-{out_l}-*-{tag}-alt-blf.fits",
            )
        )
        predicate: Callable[[Path], bool] = _is_alt_blf_name
    else:
        if parts is not None and len(parts) >= 4:
            patterns.extend(
                (
                    f"tess*-{sec_s}-{kind}-{spec.camera}-{spec.ccd}-"
                    f"{out_l}-{parts[3]}-{tag}-blf.fits",
                    f"{parts[0]}-{parts[1]}-{kind}-{spec.camera}-{spec.ccd}-"
                    f"{out_l}-{parts[3]}-{tag}-blf.fits",
                )
            )
        patterns.extend(
            (
                f"tess*-{sec_s}-{kind}-{spec.camera}-{spec.ccd}-{out_l}-*-{tag}-blf.fits",
                f"tess-{sec_s}-{kind}-{spec.camera}-{spec.ccd}-{out_l}-{tag}-blf.fits",
            )
        )
        predicate = lambda p: not _is_alt_blf_name(p)  # noqa: E731

    return _find_all_local_matches(search_dirs, patterns, predicate=predicate)


def resolve_collateral_file(
    spec: CollateralSpec,
    *,
    kind: Literal["smrow", "tvcol"],
    ddir: str | os.PathLike[str] | None = None,
    cache_dir: Path | None = None,
    offline: bool = False,
    stage_dir: Path | None = None,
    keep_raw_in_stage_only: bool = False,
    s3_fs: Any = None,
    verbose: bool = False,
) -> Path:
    """Resolve a raw ``_col.fits`` / ``_fast-col.fits`` file locally or from MAST S3."""
    out_l = spec.output_lower
    tag = spec.col_tag
    sec_s = spec.sector_subdir
    parts = spec.tpf_parts
    search_dirs = _candidate_search_dirs(ddir, cache_dir, subdir=sec_s)
    if stage_dir is not None:
        search_dirs = [stage_dir, *search_dirs]

    patterns: list[str] = []
    if parts is not None and len(parts) >= 4:
        patterns.extend(
            (
                f"tess*-{sec_s}-{kind}-{spec.camera}-{spec.ccd}-{out_l}-{parts[3]}-{tag}.fits",
                f"{parts[0]}-{parts[1]}-{kind}-{spec.camera}-{spec.ccd}-{out_l}-{parts[3]}-{tag}.fits",
            )
        )
    patterns.append(f"tess*-{sec_s}-{kind}-{spec.camera}-{spec.ccd}-{out_l}-*-{tag}.fits")

    found = _find_local_match(search_dirs, patterns)
    if found is not None:
        return found

    _guard_offline(
        offline,
        f"{kind} collateral ({sec_s}, cam={spec.camera}, ccd={spec.ccd}, out={spec.output}, {tag})",
        search_dirs,
    )

    fs = _open_anon_s3_filesystem(s3_fs=s3_fs)
    glob_pat = (
        f"{S3_TESS_PUBLIC_ROOT}/ffi/{sec_s}/*/*/{spec.camera}-{spec.ccd}/"
        f"tess*-{sec_s}-{kind}-{spec.camera}-{spec.ccd}-{out_l}-*-{tag}.fits"
    )
    matches = sorted(fs.glob(glob_pat))
    if not matches:
        raise FileNotFoundError(
            f"Collateral file not found on MAST S3 matching s3://{glob_pat}."
        )
    selected = matches
    if parts is not None and len(parts) >= 4:
        run_suffix = f"-{parts[3]}-{tag}.fits"
        run_matches = [k for k in matches if Path(k).name.endswith(run_suffix)]
        if run_matches:
            selected = run_matches
    s3_key = selected[-1]
    filename = Path(s3_key).name

    if keep_raw_in_stage_only and stage_dir is not None:
        dst = stage_dir / filename
    elif cache_dir is not None:
        dst = cache_dir / sec_s / filename
    elif ddir is not None:
        dst = Path(ddir).expanduser() / filename
    else:
        dst = resolve_cache_dir(env_var="MUNDEY_CACHE_DIR", create=True) / sec_s / filename

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
    """Stream ``TSCALn``/``TZEROn`` collateral tables in ``chunk_rows`` slices (< 2 GB RSS)."""
    if chunk_rows <= 0:
        raise ValueError(f"chunk_rows must be positive, got {chunk_rows}.")

    with fits.open(col_path, memmap=True, uint=False) as hdul:
        hdu1 = hdul[1]
        data = hdu1.data
        n_total = len(data)

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

        needed = (
            np.unique(np.asarray(row_indices, dtype=np.int64))
            if row_indices is not None
            else None
        )

        for t0 in range(0, n_total, chunk_rows):
            t1 = min(t0 + chunk_rows, n_total)
            if needed is not None:
                pos0 = int(np.searchsorted(needed, t0, side="left"))
                if pos0 >= len(needed) or int(needed[pos0]) >= t1:
                    continue

            chunk_rec = data[t0:t1]
            time_chunk = np.asarray(chunk_rec["TIME"], dtype=np.float64)
            cad_chunk = np.asarray(chunk_rec["CADENCENO"], dtype=np.int64)
            raw_f64 = np.asarray(chunk_rec[raw_colname], dtype=np.float64)

            bad = (raw_f64 < 0.0) | (~np.isfinite(raw_f64))
            if scaled_null is not None:
                bad |= raw_f64 == scaled_null
            if np.any(bad):
                raw_f64 = raw_f64.copy()
                raw_f64[bad] = np.nan

            yield t0, t1, time_chunk, cad_chunk, raw_f64


def _time_cadence_columns(
    time: np.ndarray,
    cadenceno: np.ndarray,
    time_unit: str,
) -> list[fits.Column]:
    return [
        fits.Column(
            name="TIME",
            format="D",
            unit=time_unit,
            disp="D14.7",
            array=np.asarray(time, dtype=np.float64),
        ),
        fits.Column(
            name="CADENCENO",
            format="J",
            disp="I10",
            array=np.asarray(cadenceno, dtype=np.int32),
        ),
    ]


def _write_blf_hdul(
    dst_path: str | os.PathLike[str],
    hdul: fits.HDUList,
    *,
    stage_dir: Path | None = None,
) -> Path:
    """Atomically write ``hdul`` to ``dst_path`` (optionally staging in ``stage_dir``)."""
    dst = Path(dst_path)
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
    b1d_arr = np.asarray(black1d, dtype=np.float64)
    nrows = b1d_arr.shape[1]
    cols = _time_cadence_columns(time, cadenceno, time_unit)
    cols.append(
        fits.Column(
            name="1DBLACK",
            format=f"{nrows}D",
            unit="count",
            disp="F14.7",
            dim=f"({nrows})",
            array=b1d_arr,
        )
    )
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
    return _write_blf_hdul(dst_path, fits.HDUList([fits.PrimaryHDU(), hdu1]), stage_dir=stage_dir)


def write_smear_blf(
    dst_path: str | os.PathLike[str],
    *,
    time: np.ndarray,
    cadenceno: np.ndarray,
    smear: np.ndarray,
    smear_err: np.ndarray | None = None,
    smear_sat: np.ndarray | None = None,
    unresolved_cols: np.ndarray | None = None,
    col_matched: np.ndarray | None = None,
    smear_mode: str = "standard",
    cutout_bounds: tuple[int, int, int, int] | None = None,
    cadence_match: str = "cadenceno",
    source_mask_digest: str = "NONE",
    tv_all_matched: bool | None = None,
    int_time: float | None = None,
    time_unit: str = "BJD - 2457000, days",
    stage_dir: Path | None = None,
) -> Path:
    """Atomically write a ``smrow-...-blf.fits`` cache with ``CADENCENO`` (``'J'``)."""
    sm_arr = np.asarray(smear, dtype=np.float64)
    ncols = sm_arr.shape[1]
    cols = _time_cadence_columns(time, cadenceno, time_unit)
    cols.append(
        fits.Column(
            name="SMEAR",
            format=f"{ncols}D",
            unit="e-/frame",
            disp="F14.7",
            dim=f"({ncols})",
            array=sm_arr,
        )
    )
    if smear_err is not None:
        cols.append(
            fits.Column(
                name="SMEAR_ERR",
                format=f"{ncols}E",
                unit="e-/frame",
                disp="E14.7",
                dim=f"({ncols})",
                array=np.asarray(smear_err, dtype=np.float32),
            )
        )
    if col_matched is not None:
        cols.append(
            fits.Column(
                name="COLMATCH",
                format="L",
                array=np.asarray(col_matched, dtype=bool),
            )
        )

    hdu1 = fits.BinTableHDU.from_columns(fits.ColDefs(cols), name="SMEAR")
    hdu1.header["BLFVER"] = (2, "Mundey -blf cache schema version")
    hdu1.header["SMEARMOD"] = (smear_mode, "Mundey smear estimation mode")
    hdu1.header["CADMATCH"] = (str(cadence_match), "Collateral cadence pairing mode")
    hdu1.header["SRCMASK"] = (
        str(source_mask_digest),
        "SHA-256 prefix of SOURCE_MASK (or NONE)",
    )
    if tv_all_matched is not None:
        hdu1.header["TVMATCH"] = (
            bool(tv_all_matched),
            "Whether all smrow cadences matched tvcol 1:1",
        )
    if int_time is not None:
        hdu1.header["INT_TIME"] = (float(int_time), "Single-frame integration time (s)")
    if cutout_bounds is not None:
        crval1p, cutnx, crval2p, cutny = cutout_bounds
        hdu1.header["CRVAL1P"] = (int(crval1p), "Cutout start CCD column")
        hdu1.header["CUTNX"] = (int(cutnx), "Cutout column width")
        hdu1.header["CRVAL2P"] = (int(crval2p), "Cutout start CCD row")
        hdu1.header["CUTNY"] = (int(cutny), "Cutout row height")

    hdu_list: list[fits.PrimaryHDU | fits.BinTableHDU | fits.ImageHDU] = [
        fits.PrimaryHDU(),
        hdu1,
    ]
    if smear_sat is not None:
        hdu_list.append(
            fits.ImageHDU(data=np.asarray(smear_sat, dtype=np.uint8), name="SMROW_SAT")
        )
    if unresolved_cols is not None:
        hdu_list.append(
            fits.ImageHDU(data=np.asarray(unresolved_cols, dtype=np.uint8), name="ALT_UNRES")
        )
    return _write_blf_hdul(dst_path, fits.HDUList(hdu_list), stage_dir=stage_dir)
