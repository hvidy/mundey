"""JAX ``@in_float64`` TESS bright-star calibration engine (``mundey`` v2.0).

Implements physical TESS readout calibration in locally scoped 64-bit JAX:
1. Frame-scaled 2D black & output mean-black subtraction with sentinel masking.
2. Cancellation-free discrete orthonormal Legendre 1D-black fitting on clean
   trailing virtual columns (``cols 3..9``) with AIC degree selection.
3. Causal serial LDE undershoot inversion via ``jax.lax.associative_scan`` in
   ADU *before* per-frame ``[0, maxDomain]``-clamped nonlinearity and gain, with
   steady-state restart across ``NaN`` gaps.
4. Exact ``CADENCENO`` collateral pairing across unified 2048-column multi-output
   (``'A'``, ``'B'``, ``'C'``, ``'D'``) CCD cutouts in 4,096-cadence chunks.
5. Concavity-guarded collateral-anchored ``smear='alternate'`` and FFI cutout
   desmearing (:func:`desmear_calibrated_cutout`) under ``HALOCVER = 1``.
"""

from __future__ import annotations

import contextlib
import errno
import hashlib
import importlib.resources
import os
import warnings
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import jax
import jax.numpy as jnp
import lightkurve as lk
import numpy as np
from astropy.io import fits
from scipy import linalg as sp_linalg

from ._contract import (
    PIXEL_QUALITY_EXTNAME,
    PixelQuality,
    _is_fred_path,
    append_smear_err_column,
    resolve_cache_dir,
    write_mundey_provenance,
    write_pixel_quality,
)
from ._precision import in_float64
from .archive import (
    DEFAULT_CHUNK_ROWS,
    CadenceKind,
    CollateralSpec,
    StaticModelPaths,
    collateral_blf_cache_path,
    collateral_suffix,
    find_all_cached_blf,
    infer_cadence_kind,
    iter_collateral_chunks,
    read_collateral_headers_and_index,
    resolve_collateral_file,
    resolve_static_models,
    write_onedblack_blf,
    write_smear_blf,
)
from .smear import (
    TESS_PARALLEL_SHIFTS,
    AlternateSmearConfig,
    _build_clean_rows_2d,
    _detect_bloom_columns,  # noqa: F401
    _dilate_vertical_bool,
    _estimate_alternate_smear_stamp,
    _eval_bloom_from_stats,
    _select_bucket,
    _solve_huber_smear_location,  # noqa: F401
    _unit_normal_median_variance,
)

_MISSING = object()

CadenceMatchMode = Literal["cadenceno", "positional"]
MissingCadenceMode = Literal["nan", "error"]
SmearMode = Literal["automatic", "standard", "alternate"]

OUTPUT_ORDER: tuple[str, str, str, str] = ("A", "B", "C", "D")
OUTPUT_BOUNDS_0BASED: dict[str, tuple[int, int, int, int]] = {
    "A": (44, 556, 2092, 2103),
    "B": (556, 1068, 2103, 2114),
    "C": (1068, 1580, 2114, 2125),
    "D": (1580, 2092, 2125, 2136),
}
OUTPUT_REVERSE_READOUT: dict[str, bool] = {"A": False, "B": True, "C": False, "D": True}

DEFAULT_2MIN_EXPOSURE_SECONDS: float = 95.04  # 60 frames * 0.8 CRMITEN * 1.98 s
CRMITEN_BLOCK10_VAR_FACTOR: float = 0.84285
_ONEDBLACK_CAD_BUCKETS: tuple[int, ...] = (32, 256, 4096)
_PRESMEAR_WIDTH_BUCKETS: tuple[int, ...] = (16, 64, 512, 1056, 2048)


@dataclass(frozen=True)
class CalibrationConfig:
    """User-facing configuration for :func:`recalibrate_tpf`."""

    smear: SmearMode = "automatic"
    cadence_match: CadenceMatchMode = "cadenceno"
    missing_cadence: MissingCadenceMode = "nan"
    cache_dir: str | os.PathLike[str] | None = None
    ddir: str | os.PathLike[str] | None = None
    offline: bool = False
    strict_hpc: bool = False
    chunk_rows: int = DEFAULT_CHUNK_ROWS
    verbose: bool = False


@dataclass(frozen=True)
class LinearityChannel:
    """Parsed TESS ``linearity.xml`` polynomial specification for a single CCD output."""

    camera: int
    ccd: int
    output: str
    originx: float
    scalex: float
    offsetx: float
    max_domain: float
    coeffs: tuple[float, ...]


def _build_ccid80_strap_cols_1based_imaging() -> tuple[int, ...]:
    left: list[int] = []
    for base in (77, 205, 333, 461, 589, 717, 845, 973):
        left.extend((base, base + 1, base + 4, base + 5, base + 8, base + 9))
    left.extend((516, 517, 520, 521, 524, 525, 673, 674, 675, 676))
    left_sorted = sorted(left)
    return tuple(left_sorted + sorted(2049 - c for c in left_sorted))


CCID80_STRAP_COLS_IMAGING_1BASED: tuple[int, ...] = _build_ccid80_strap_cols_1based_imaging()
CCID80_STRAP_COLS_RAWX_1BASED: np.ndarray = (
    np.asarray(CCID80_STRAP_COLS_IMAGING_1BASED, dtype=np.int64) + 44
)


def strap_mask_for_rawx_columns(cols_1based: np.ndarray, *, buffer_px: int = 2) -> np.ndarray:
    """Return boolean strap mask for 1-based FITS ``RAWX`` columns with ``±buffer_px``."""
    cols = np.asarray(cols_1based, dtype=np.int64)[:, np.newaxis]
    dist = np.min(np.abs(cols - CCID80_STRAP_COLS_RAWX_1BASED[np.newaxis, :]), axis=1)
    return np.asarray(dist <= int(buffer_px), dtype=bool)


def _build_orthonormal_legendre_basis(n_rows: int, max_degree: int) -> np.ndarray:
    x = np.linspace(-1.0, 1.0, n_rows, dtype=np.float64)
    q, r = np.linalg.qr(np.polynomial.legendre.legvander(x, max_degree), mode="reduced")
    return np.asarray(q * np.where(np.diag(r) >= 0.0, 1.0, -1.0)[np.newaxis, :], dtype=np.float64)


_LEGENDRE_Q_SCI_NP: np.ndarray = _build_orthonormal_legendre_basis(2048, 49)
_LEGENDRE_Q_OV_NP: np.ndarray = _build_orthonormal_legendre_basis(30, 2)


def load_linearity_channel(
    xml_path: str | os.PathLike[str], camera: int, ccd: int, output: str
) -> LinearityChannel:
    """Parse ``linearityPoly`` parameters (including ``maxDomain``) from ``linearity.xml``."""
    root = ET.parse(Path(xml_path)).getroot()
    xpath = (
        f".//*[@cameraNumber='{int(camera)}']"
        f"[@ccdNumber='{int(ccd)}']"
        f"[@ccdOutput='{output.upper()}']/linearityPoly"
    )
    poly = root.find(xpath)
    if poly is None:
        raise KeyError(
            f"Linearity polynomial not found in {xml_path} for "
            f"camera={camera}, ccd={ccd}, output={output}."
        )
    coeffs = tuple(float(elem.attrib["value"]) for elem in poly.findall("coeffs/coeff"))
    if not coeffs:
        raise ValueError(
            f"Empty linearity coefficient list for camera={camera}, ccd={ccd}, output={output}."
        )
    return LinearityChannel(
        camera=int(camera),
        ccd=int(ccd),
        output=output.upper(),
        originx=float(poly.attrib["originx"]),
        scalex=float(poly.attrib["scalex"]),
        offsetx=float(poly.attrib.get("offsetx", "0.0")),
        max_domain=float(poly.attrib["maxDomain"]),
        coeffs=coeffs,
    )


def linearity_max_electrons_per_frame(channel: LinearityChannel, gain: float) -> float:
    """Return exact calibrated ``e-/frame`` evaluated at ``dn_per_frame == channel.max_domain``."""
    s_max = channel.offsetx + (channel.max_domain - channel.originx) * channel.scalex
    p_max = float(np.polyval(np.asarray(channel.coeffs, dtype=np.float64)[::-1], s_max))
    return float(channel.max_domain * p_max * float(gain))


def load_undershoot_coeff(camera: int, ccd: int, output: str) -> float:
    """Load the dimensionless LDE undershoot coefficient from bundled ``tess-undershoot.xml``."""
    resource = importlib.resources.files("mundey").joinpath("data/tess-undershoot.xml")
    with resource.open("rb") as fh:
        root = ET.parse(fh).getroot()
    xpath = (
        f".//*[@cameraNumber='{int(camera)}']"
        f"[@ccdNumber='{int(ccd)}']"
        f"[@ccdOutput='{output.upper()}']"
    )
    node = root.find(xpath)
    if node is None:
        raise KeyError(
            f"Undershoot factor not found for camera={camera}, ccd={ccd}, output={output}."
        )
    return float(node.attrib["undershoot"])


def _hdr_lookup(key: str, *headers: fits.Header | None) -> Any:
    for hdr in headers:
        if hdr is not None and key in hdr:
            return hdr[key]
    return None


def effective_frame_count(
    hdr0: fits.Header | None,
    hdr1: fits.Header | None,
    *,
    fallback_hdr0: fits.Header | None = None,
    fallback_hdr1: fits.Header | None = None,
) -> float:
    """Compute effective co-added frame count ``N_eff = factor * NUM_FRM`` from FITS headers."""
    hdrs = (hdr0, hdr1, fallback_hdr0, fallback_hdr1)
    crmiten = bool(_hdr_lookup("CRMITEN", *hdrs) or False)
    if crmiten:
        crblksz_raw = _hdr_lookup("CRBLKSZ", *hdrs)
        crblksz = float(crblksz_raw) if crblksz_raw is not None else 10.0
        if crblksz <= 2.0:
            raise ValueError(f"Invalid CRBLKSZ={crblksz} with CRMITEN=True.")
        factor = (crblksz - 2.0) / crblksz
    else:
        factor = 1.0
    num_frm_raw = _hdr_lookup("NUM_FRM", *hdrs)
    if num_frm_raw is None:
        raise KeyError("Missing NUM_FRM keyword in FITS headers.")
    num_frm = float(num_frm_raw)
    if num_frm <= 0.0:
        raise ValueError(f"NUM_FRM must be positive, got {num_frm}.")
    return factor * num_frm


def _crmiten_variance_factor(
    hdr0: fits.Header | None,
    hdr1: fits.Header | None,
    *,
    fallback_hdr0: fits.Header | None = None,
    fallback_hdr1: fits.Header | None = None,
) -> float:
    hdrs = (hdr0, hdr1, fallback_hdr0, fallback_hdr1)
    if not bool(_hdr_lookup("CRMITEN", *hdrs) or False):
        return 1.0
    crblksz_raw = _hdr_lookup("CRBLKSZ", *hdrs)
    crblksz = int(round(float(crblksz_raw))) if crblksz_raw is not None else 10
    return CRMITEN_BLOCK10_VAR_FACTOR if crblksz == 10 else 1.0


def _frame_transfer_row_weight(
    hdr0: fits.Header | None = None,
    hdr1: fits.Header | None = None,
) -> float:
    """Return dimensionless shutterless row-smear ratio ``t_row / int_time``."""
    readtime, int_time = None, None
    for hdr in (hdr1, hdr0):
        if hdr is None:
            continue
        if readtime is None and "READTIME" in hdr:
            with contextlib.suppress(TypeError, ValueError):
                v = float(hdr["READTIME"])
                if np.isfinite(v) and v > 0.0:
                    readtime = v
        if int_time is None and "INT_TIME" in hdr:
            with contextlib.suppress(TypeError, ValueError):
                v = float(hdr["INT_TIME"])
                if np.isfinite(v) and v > 0.0:
                    int_time = v
    if readtime is not None and int_time is not None:
        return float(readtime / (int_time * float(TESS_PARALLEL_SHIFTS)))

    cr_factor = 1.0
    if bool(_hdr_lookup("CRMITEN", hdr0, hdr1) or False):
        raw_blk = _hdr_lookup("CRBLKSZ", hdr0, hdr1)
        crblksz = 10.0
        if raw_blk is not None:
            with contextlib.suppress(TypeError, ValueError):
                v = float(raw_blk)
                if np.isfinite(v) and v > 2.0:
                    crblksz = v
        cr_factor = (crblksz - 2.0) / crblksz

    deadc = 0.99
    for hdr in (hdr1, hdr0):
        if hdr is not None and "DEADC" in hdr:
            with contextlib.suppress(TypeError, ValueError):
                v = float(hdr["DEADC"])
                if np.isfinite(v) and 0.0 < v < 1.0:
                    adj = v / cr_factor
                    if 0.5 < adj < 1.0:
                        deadc = adj
                        break
                    if 0.5 < v < 1.0:
                        deadc = v
                        break
    return float((1.0 - deadc) / (deadc * float(TESS_PARALLEL_SHIFTS)))


def _extract_positive_int_card(hdr: fits.Header | None, keys: tuple[str, ...]) -> int | None:
    if hdr is None:
        return None
    for key in keys:
        if key in hdr and not isinstance(hdr[key], bool | np.bool_):
            with contextlib.suppress(TypeError, ValueError):
                val = float(hdr[key])
                if np.isfinite(val) and val > 0.0:
                    return int(round(val))
    return None


def _stamp_ccd_origin_1based(
    hdr1: fits.Header | None,
    hdr2: fits.Header | None = None,
    *,
    default: tuple[int, int] | None = None,
) -> tuple[int, int]:
    """Resolve 1-based CCD ``(col0, row0)``, preferring ``HDU[1]`` ``1CRV5P``/``2CRV5P``."""
    col0 = _extract_positive_int_card(hdr1, ("1CRV5P", "1CRV4P")) or _extract_positive_int_card(
        hdr2, ("CRVAL1P",)
    )
    row0 = _extract_positive_int_card(hdr1, ("2CRV5P", "2CRV4P")) or _extract_positive_int_card(
        hdr2, ("CRVAL2P",)
    )
    if col0 is None:
        if default is None:
            raise KeyError("Missing CCD column origin ('1CRV5P', '1CRV4P', or 'CRVAL1P').")
        col0 = int(default[0])
    if row0 is None:
        if default is None:
            raise KeyError("Missing CCD row origin ('2CRV5P', '2CRV4P', or 'CRVAL2P').")
        row0 = int(default[1])
    return col0, row0


def _sync_stamp_ccd_origin(
    hdr1: fits.Header | None, hdr2: fits.Header | None, col0_1based: int, row0_1based: int
) -> None:
    col_i, row_i = int(col0_1based), int(row0_1based)
    if hdr1 is not None:
        for key in list(hdr1.keys()):
            if len(key) >= 6 and key.endswith("P") and key[1:4] == "CRV" and key[4:-1].isdigit():
                if key[0] == "1" and hdr1[key] != col_i:
                    hdr1[key] = col_i
                elif key[0] == "2" and hdr1[key] != row_i:
                    hdr1[key] = row_i
    if hdr2 is not None:
        hdr2["CRVAL1P"] = col_i
        hdr2["CRVAL2P"] = row_i


def _source_mask_digest(source_mask_2d: np.ndarray | None) -> str:
    """Return a deterministic 16-hex-char SHA-256 digest of ``source_mask_2d``, or ``'NONE'``."""
    if source_mask_2d is None:
        return "NONE"
    mask_bool = np.asarray(source_mask_2d, dtype=bool)
    if not np.any(mask_bool):
        return "NONE"
    h = hashlib.sha256()
    h.update(f"{mask_bool.shape[0]}x{mask_bool.shape[1]}:".encode("ascii"))
    h.update(np.ascontiguousarray(mask_bool, dtype=np.uint8).tobytes())
    return h.hexdigest()[:16]


def _effective_exposure_seconds(
    hdr0: fits.Header | None, hdr1: fits.Header | None
) -> float:
    int_time = _hdr_lookup("INT_TIME", hdr1, hdr0)
    if int_time is not None and _hdr_lookup("NUM_FRM", hdr1, hdr0) is not None:
        with contextlib.suppress(KeyError, TypeError, ValueError):
            val = float(int_time)
            if val > 0.0:
                return float(effective_frame_count(hdr0, hdr1) * val)
    deadc = 0.99
    deadc_raw = _hdr_lookup("DEADC", hdr1, hdr0)
    if deadc_raw is not None:
        with contextlib.suppress(TypeError, ValueError):
            v = float(deadc_raw)
            if 0.0 < v <= 1.0:
                deadc = v
    timedel_raw = _hdr_lookup("TIMEDEL", hdr1, hdr0)
    if timedel_raw is not None:
        with contextlib.suppress(TypeError, ValueError):
            dt_days = float(timedel_raw)
            if dt_days > 0.0:
                return float(dt_days * 86400.0 * deadc)
    return DEFAULT_2MIN_EXPOSURE_SECONDS


def output_for_rawx_col(col_1based: int) -> str:
    """Map a 1-based FITS ``RAWX`` column (``45..2092``) to TESS output ``'A'..'D'``."""
    if 45 <= col_1based <= 556:
        return "A"
    if 557 <= col_1based <= 1068:
        return "B"
    if 1069 <= col_1based <= 1580:
        return "C"
    if 1581 <= col_1based <= 2092:
        return "D"
    raise ValueError(
        f"TPF column {int(col_1based)} falls outside TESS science outputs 45..2092."
    )


def _compute_quality_mask(
    hdu1: fits.BinTableHDU,
    *,
    quality_bitmask: Any = "default",
    target_obj: Any = None,
) -> np.ndarray:
    """Return boolean ``(N_raw,)`` mask of cadences passing ``QUALITY`` bitmask checks."""
    n_raw = len(hdu1.data)
    if target_obj is not None:
        qmask_attr = getattr(target_obj, "quality_mask", None)
        if qmask_attr is not None:
            arr = np.asarray(qmask_attr, dtype=bool)
            if arr.shape == (n_raw,):
                return arr

    if "QUALITY" not in (hdu1.data.names or ()) or quality_bitmask in (None, 0, "none"):
        return np.ones(n_raw, dtype=bool)

    qual = np.asarray(hdu1.data["QUALITY"], dtype=np.int64)
    try:
        return np.asarray(
            lk.utils.TessQualityFlags.create_quality_mask(qual, bitmask=quality_bitmask),
            dtype=bool,
        )
    except Exception:
        return qual == 0


def match_cadenceno(
    tpf_cadenceno: np.ndarray,
    col_cadenceno: np.ndarray,
    *,
    cadence_match: CadenceMatchMode = "cadenceno",
    label: str = "collateral",
) -> tuple[np.ndarray, np.ndarray]:
    """Align TPF cadences to a collateral or ``-blf.fits`` table."""
    tpf_cad = np.asarray(tpf_cadenceno, dtype=np.int64)
    col_cad = np.asarray(col_cadenceno, dtype=np.int64)
    n_tpf, n_col = len(tpf_cad), len(col_cad)
    if n_tpf == 0:
        raise ValueError("TPF contains zero cadences.")
    if n_col == 0:
        raise ValueError(f"{label} contains zero cadences.")

    if cadence_match == "positional":
        if n_col < n_tpf:
            raise ValueError(
                f"Positional cadence pairing failed for {label}: table has {n_col} rows, "
                f"fewer than TPF's {n_tpf} rows."
            )
        return np.ones(n_tpf, dtype=bool), np.arange(n_tpf, dtype=np.int64)

    if cadence_match != "cadenceno":
        raise ValueError(
            f"Invalid cadence_match={cadence_match!r}; expected 'cadenceno' or 'positional'."
        )

    order = np.argsort(col_cad, kind="mergesort")
    sorted_col = col_cad[order]
    pos = np.searchsorted(sorted_col, tpf_cad, side="left")
    in_bounds = (pos >= 0) & (pos < n_col)
    safe_pos = np.where(in_bounds, pos, 0)
    matched = in_bounds & (sorted_col[safe_pos] == tpf_cad)
    if not np.any(matched):
        t_rng = f"{int(tpf_cad.min())}..{int(tpf_cad.max())}"
        c_rng = f"{int(col_cad.min())}..{int(col_cad.max())}"
        raise ValueError(
            f"Zero matching CADENCENO rows between TPF ({t_rng}) and {label} ({c_rng})."
        )
    return matched, np.where(matched, order[safe_pos], -1).astype(np.int64)


def _enforce_missing_cadence_policy(
    tpf_cadenceno: np.ndarray,
    combined_matched: np.ndarray,
    *,
    missing_cadence: MissingCadenceMode,
) -> int:
    if missing_cadence not in ("nan", "error"):
        raise ValueError(
            f"Invalid missing_cadence={missing_cadence!r}; expected 'nan' or 'error'."
        )
    n_tpf, n_matched = len(tpf_cadenceno), int(np.sum(combined_matched))
    if n_matched == 0:
        raise ValueError(
            "Zero TPF cadences have complete collateral CADENCENO coverage across outputs."
        )
    n_missing = n_tpf - n_matched
    if n_missing == 0:
        return 0
    missing_cads = np.asarray(tpf_cadenceno, dtype=np.int64)[~combined_matched]
    msg = (
        f"{n_missing} / {n_tpf} TPF cadences lack matching collateral CADENCENO "
        f"(missing CADENCENO range {int(missing_cads.min())}..{int(missing_cads.max())})"
    )
    if missing_cadence == "error":
        raise ValueError(f"{msg}.")
    warnings.warn(f"{msg}; masking unmatched cadences to NaN.", UserWarning, stacklevel=3)
    return n_missing


# ---------------------------------------------------------------------------
# 1D-black Legendre fitting & fused causal LDE undershoot + linearity
# ---------------------------------------------------------------------------


def _repair_partial_onedblack_cadence(
    y_cad: np.ndarray, q_sci_np: np.ndarray, q_ov_np: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Exact masked Cholesky-Legendre least-squares fit for a cadence with ``NaN`` row gaps."""
    y_sci = y_cad[:2048]
    m_sci = np.isfinite(y_sci)
    n_fin = int(np.sum(m_sci))
    if n_fin < 64:
        return np.full(2078, np.nan, dtype=np.float64), np.full(2078, np.nan, dtype=np.float64)

    q_m, y_m = q_sci_np[m_sci], y_sci[m_sci]
    chol_l = sp_linalg.cholesky(
        q_m.T @ q_m + 1e-14 * np.eye(50, dtype=np.float64), lower=True, check_finite=False
    )
    z_orth = sp_linalg.solve_triangular(chol_l, q_m.T @ y_m, lower=True, check_finite=False)
    c_full = sp_linalg.solve_triangular(chol_l.T, z_orth, lower=False, check_finite=False)
    res_full = y_m - q_m @ c_full
    norm_perp_sq = float(np.dot(res_full, res_full))

    suffix_excl = np.cumsum((z_orth * z_orth)[::-1])[::-1]
    tail_all = np.concatenate([suffix_excl[1:], np.zeros(1, dtype=np.float64)])
    rss_20_49 = norm_perp_sq + tail_all[20:50]

    n_f = float(n_fin)
    max_deg = min(49, max(20, (n_fin // 2) - 1))
    ks = np.arange(20, 50, dtype=np.float64)
    aics = np.where(
        ks <= float(max_deg),
        2.0 * (ks + 1.0) + n_f * np.log(np.maximum(rss_20_49 / n_f, 1e-30)),
        np.inf,
    )
    best_rel = int(np.argmin(aics))
    best_k = best_rel + 20
    c_sel = sp_linalg.solve_triangular(
        chol_l.T, np.where(np.arange(50) <= best_k, z_orth, 0.0), lower=False, check_finite=False
    )
    b1d_sci = q_sci_np @ c_sel
    sigma2_sci = float(rss_20_49[best_rel]) / max(n_f - float(best_k + 1), 1.0)
    v_t = sp_linalg.solve_triangular(
        chol_l[: best_k + 1, : best_k + 1],
        q_sci_np[:, : best_k + 1].T,
        lower=True,
        check_finite=False,
    )
    err_sci = np.sqrt(np.maximum(sigma2_sci * np.sum(v_t * v_t, axis=0), 0.0))

    y_ov = y_cad[2048:2078]
    m_ov = np.isfinite(y_ov)
    n_ov = int(np.sum(m_ov))
    if n_ov >= 6:
        q_ov_m, y_ov_m = q_ov_np[m_ov], y_ov[m_ov]
        chol_ov = sp_linalg.cholesky(
            q_ov_m.T @ q_ov_m + 1e-14 * np.eye(3, dtype=np.float64),
            lower=True,
            check_finite=False,
        )
        z_ov = sp_linalg.solve_triangular(
            chol_ov, q_ov_m.T @ y_ov_m, lower=True, check_finite=False
        )
        c_ov = sp_linalg.solve_triangular(chol_ov.T, z_ov, lower=False, check_finite=False)
        b1d_ov = q_ov_np @ c_ov
        res_ov = y_ov_m - q_ov_m @ c_ov
        sigma2_ov = float(np.dot(res_ov, res_ov)) / max(float(n_ov - 3), 1.0)
        v_ov_t = sp_linalg.solve_triangular(chol_ov, q_ov_np.T, lower=True, check_finite=False)
        err_ov = np.sqrt(np.maximum(sigma2_ov * np.sum(v_ov_t * v_ov_t, axis=0), 0.0))
    else:
        b1d_ov = np.full(30, b1d_sci[-1], dtype=np.float64)
        err_ov = np.full(30, err_sci[-1], dtype=np.float64)
    return np.concatenate([b1d_sci, b1d_ov]), np.concatenate([err_sci, err_ov])


@jax.jit
def _fit_onedblack_legendre_jit(
    tv_debiased_adu: jax.Array, q_sci: jax.Array, q_ov: jax.Array
) -> tuple[jax.Array, jax.Array, jax.Array, jax.Array]:
    """JIT-compiled 64-bit orthonormal Legendre 1D-black solver for ``(N_c, 2078, 11)``."""
    v_clean = tv_debiased_adu[:, :, 3:10]
    abs_dev = jnp.abs(v_clean - jnp.nanmedian(v_clean, axis=2, keepdims=True))
    row_mad = 1.4826 * jnp.nanmedian(abs_dev, axis=2, keepdims=True)
    global_mad = 1.4826 * jnp.nanmedian(abs_dev, axis=(1, 2), keepdims=True)
    sigma_row = jnp.maximum(row_mad, jnp.maximum(global_mad, 1e-6))

    valid_pix = jnp.isfinite(v_clean) & (abs_dev <= 5.0 * sigma_row)
    use_pix = jnp.where(jnp.any(valid_pix, axis=2, keepdims=True), valid_pix, jnp.isfinite(v_clean))
    count_pix = jnp.sum(use_pix, axis=2)
    y_all = jnp.where(
        count_pix > 0,
        jnp.sum(jnp.where(use_pix, v_clean, 0.0), axis=2) / jnp.maximum(count_pix, 1),
        jnp.nan,
    )

    y_sci_raw, y_ov_raw = y_all[:, :2048], y_all[:, 2048:2078]
    m_sci, m_ov = jnp.isfinite(y_sci_raw), jnp.isfinite(y_ov_raw)
    n_fin_sci = jnp.sum(m_sci, axis=1)
    cad_valid = n_fin_sci >= 64
    needs_repair = cad_valid & ((n_fin_sci < 2048) | (jnp.sum(m_ov, axis=1) < 30))

    y_sci = jnp.where(m_sci, y_sci_raw, 0.0)
    coeffs = y_sci @ q_sci
    y_perp = y_sci - coeffs @ q_sci.T
    norm_perp_sq = jnp.sum(y_perp * y_perp, axis=1)

    suffix_excl = jnp.cumsum((coeffs * coeffs)[:, ::-1], axis=1)[:, ::-1]
    tail_all = jnp.concatenate(
        [suffix_excl[:, 1:], jnp.zeros((coeffs.shape[0], 1), dtype=coeffs.dtype)], axis=1
    )
    rss_20_49 = norm_perp_sq[:, None] + tail_all[:, 20:50]

    ks = jnp.arange(20, 50, dtype=coeffs.dtype)[None, :]
    n_rows_f = jnp.asarray(2048.0, dtype=coeffs.dtype)
    aics = 2.0 * (ks + 1.0) + n_rows_f * jnp.log(jnp.maximum(rss_20_49 / n_rows_f, 1e-30))
    best_rel = jnp.argmin(aics, axis=1)
    best_k = best_rel + 20

    b1d_sci = jnp.where(jnp.arange(50)[None, :] <= best_k[:, None], coeffs, 0.0) @ q_sci.T
    best_rss = jnp.take_along_axis(rss_20_49, best_rel[:, None], axis=1)[:, 0]
    sigma2_sci = best_rss / jnp.maximum(n_rows_f - (best_k.astype(coeffs.dtype) + 1.0), 1.0)
    lev_sci = jnp.take(jnp.cumsum(q_sci * q_sci, axis=1), best_k, axis=1).T
    err_sci = jnp.sqrt(jnp.maximum(sigma2_sci[:, None] * lev_sci, 0.0))

    y_ov = jnp.where(m_ov, y_ov_raw, 0.0)
    b1d_ov = (y_ov @ q_ov) @ q_ov.T
    res_ov = y_ov - b1d_ov
    sigma2_ov = jnp.sum(res_ov * res_ov, axis=1) / 27.0
    err_ov = jnp.sqrt(jnp.maximum(sigma2_ov[:, None] * jnp.sum(q_ov * q_ov, axis=1)[None, :], 0.0))

    return (
        jnp.where(cad_valid[:, None], jnp.concatenate([b1d_sci, b1d_ov], axis=1), jnp.nan),
        jnp.where(cad_valid[:, None], jnp.concatenate([err_sci, err_ov], axis=1), jnp.nan),
        y_all,
        needs_repair,
    )


@in_float64
def _fit_onedblack_static_chunk(tv_padded_adu: Any) -> tuple[np.ndarray, np.ndarray]:
    b1d, err, y_all, needs_repair = _fit_onedblack_legendre_jit(
        tv_padded_adu,
        jnp.asarray(_LEGENDRE_Q_SCI_NP, dtype=jnp.float64),
        jnp.asarray(_LEGENDRE_Q_OV_NP, dtype=jnp.float64),
    )
    b1d_np = np.asarray(b1d, dtype=np.float64)
    err_np = np.asarray(err, dtype=np.float64)
    repair_mask = np.asarray(needs_repair, dtype=bool)
    if np.any(repair_mask):
        b1d_np, err_np = b1d_np.copy(), err_np.copy()
        y_all_np = np.asarray(y_all, dtype=np.float64)
        for idx in np.flatnonzero(repair_mask):
            b1d_np[int(idx)], err_np[int(idx)] = _repair_partial_onedblack_cadence(
                y_all_np[int(idx)], _LEGENDRE_Q_SCI_NP, _LEGENDRE_Q_OV_NP
            )
    return b1d_np, err_np


def fit_onedblack_legendre(tv_debiased_adu: Any) -> tuple[np.ndarray, np.ndarray]:
    """Fit cancellation-free orthonormal Legendre 1D black in locally scoped ``float64``."""
    arr = np.asarray(tv_debiased_adu, dtype=np.float64)
    n_c = arr.shape[0]
    if n_c == 0:
        return np.empty((0, 2078), dtype=np.float64), np.empty((0, 2078), dtype=np.float64)
    if n_c <= _ONEDBLACK_CAD_BUCKETS[-1]:
        bucket_n = _select_bucket(n_c, _ONEDBLACK_CAD_BUCKETS)
        if bucket_n == n_c:
            return _fit_onedblack_static_chunk(arr)
        buf = np.full((bucket_n, 2078, 11), np.nan, dtype=np.float64)
        buf[:n_c] = arr
        b1d, err = _fit_onedblack_static_chunk(buf)
        return b1d[:n_c], err[:n_c]

    b1d_out = np.empty((n_c, 2078), dtype=np.float64)
    err_out = np.empty((n_c, 2078), dtype=np.float64)
    step = _ONEDBLACK_CAD_BUCKETS[-1]
    buf = np.empty((step, 2078, 11), dtype=np.float64)
    for t0 in range(0, n_c, step):
        t1 = min(t0 + step, n_c)
        m = t1 - t0
        if m == step:
            b_c, e_c = _fit_onedblack_static_chunk(arr[t0:t1])
        else:
            buf[:m], buf[m:] = arr[t0:t1], np.nan
            b_c, e_c = _fit_onedblack_static_chunk(buf)
        b1d_out[t0:t1], err_out[t0:t1] = b_c[:m], e_c[:m]
    return b1d_out, err_out


def _undershoot_combine(
    prev: tuple[jax.Array, jax.Array], curr: tuple[jax.Array, jax.Array]
) -> tuple[jax.Array, jax.Array]:
    return curr[0] * prev[0], curr[0] * prev[1] + curr[1]


def _invert_lde_undershoot_forward_impl(w_adu: jax.Array, u_coeff: jax.Array) -> jax.Array:
    """Invert causal serial LDE undershoot with steady-state restart after any ``NaN`` pixel."""
    if w_adu.shape[-1] == 0:
        return w_adu
    finite = jnp.isfinite(w_adu)
    w_clean = jnp.where(finite, w_adu, 0.0)
    prev_finite = jnp.concatenate([jnp.zeros_like(finite[..., :1]), finite[..., :-1]], axis=-1)
    is_restart = finite & (~prev_finite)
    is_continue = finite & prev_finite
    a = jnp.where(is_continue, u_coeff, 0.0)
    b = jnp.where(
        is_restart, w_clean / (1.0 - u_coeff), jnp.where(is_continue, w_clean, 0.0)
    )
    x_scan = (
        jax.lax.associative_scan(_undershoot_combine, (a, b), axis=-1)[1]
        if w_adu.shape[-1] > 1
        else b
    )
    return jnp.where(finite, x_scan, jnp.nan)


@jax.jit
def _calibrate_pre_smear_2d_jit(
    w_adu_2d: jax.Array,
    u_coeff: jax.Array,
    n_eff: jax.Array,
    originx: jax.Array,
    scalex: jax.Array,
    offsetx: jax.Array,
    max_domain: jax.Array,
    coeffs_desc: jax.Array,
    deriv_coeffs_desc: jax.Array,
    gain: jax.Array,
) -> tuple[jax.Array, jax.Array, jax.Array]:
    """Fused 2D ``(B, W)`` 64-bit undershoot -> per-frame DN -> clamped linearity & gain."""
    dn_per_frame = _invert_lde_undershoot_forward_impl(w_adu_2d, u_coeff) / n_eff
    sat_mask = jnp.isfinite(dn_per_frame) & (dn_per_frame > max_domain)
    d_clamped = jnp.clip(dn_per_frame, 0.0, max_domain)
    s = offsetx + (d_clamped - originx) * scalex
    poly_val = jnp.polyval(coeffs_desc, s)
    e_per_frame = dn_per_frame * poly_val * gain
    in_domain = (dn_per_frame >= 0.0) & (dn_per_frame <= max_domain)
    local_gain = (
        poly_val + jnp.where(in_domain, d_clamped * jnp.polyval(deriv_coeffs_desc, s) * scalex, 0.0)
    ) * gain
    return e_per_frame, sat_mask, local_gain


def _prepare_poly_coeffs_desc_np(channel: LinearityChannel) -> tuple[np.ndarray, np.ndarray]:
    raw_asc = np.asarray(channel.coeffs, dtype=np.float64)
    pad_len = max(8, len(raw_asc))
    coeffs_asc = np.zeros(pad_len, dtype=np.float64)
    coeffs_asc[: len(raw_asc)] = raw_asc
    deriv_asc = coeffs_asc[1:] * np.arange(1, pad_len, dtype=np.float64)
    return coeffs_asc[::-1].copy(), deriv_asc[::-1].copy()


@in_float64
def _calibrate_pre_smear_static_chunk(
    w_buf_2d: Any,
    *,
    u_coeff: float,
    n_eff: float,
    channel: LinearityChannel,
    coeffs_desc_np: np.ndarray,
    deriv_desc_np: np.ndarray,
    gain: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    e_out, sat_out, g_out = _calibrate_pre_smear_2d_jit(
        w_buf_2d,
        jnp.asarray(u_coeff, dtype=jnp.float64),
        jnp.asarray(n_eff, dtype=jnp.float64),
        jnp.asarray(channel.originx, dtype=jnp.float64),
        jnp.asarray(channel.scalex, dtype=jnp.float64),
        jnp.asarray(channel.offsetx, dtype=jnp.float64),
        jnp.asarray(channel.max_domain, dtype=jnp.float64),
        jnp.asarray(coeffs_desc_np, dtype=jnp.float64),
        jnp.asarray(deriv_desc_np, dtype=jnp.float64),
        jnp.asarray(gain, dtype=jnp.float64),
    )
    return (
        np.asarray(e_out, dtype=np.float64),
        np.asarray(sat_out, dtype=bool),
        np.asarray(g_out, dtype=np.float64),
    )


def calibrate_pre_smear_e_per_frame(
    y_adu: Any,
    *,
    u_coeff: float,
    n_eff: float,
    channel: LinearityChannel,
    gain: float,
    reverse_readout: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run shape-bucketed ``@in_float64`` undershoot + ``/ N_eff`` + clamped linearity & gain."""
    arr = np.asarray(y_adu, dtype=np.float64)
    nx = arr.shape[-1]
    if nx == 0 or arr.size == 0:
        return (
            np.empty_like(arr, dtype=np.float64),
            np.zeros(arr.shape, dtype=bool),
            np.empty_like(arr, dtype=np.float64),
        )
    w_2d = (arr[..., ::-1] if reverse_readout else arr).reshape(-1, nx)
    b_tot = w_2d.shape[0]
    w_bucket = _select_bucket(nx, _PRESMEAR_WIDTH_BUCKETS)
    b_step = (
        512
        if b_tot <= 512
        else (131072 if (8192 < b_tot <= 131072 and w_bucket <= 64) else 8192)
    )
    coeffs_desc_np, deriv_desc_np = _prepare_poly_coeffs_desc_np(channel)

    e_2d = np.empty((b_tot, nx), dtype=np.float64)
    sat_2d = np.empty((b_tot, nx), dtype=bool)
    g_2d = np.empty((b_tot, nx), dtype=np.float64)
    buf = np.full((b_step, w_bucket), np.nan, dtype=np.float64)
    for b0 in range(0, b_tot, b_step):
        b1 = min(b0 + b_step, b_tot)
        m = b1 - b0
        buf[:m, :nx] = w_2d[b0:b1]
        if m < b_step:
            buf[m:, :nx] = np.nan
        ec, sc, gc = _calibrate_pre_smear_static_chunk(
            buf,
            u_coeff=u_coeff,
            n_eff=n_eff,
            channel=channel,
            coeffs_desc_np=coeffs_desc_np,
            deriv_desc_np=deriv_desc_np,
            gain=gain,
        )
        e_2d[b0:b1], sat_2d[b0:b1], g_2d[b0:b1] = ec[:m, :nx], sc[:m, :nx], gc[:m, :nx]

    e_full, sat_full, g_full = (
        e_2d.reshape(arr.shape),
        sat_2d.reshape(arr.shape),
        g_2d.reshape(arr.shape),
    )
    if reverse_readout:
        return e_full[..., ::-1].copy(), sat_full[..., ::-1].copy(), g_full[..., ::-1].copy()
    return e_full, sat_full, g_full


def _build_sat_bleed_mask_from_raw_sat(
    raw_sat_2d: np.ndarray, col_outputs: tuple[str, ...]
) -> np.ndarray:
    """Build 2D ``PixelQuality.SAT_BLEED`` mask including downstream undershoot neighbours."""
    ny, nx = raw_sat_2d.shape
    dilated = _dilate_vertical_bool(raw_sat_2d, radius=1)
    downstream = np.zeros((ny, nx), dtype=bool)
    for c in range(nx):
        if not np.any(raw_sat_2d[:, c]):
            continue
        out = col_outputs[c]
        c_down = (c - 1) if OUTPUT_REVERSE_READOUT[out] else (c + 1)
        if 0 <= c_down < nx and col_outputs[c_down] == out:
            downstream[:, c_down] |= raw_sat_2d[:, c]
    return dilated | downstream


# ---------------------------------------------------------------------------
# Collateral streaming calibration & -blf.fits cache helpers
# ---------------------------------------------------------------------------


def _calibrate_tvcol_raw_to_blf(
    tv_col_path: Path,
    blf_dst: Path,
    *,
    output: str,
    black2d_img: np.ndarray,
    sci_hdr0: fits.Header | None = None,
    sci_hdr1: fits.Header | None = None,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
    stage_dir: Path | None = None,
) -> Path:
    tv_hdr0, tv_hdr1, tv_time, tv_cad = read_collateral_headers_and_index(tv_col_path)
    n_total = len(tv_cad)
    n_eff_tv = effective_frame_count(
        tv_hdr0, tv_hdr1, fallback_hdr0=sci_hdr0, fallback_hdr1=sci_hdr1
    )
    fxdoff = float(_hdr_lookup("FXDOFF", tv_hdr1, sci_hdr1) or 209700.0)
    meanblc_val = _hdr_lookup(f"MEANBLC{output}", tv_hdr1, sci_hdr1)
    if meanblc_val is None:
        raise KeyError(f"Missing MEANBLC{output} in collateral and science FITS headers.")

    _, _, tx1, tx2 = OUTPUT_BOUNDS_0BASED[output]
    b2d_tv = np.asarray(black2d_img[:, tx1:tx2], dtype=np.float64)
    bias_shift = -fxdoff + n_eff_tv * (float(meanblc_val) - b2d_tv)[np.newaxis, :, :]
    b1d_all = np.full((n_total, 2078), np.nan, dtype=np.float64)
    b1d_err_all = np.full((n_total, 2078), np.nan, dtype=np.float32)

    for t0, t1, _tc, _cc, raw_chunk in iter_collateral_chunks(
        tv_col_path, "TVCOL_RAW", chunk_rows=chunk_rows
    ):
        b1d_c, err_c = fit_onedblack_legendre(raw_chunk + bias_shift)
        b1d_all[t0:t1], b1d_err_all[t0:t1] = b1d_c, err_c.astype(np.float32)

    return write_onedblack_blf(
        blf_dst,
        time=tv_time,
        cadenceno=tv_cad,
        black1d=b1d_all,
        black1d_err=b1d_err_all,
        n_eff=n_eff_tv,
        time_unit=str(tv_hdr1.get("TUNIT1", "BJD - 2457000, days")),
        stage_dir=stage_dir,
    )


def _calibrate_smrow_chunk_e_per_frame(
    sm_raw_chunk: np.ndarray,
    *,
    output: str,
    fxdoff: float,
    meanblc: float,
    gain: float,
    read_noise_e: float,
    n_eff_sm: float,
    n_eff_tv: float,
    b2d_sm: np.ndarray,
    b1d_sm: np.ndarray,
    b1d_err_sm: np.ndarray | None,
    u_coeff: float,
    lin_channel: LinearityChannel,
    crmiten_var_factor: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    debiased = (
        sm_raw_chunk
        - fxdoff
        + n_eff_sm * (meanblc - b2d_sm)[np.newaxis, :, :]
        - (n_eff_sm / n_eff_tv) * b1d_sm[:, :, np.newaxis]
    )
    sm_e_np, sm_sat_np, sm_gain_np = calibrate_pre_smear_e_per_frame(
        debiased,
        u_coeff=u_coeff,
        n_eff=n_eff_sm,
        channel=lin_channel,
        gain=gain,
        reverse_readout=OUTPUT_REVERSE_READOUT[output],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        smear_med = np.nanmedian(sm_e_np, axis=1)
        mad_rows = 1.4826 * np.nanmedian(np.abs(sm_e_np - smear_med[:, np.newaxis, :]), axis=1)

    n_valid = np.sum(np.isfinite(sm_e_np), axis=1)
    med_var = _unit_normal_median_variance(n_valid)
    emp_var = med_var * (mad_rows * mad_rows)
    theo_var = med_var * (
        float(crmiten_var_factor)
        * (np.maximum(smear_med, 0.0) + read_noise_e * read_noise_e)
        / n_eff_sm
    )
    if b1d_err_sm is not None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean_g = np.nanmedian(sm_gain_np, axis=1)
            b1d_sm_med = np.nanmean(b1d_err_sm, axis=1)[:, np.newaxis]
        mean_g = np.where(np.isfinite(mean_g), mean_g, gain)
        b1d_var = np.where(np.isfinite(b1d_sm_med), ((mean_g * b1d_sm_med) / n_eff_tv) ** 2, 0.0)
    else:
        b1d_var = 0.0

    gated_var = np.where(emp_var > 2.5 * theo_var, emp_var, theo_var)
    smear_err = np.where(n_valid > 0, np.sqrt(gated_var + b1d_var), np.nan)
    return smear_med, smear_err, np.any(sm_sat_np, axis=1)


def _stream_calibrate_smrow_rows(
    sm_col_path: Path,
    *,
    output: str,
    black2d_img: np.ndarray,
    lin_channel: LinearityChannel,
    u_coeff: float,
    b1d_sm_rows: np.ndarray,
    b1d_err_sm_rows: np.ndarray | None,
    n_eff_tv: float,
    target_matched: np.ndarray,
    target_to_sm_idx: np.ndarray,
    sci_hdr0: fits.Header | None = None,
    sci_hdr1: fits.Header | None = None,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sm_hdr0, sm_hdr1, _, _ = read_collateral_headers_and_index(sm_col_path)
    n_eff_sm = effective_frame_count(
        sm_hdr0, sm_hdr1, fallback_hdr0=sci_hdr0, fallback_hdr1=sci_hdr1
    )
    crmiten_var_sm = _crmiten_variance_factor(
        sm_hdr0, sm_hdr1, fallback_hdr0=sci_hdr0, fallback_hdr1=sci_hdr1
    )

    def _hf(key: str, default: float | None = None) -> float:
        val = _hdr_lookup(key, sm_hdr1, sci_hdr1)
        if val is not None:
            return float(val)
        if default is not None:
            return default
        raise KeyError(f"Missing {key} in collateral and science FITS headers.")

    fxdoff, meanblc = _hf("FXDOFF", 209700.0), _hf(f"MEANBLC{output}")
    gain, read_noise_e = _hf(f"GAIN{output}"), _hf(f"READNOI{output}", 8.5)
    sx1, sx2, _, _ = OUTPUT_BOUNDS_0BASED[output]
    b2d_sm = np.asarray(black2d_img[2058:2068, sx1:sx2], dtype=np.float64)

    n_target = len(target_matched)
    smear_out = np.full((n_target, 512), np.nan, dtype=np.float64)
    smear_err_out = np.full((n_target, 512), np.nan, dtype=np.float64)
    smear_sat_out = np.zeros((n_target, 512), dtype=bool)
    for t0, t1, _tc, _cc, raw_chunk in iter_collateral_chunks(
        sm_col_path,
        "SMROW_RAW",
        chunk_rows=chunk_rows,
        row_indices=target_to_sm_idx[target_matched],
    ):
        in_chunk = (target_to_sm_idx >= t0) & (target_to_sm_idx < t1) & target_matched
        if not np.any(in_chunk):
            continue
        dst_sel = np.flatnonzero(in_chunk)
        sm_med, sm_err, sm_sat_2d = _calibrate_smrow_chunk_e_per_frame(
            raw_chunk[target_to_sm_idx[dst_sel] - t0],
            output=output,
            fxdoff=fxdoff,
            meanblc=meanblc,
            gain=gain,
            read_noise_e=read_noise_e,
            n_eff_sm=n_eff_sm,
            n_eff_tv=n_eff_tv,
            b2d_sm=b2d_sm,
            b1d_sm=b1d_sm_rows[dst_sel],
            b1d_err_sm=b1d_err_sm_rows[dst_sel] if b1d_err_sm_rows is not None else None,
            u_coeff=u_coeff,
            lin_channel=lin_channel,
            crmiten_var_factor=crmiten_var_sm,
        )
        smear_out[dst_sel], smear_err_out[dst_sel], smear_sat_out[dst_sel] = (
            sm_med,
            sm_err,
            sm_sat_2d,
        )
    return smear_out, smear_err_out, smear_sat_out


def _calibrate_smrow_raw_to_standard_blf(
    sm_col_path: Path,
    tv_blf_path: Path,
    blf_dst: Path,
    *,
    output: str,
    black2d_img: np.ndarray,
    lin_channel: LinearityChannel,
    u_coeff: float,
    sci_hdr0: fits.Header | None = None,
    sci_hdr1: fits.Header | None = None,
    cadence_match: CadenceMatchMode = "cadenceno",
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
    stage_dir: Path | None = None,
) -> Path:
    sm_hdr0, sm_hdr1, sm_time, sm_cad = read_collateral_headers_and_index(sm_col_path)
    n_eff_sm = effective_frame_count(
        sm_hdr0, sm_hdr1, fallback_hdr0=sci_hdr0, fallback_hdr1=sci_hdr1
    )
    with fits.open(tv_blf_path, memmap=True) as tv_hdul:
        tv_hdu1 = tv_hdul[1]
        tv_cad = np.asarray(tv_hdu1.data["CADENCENO"], dtype=np.int64)
        tv_b1d_sm = np.asarray(tv_hdu1.data["1DBLACK"][:, 2058:2068], dtype=np.float64)
        tv_err_sm = (
            np.asarray(tv_hdu1.data["1DBLACK_ERR"][:, 2058:2068], dtype=np.float64)
            if "1DBLACK_ERR" in (tv_hdu1.data.names or ())
            else None
        )
        n_eff_tv = float(tv_hdu1.header.get("NEFF", n_eff_sm))

    matched_tv, tv_idx = match_cadenceno(
        sm_cad, tv_cad, cadence_match=cadence_match, label=f"tvcol-blf ({output})"
    )
    safe_tv = np.where(matched_tv, tv_idx, 0)
    b1d_al = np.where(matched_tv[:, None], tv_b1d_sm[safe_tv], np.nan)
    err_al = (
        np.where(matched_tv[:, None], tv_err_sm[safe_tv], np.nan) if tv_err_sm is not None else None
    )
    smear_all, smear_err_all, smear_sat_all = _stream_calibrate_smrow_rows(
        sm_col_path,
        output=output,
        black2d_img=black2d_img,
        lin_channel=lin_channel,
        u_coeff=u_coeff,
        b1d_sm_rows=b1d_al,
        b1d_err_sm_rows=err_al,
        n_eff_tv=n_eff_tv,
        target_matched=matched_tv,
        target_to_sm_idx=np.arange(len(sm_cad), dtype=np.int64),
        sci_hdr0=sci_hdr0,
        sci_hdr1=sci_hdr1,
        chunk_rows=chunk_rows,
    )
    tv_all_matched = bool(
        np.all(matched_tv) and np.array_equal(tv_idx, np.arange(len(sm_cad), dtype=np.int64))
    )
    return write_smear_blf(
        blf_dst,
        time=sm_time,
        cadenceno=sm_cad,
        smear=smear_all,
        smear_err=smear_err_all.astype(np.float32),
        smear_sat=smear_sat_all,
        smear_mode="standard",
        cadence_match=cadence_match,
        source_mask_digest="NONE",
        tv_all_matched=tv_all_matched,
        int_time=float(_hdr_lookup("INT_TIME", sm_hdr1, sci_hdr1) or 1.98),
        time_unit=str(sm_hdr1.get("TUNIT1", "BJD - 2457000, days")),
        stage_dir=stage_dir,
    )


def _smear_cache_matches_cutout(
    hdr: fits.Header,
    smear_mode: str,
    cutout_bounds: tuple[int, int, int, int] | None = None,
    *,
    cadence_match: str = "cadenceno",
    source_mask_digest: str = "NONE",
) -> bool:
    if "SMEARMOD" in hdr and str(hdr["SMEARMOD"]).strip() != smear_mode:
        return False
    if str(hdr.get("CADMATCH", "")).strip() != str(cadence_match):
        if not (smear_mode == "standard" and hdr.get("TVMATCH") is True):
            return False
    if str(hdr.get("SRCMASK", "")).strip() != str(source_mask_digest):
        return False
    if smear_mode == "standard":
        return "CRVAL1P" not in hdr and "CUTNX" not in hdr
    if smear_mode == "alternate":
        if str(hdr.get("SMEARMOD", "")).strip() != "alternate" or cutout_bounds is None:
            return False
        for key, val in zip(("CRVAL1P", "CUTNX", "CRVAL2P", "CUTNY"), cutout_bounds, strict=True):
            if key not in hdr or int(hdr[key]) != int(val):
                return False
    return True


def _find_valid_tvcol_blf(
    spec: CollateralSpec, *, ddir: str | os.PathLike[str] | None, cache_dir: Path | None
) -> Path | None:
    for cand in find_all_cached_blf(
        spec, kind="tvcol", smear_mode="standard", ddir=ddir, cache_dir=cache_dir
    ):
        with contextlib.suppress(OSError), fits.open(cand, memmap=True) as hdul:
            if len(hdul) >= 2 and "1DBLACK" in (hdul[1].columns.names or ()):
                return cand
    return None


def _find_valid_smrow_blf(
    spec: CollateralSpec,
    *,
    smear_mode: str,
    ddir: str | os.PathLike[str] | None,
    cache_dir: Path | None,
    cutout_bounds: tuple[int, int, int, int] | None = None,
    require_alt_unres: bool = False,
    required_cadenceno: np.ndarray | None = None,
    cadence_match: CadenceMatchMode = "cadenceno",
    source_mask_digest: str = "NONE",
) -> Path | None:
    req_cad = (
        np.asarray(required_cadenceno, dtype=np.int64) if required_cadenceno is not None else None
    )
    for cand in find_all_cached_blf(
        spec, kind="smrow", smear_mode=smear_mode, ddir=ddir, cache_dir=cache_dir
    ):
        with contextlib.suppress(OSError), fits.open(cand, memmap=True) as hdul:
            if len(hdul) < 2 or "SMEAR" not in (hdul[1].columns.names or ()):
                continue
            if not _smear_cache_matches_cutout(
                hdul[1].header,
                smear_mode,
                cutout_bounds,
                cadence_match=cadence_match,
                source_mask_digest=source_mask_digest,
            ):
                continue
            if require_alt_unres and "ALT_UNRES" not in hdul:
                continue
            if req_cad is not None:
                if "CADENCENO" not in (hdul[1].columns.names or ()):
                    continue
                cached_cad = np.asarray(hdul[1].data["CADENCENO"], dtype=np.int64)
                if cadence_match == "positional":
                    if len(cached_cad) < len(req_cad) or not np.array_equal(
                        cached_cad[: len(req_cad)], req_cad
                    ):
                        continue
                elif not np.all(np.isin(req_cad, cached_cad)):
                    continue
            return cand
    return None


def _alt_blf_has_more_cadences(
    path: Path,
    cutout_bounds: tuple[int, int, int, int],
    n_new: int,
    *,
    cadence_match: CadenceMatchMode = "cadenceno",
    source_mask_digest: str = "NONE",
) -> bool:
    if not path.is_file():
        return False
    with contextlib.suppress(OSError), fits.open(path, memmap=True) as hdul:
        if (
            len(hdul) >= 2
            and "SMEAR" in (hdul[1].columns.names or ())
            and "CADENCENO" in (hdul[1].columns.names or ())
            and "ALT_UNRES" in hdul
            and _smear_cache_matches_cutout(
                hdul[1].header,
                "alternate",
                cutout_bounds,
                cadence_match=cadence_match,
                source_mask_digest=source_mask_digest,
            )
        ):
            return len(hdul[1].data["CADENCENO"]) > n_new
    return False


def build_sector_output_blf_caches(
    *,
    sector: int,
    camera: int,
    ccd: int,
    output: str,
    cadence: CadenceKind = "2min",
    models: StaticModelPaths | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
    ddir: str | os.PathLike[str] | None = None,
    stage_dir: Path | None = None,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
    offline: bool = False,
    strict_hpc: bool = False,
    cleanup_staged_raw: bool = True,
    s3_fs: Any = None,
    verbose: bool = False,
) -> tuple[Path, Path]:
    """Build or reuse ``tvcol`` and ``smrow`` ``-blf.fits`` caches for one CCD output."""
    out_upper = output.upper()
    if out_upper not in OUTPUT_ORDER:
        raise ValueError(f"Invalid CCD output {output!r}; expected one of {OUTPUT_ORDER}.")

    resolved_cache = resolve_cache_dir(
        cache_dir, env_var="MUNDEY_CACHE_DIR", strict_hpc=strict_hpc, create=True
    )
    spec = CollateralSpec(sector, camera, ccd, out_upper, cadence=cadence)
    sec_dir = resolved_cache / spec.sector_subdir

    tv_cached = _find_valid_tvcol_blf(spec, ddir=ddir, cache_dir=resolved_cache)
    sm_cached = _find_valid_smrow_blf(
        spec, smear_mode="standard", ddir=ddir, cache_dir=resolved_cache
    )
    if tv_cached is not None and sm_cached is not None:
        return tv_cached, sm_cached

    if models is None:
        models = resolve_static_models(
            camera,
            ccd,
            ddir=ddir,
            cache_dir=resolved_cache,
            offline=offline,
            stage_dir=stage_dir,
            s3_fs=s3_fs,
            verbose=verbose,
        )
    with fits.open(models.black2d_path, memmap=True) as b2d_hdul:
        black2d_img = np.asarray(b2d_hdul[1].data, dtype=np.float64)
    lin_channel = load_linearity_channel(models.linearity_path, camera, ccd, out_upper)
    u_coeff = load_undershoot_coeff(camera, ccd, out_upper)
    stage_root = stage_dir.resolve(strict=False) if stage_dir is not None else None

    if tv_cached is None:
        tv_raw = resolve_collateral_file(
            spec,
            kind="tvcol",
            ddir=ddir,
            cache_dir=resolved_cache,
            offline=offline,
            stage_dir=stage_dir,
            keep_raw_in_stage_only=cleanup_staged_raw,
            s3_fs=s3_fs,
            verbose=verbose,
        )
        tv_dst = collateral_blf_cache_path(sec_dir, spec, kind="tvcol", col_stem=tv_raw.stem)
        tv_cached = _calibrate_tvcol_raw_to_blf(
            tv_raw,
            tv_dst,
            output=out_upper,
            black2d_img=black2d_img,
            chunk_rows=chunk_rows,
            stage_dir=stage_dir,
        )
        if (
            cleanup_staged_raw
            and stage_root is not None
            and tv_raw.resolve(strict=False).parent == stage_root
        ):
            with contextlib.suppress(OSError):
                tv_raw.unlink()

    if sm_cached is None:
        sm_raw = resolve_collateral_file(
            spec,
            kind="smrow",
            ddir=ddir,
            cache_dir=resolved_cache,
            offline=offline,
            stage_dir=stage_dir,
            keep_raw_in_stage_only=cleanup_staged_raw,
            s3_fs=s3_fs,
            verbose=verbose,
        )
        sm_dst = collateral_blf_cache_path(
            sec_dir, spec, kind="smrow", smear_mode="standard", col_stem=sm_raw.stem
        )
        sm_cached = _calibrate_smrow_raw_to_standard_blf(
            sm_raw,
            tv_cached,
            sm_dst,
            output=out_upper,
            black2d_img=black2d_img,
            lin_channel=lin_channel,
            u_coeff=u_coeff,
            chunk_rows=chunk_rows,
            stage_dir=stage_dir,
        )
        if (
            cleanup_staged_raw
            and stage_root is not None
            and sm_raw.resolve(strict=False).parent == stage_root
        ):
            with contextlib.suppress(OSError):
                sm_raw.unlink()

    return tv_cached, sm_cached


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
    """Prefetch static models and build ``-blf.fits`` caches in ``chunk_rows`` slices."""
    resolved_cache = resolve_cache_dir(
        cache_dir, env_var="MUNDEY_CACHE_DIR", strict_hpc=strict_hpc, create=True
    )
    if stage_dir is not None:
        resolved_stage: Path | None = Path(stage_dir).expanduser()
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
        if out not in OUTPUT_ORDER:
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
                strict_hpc=strict_hpc,
                cleanup_staged_raw=cleanup_staged_raw,
                s3_fs=s3_fs,
                verbose=verbose,
            )
            published.extend([tv_blf, sm_blf])
    return published


# ---------------------------------------------------------------------------
# Decomposed chunked TPF recalibration pipeline
# ---------------------------------------------------------------------------


@dataclass
class _OutputSliceParams:
    output: str
    c0: int
    c1: int
    sm_c0: int
    sm_c1: int
    lin_channel: LinearityChannel
    u_coeff: float
    gain_sci: float
    meanblc_sci: float
    readnoi_sci: float
    b2d_sci: np.ndarray
    b1d_sci_rows: np.ndarray
    b1d_err_sci_rows: np.ndarray
    n_eff_tv: float


@dataclass
class _CalibrationContext:
    hdul: fits.HDUList
    sector: int
    camera: int
    ccd: int
    ticid: int | None
    n_tpf: int
    ny: int
    nx: int
    col0_1based: int
    row0_1based: int
    x1_ccd0: int
    y1_ccd0: int
    cutout_bounds: tuple[int, int, int, int]
    cols_1based: np.ndarray
    col_outputs: tuple[str, ...]
    unique_outputs: tuple[str, ...]
    cadence_kind: CadenceKind
    tpf_parts: tuple[str, ...] | None
    resolved_cache: Path | None
    blf_write_base: Path
    warn_fred_on_write: bool
    black2d_img: np.ndarray
    flat_cutout: np.ndarray
    models: StaticModelPaths
    n_eff_sci: float
    crmiten_var_sci: float
    fxdoff_sci: float
    int_time: float
    tpf_time: np.ndarray
    tpf_cad: np.ndarray
    time_unit: str
    raw_null_val: int
    source_mask_2d: np.ndarray
    srcmask_digest: str
    strap_col_mask: np.ndarray
    quality_mask: np.ndarray


def _raw_cnts_null_val(hdu1: fits.BinTableHDU) -> int:
    for idx, name in enumerate(hdu1.columns.names or (), start=1):
        if name == "RAW_CNTS":
            tnull_key = f"TNULL{idx}"
            if tnull_key in hdu1.header:
                with contextlib.suppress(TypeError, ValueError):
                    return int(hdu1.header[tnull_key])
            break
    return -1


def _extract_raw_cnts_chunk_f64(
    hdu1: fits.BinTableHDU, t0: int, t1: int, null_val: int
) -> np.ndarray:
    raw_arr = hdu1.data["RAW_CNTS"][t0:t1]
    raw_f64 = np.asarray(raw_arr, dtype=np.float64)
    bad = (raw_f64 < 0.0) | (~np.isfinite(raw_f64))
    if np.issubdtype(raw_arr.dtype, np.integer):
        bad |= (raw_arr == null_val) | (raw_arr == -1)
    if np.any(bad):
        raw_f64 = raw_f64.copy()
        raw_f64[bad] = np.nan
    return raw_f64


def _resolve_cache_paths(
    *,
    sector: int,
    cache_dir: str | os.PathLike[str] | None,
    ddir: str | os.PathLike[str] | None,
    strict_hpc: bool,
) -> tuple[Path | None, Path, bool]:
    warn_fred = False
    if cache_dir is not None:
        resolved = resolve_cache_dir(
            cache_dir, env_var="MUNDEY_CACHE_DIR", strict_hpc=strict_hpc, create=True
        )
    elif ddir is not None:
        ddir_path = Path(ddir).expanduser()
        if _is_fred_path(ddir_path):
            env_vars = ("MUNDEY_CACHE_DIR", "TESS_BRIGHT_CACHE_DIR", "SCRATCH", "XDG_CACHE_HOME")
            if any(os.environ.get(v, "").strip() for v in env_vars):
                resolved = resolve_cache_dir(
                    None, env_var="MUNDEY_CACHE_DIR", strict_hpc=strict_hpc, create=True
                )
            elif strict_hpc:
                resolve_cache_dir(
                    ddir_path, env_var="MUNDEY_CACHE_DIR", strict_hpc=True, create=False
                )
                resolved = None
            else:
                resolved, warn_fred = None, True
        else:
            resolved = None
    else:
        resolved = resolve_cache_dir(
            None, env_var="MUNDEY_CACHE_DIR", strict_hpc=strict_hpc, create=True
        )

    blf_base = (
        (resolved / f"s{sector:04d}")
        if resolved is not None
        else Path(ddir).expanduser()  # type: ignore[arg-type]
    )
    return resolved, blf_base, warn_fred


def _resolve_calibration_context(
    hdul: fits.HDUList,
    *,
    source_path: Path | None,
    cache_dir: str | os.PathLike[str] | None,
    ddir: str | os.PathLike[str] | None,
    offline: bool,
    strict_hpc: bool,
    s3_fs: Any,
    verbose: bool,
    quality_bitmask: Any,
    target_obj: Any,
) -> _CalibrationContext:
    hdr0, hdu1, hdr2 = hdul[0].header, hdul[1], hdul[2].header
    hdr1 = hdu1.header
    if "RAW_CNTS" not in (hdu1.data.names or ()):
        raise ValueError(
            "TPF HDU[1] lacks a 'RAW_CNTS' column; use mundey.desmear_calibrated_cutout(tpf) "
            "for TessCut / FFI cutouts without collateral data."
        )

    sector, camera, ccd = int(hdr0["SECTOR"]), int(hdr0["CAMERA"]), int(hdr0["CCD"])
    ticid = int(hdr0["TICID"]) if "TICID" in hdr0 else None
    n_tpf, ny, nx = hdu1.data["RAW_CNTS"].shape

    col0_1based, row0_1based = _stamp_ccd_origin_1based(hdr1, hdr2)
    x1_ccd0, y1_ccd0 = col0_1based - 1, row0_1based - 1
    if row0_1based < 1 or y1_ccd0 + ny > 2048:
        raise ValueError(
            f"TPF rows {row0_1based}..{row0_1based + ny - 1} fall outside "
            "TESS imaging rows 1..2048."
        )
    _sync_stamp_ccd_origin(hdr1, hdr2, col0_1based, row0_1based)

    cols_1based = np.arange(col0_1based, col0_1based + nx, dtype=np.int64)
    col_outputs = tuple(output_for_rawx_col(int(c)) for c in cols_1based)
    unique_outputs = tuple(out for out in OUTPUT_ORDER if out in col_outputs)
    cadence_kind = infer_cadence_kind(tpf_path=source_path, hdu0_header=hdr0, hdu1_header=hdr1)
    parts_list = source_path.name.split("-") if source_path is not None else None
    tpf_parts = tuple(parts_list) if (parts_list is not None and len(parts_list) >= 4) else None

    resolved_cache, blf_write_base, warn_fred = _resolve_cache_paths(
        sector=sector, cache_dir=cache_dir, ddir=ddir, strict_hpc=strict_hpc
    )
    models = resolve_static_models(
        camera,
        ccd,
        ddir=ddir,
        cache_dir=resolved_cache,
        offline=offline,
        s3_fs=s3_fs,
        verbose=verbose,
    )
    with fits.open(models.black2d_path, memmap=True) as b2d_hdul:
        black2d_img = np.asarray(b2d_hdul[1].data, dtype=np.float64)
    with fits.open(models.flat_path, memmap=True) as flat_hdul:
        flat_cutout = np.asarray(
            flat_hdul[1].data[y1_ccd0 : y1_ccd0 + ny, x1_ccd0 : x1_ccd0 + nx], dtype=np.float64
        )

    existing_pq = (
        np.asarray(hdul[PIXEL_QUALITY_EXTNAME].data, dtype=np.uint16)
        if PIXEL_QUALITY_EXTNAME in hdul
        else np.zeros((ny, nx), dtype=np.uint16)
    )
    source_mask_2d = (existing_pq & int(PixelQuality.SOURCE_MASK)) != 0
    return _CalibrationContext(
        hdul=hdul,
        sector=sector,
        camera=camera,
        ccd=ccd,
        ticid=ticid,
        n_tpf=n_tpf,
        ny=ny,
        nx=nx,
        col0_1based=col0_1based,
        row0_1based=row0_1based,
        x1_ccd0=x1_ccd0,
        y1_ccd0=y1_ccd0,
        cutout_bounds=(col0_1based, nx, row0_1based, ny),
        cols_1based=cols_1based,
        col_outputs=col_outputs,
        unique_outputs=unique_outputs,
        cadence_kind=cadence_kind,
        tpf_parts=tpf_parts,
        resolved_cache=resolved_cache,
        blf_write_base=blf_write_base,
        warn_fred_on_write=warn_fred,
        black2d_img=black2d_img,
        flat_cutout=flat_cutout,
        models=models,
        n_eff_sci=effective_frame_count(hdr0, hdr1),
        crmiten_var_sci=_crmiten_variance_factor(hdr0, hdr1),
        fxdoff_sci=float(hdr1.get("FXDOFF", 209700.0)),
        int_time=float(hdr1["INT_TIME"]),
        tpf_time=np.asarray(hdu1.data["TIME"], dtype=np.float64),
        tpf_cad=np.asarray(hdu1.data["CADENCENO"], dtype=np.int64),
        time_unit=str(hdr1.get("TUNIT1", "BJD - 2457000, days")),
        raw_null_val=_raw_cnts_null_val(hdu1),
        source_mask_2d=source_mask_2d,
        srcmask_digest=_source_mask_digest(source_mask_2d),
        strap_col_mask=strap_mask_for_rawx_columns(cols_1based, buffer_px=2),
        quality_mask=_compute_quality_mask(
            hdu1, quality_bitmask=quality_bitmask, target_obj=target_obj
        ),
    )


def _output_slice_bounds(ctx: _CalibrationContext, out: str) -> tuple[int, int, int, int]:
    local_cols = np.flatnonzero(np.asarray([o == out for o in ctx.col_outputs], dtype=bool))
    c0_loc, c1_loc = int(local_cols[0]), int(local_cols[-1]) + 1
    sm_c0 = int(ctx.cols_1based[c0_loc] - 1 - OUTPUT_BOUNDS_0BASED[out][0])
    return c0_loc, c1_loc, sm_c0, sm_c0 + (c1_loc - c0_loc)


def _try_load_alt_caches(
    ctx: _CalibrationContext,
    *,
    ddir: str | os.PathLike[str] | None,
    cadence_match: CadenceMatchMode,
) -> tuple[bool, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    alt_paths: dict[str, Path] = {}
    for out in ctx.unique_outputs:
        spec = CollateralSpec(
            ctx.sector, ctx.camera, ctx.ccd, out, ctx.cadence_kind, ctx.tpf_parts, ctx.ticid
        )
        found = _find_valid_smrow_blf(
            spec,
            smear_mode="alternate",
            ddir=ddir,
            cache_dir=ctx.resolved_cache,
            cutout_bounds=ctx.cutout_bounds,
            require_alt_unres=True,
            required_cadenceno=ctx.tpf_cad,
            cadence_match=cadence_match,
            source_mask_digest=ctx.srcmask_digest,
        )
        if found is None:
            empty = np.empty((0, 0), dtype=np.float64)
            return False, empty, empty, np.zeros(ctx.nx, dtype=bool), np.ones(ctx.n_tpf, dtype=bool)
        alt_paths[out] = found

    hdr1 = ctx.hdul[1].header
    sm_stamp = np.full((ctx.n_tpf, ctx.nx), np.nan, dtype=np.float64)
    err_stamp = np.full((ctx.n_tpf, ctx.nx), np.nan, dtype=np.float64)
    unres_cols = np.zeros(ctx.nx, dtype=bool)
    alt_all_matched = np.ones(ctx.n_tpf, dtype=bool)

    for out in ctx.unique_outputs:
        c0, c1, sm_c0, sm_c1 = _output_slice_bounds(ctx, out)
        readnoi_sci = float(hdr1.get(f"READNOI{out}", 8.5))
        with fits.open(alt_paths[out], memmap=True) as alt_hdul:
            alt_hdu1 = alt_hdul[1]
            alt_cad = np.asarray(alt_hdu1.data["CADENCENO"], dtype=np.int64)
            alt_matched, alt_idx = match_cadenceno(
                ctx.tpf_cad, alt_cad, cadence_match=cadence_match, label=f"smrow-alt-blf ({out})"
            )
            safe_idx = np.where(alt_matched, alt_idx, 0)
            if "COLMATCH" in (alt_hdu1.data.names or ()):
                alt_matched &= np.asarray(alt_hdu1.data["COLMATCH"], dtype=bool)[safe_idx]
            alt_all_matched &= alt_matched

            alt_sm_512 = np.asarray(alt_hdu1.data["SMEAR"], dtype=np.float64)[safe_idx]
            if "SMEAR_ERR" in (alt_hdu1.data.names or ()):
                alt_err_512 = np.asarray(alt_hdu1.data["SMEAR_ERR"], dtype=np.float64)[safe_idx]
            else:
                alt_err_512 = np.sqrt(
                    float(_unit_normal_median_variance(10))
                    * ctx.crmiten_var_sci
                    * (np.maximum(alt_sm_512, 0.0) + readnoi_sci**2)
                    / ctx.n_eff_sci
                )
            sm_stamp[:, c0:c1] = np.where(alt_matched[:, None], alt_sm_512[:, sm_c0:sm_c1], np.nan)
            err_stamp[:, c0:c1] = np.where(
                alt_matched[:, None], alt_err_512[:, sm_c0:sm_c1], np.nan
            )
            unres_cols[c0:c1] = np.asarray(alt_hdul["ALT_UNRES"].data, dtype=bool)[sm_c0:sm_c1]

    return True, sm_stamp, err_stamp, unres_cols, alt_all_matched


@dataclass
class _CollateralData:
    out_params: list[_OutputSliceParams]
    combined_matched: np.ndarray
    collat_smear_stamp: np.ndarray
    collat_smear_err_stamp: np.ndarray
    collat_sm_sat_2d: np.ndarray
    per_out_collat_512: dict[str, np.ndarray]
    per_out_collat_err_512: dict[str, np.ndarray]
    off_target_smears: list[np.ndarray]
    off_target_errs: list[np.ndarray]
    max_e_per_frame_list: list[float]


def _load_output_smrow(
    ctx: _CalibrationContext,
    spec: CollateralSpec,
    *,
    out: str,
    eff_smear: Literal["standard", "alternate"],
    tv_blf: Path,
    b1d_all: np.ndarray,
    b1d_err_all: np.ndarray,
    n_eff_tv: float,
    lin_ch: LinearityChannel,
    u_coeff: float,
    readnoi_sci: float,
    max_sm_e: float,
    ddir: str | os.PathLike[str] | None,
    cadence_match: CadenceMatchMode,
    offline: bool,
    chunk_rows: int,
    s3_fs: Any,
    verbose: bool,
    warn_fred_once: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    hdr0, hdr1 = ctx.hdul[0].header, ctx.hdul[1].header
    std_sm_blf = _find_valid_smrow_blf(
        spec,
        smear_mode="standard",
        ddir=ddir,
        cache_dir=ctx.resolved_cache,
        cutout_bounds=ctx.cutout_bounds,
        cadence_match=cadence_match,
        source_mask_digest="NONE",
    )
    if std_sm_blf is None:
        if eff_smear == "standard":
            warn_fred_once()
        sm_raw = resolve_collateral_file(
            spec,
            kind="smrow",
            ddir=ddir,
            cache_dir=ctx.resolved_cache,
            offline=offline,
            s3_fs=s3_fs,
            verbose=verbose,
        )
        if eff_smear == "standard":
            sm_dst = collateral_blf_cache_path(
                ctx.blf_write_base,
                spec,
                kind="smrow",
                smear_mode="standard",
                col_stem=sm_raw.stem if ctx.tpf_parts is None else None,
            )
            std_sm_blf = _calibrate_smrow_raw_to_standard_blf(
                sm_raw,
                tv_blf,
                sm_dst,
                output=out,
                black2d_img=ctx.black2d_img,
                lin_channel=lin_ch,
                u_coeff=u_coeff,
                sci_hdr0=hdr0,
                sci_hdr1=hdr1,
                cadence_match=cadence_match,
                chunk_rows=chunk_rows,
            )
        else:
            _, _, _, sm_cad = read_collateral_headers_and_index(sm_raw)
            sm_matched, sm_idx = match_cadenceno(
                ctx.tpf_cad, sm_cad, cadence_match=cadence_match, label=f"smrow ({out})"
            )
            sm_512, sm_err_512, sm_sat_cad = _stream_calibrate_smrow_rows(
                sm_raw,
                output=out,
                black2d_img=ctx.black2d_img,
                lin_channel=lin_ch,
                u_coeff=u_coeff,
                b1d_sm_rows=b1d_all[:, 2058:2068],
                b1d_err_sm_rows=b1d_err_all[:, 2058:2068],
                n_eff_tv=n_eff_tv,
                target_matched=sm_matched,
                target_to_sm_idx=sm_idx,
                sci_hdr0=hdr0,
                sci_hdr1=hdr1,
                chunk_rows=chunk_rows,
            )
            return sm_matched, sm_512, sm_err_512, sm_sat_cad

    with fits.open(std_sm_blf, memmap=True) as sm_hdul:
        sm_hdu1 = sm_hdul[1]
        sm_cad = np.asarray(sm_hdu1.data["CADENCENO"], dtype=np.int64)
        sm_matched, sm_idx = match_cadenceno(
            ctx.tpf_cad, sm_cad, cadence_match=cadence_match, label=f"smrow ({out})"
        )
        safe_sm = np.where(sm_matched, sm_idx, 0)
        sm_512 = np.asarray(sm_hdu1.data["SMEAR"], dtype=np.float64)[safe_sm]
        if "SMEAR_ERR" in (sm_hdu1.data.names or ()):
            sm_err_512 = np.asarray(sm_hdu1.data["SMEAR_ERR"], dtype=np.float64)[safe_sm]
        else:
            sm_err_512 = np.sqrt(
                float(_unit_normal_median_variance(10))
                * ctx.crmiten_var_sci
                * (np.maximum(sm_512, 0.0) + readnoi_sci**2)
                / ctx.n_eff_sci
            )
        sm_512 = np.where(sm_matched[:, None], sm_512, np.nan)
        sm_err_512 = np.where(sm_matched[:, None], sm_err_512, np.nan)
        if "SMROW_SAT" in sm_hdul:
            sm_sat_cad = np.where(
                sm_matched[:, None],
                np.asarray(sm_hdul["SMROW_SAT"].data, dtype=bool)[safe_sm],
                False,
            )
        else:
            sm_sat_cad = np.isfinite(sm_512) & (sm_512 >= max_sm_e * (1.0 - 1e-12))
    return sm_matched, sm_512, sm_err_512, sm_sat_cad


def _load_output_collateral(
    ctx: _CalibrationContext,
    *,
    eff_smear: Literal["standard", "alternate"],
    all_outputs_have_alt_cache: bool,
    ddir: str | os.PathLike[str] | None,
    cadence_match: CadenceMatchMode,
    offline: bool,
    chunk_rows: int,
    s3_fs: Any,
    verbose: bool,
) -> _CollateralData:
    hdr0, hdr1 = ctx.hdul[0].header, ctx.hdul[1].header
    combined_matched = np.ones(ctx.n_tpf, dtype=bool)
    collat_smear_stamp = np.full((ctx.n_tpf, ctx.nx), np.nan, dtype=np.float64)
    collat_smear_err_stamp = np.full((ctx.n_tpf, ctx.nx), np.nan, dtype=np.float64)
    collat_sm_sat_2d = np.zeros((ctx.n_tpf, ctx.nx), dtype=bool)

    out_params: list[_OutputSliceParams] = []
    per_out_collat_512: dict[str, np.ndarray] = {}
    per_out_collat_err_512: dict[str, np.ndarray] = {}
    off_target_smears: list[np.ndarray] = []
    off_target_errs: list[np.ndarray] = []
    max_e_per_frame_list: list[float] = []
    warned_fred = False

    def _warn_fred_once() -> None:
        nonlocal warned_fred
        if ctx.warn_fred_on_write and not warned_fred:
            resolve_cache_dir(
                ctx.blf_write_base, env_var="MUNDEY_CACHE_DIR", strict_hpc=False, create=False
            )
            warned_fred = True

    for out in ctx.unique_outputs:
        c0, c1, sm_c0, sm_c1 = _output_slice_bounds(ctx, out)
        spec = CollateralSpec(
            ctx.sector, ctx.camera, ctx.ccd, out, ctx.cadence_kind, ctx.tpf_parts, ctx.ticid
        )
        lin_ch = load_linearity_channel(ctx.models.linearity_path, ctx.camera, ctx.ccd, out)
        u_coeff = load_undershoot_coeff(ctx.camera, ctx.ccd, out)
        gain_sci = float(hdr1[f"GAIN{out}"])
        meanblc_sci = float(hdr1[f"MEANBLC{out}"])
        readnoi_sci = float(hdr1.get(f"READNOI{out}", 8.5))
        max_sm_e = linearity_max_electrons_per_frame(lin_ch, gain_sci)
        max_e_per_frame_list.append(max_sm_e)

        tv_blf = _find_valid_tvcol_blf(spec, ddir=ddir, cache_dir=ctx.resolved_cache)
        if tv_blf is None:
            _warn_fred_once()
            tv_raw = resolve_collateral_file(
                spec,
                kind="tvcol",
                ddir=ddir,
                cache_dir=ctx.resolved_cache,
                offline=offline,
                s3_fs=s3_fs,
                verbose=verbose,
            )
            tv_dst = collateral_blf_cache_path(
                ctx.blf_write_base,
                spec,
                kind="tvcol",
                col_stem=tv_raw.stem if ctx.tpf_parts is None else None,
            )
            tv_blf = _calibrate_tvcol_raw_to_blf(
                tv_raw,
                tv_dst,
                output=out,
                black2d_img=ctx.black2d_img,
                sci_hdr0=hdr0,
                sci_hdr1=hdr1,
                chunk_rows=chunk_rows,
            )

        with fits.open(tv_blf, memmap=True) as tv_hdul:
            tv_hdu1 = tv_hdul[1]
            tv_cad = np.asarray(tv_hdu1.data["CADENCENO"], dtype=np.int64)
            tv_matched, tv_idx = match_cadenceno(
                ctx.tpf_cad, tv_cad, cadence_match=cadence_match, label=f"tvcol ({out})"
            )
            combined_matched &= tv_matched
            safe_tv = np.where(tv_matched, tv_idx, 0)
            b1d_all = np.asarray(tv_hdu1.data["1DBLACK"], dtype=np.float64)[safe_tv]
            b1d_err_all = (
                np.asarray(tv_hdu1.data["1DBLACK_ERR"], dtype=np.float64)[safe_tv]
                if "1DBLACK_ERR" in (tv_hdu1.data.names or ())
                else np.full_like(b1d_all, 0.5)
            )
            n_eff_tv = float(tv_hdu1.header.get("NEFF", ctx.n_eff_sci))

        out_params.append(
            _OutputSliceParams(
                output=out,
                c0=c0,
                c1=c1,
                sm_c0=sm_c0,
                sm_c1=sm_c1,
                lin_channel=lin_ch,
                u_coeff=u_coeff,
                gain_sci=gain_sci,
                meanblc_sci=meanblc_sci,
                readnoi_sci=readnoi_sci,
                b2d_sci=ctx.black2d_img[
                    ctx.y1_ccd0 : ctx.y1_ccd0 + ctx.ny, ctx.x1_ccd0 + c0 : ctx.x1_ccd0 + c1
                ],
                b1d_sci_rows=np.where(
                    tv_matched[:, None], b1d_all[:, ctx.y1_ccd0 : ctx.y1_ccd0 + ctx.ny], np.nan
                ),
                b1d_err_sci_rows=np.where(
                    tv_matched[:, None], b1d_err_all[:, ctx.y1_ccd0 : ctx.y1_ccd0 + ctx.ny], np.nan
                ),
                n_eff_tv=n_eff_tv,
            )
        )
        if all_outputs_have_alt_cache:
            continue

        sm_matched, sm_512, sm_err_512, sm_sat_cad = _load_output_smrow(
            ctx,
            spec,
            out=out,
            eff_smear=eff_smear,
            tv_blf=tv_blf,
            b1d_all=b1d_all,
            b1d_err_all=b1d_err_all,
            n_eff_tv=n_eff_tv,
            lin_ch=lin_ch,
            u_coeff=u_coeff,
            readnoi_sci=readnoi_sci,
            max_sm_e=max_sm_e,
            ddir=ddir,
            cadence_match=cadence_match,
            offline=offline,
            chunk_rows=chunk_rows,
            s3_fs=s3_fs,
            verbose=verbose,
            warn_fred_once=_warn_fred_once,
        )
        combined_matched &= sm_matched
        per_out_collat_512[out], per_out_collat_err_512[out] = sm_512, sm_err_512
        collat_smear_stamp[:, c0:c1] = sm_512[:, sm_c0:sm_c1]
        collat_smear_err_stamp[:, c0:c1] = sm_err_512[:, sm_c0:sm_c1]
        collat_sm_sat_2d[:, c0:c1] = sm_sat_cad[:, sm_c0:sm_c1]

        good_sm = sm_matched & ctx.quality_mask
        sm_sat_512 = np.any(sm_sat_cad[good_sm] if np.any(good_sm) else sm_sat_cad, axis=0)
        cols_512 = np.arange(512)
        off_mask = ((cols_512 < sm_c0 - 10) | (cols_512 >= sm_c1 + 10)) & (~sm_sat_512)
        n_off = int(np.sum(off_mask))
        if n_off >= 8:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                off_med = np.nanmedian(sm_512[:, off_mask], axis=1)
                off_mad = 1.4826 * np.nanmedian(
                    np.abs(sm_512[:, off_mask] - off_med[:, None]), axis=1
                )
                se_theo = np.sqrt(float(_unit_normal_median_variance(n_off))) * np.nanmedian(
                    sm_err_512[:, off_mask], axis=1
                )
            off_target_smears.append(off_med)
            safe_se = np.where(np.isfinite(se_theo), se_theo, 0.0)
            off_target_errs.append(np.maximum(off_mad, safe_se))

    return _CollateralData(
        out_params=out_params,
        combined_matched=combined_matched,
        collat_smear_stamp=collat_smear_stamp,
        collat_smear_err_stamp=collat_smear_err_stamp,
        collat_sm_sat_2d=collat_sm_sat_2d,
        per_out_collat_512=per_out_collat_512,
        per_out_collat_err_512=per_out_collat_err_512,
        off_target_smears=off_target_smears,
        off_target_errs=off_target_errs,
        max_e_per_frame_list=max_e_per_frame_list,
    )


def _calibrate_presmear_chunk(
    ctx: _CalibrationContext,
    out_params: list[_OutputSliceParams],
    combined_matched: np.ndarray,
    t0: int,
    t1: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m = t1 - t0
    raw_chunk = _extract_raw_cnts_chunk_f64(ctx.hdul[1], t0, t1, ctx.raw_null_val)
    e_chunk = np.full((m, ctx.ny, ctx.nx), np.nan, dtype=np.float64)
    sig_chunk = np.full((m, ctx.ny, ctx.nx), np.nan, dtype=np.float64)
    sat_chunk = np.zeros((m, ctx.ny, ctx.nx), dtype=bool)
    matched_c = combined_matched[t0:t1]

    for op in out_params:
        sci_debiased = (
            raw_chunk[:, :, op.c0 : op.c1]
            - ctx.fxdoff_sci
            + ctx.n_eff_sci * (op.meanblc_sci - op.b2d_sci)[np.newaxis, :, :]
            - (ctx.n_eff_sci / op.n_eff_tv) * op.b1d_sci_rows[t0:t1, :, np.newaxis]
        )
        e_out, sat_out, g_out = calibrate_pre_smear_e_per_frame(
            sci_debiased,
            u_coeff=op.u_coeff,
            n_eff=ctx.n_eff_sci,
            channel=op.lin_channel,
            gain=op.gain_sci,
            reverse_readout=OUTPUT_REVERSE_READOUT[op.output],
        )
        n_loc = op.c1 - op.c0
        u_mult_pix = np.full(n_loc, 1.0 / (1.0 - op.u_coeff * op.u_coeff), dtype=np.float64)
        u_mult_pix[(n_loc - 1) if OUTPUT_REVERSE_READOUT[op.output] else 0] = 1.0 / (
            (1.0 - op.u_coeff) ** 2
        )
        u_mult_b1d = 1.0 / ((1.0 - op.u_coeff) ** 2)
        shot_read_var = (
            ctx.crmiten_var_sci
            * (np.maximum(e_out, 0.0) + op.readnoi_sci * op.readnoi_sci)
            / ctx.n_eff_sci
        )
        b1d_var = ((g_out * op.b1d_err_sci_rows[t0:t1, :, np.newaxis]) / op.n_eff_tv) ** 2
        b1d_var = np.where(np.isfinite(b1d_var), b1d_var, 0.0)
        sig_out = np.sqrt(u_mult_pix[None, None, :] * shot_read_var + u_mult_b1d * b1d_var)
        e_chunk[:, :, op.c0 : op.c1] = e_out
        sig_chunk[:, :, op.c0 : op.c1] = sig_out
        sat_chunk[:, :, op.c0 : op.c1] = sat_out

    if not np.all(matched_c):
        e_chunk[~matched_c], sig_chunk[~matched_c], sat_chunk[~matched_c] = np.nan, np.nan, False
    return e_chunk, sig_chunk, sat_chunk


def _persist_alt_caches(
    ctx: _CalibrationContext,
    *,
    per_out_collat_512: dict[str, np.ndarray],
    per_out_collat_err_512: dict[str, np.ndarray],
    m_hat_stamp: np.ndarray,
    m_err_stamp: np.ndarray,
    unresolved_cols: np.ndarray,
    combined_matched: np.ndarray,
    cadence_match: CadenceMatchMode,
) -> None:
    if ctx.warn_fred_on_write:
        resolve_cache_dir(
            ctx.blf_write_base, env_var="MUNDEY_CACHE_DIR", strict_hpc=False, create=False
        )
    for out in ctx.unique_outputs:
        if out not in per_out_collat_512:
            continue
        c0, c1, sm_c0, sm_c1 = _output_slice_bounds(ctx, out)
        alt_sm_512 = per_out_collat_512[out].copy()
        alt_err_512 = per_out_collat_err_512[out].copy()
        alt_sm_512[:, sm_c0:sm_c1] = m_hat_stamp[:, c0:c1]
        alt_err_512[:, sm_c0:sm_c1] = m_err_stamp[:, c0:c1]
        alt_unres_512 = np.zeros(512, dtype=bool)
        alt_unres_512[sm_c0:sm_c1] = unresolved_cols[c0:c1]

        spec = CollateralSpec(
            ctx.sector, ctx.camera, ctx.ccd, out, ctx.cadence_kind, ctx.tpf_parts, ctx.ticid
        )
        alt_dst = collateral_blf_cache_path(
            ctx.blf_write_base,
            spec,
            kind="smrow",
            smear_mode="alternate",
            col_stem=(
                f"tess-s{ctx.sector:04d}-smrow-{ctx.camera}-{ctx.ccd}-{out.lower()}-"
                f"{int(ctx.ticid):016d}-{collateral_suffix(ctx.cadence_kind)}"
                if ctx.ticid is not None and ctx.tpf_parts is None
                else None
            ),
        )
        if _alt_blf_has_more_cadences(
            alt_dst,
            ctx.cutout_bounds,
            len(ctx.tpf_cad),
            cadence_match=cadence_match,
            source_mask_digest=ctx.srcmask_digest,
        ):
            continue
        try:
            write_smear_blf(
                alt_dst,
                time=ctx.tpf_time,
                cadenceno=ctx.tpf_cad,
                smear=alt_sm_512,
                smear_err=alt_err_512.astype(np.float32),
                unresolved_cols=alt_unres_512,
                col_matched=combined_matched,
                smear_mode="alternate",
                cutout_bounds=ctx.cutout_bounds,
                cadence_match=cadence_match,
                source_mask_digest=ctx.srcmask_digest,
                int_time=ctx.int_time,
                time_unit=ctx.time_unit,
            )
        except OSError as exc:
            if exc.errno in (errno.EACCES, errno.EPERM, errno.EROFS):
                warnings.warn(
                    f"Skipping alternate smear cache write to read-only path {alt_dst}: {exc}",
                    UserWarning,
                    stacklevel=2,
                )
            else:
                raise


def _write_calibrated_chunk(
    hdu1: fits.BinTableHDU,
    *,
    t0: int,
    t1: int,
    e_chunk: np.ndarray,
    sig_chunk: np.ndarray,
    m_hat_c: np.ndarray,
    cal_smear_err_c: np.ndarray,
    denom_rate: np.ndarray,
    matched_c: np.ndarray,
) -> None:
    cal_flux = (e_chunk - m_hat_c[:, np.newaxis, :]) / denom_rate
    cal_flux_err = np.hypot(sig_chunk / denom_rate, cal_smear_err_c[:, np.newaxis, :])
    if not np.all(matched_c):
        cal_flux[~matched_c], cal_flux_err[~matched_c] = np.nan, np.nan
    hdu1.data["FLUX"][t0:t1] = cal_flux.astype(np.float32)
    if "FLUX_ERR" in (hdu1.data.names or ()):
        hdu1.data["FLUX_ERR"][t0:t1] = cal_flux_err.astype(np.float32)


def _finalize_calibrated_hdulist(
    ctx: _CalibrationContext,
    *,
    cal_smear_err: np.ndarray,
    pq_out: np.ndarray,
    eff_smear: str,
    cadence_match: CadenceMatchMode,
    n_missing: int,
    target_obj: Any,
) -> None:
    hdul = ctx.hdul
    col_names = hdul[1].data.names or ()
    if "FLUX_BKG" in col_names:
        hdul[1].data["FLUX_BKG"] = np.float32(0.0)
    if "FLUX_BKG_ERR" in col_names:
        hdul[1].data["FLUX_BKG_ERR"] = np.float32(0.0)

    append_smear_err_column(hdul, cal_smear_err.astype(np.float32))
    write_pixel_quality(hdul, pq_out)
    write_mundey_provenance(
        target_obj if target_obj is not None else hdul,
        bkgsub=False,
        smear_mode=eff_smear,
    )
    hdul[0].header["CADMATCH"] = (cadence_match, "Collateral cadence pairing mode")
    hdul[0].header["CADMISS"] = (int(n_missing), "Unmatched TPF cadences masked to NaN")
    if target_obj is not None:
        meta_dict = getattr(getattr(target_obj, "meta", None), "_dict", None)
        if isinstance(meta_dict, dict):
            meta_dict["CADMATCH"] = cadence_match
            meta_dict["CADMISS"] = int(n_missing)


def _scan_pass1_bloom(
    ctx: _CalibrationContext,
    out_params: list[_OutputSliceParams],
    combined_matched: np.ndarray,
    good_cad: np.ndarray,
    collat_smear_stamp: np.ndarray,
    chunk_rows: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[np.ndarray, np.ndarray] | None]:
    """Stream chunks to accumulate quality-filtered ``raw_sat_2d`` & exact 2D bloom stats."""
    raw_sat_2d = np.zeros((ctx.ny, ctx.nx), dtype=bool)
    w_sum = np.zeros((ctx.ny, ctx.nx), dtype=np.float64)
    w_num = np.zeros((ctx.ny, ctx.nx), dtype=np.float64)
    n_pix_2d = np.zeros((ctx.ny, ctx.nx), dtype=np.int64)
    min_z_2d = np.full((ctx.ny, ctx.nx), np.inf, dtype=np.float64)
    single_chunk_cache: tuple[np.ndarray, np.ndarray] | None = None

    for t0 in range(0, ctx.n_tpf, chunk_rows):
        t1 = min(t0 + chunk_rows, ctx.n_tpf)
        e_c, sig_c, sat_c = _calibrate_presmear_chunk(ctx, out_params, combined_matched, t0, t1)
        if ctx.n_tpf <= chunk_rows:
            single_chunk_cache = (e_c, sig_c)
        gc = good_cad[t0:t1]
        if np.any(gc):
            raw_sat_2d |= np.any(sat_c[gc], axis=0)
            sub_g = e_c[gc] - collat_smear_stamp[t0:t1][gc, None, :]
            sig_g = sig_c[gc]
            val_g = np.isfinite(sub_g) & np.isfinite(sig_g)
            sig_safe = np.maximum(np.where(val_g, sig_g, 1.0), 1e-6)
            wg = np.where(val_g, 1.0 / (sig_safe * sig_safe), 0.0)
            w_sum += np.sum(wg, axis=0)
            w_num += np.sum(np.where(val_g, wg * sub_g, 0.0), axis=0)
            n_pix_2d += np.sum(val_g, axis=0)
            z_g = np.where(val_g, sub_g / sig_safe, np.inf)
            min_z_2d = np.minimum(min_z_2d, np.min(z_g, axis=0))

    sat_bleed_2d = _build_sat_bleed_mask_from_raw_sat(raw_sat_2d, ctx.col_outputs)
    clean_rows_2d = _build_clean_rows_2d(raw_sat_2d, sat_bleed_2d, ctx.source_mask_2d)
    col_has_bloom = _eval_bloom_from_stats(w_num, w_sum, n_pix_2d, min_z_2d, clean_rows_2d)
    return raw_sat_2d, sat_bleed_2d, col_has_bloom, single_chunk_cache


def _run_cold_alternate_pass(
    ctx: _CalibrationContext,
    *,
    out_params: list[_OutputSliceParams],
    combined_matched: np.ndarray,
    good_cad: np.ndarray,
    collat_smear_stamp: np.ndarray,
    collat_smear_err_stamp: np.ndarray,
    collat_sm_sat_cols: np.ndarray,
    off_target_smears: list[np.ndarray],
    off_target_errs: list[np.ndarray],
    max_e_per_frame_list: list[float],
    raw_sat_2d: np.ndarray,
    sat_bleed_2d: np.ndarray,
    col_has_bloom: np.ndarray,
    single_chunk_cache: tuple[np.ndarray, np.ndarray] | None,
    denom_rate: np.ndarray,
    smear_rate_denom: np.ndarray,
    chunk_rows: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if off_target_smears:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            off_sm_vec: np.ndarray | None = np.nanmean(np.vstack(off_target_smears), axis=0)
            off_err_vec: np.ndarray | None = np.nanmean(np.vstack(off_target_errs), axis=0)
    else:
        off_sm_vec, off_err_vec = None, None

    alt_cfg = AlternateSmearConfig(
        row0_1based=ctx.row0_1based,
        row_smear_weight=_frame_transfer_row_weight(ctx.hdul[0].header, ctx.hdul[1].header),
        saturation_limit=float(np.mean(max_e_per_frame_list)) if max_e_per_frame_list else None,
    )
    c_bleed_mask = np.any(raw_sat_2d, axis=0) | collat_sm_sat_cols | col_has_bloom

    if single_chunk_cache is not None:
        e_all, sig_all = single_chunk_cache
        m_hat_stamp, m_err_stamp, unresolved_cols, _ = _estimate_alternate_smear_stamp(
            sci_pre_smear=e_all,
            sci_sigma_pre=sig_all,
            flat_2d=ctx.flat_cutout,
            collat_smear_stamp=collat_smear_stamp,
            collat_smear_err_stamp=collat_smear_err_stamp,
            collat_sm_sat_cols=collat_sm_sat_cols,
            off_target_smear=off_sm_vec,
            off_target_smear_err=off_err_vec,
            raw_sat_2d=raw_sat_2d,
            sat_bleed_2d=sat_bleed_2d,
            strap_col_mask=ctx.strap_col_mask,
            source_mask_2d=ctx.source_mask_2d,
            config=alt_cfg,
            good_cad_mask=good_cad,
            precomputed_bleed_mask=c_bleed_mask,
        )
        cal_smear_err = m_err_stamp / smear_rate_denom
        if not np.all(combined_matched):
            cal_smear_err[~combined_matched] = np.nan
        _write_calibrated_chunk(
            ctx.hdul[1],
            t0=0,
            t1=ctx.n_tpf,
            e_chunk=e_all,
            sig_chunk=sig_all,
            m_hat_c=m_hat_stamp,
            cal_smear_err_c=cal_smear_err,
            denom_rate=denom_rate,
            matched_c=combined_matched,
        )
        return m_hat_stamp, m_err_stamp, cal_smear_err, unresolved_cols

    m_hat_stamp = np.full((ctx.n_tpf, ctx.nx), np.nan, dtype=np.float64)
    m_err_stamp = np.full((ctx.n_tpf, ctx.nx), np.nan, dtype=np.float64)
    cal_smear_err = np.full((ctx.n_tpf, ctx.nx), np.nan, dtype=np.float64)
    unresolved_cols = np.zeros(ctx.nx, dtype=bool)
    for t0 in range(0, ctx.n_tpf, chunk_rows):
        t1 = min(t0 + chunk_rows, ctx.n_tpf)
        e_c, sig_c, _ = _calibrate_presmear_chunk(ctx, out_params, combined_matched, t0, t1)
        mh_c, me_c, unres_c, _ = _estimate_alternate_smear_stamp(
            sci_pre_smear=e_c,
            sci_sigma_pre=sig_c,
            flat_2d=ctx.flat_cutout,
            collat_smear_stamp=collat_smear_stamp[t0:t1],
            collat_smear_err_stamp=collat_smear_err_stamp[t0:t1],
            collat_sm_sat_cols=collat_sm_sat_cols,
            off_target_smear=off_sm_vec[t0:t1] if off_sm_vec is not None else None,
            off_target_smear_err=off_err_vec[t0:t1] if off_err_vec is not None else None,
            raw_sat_2d=raw_sat_2d,
            sat_bleed_2d=sat_bleed_2d,
            strap_col_mask=ctx.strap_col_mask,
            source_mask_2d=ctx.source_mask_2d,
            config=alt_cfg,
            good_cad_mask=good_cad[t0:t1],
            precomputed_bleed_mask=c_bleed_mask,
        )
        m_hat_stamp[t0:t1], m_err_stamp[t0:t1] = mh_c, me_c
        unresolved_cols |= unres_c
        cse_c = me_c / smear_rate_denom
        if not np.all(combined_matched[t0:t1]):
            cse_c[~combined_matched[t0:t1]] = np.nan
        cal_smear_err[t0:t1] = cse_c
        _write_calibrated_chunk(
            ctx.hdul[1],
            t0=t0,
            t1=t1,
            e_chunk=e_c,
            sig_chunk=sig_c,
            m_hat_c=mh_c,
            cal_smear_err_c=cse_c,
            denom_rate=denom_rate,
            matched_c=combined_matched[t0:t1],
        )
    return m_hat_stamp, m_err_stamp, cal_smear_err, unresolved_cols


def _apply_fixed_smear_chunks(
    ctx: _CalibrationContext,
    out_params: list[_OutputSliceParams],
    combined_matched: np.ndarray,
    m_hat_stamp: np.ndarray,
    m_err_stamp: np.ndarray,
    denom_rate: np.ndarray,
    smear_rate_denom: np.ndarray,
    chunk_rows: int,
    *,
    good_cad: np.ndarray | None = None,
    single_chunk_cache: tuple[np.ndarray, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    cal_smear_err = m_err_stamp / smear_rate_denom
    if not np.all(combined_matched):
        cal_smear_err[~combined_matched] = np.nan
    raw_sat_2d = np.zeros((ctx.ny, ctx.nx), dtype=bool) if good_cad is not None else None
    for t0 in range(0, ctx.n_tpf, chunk_rows):
        t1 = min(t0 + chunk_rows, ctx.n_tpf)
        if single_chunk_cache is not None:
            e_c, sig_c = single_chunk_cache
        else:
            e_c, sig_c, sat_c = _calibrate_presmear_chunk(ctx, out_params, combined_matched, t0, t1)
            if raw_sat_2d is not None and good_cad is not None and np.any(good_cad[t0:t1]):
                raw_sat_2d |= np.any(sat_c[good_cad[t0:t1]], axis=0)
        _write_calibrated_chunk(
            ctx.hdul[1],
            t0=t0,
            t1=t1,
            e_chunk=e_c,
            sig_chunk=sig_c,
            m_hat_c=m_hat_stamp[t0:t1],
            cal_smear_err_c=cal_smear_err[t0:t1],
            denom_rate=denom_rate,
            matched_c=combined_matched[t0:t1],
        )
    return cal_smear_err, raw_sat_2d


def _recalibrate_hdulist_inplace(
    hdul: fits.HDUList,
    *,
    source_path: Path | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
    ddir: str | os.PathLike[str] | None = None,
    smear: SmearMode = "automatic",
    cadence_match: CadenceMatchMode = "cadenceno",
    missing_cadence: MissingCadenceMode = "nan",
    offline: bool = False,
    strict_hpc: bool = False,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
    s3_fs: Any = None,
    verbose: bool = False,
    quality_bitmask: Any = "default",
    target_obj: Any = None,
) -> None:
    """Execute the chunked v2.0 calibration pipeline on ``hdul`` in place."""
    if smear not in ("automatic", "standard", "alternate"):
        raise ValueError(
            f"Invalid smear mode {smear!r}; expected 'automatic', 'standard', or 'alternate'."
        )
    if cadence_match not in ("cadenceno", "positional"):
        raise ValueError(
            f"Invalid cadence_match={cadence_match!r}; expected 'cadenceno' or 'positional'."
        )
    if missing_cadence not in ("nan", "error"):
        raise ValueError(
            f"Invalid missing_cadence={missing_cadence!r}; expected 'nan' or 'error'."
        )

    ctx = _resolve_calibration_context(
        hdul,
        source_path=source_path,
        cache_dir=cache_dir,
        ddir=ddir,
        offline=offline,
        strict_hpc=strict_hpc,
        s3_fs=s3_fs,
        verbose=verbose,
        quality_bitmask=quality_bitmask,
        target_obj=target_obj,
    )
    defer_auto_smear = (smear == "automatic") and (ctx.y1_ccd0 + ctx.ny < 2048)
    eff_smear: Literal["standard", "alternate"] = (
        ("alternate" if (ctx.y1_ccd0 + ctx.ny >= 2048) else "standard")
        if smear == "automatic"
        else smear  # type: ignore[assignment]
    )
    if verbose:
        print(
            f"Recalibrating Sector {ctx.sector} Camera {ctx.camera} CCD {ctx.ccd} "
            f"Outputs {list(ctx.unique_outputs)} "
            f"(smear={smear if defer_auto_smear else eff_smear!r}, cadence={ctx.cadence_kind!r})"
        )

    all_outputs_have_alt_cache = False
    alt_sm_stamp = np.empty((0, 0), dtype=np.float64)
    alt_err_stamp = np.empty((0, 0), dtype=np.float64)
    alt_unres_cols = np.zeros(ctx.nx, dtype=bool)
    alt_matched = np.ones(ctx.n_tpf, dtype=bool)
    if eff_smear == "alternate":
        all_outputs_have_alt_cache, alt_sm_stamp, alt_err_stamp, alt_unres_cols, alt_matched = (
            _try_load_alt_caches(ctx, ddir=ddir, cadence_match=cadence_match)
        )

    collat = _load_output_collateral(
        ctx,
        eff_smear=eff_smear,
        all_outputs_have_alt_cache=all_outputs_have_alt_cache,
        ddir=ddir,
        cadence_match=cadence_match,
        offline=offline,
        chunk_rows=chunk_rows,
        s3_fs=s3_fs,
        verbose=verbose,
    )
    combined_matched = (
        collat.combined_matched & alt_matched
        if all_outputs_have_alt_cache
        else collat.combined_matched
    )

    n_missing = _enforce_missing_cadence_policy(
        ctx.tpf_cad, combined_matched, missing_cadence=missing_cadence
    )
    good_cad = combined_matched & ctx.quality_mask
    if not np.any(good_cad):
        good_cad = combined_matched
    collat_sm_sat_cols = np.any(collat.collat_sm_sat_2d[good_cad], axis=0)

    safe_flat = np.where(
        np.isfinite(ctx.flat_cutout) & (ctx.flat_cutout > 0.0), ctx.flat_cutout, np.nan
    )
    denom_rate = safe_flat[np.newaxis, :, :] * ctx.int_time
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        col_mean_flat = np.nanmean(safe_flat, axis=0)
    col_mean_flat = np.where(np.isfinite(col_mean_flat) & (col_mean_flat > 0.0), col_mean_flat, 1.0)
    smear_rate_denom = col_mean_flat[np.newaxis, :] * ctx.int_time

    pq_out = np.zeros((ctx.ny, ctx.nx), dtype=np.uint16)
    pq_out[:, ctx.strap_col_mask] |= np.uint16(PixelQuality.STRAP_COL)
    pq_out[ctx.source_mask_2d] |= np.uint16(PixelQuality.SOURCE_MASK)

    if (eff_smear == "standard" and not defer_auto_smear) or all_outputs_have_alt_cache:
        m_hat_stamp = alt_sm_stamp if all_outputs_have_alt_cache else collat.collat_smear_stamp
        m_err_stamp = alt_err_stamp if all_outputs_have_alt_cache else collat.collat_smear_err_stamp
        cal_smear_err, raw_sat_2d = _apply_fixed_smear_chunks(
            ctx,
            collat.out_params,
            combined_matched,
            m_hat_stamp,
            m_err_stamp,
            denom_rate,
            smear_rate_denom,
            chunk_rows,
            good_cad=good_cad,
        )
        assert raw_sat_2d is not None
        pq_out[_build_sat_bleed_mask_from_raw_sat(raw_sat_2d, ctx.col_outputs)] |= np.uint16(
            PixelQuality.SAT_BLEED
        )
        if all_outputs_have_alt_cache:
            pq_out[:, alt_unres_cols] |= np.uint16(PixelQuality.SMEAR_UNRESOLVED)
        _finalize_calibrated_hdulist(
            ctx,
            cal_smear_err=cal_smear_err,
            pq_out=pq_out,
            eff_smear=eff_smear,
            cadence_match=cadence_match,
            n_missing=n_missing,
            target_obj=target_obj,
        )
        return

    raw_sat_2d, sat_bleed_2d, col_has_bloom, single_chunk_cache = _scan_pass1_bloom(
        ctx, collat.out_params, combined_matched, good_cad, collat.collat_smear_stamp, chunk_rows
    )
    pq_out[sat_bleed_2d] |= np.uint16(PixelQuality.SAT_BLEED)

    # Deferred smear="automatic": NEVER trigger on ordinary science raw_sat_2d!
    if defer_auto_smear and (
        (ctx.y1_ccd0 + ctx.ny >= 2048)
        or bool(np.any(collat_sm_sat_cols))
        or bool(np.any(col_has_bloom))
    ):
        eff_smear = "alternate"
        hit, alt_sm_stamp, alt_err_stamp, alt_unres_cols, alt_matched = _try_load_alt_caches(
            ctx, ddir=ddir, cadence_match=cadence_match
        )
        if hit:
            all_outputs_have_alt_cache = True
            new_combined = combined_matched & alt_matched
            if not np.array_equal(new_combined, combined_matched):
                combined_matched = new_combined
                n_missing = _enforce_missing_cadence_policy(
                    ctx.tpf_cad, combined_matched, missing_cadence=missing_cadence
                )

    if eff_smear == "standard" or all_outputs_have_alt_cache:
        m_hat_stamp = alt_sm_stamp if all_outputs_have_alt_cache else collat.collat_smear_stamp
        m_err_stamp = alt_err_stamp if all_outputs_have_alt_cache else collat.collat_smear_err_stamp
        if all_outputs_have_alt_cache:
            pq_out[:, alt_unres_cols] |= np.uint16(PixelQuality.SMEAR_UNRESOLVED)
        cal_smear_err, _ = _apply_fixed_smear_chunks(
            ctx,
            collat.out_params,
            combined_matched,
            m_hat_stamp,
            m_err_stamp,
            denom_rate,
            smear_rate_denom,
            chunk_rows,
            single_chunk_cache=single_chunk_cache,
        )
    else:
        m_hat_stamp, m_err_stamp, cal_smear_err, unresolved_cols = _run_cold_alternate_pass(
            ctx,
            out_params=collat.out_params,
            combined_matched=combined_matched,
            good_cad=good_cad,
            collat_smear_stamp=collat.collat_smear_stamp,
            collat_smear_err_stamp=collat.collat_smear_err_stamp,
            collat_sm_sat_cols=collat_sm_sat_cols,
            off_target_smears=collat.off_target_smears,
            off_target_errs=collat.off_target_errs,
            max_e_per_frame_list=collat.max_e_per_frame_list,
            raw_sat_2d=raw_sat_2d,
            sat_bleed_2d=sat_bleed_2d,
            col_has_bloom=col_has_bloom,
            single_chunk_cache=single_chunk_cache,
            denom_rate=denom_rate,
            smear_rate_denom=smear_rate_denom,
            chunk_rows=chunk_rows,
        )
        pq_out[:, unresolved_cols] |= np.uint16(PixelQuality.SMEAR_UNRESOLVED)
        _persist_alt_caches(
            ctx,
            per_out_collat_512=collat.per_out_collat_512,
            per_out_collat_err_512=collat.per_out_collat_err_512,
            m_hat_stamp=m_hat_stamp,
            m_err_stamp=m_err_stamp,
            unresolved_cols=unresolved_cols,
            combined_matched=combined_matched,
            cadence_match=cadence_match,
        )

    _finalize_calibrated_hdulist(
        ctx,
        cal_smear_err=cal_smear_err,
        pq_out=pq_out,
        eff_smear=eff_smear,
        cadence_match=cadence_match,
        n_missing=n_missing,
        target_obj=target_obj,
    )


def _clone_input_tpf(
    tpf: str | os.PathLike[str] | lk.TargetPixelFile,
) -> tuple[fits.HDUList, Path | None, Any]:
    if isinstance(tpf, str | os.PathLike):
        src_path = Path(tpf).expanduser()
        with fits.open(src_path, mode="readonly", memmap=False) as src_hdul:
            return fits.HDUList([h.copy() for h in src_hdul]), src_path, "default"
    if hasattr(tpf, "hdu") and isinstance(tpf.hdu, fits.HDUList):
        hdul = fits.HDUList([h.copy() for h in tpf.hdu])
        raw_path = getattr(tpf, "path", None)
        src_path = (
            Path(raw_path).expanduser()
            if isinstance(raw_path, str | os.PathLike) and os.fspath(raw_path)
            else None
        )
        qmask = getattr(tpf, "quality_bitmask", _MISSING)
        return hdul, src_path, ("default" if qmask is _MISSING else qmask)
    raise TypeError(f"Expected file path or lightkurve.TargetPixelFile, got {type(tpf)!r}.")


def recalibrate_tpf(
    tpf: str | os.PathLike[str] | lk.TargetPixelFile,
    *,
    config: CalibrationConfig | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
    ddir: str | os.PathLike[str] | None = None,
    smear: SmearMode | None = None,
    cadence_match: CadenceMatchMode | None = None,
    missing_cadence: MissingCadenceMode | None = None,
    offline: bool | None = None,
    strict_hpc: bool | None = None,
    chunk_rows: int | None = None,
    verbose: bool | None = None,
) -> lk.TargetPixelFile:
    """Recalibrate a TESS Target Pixel File from ``RAW_CNTS`` under ``HALOCVER = 1``."""
    cfg = config if config is not None else CalibrationConfig()
    hdul, src_path, quality_bitmask = _clone_input_tpf(tpf)
    _recalibrate_hdulist_inplace(
        hdul,
        source_path=src_path,
        cache_dir=cache_dir if cache_dir is not None else cfg.cache_dir,
        ddir=ddir if ddir is not None else cfg.ddir,
        smear=smear if smear is not None else cfg.smear,
        cadence_match=cadence_match if cadence_match is not None else cfg.cadence_match,
        missing_cadence=missing_cadence if missing_cadence is not None else cfg.missing_cadence,
        offline=offline if offline is not None else cfg.offline,
        strict_hpc=strict_hpc if strict_hpc is not None else cfg.strict_hpc,
        chunk_rows=chunk_rows if chunk_rows is not None else cfg.chunk_rows,
        verbose=verbose if verbose is not None else cfg.verbose,
        quality_bitmask=quality_bitmask,
    )
    out_tpf = lk.TessTargetPixelFile(hdul, quality_bitmask=quality_bitmask)
    if src_path is not None and hasattr(out_tpf, "path"):
        out_tpf.path = str(src_path)
    return out_tpf


def desmear_calibrated_cutout(
    tpf: str | os.PathLike[str] | lk.TargetPixelFile,
    *,
    saturation_flux_threshold: float = 1.0e5,
    verbose: bool = False,
) -> lk.TargetPixelFile:
    """Repair collateral-bleed smear corruption on calibrated FFI cutouts (e.g. ``TessCut``)."""
    hdul, src_path, quality_bitmask = _clone_input_tpf(tpf)
    hdr0, hdu1, hdr2 = hdul[0].header, hdul[1], hdul[2].header
    hdr1 = hdu1.header

    names = hdu1.data.names or ()
    if "FLUX" not in names:
        raise ValueError("Cutout HDU[1] lacks a 'FLUX' column.")

    flux_in = np.asarray(hdu1.data["FLUX"], dtype=np.float64)
    n_tpf, ny, nx = flux_in.shape
    quality_mask = _compute_quality_mask(hdu1, quality_bitmask=quality_bitmask)

    bkgsub_hdr = hdr0.get("BKGSUB", hdr1.get("BKGSUB", None))
    backapp_hdr = bool(hdr1.get("BACKAPP", hdr0.get("BACKAPP", False)))
    has_bkg_col = "FLUX_BKG" in names
    if has_bkg_col:
        bkg_arr = np.asarray(hdu1.data["FLUX_BKG"], dtype=np.float64)
        bkg_has_finite = bool(np.any(np.isfinite(bkg_arr)))
        bkg_nonzero = bool(np.any(np.isfinite(bkg_arr) & (bkg_arr != 0.0)))
    else:
        bkg_arr = np.zeros_like(flux_in)
        bkg_has_finite, bkg_nonzero = False, False

    bkg_unrestorable = (not has_bkg_col) or (has_bkg_col and not bkg_has_finite)
    should_restore_bkg = bkg_nonzero and (
        (bkgsub_hdr is True) or (bkgsub_hdr is None and backapp_hdr)
    )
    preserve_bkgsub = (
        not should_restore_bkg
        and bkg_unrestorable
        and ((bkgsub_hdr is True) or (bkgsub_hdr is None and backapp_hdr))
    )
    if should_restore_bkg:
        u_in = flux_in + np.where(np.isfinite(bkg_arr), bkg_arr, 0.0)
    else:
        u_in = flux_in.copy()

    if "FLUX_ERR" in names:
        err_in = np.asarray(hdu1.data["FLUX_ERR"], dtype=np.float64)
        if should_restore_bkg and "FLUX_BKG_ERR" in names:
            berr = np.where(
                np.isfinite(hdu1.data["FLUX_BKG_ERR"]),
                np.asarray(hdu1.data["FLUX_BKG_ERR"], dtype=np.float64),
                0.0,
            )
            sigma_in = np.sqrt(np.maximum(err_in * err_in - berr * berr, 1e-6))
        else:
            sigma_in = np.where(np.isfinite(err_in) & (err_in > 0.0), err_in, 1.0)
    else:
        sigma_in = np.sqrt(
            np.maximum(np.abs(u_in), 1.0) / _effective_exposure_seconds(hdr0, hdr1)
        )

    col0_1based, row0_1based = _stamp_ccd_origin_1based(hdr1, hdr2, default=(45, 1))
    _sync_stamp_ccd_origin(hdr1, hdr2, col0_1based, row0_1based)
    cols_1based = np.arange(col0_1based, col0_1based + nx, dtype=np.int64)
    if np.all((cols_1based >= 45) & (cols_1based <= 2092)):
        col_outputs = tuple(output_for_rawx_col(int(c)) for c in cols_1based)
        strap_col_mask = strap_mask_for_rawx_columns(cols_1based, buffer_px=2)
    else:
        col_outputs = tuple("A" for _ in range(nx))
        strap_col_mask = np.zeros(nx, dtype=bool)

    existing_pq = (
        np.asarray(hdul[PIXEL_QUALITY_EXTNAME].data, dtype=np.uint16)
        if PIXEL_QUALITY_EXTNAME in hdul
        else np.zeros((ny, nx), dtype=np.uint16)
    )
    source_mask_2d = (existing_pq & int(PixelQuality.SOURCE_MASK)) != 0

    good_cad = quality_mask if np.any(quality_mask) else np.ones(n_tpf, dtype=bool)
    sci_sat_3d = np.isfinite(u_in) & (u_in > float(saturation_flux_threshold))
    raw_sat_2d = np.any(sci_sat_3d[good_cad], axis=0) | (
        (existing_pq & int(PixelQuality.SAT_BLEED)) != 0
    )
    sat_bleed_2d = _build_sat_bleed_mask_from_raw_sat(raw_sat_2d, col_outputs) | raw_sat_2d

    if verbose:
        print(f"Desmearing calibrated cutout of shape {(n_tpf, ny, nx)}")

    delta_m, delta_m_err, unresolved_cols, zero_floor_clamped = _estimate_alternate_smear_stamp(
        sci_pre_smear=u_in,
        sci_sigma_pre=sigma_in,
        flat_2d=np.ones((ny, nx), dtype=np.float64),
        collat_smear_stamp=None,
        collat_smear_err_stamp=None,
        collat_sm_sat_cols=None,
        off_target_smear=None,
        off_target_smear_err=None,
        raw_sat_2d=raw_sat_2d,
        sat_bleed_2d=sat_bleed_2d,
        strap_col_mask=strap_col_mask,
        source_mask_2d=source_mask_2d,
        config=AlternateSmearConfig(
            row0_1based=row0_1based,
            row_smear_weight=_frame_transfer_row_weight(hdr0, hdr1),
            saturation_limit=float(saturation_flux_threshold),
        ),
        bloom_detect_cube=flux_in if should_restore_bkg else None,
        good_cad_mask=good_cad,
    )

    effective_bkgsub = preserve_bkgsub
    if zero_floor_clamped and not preserve_bkgsub:
        warnings.warn(
            "All cutout columns are bloomed with no collateral data or FLUX_BKG reference; "
            "a uniform sky pedestal is degenerate with uniform smear bloom, so the output "
            "floor has been zero-referenced (stamping BKGSUB=True).",
            UserWarning,
            stacklevel=2,
        )
        effective_bkgsub = True

    pq_out = np.zeros((ny, nx), dtype=np.uint16)
    pq_out[:, strap_col_mask] |= np.uint16(PixelQuality.STRAP_COL)
    pq_out[sat_bleed_2d] |= np.uint16(PixelQuality.SAT_BLEED)
    pq_out[:, unresolved_cols] |= np.uint16(PixelQuality.SMEAR_UNRESOLVED)
    pq_out[source_mask_2d] |= np.uint16(PixelQuality.SOURCE_MASK)

    flux_out = u_in - delta_m[:, np.newaxis, :]
    flux_err_out = np.hypot(sigma_in, delta_m_err[:, np.newaxis, :])
    finite_cad = np.any(np.isfinite(u_in), axis=(1, 2))
    if not np.all(finite_cad):
        flux_out[~finite_cad], flux_err_out[~finite_cad], delta_m_err[~finite_cad] = (
            np.nan,
            np.nan,
            np.nan,
        )

    hdul[1].data["FLUX"] = flux_out.astype(np.float32)
    if "FLUX_ERR" in names:
        hdul[1].data["FLUX_ERR"] = flux_err_out.astype(np.float32)
    if not effective_bkgsub:
        if "FLUX_BKG" in names:
            hdul[1].data["FLUX_BKG"] = np.float32(0.0)
        if "FLUX_BKG_ERR" in names:
            hdul[1].data["FLUX_BKG_ERR"] = np.float32(0.0)

    append_smear_err_column(hdul, delta_m_err.astype(np.float32))
    write_pixel_quality(hdul, pq_out)
    write_mundey_provenance(hdul, bkgsub=effective_bkgsub, smear_mode="cutout")
    hdul[0].header["CADMATCH"] = ("none", "Collateral cadence pairing mode")
    hdul[0].header["CADMISS"] = (0, "Unmatched TPF cadences masked to NaN")

    out_tpf = lk.TessTargetPixelFile(hdul, quality_bitmask=quality_bitmask)
    if src_path is not None and hasattr(out_tpf, "path"):
        out_tpf.path = str(src_path)
    return out_tpf
