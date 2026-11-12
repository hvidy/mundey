"""JAX ``@in_float64`` TESS bright-star calibration engine (``mundey`` v2.0).

Implements physical TESS readout calibration in locally scoped 64-bit JAX:
1. Frame-scaled 2D black & output mean-black subtraction with sentinel masking.
2. Cancellation-free discrete orthonormal Legendre 1D-black fitting on clean
   trailing virtual columns (``cols 3..9``) with AIC degree selection.
3. Causal serial LDE undershoot inversion via ``jax.lax.associative_scan`` in
   ADU *before* per-frame ``[0, maxDomain]``-clamped nonlinearity and gain.
4. Exact ``CADENCENO`` collateral pairing across unified 2048-column multi-output
   (``'A'``, ``'B'``, ``'C'``, ``'D'``) CCD cutouts.
5. Concavity-guarded collateral-anchored ``smear='alternate'`` and FFI cutout
   desmearing (:func:`desmear_calibrated_cutout`) under ``HALOCVER = 1``.
"""

from __future__ import annotations

import contextlib
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

from ._contract import (
    PIXEL_QUALITY_EXTNAME,
    PixelQuality,
    append_smear_err_column,
    resolve_cache_dir,
    write_mundey_provenance,
    write_pixel_quality,
)
from ._precision import in_float64
from .archive import (
    DEFAULT_CHUNK_ROWS,
    CadenceKind,
    StaticModelPaths,
    collateral_blf_cache_path,
    find_cached_blf,
    infer_cadence_kind,
    iter_collateral_chunks,
    read_collateral_headers_and_index,
    resolve_collateral_file,
    resolve_static_models,
    write_onedblack_blf,
    write_smear_blf,
)

CadenceMatchMode = Literal["cadenceno", "positional"]
MissingCadenceMode = Literal["nan", "error"]
SmearMode = Literal["automatic", "standard", "alternate"]

# Output geometry on full 2136-column TESS CCD frame (0-based slice bounds):
# Science imaging columns occupy 0-based [44, 2092) -> 1-based RAWX 45..2092.
OUTPUT_ORDER: tuple[str, str, str, str] = ("A", "B", "C", "D")
OUTPUT_BOUNDS_0BASED: dict[str, tuple[int, int, int, int]] = {
    "A": (44, 556, 2092, 2103),
    "B": (556, 1068, 2103, 2114),
    "C": (1068, 1580, 2114, 2125),
    "D": (1580, 2092, 2125, 2136),
}
OUTPUT_SCI0_OFFSET: dict[str, int] = {
    "A": 0,
    "B": 512,
    "C": 1024,
    "D": 1536,
}
OUTPUT_REVERSE_READOUT: dict[str, bool] = {
    "A": False,
    "B": True,
    "C": False,
    "D": True,
}


def _build_ccid80_strap_cols_1based_imaging() -> tuple[int, ...]:
    """Return the 116 TESS CCID-80 routing-strap columns (1-based imaging 1..2048).

    Constructed from TESS Instrument Handbook Table 4.2 / Table A.3 with exact
    bilateral CCID-80 symmetry ``c <-> 2049 - c``. Add 44 for 1-based FITS
    ``RAWX`` coordinates (``45..2092``).
    """
    left: list[int] = []
    for base in (77, 205, 333, 461, 589, 717, 845, 973):
        left.extend((base, base + 1, base + 4, base + 5, base + 8, base + 9))
    left.extend((516, 517, 520, 521, 524, 525))
    left.extend((673, 674, 675, 676))
    left_sorted = sorted(left)
    right_sorted = sorted(2049 - c for c in left_sorted)
    return tuple(left_sorted + right_sorted)


CCID80_STRAP_COLS_IMAGING_1BASED: tuple[int, ...] = _build_ccid80_strap_cols_1based_imaging()
CCID80_STRAP_COLS_RAWX_1BASED: np.ndarray = (
    np.asarray(CCID80_STRAP_COLS_IMAGING_1BASED, dtype=np.int64) + 44
)


def strap_mask_for_rawx_columns(
    cols_1based: np.ndarray,
    *,
    buffer_px: int = 2,
) -> np.ndarray:
    """Return boolean strap mask for 1-based FITS ``RAWX`` columns with ``±buffer_px``."""
    cols = np.asarray(cols_1based, dtype=np.int64)[:, np.newaxis]
    dist = np.min(np.abs(cols - CCID80_STRAP_COLS_RAWX_1BASED[np.newaxis, :]), axis=1)
    return np.asarray(dist <= int(buffer_px), dtype=bool)


def _build_orthonormal_legendre_basis(n_rows: int, max_degree: int) -> np.ndarray:
    """Construct discrete orthonormal Legendre basis ``Q`` of shape ``(n_rows, max_degree + 1)``."""
    x = np.linspace(-1.0, 1.0, n_rows, dtype=np.float64)
    v = np.polynomial.legendre.legvander(x, max_degree)
    q, r = np.linalg.qr(v, mode="reduced")
    signs = np.where(np.diag(r) >= 0.0, 1.0, -1.0)
    return np.asarray(q * signs[np.newaxis, :], dtype=np.float64)


# Precomputed pure-NumPy float64 bases (zero import-time JAX x64 state change):
_LEGENDRE_Q_SCI_NP: np.ndarray = _build_orthonormal_legendre_basis(2048, 49)  # (2048, 50)
_LEGENDRE_Q_OV_NP: np.ndarray = _build_orthonormal_legendre_basis(30, 2)  # (30, 3)


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


def load_linearity_channel(
    xml_path: str | os.PathLike[str],
    camera: int,
    ccd: int,
    output: str,
) -> LinearityChannel:
    """Parse ``linearityPoly`` parameters (including ``maxDomain``) from ``linearity.xml``."""
    tree = ET.parse(Path(xml_path))
    root = tree.getroot()
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


def effective_frame_count(
    hdr0: fits.Header | None,
    hdr1: fits.Header | None,
    *,
    fallback_hdr0: fits.Header | None = None,
    fallback_hdr1: fits.Header | None = None,
) -> float:
    """Compute effective co-added frame count ``N_eff = factor * NUM_FRM`` from FITS headers.

    Supports both 2-minute cosmic-ray-mitigated data (``CRMITEN=True``, trimming 2
    extrema per ``CRBLKSZ``-frame block) and 20-second fast data (``CRMITEN=False``,
    ``factor = 1.0``). Reads independently from science or collateral headers, using
    ``fallback_hdr*`` only when minimal synthetic collateral tables omit keywords.
    """

    def _lookup(key: str) -> Any:
        for hdr in (hdr0, hdr1, fallback_hdr0, fallback_hdr1):
            if hdr is not None and key in hdr:
                return hdr[key]
        return None

    crmiten_raw = _lookup("CRMITEN")
    crmiten = bool(crmiten_raw) if crmiten_raw is not None else False
    if crmiten:
        crblksz_raw = _lookup("CRBLKSZ")
        crblksz = float(crblksz_raw) if crblksz_raw is not None else 10.0
        if crblksz <= 2.0:
            raise ValueError(f"Invalid CRBLKSZ={crblksz} with CRMITEN=True.")
        factor = (crblksz - 2.0) / crblksz
    else:
        factor = 1.0

    num_frm_raw = _lookup("NUM_FRM")
    if num_frm_raw is None:
        raise KeyError("Missing NUM_FRM keyword in FITS headers.")
    num_frm = float(num_frm_raw)
    if num_frm <= 0.0:
        raise ValueError(f"NUM_FRM must be positive, got {num_frm}.")
    return factor * num_frm


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


def match_cadenceno(
    tpf_cadenceno: np.ndarray,
    col_cadenceno: np.ndarray,
    *,
    cadence_match: CadenceMatchMode = "cadenceno",
    label: str = "collateral",
) -> tuple[np.ndarray, np.ndarray]:
    """Align TPF cadences to a collateral or ``-blf.fits`` table.

    Parameters
    ----------
    tpf_cadenceno : np.ndarray
        1D integer ``CADENCENO`` array from ``tpf.hdu[1].data["CADENCENO"]``.
    col_cadenceno : np.ndarray
        1D integer ``CADENCENO`` array from collateral or ``-blf.fits`` table.
    cadence_match : {"cadenceno", "positional"}, default "cadenceno"
        Exact integer ``CADENCENO`` lookup (``"cadenceno"``) or legacy row-index
        pairing (``"positional"``).
    label : str
        Human-readable table description for error messages.

    Returns
    -------
    matched : np.ndarray
        Boolean mask ``(N_tpf,)`` of matched TPF rows.
    col_idx : np.ndarray
        ``int64`` index ``(N_tpf,)`` into ``col_cadenceno`` (``-1`` where unmatched).
    """
    tpf_cad = np.asarray(tpf_cadenceno, dtype=np.int64)
    col_cad = np.asarray(col_cadenceno, dtype=np.int64)
    n_tpf = len(tpf_cad)
    n_col = len(col_cad)

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

    col_idx = np.where(matched, order[safe_pos], -1).astype(np.int64)
    return matched, col_idx


def _enforce_missing_cadence_policy(
    tpf_cadenceno: np.ndarray,
    combined_matched: np.ndarray,
    *,
    missing_cadence: MissingCadenceMode,
) -> int:
    """Validate or warn on unmatched TPF cadences and return ``n_missing``."""
    if missing_cadence not in ("nan", "error"):
        raise ValueError(
            f"Invalid missing_cadence={missing_cadence!r}; expected 'nan' or 'error'."
        )

    n_tpf = len(tpf_cadenceno)
    n_matched = int(np.sum(combined_matched))
    if n_matched == 0:
        raise ValueError(
            "Zero TPF cadences have complete collateral CADENCENO coverage across outputs."
        )

    n_missing = n_tpf - n_matched
    if n_missing == 0:
        return 0

    missing_cads = np.asarray(tpf_cadenceno, dtype=np.int64)[~combined_matched]
    c_min = int(missing_cads.min())
    c_max = int(missing_cads.max())
    msg = (
        f"{n_missing} / {n_tpf} TPF cadences lack matching collateral CADENCENO "
        f"(missing CADENCENO range {c_min}..{c_max})"
    )
    if missing_cadence == "error":
        raise ValueError(f"{msg}.")

    warnings.warn(
        f"{msg}; masking unmatched cadences to NaN.",
        UserWarning,
        stacklevel=3,
    )
    return n_missing


@jax.jit
def _fit_onedblack_legendre_jit(
    tv_debiased_adu: jax.Array,
    q_sci: jax.Array,
    q_ov: jax.Array,
) -> tuple[jax.Array, jax.Array]:
    """JIT-compiled 64-bit orthonormal Legendre 1D-black solver for ``(N_c, 2078, 11)``."""
    # Exclude leading settling virtual columns 0..2 and trailing column 10 -> keep 3..9
    v_clean = tv_debiased_adu[:, :, 3:10]  # (N_c, 2078, 7)

    # Row-wise 5-sigma cosmic-ray rejection across the 7 clean virtual columns
    row_med = jnp.nanmedian(v_clean, axis=2, keepdims=True)
    abs_dev = jnp.abs(v_clean - row_med)
    row_mad = 1.4826 * jnp.nanmedian(abs_dev, axis=2, keepdims=True)
    global_mad = 1.4826 * jnp.nanmedian(abs_dev, axis=(1, 2), keepdims=True)
    sigma_floor = jnp.maximum(global_mad, 1e-6)
    sigma_row = jnp.maximum(row_mad, sigma_floor)

    valid_pix = jnp.isfinite(v_clean) & (abs_dev <= 5.0 * sigma_row)
    # Fallback to all finite pixels in a row if every column was rejected
    any_valid = jnp.any(valid_pix, axis=2, keepdims=True)
    use_pix = jnp.where(any_valid, valid_pix, jnp.isfinite(v_clean))
    count_pix = jnp.sum(use_pix, axis=2)
    y_all = jnp.sum(jnp.where(use_pix, v_clean, 0.0), axis=2) / jnp.maximum(count_pix, 1)
    y_all = jnp.where(count_pix > 0, y_all, jnp.nan)  # (N_c, 2078)

    y_sci_raw = y_all[:, :2048]
    cad_valid = jnp.sum(jnp.isfinite(y_sci_raw), axis=1) >= 64  # (N_c,)
    cad_med_sci = jnp.where(cad_valid, jnp.nanmedian(y_sci_raw, axis=1), 0.0)
    y_sci = jnp.where(jnp.isfinite(y_sci_raw), y_sci_raw, cad_med_sci[:, None])

    # Orthonormal projection C = Y @ Q, Y_perp = Y - C @ Q.T
    coeffs = y_sci @ q_sci  # (N_c, 50)
    y_perp = y_sci - coeffs @ q_sci.T
    norm_perp_sq = jnp.sum(y_perp * y_perp, axis=1)  # (N_c,)

    # Strictly non-negative orthogonal-complement RSS_{t, k} for k in 20..49
    # tail_sq[t, k] = sum_{j = k+1}^{49} C_{t, j}^2
    c_sq = coeffs * coeffs
    suffix_excl = jnp.cumsum(c_sq[:, ::-1], axis=1)[:, ::-1]
    tail_all = jnp.concatenate(
        [suffix_excl[:, 1:], jnp.zeros((coeffs.shape[0], 1), dtype=coeffs.dtype)],
        axis=1,
    )  # (N_c, 50) for k = 0..49
    rss_20_49 = norm_perp_sq[:, None] + tail_all[:, 20:50]  # (N_c, 30)

    ks = jnp.arange(20, 50, dtype=coeffs.dtype)[None, :]  # degrees 20..49
    n_rows_f = jnp.asarray(2048.0, dtype=coeffs.dtype)
    aics = 2.0 * (ks + 1.0) + n_rows_f * jnp.log(jnp.maximum(rss_20_49, 1e-30))
    best_rel = jnp.argmin(aics, axis=1)  # (N_c,) in 0..29
    best_k = best_rel + 20  # (N_c,) in 20..49

    deg_idx = jnp.arange(50)[None, :]
    keep_deg = deg_idx <= best_k[:, None]
    coeffs_sel = jnp.where(keep_deg, coeffs, 0.0)
    b1d_sci = coeffs_sel @ q_sci.T  # (N_c, 2048)

    best_rss = jnp.take_along_axis(rss_20_49, best_rel[:, None], axis=1)[:, 0]
    dof_sci = n_rows_f - (best_k.astype(coeffs.dtype) + 1.0)
    sigma2_sci = best_rss / jnp.maximum(dof_sci, 1.0)
    # Row leverage diag(Q_{:, :k+1} Q_{:, :k+1}^T) = sum_{j <= k} Q_{r, j}^2
    q_sci_sq_cum = jnp.cumsum(q_sci * q_sci, axis=1)  # (2048, 50)
    lev_sci = jnp.take(q_sci_sq_cum, best_k, axis=1).T  # (N_c, 2048)
    err_sci = jnp.sqrt(jnp.maximum(sigma2_sci[:, None] * lev_sci, 0.0))

    # Collateral rows 2048..2077 (30 rows): separate degree-2 orthonormal Legendre fit
    y_ov_raw = y_all[:, 2048:2078]
    ov_valid = jnp.sum(jnp.isfinite(y_ov_raw), axis=1) >= 6
    cad_med_ov = jnp.where(ov_valid, jnp.nanmedian(y_ov_raw, axis=1), cad_med_sci)
    y_ov = jnp.where(jnp.isfinite(y_ov_raw), y_ov_raw, cad_med_ov[:, None])

    c_ov = y_ov @ q_ov  # (N_c, 3)
    b1d_ov = c_ov @ q_ov.T  # (N_c, 30)
    res_ov = y_ov - b1d_ov
    sigma2_ov = jnp.sum(res_ov * res_ov, axis=1) / 27.0
    lev_ov = jnp.sum(q_ov * q_ov, axis=1)[None, :]  # (1, 30)
    err_ov = jnp.sqrt(jnp.maximum(sigma2_ov[:, None] * lev_ov, 0.0))

    b1d_full = jnp.concatenate([b1d_sci, b1d_ov], axis=1)
    err_full = jnp.concatenate([err_sci, err_ov], axis=1)
    b1d_out = jnp.where(cad_valid[:, None], b1d_full, jnp.nan)
    err_out = jnp.where(cad_valid[:, None], err_full, jnp.nan)
    return b1d_out, err_out


_ONEDBLACK_CAD_BUCKETS: tuple[int, ...] = (32, 256, 4096)
_PRESMEAR_WIDTH_BUCKETS: tuple[int, ...] = (16, 64, 512, 1056, 2048)
_HUBER_ROW_BUCKETS: tuple[int, ...] = (64, 256, 2048)


def _select_bucket(size: int, buckets: tuple[int, ...]) -> int:
    for b in buckets:
        if size <= b:
            return b
    return size


@in_float64
def _fit_onedblack_static_chunk(tv_padded_adu: Any) -> tuple[np.ndarray, np.ndarray]:
    q_sci = jnp.asarray(_LEGENDRE_Q_SCI_NP, dtype=jnp.float64)
    q_ov = jnp.asarray(_LEGENDRE_Q_OV_NP, dtype=jnp.float64)
    b1d, err = _fit_onedblack_legendre_jit(tv_padded_adu, q_sci, q_ov)
    return np.asarray(b1d, dtype=np.float64), np.asarray(err, dtype=np.float64)


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
            buf[:m] = arr[t0:t1]
            buf[m:] = np.nan
            b_c, e_c = _fit_onedblack_static_chunk(buf)
        b1d_out[t0:t1] = b_c[:m]
        err_out[t0:t1] = e_c[:m]
    return b1d_out, err_out


def _undershoot_combine(
    prev: tuple[jax.Array, jax.Array],
    curr: tuple[jax.Array, jax.Array],
) -> tuple[jax.Array, jax.Array]:
    """Associative affine composition ``q(p(x)) = (q_a * p_a) x + (q_a * p_b + q_b)``."""
    p_a, p_b = prev
    q_a, q_b = curr
    return q_a * p_a, q_a * p_b + q_b


def _invert_lde_undershoot_forward_impl(
    w_adu: jax.Array,
    u_coeff: jax.Array,
) -> jax.Array:
    """Invert causal serial LDE undershoot ``w_m = x_m - u * x_{m-1}`` along ``axis=-1``."""
    if w_adu.shape[-1] == 0:
        return w_adu

    finite = jnp.isfinite(w_adu)
    w_clean = jnp.where(finite, w_adu, 0.0)

    # Steady-state boundary x_0 = w_0 / (1 - u); causal recursion x_m = w_m + u * x_{m-1}
    b0 = jnp.where(finite[..., :1], w_clean[..., :1] / (1.0 - u_coeff), 0.0)
    a0 = jnp.zeros_like(b0)

    if w_adu.shape[-1] > 1:
        # Non-finite interior/padded pixels pass previous state through cleanly
        a_rest = jnp.where(finite[..., 1:], u_coeff, 1.0)
        b_rest = jnp.where(finite[..., 1:], w_clean[..., 1:], 0.0)
        a = jnp.concatenate([a0, a_rest], axis=-1)
        b = jnp.concatenate([b0, b_rest], axis=-1)
        _, x_scan = jax.lax.associative_scan(_undershoot_combine, (a, b), axis=-1)
    else:
        x_scan = b0

    return jnp.where(finite, x_scan, jnp.nan)


_invert_lde_undershoot_forward_jit = jax.jit(_invert_lde_undershoot_forward_impl)


@in_float64
def _invert_lde_undershoot_static(w_2d: Any, u_coeff: float) -> np.ndarray:
    res = _invert_lde_undershoot_forward_jit(w_2d, jnp.asarray(u_coeff, dtype=jnp.float64))
    return np.asarray(res, dtype=np.float64)


def invert_lde_undershoot_adu(
    y_adu: Any,
    u_coeff: float,
    *,
    reverse_readout: bool = False,
) -> np.ndarray:
    """Invert causal LDE undershoot in ADU along physical readout direction."""
    arr = np.asarray(y_adu, dtype=np.float64)
    nx = arr.shape[-1]
    if nx == 0 or arr.size == 0:
        return arr.copy()
    w = arr[..., ::-1] if reverse_readout else arr
    w_2d = w.reshape(-1, nx)
    b_tot = w_2d.shape[0]
    w_bucket = _select_bucket(nx, _PRESMEAR_WIDTH_BUCKETS)
    b_step = 512 if b_tot <= 512 else 8192

    out_2d = np.empty((b_tot, nx), dtype=np.float64)
    buf = np.full((b_step, w_bucket), np.nan, dtype=np.float64)
    for b0 in range(0, b_tot, b_step):
        b1 = min(b0 + b_step, b_tot)
        m = b1 - b0
        buf[:m, :nx] = w_2d[b0:b1]
        if m < b_step:
            buf[m:, :nx] = np.nan
        res_chunk = _invert_lde_undershoot_static(buf, u_coeff)
        out_2d[b0:b1] = res_chunk[:m, :nx]

    out = out_2d.reshape(arr.shape)
    return out[..., ::-1].copy() if reverse_readout else out


def _apply_clamped_linearity_impl(
    dn_per_frame: jax.Array,
    originx: jax.Array,
    scalex: jax.Array,
    offsetx: jax.Array,
    max_domain: jax.Array,
    coeffs_desc: jax.Array,
    deriv_coeffs_desc: jax.Array,
    gain: jax.Array,
) -> tuple[jax.Array, jax.Array, jax.Array]:
    """Evaluate ``[0, maxDomain]``-clamped nonlinearity polynomial, saturation mask, and slope."""
    sat_mask = jnp.isfinite(dn_per_frame) & (dn_per_frame > max_domain)
    d_clamped = jnp.clip(dn_per_frame, 0.0, max_domain)
    s = offsetx + (d_clamped - originx) * scalex
    poly_val = jnp.polyval(coeffs_desc, s)
    e_per_frame = dn_per_frame * poly_val * gain

    in_domain = (dn_per_frame >= 0.0) & (dn_per_frame <= max_domain)
    dpoly_ds = jnp.polyval(deriv_coeffs_desc, s)
    local_gain = (
        poly_val + jnp.where(in_domain, d_clamped * dpoly_ds * scalex, 0.0)
    ) * gain
    return e_per_frame, sat_mask, local_gain


_apply_clamped_linearity_jit = jax.jit(_apply_clamped_linearity_impl)


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
    x_adu = _invert_lde_undershoot_forward_impl(w_adu_2d, u_coeff)
    dn_per_frame = x_adu / n_eff
    return _apply_clamped_linearity_impl(
        dn_per_frame,
        originx,
        scalex,
        offsetx,
        max_domain,
        coeffs_desc,
        deriv_coeffs_desc,
        gain,
    )


def _prepare_poly_coeffs_desc_np(channel: LinearityChannel) -> tuple[np.ndarray, np.ndarray]:
    """Return fixed-length descending polynomial & derivative coefficients in ``float64``."""
    raw_asc = np.asarray(channel.coeffs, dtype=np.float64)
    pad_len = max(8, len(raw_asc))
    coeffs_asc = np.zeros(pad_len, dtype=np.float64)
    coeffs_asc[: len(raw_asc)] = raw_asc
    coeffs_desc = coeffs_asc[::-1].copy()
    deriv_asc = coeffs_asc[1:] * np.arange(1, pad_len, dtype=np.float64)
    deriv_desc = deriv_asc[::-1].copy()
    return coeffs_desc, deriv_desc


@in_float64
def apply_clamped_linearity_and_gain(
    dn_per_frame: Any,
    channel: LinearityChannel,
    gain: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert per-frame DN to ``e-/frame`` with ``linearity.xml`` clamped to ``[0, maxDomain]``."""
    coeffs_desc_np, deriv_desc_np = _prepare_poly_coeffs_desc_np(channel)
    e_out, sat_out, g_out = _apply_clamped_linearity_jit(
        dn_per_frame,
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

    w = arr[..., ::-1] if reverse_readout else arr
    w_2d = w.reshape(-1, nx)
    b_tot = w_2d.shape[0]
    w_bucket = _select_bucket(nx, _PRESMEAR_WIDTH_BUCKETS)
    if b_tot <= 512:
        b_step = 512
    elif b_tot <= 8192:
        b_step = 8192
    elif b_tot <= 131072 and w_bucket <= 64:
        b_step = 131072
    else:
        b_step = 8192
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
        e_2d[b0:b1] = ec[:m, :nx]
        sat_2d[b0:b1] = sc[:m, :nx]
        g_2d[b0:b1] = gc[:m, :nx]

    e_full = e_2d.reshape(arr.shape)
    sat_full = sat_2d.reshape(arr.shape)
    g_full = g_2d.reshape(arr.shape)
    if reverse_readout:
        return (
            e_full[..., ::-1].copy(),
            sat_full[..., ::-1].copy(),
            g_full[..., ::-1].copy(),
        )
    return e_full, sat_full, g_full


@jax.jit
def _huber_weighted_smear_location_2d_jit(
    residuals: jax.Array,
    inv_wing_weights: jax.Array,
    sigma_pix: jax.Array,
    clean_mask: jax.Array,
) -> tuple[jax.Array, jax.Array]:
    """2D shape-bucketed MAD-scaled Huber IRLS row-location solver ``(M, H) -> (M,)``."""
    valid = clean_mask & jnp.isfinite(residuals) & jnp.isfinite(inv_wing_weights)
    w_base = jnp.where(valid, jnp.maximum(inv_wing_weights, 0.0), 0.0)

    # Identify far-wing rows (top >= 5% of inverse-wing weight) for initial median & scale
    w_max = jnp.max(w_base, axis=1, keepdims=True)
    far_mask = valid & (w_base >= 0.05 * w_max)
    far_count = jnp.sum(far_mask, axis=1, keepdims=True)
    scale_mask = jnp.where(far_count >= 5, far_mask, valid)

    # Initial weighted location over far clean rows
    w_init = jnp.where(scale_mask, w_base, 0.0)
    w_init_sum = jnp.sum(w_init, axis=1)
    loc = jnp.sum(w_init * jnp.where(valid, residuals, 0.0), axis=1) / jnp.maximum(
        w_init_sum, 1e-30
    )

    sigma_eff = jnp.where(w_base > 0.0, jax.lax.rsqrt(jnp.maximum(w_base, 1e-30)), 1.0)
    c_huber = jnp.asarray(1.345, dtype=residuals.dtype)

    # Fixed preliminary MAD scale tau_0 (canonical Huber 1981 IRLS without inner sorts)
    r0 = jnp.where(valid, residuals - loc[:, None], 0.0)
    std_r0 = jnp.where(scale_mask, jnp.abs(r0) / jnp.maximum(sigma_eff, 1e-12), jnp.nan)
    tau0 = jnp.maximum(1.0, 1.4826 * jnp.nanmedian(std_r0, axis=1))
    inv_denom_scale = 1.0 / jnp.maximum(tau0[:, None] * sigma_eff, 1e-12)
    res_valid = jnp.where(valid, residuals, 0.0)

    def _irls_step(_idx: int, m_curr: jax.Array) -> jax.Array:
        r = jnp.where(valid, residuals - m_curr[:, None], 0.0)
        u = jnp.abs(r) * inv_denom_scale
        w_rob = jnp.where(u <= c_huber, 1.0, c_huber / jnp.maximum(u, 1e-12))
        w_tot = w_base * w_rob
        denom = jnp.sum(w_tot, axis=1)
        m_next = jnp.sum(w_tot * res_valid, axis=1) / jnp.maximum(denom, 1e-30)
        return jnp.where(denom > 0.0, m_next, m_curr)

    loc_final = jax.lax.fori_loop(0, 6, _irls_step, loc)

    # Sandwich variance of Huber location
    r_fin = jnp.where(valid, residuals - loc_final[:, None], 0.0)
    u_fin = jnp.abs(r_fin) * inv_denom_scale
    w_rob_fin = jnp.where(u_fin <= c_huber, 1.0, c_huber / jnp.maximum(u_fin, 1e-12))
    w_tot_fin = w_base * w_rob_fin
    denom_fin = jnp.sum(w_tot_fin, axis=1)

    pix_var = jnp.maximum(r_fin * r_fin, sigma_pix * sigma_pix)
    num_var = jnp.sum((w_tot_fin * w_tot_fin) * jnp.where(valid, pix_var, 0.0), axis=1)
    loc_err = jnp.sqrt(jnp.maximum(num_var, 0.0)) / jnp.maximum(denom_fin, 1e-30)

    has_any = jnp.any(valid, axis=1)
    return (
        jnp.where(has_any, loc_final, jnp.nan),
        jnp.where(has_any, loc_err, jnp.nan),
    )


@in_float64
def _solve_huber_static_chunk(
    res_buf: Any,
    w_buf: Any,
    sig_buf: Any,
    mask_buf: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    loc, err = _huber_weighted_smear_location_2d_jit(
        res_buf,
        w_buf,
        sig_buf,
        jnp.asarray(mask_buf, dtype=bool),
    )
    return np.asarray(loc, dtype=np.float64), np.asarray(err, dtype=np.float64)


def _solve_huber_smear_location(
    residuals: Any,
    inv_wing_weights: Any,
    sigma_pix: Any,
    clean_row_mask: Any,
) -> tuple[np.ndarray, np.ndarray]:
    """Run shape-bucketed ``@in_float64`` MAD-scaled Huber smear solver."""
    res_np = np.asarray(residuals, dtype=np.float64)  # (N_t, ny, n_bleed)
    w_np = np.asarray(inv_wing_weights, dtype=np.float64)
    sig_np = np.asarray(sigma_pix, dtype=np.float64)
    mask_np = np.asarray(clean_row_mask, dtype=bool)  # (ny, n_bleed)

    n_t, ny, n_bleed = res_np.shape
    m_tot = n_t * n_bleed
    if m_tot == 0 or ny == 0:
        return (
            np.full((n_t, n_bleed), np.nan, dtype=np.float64),
            np.full((n_t, n_bleed), np.nan, dtype=np.float64),
        )

    res_2d = np.transpose(res_np, (0, 2, 1)).reshape(m_tot, ny)
    w_2d = np.transpose(w_np, (0, 2, 1)).reshape(m_tot, ny)
    sig_2d = np.transpose(sig_np, (0, 2, 1)).reshape(m_tot, ny)
    mask_2d = np.broadcast_to(mask_np.T[None, :, :], (n_t, n_bleed, ny)).reshape(m_tot, ny)

    h_bucket = _select_bucket(ny, _HUBER_ROW_BUCKETS)
    m_step = 64 if m_tot <= 64 else 1024

    loc_1d = np.empty(m_tot, dtype=np.float64)
    err_1d = np.empty(m_tot, dtype=np.float64)
    res_buf = np.zeros((m_step, h_bucket), dtype=np.float64)
    w_buf = np.zeros((m_step, h_bucket), dtype=np.float64)
    sig_buf = np.ones((m_step, h_bucket), dtype=np.float64)
    mask_buf = np.zeros((m_step, h_bucket), dtype=bool)

    for m0 in range(0, m_tot, m_step):
        m1 = min(m0 + m_step, m_tot)
        k = m1 - m0
        res_buf[:k, :ny] = res_2d[m0:m1]
        w_buf[:k, :ny] = w_2d[m0:m1]
        sig_buf[:k, :ny] = sig_2d[m0:m1]
        mask_buf[:k, :ny] = mask_2d[m0:m1]
        if k < m_step:
            mask_buf[k:, :] = False
        lc, ec = _solve_huber_static_chunk(res_buf, w_buf, sig_buf, mask_buf)
        loc_1d[m0:m1] = lc[:k]
        err_1d[m0:m1] = ec[:k]

    return loc_1d.reshape(n_t, n_bleed), err_1d.reshape(n_t, n_bleed)


def _dilate_vertical_bool(mask_2d: np.ndarray, radius: int) -> np.ndarray:
    """Dilate a 2D boolean mask ``(ny, nx)`` along rows (``axis=0``) by ``±radius``."""
    ny, nx = mask_2d.shape
    out = np.array(mask_2d, dtype=bool, copy=True)
    for shift in range(1, radius + 1):
        if shift < ny:
            out[shift:, :] |= mask_2d[:-shift, :]
            out[:-shift, :] |= mask_2d[shift:, :]
        else:
            has_col = np.any(mask_2d, axis=0, keepdims=True)
            out |= np.broadcast_to(has_col, (ny, nx))
            break
    return out


def _build_sat_bleed_mask_2d(
    sci_sat_3d: np.ndarray,
    col_outputs: tuple[str, ...],
) -> np.ndarray:
    """Build 2D ``PixelQuality.SAT_BLEED`` mask including downstream undershoot neighbours."""
    raw_sat_2d = np.any(sci_sat_3d, axis=0)  # (ny, nx)
    ny, nx = raw_sat_2d.shape
    dilated = _dilate_vertical_bool(raw_sat_2d, radius=1)

    # Flag 1-pixel downstream serial-readout neighbour within the same CCD output
    downstream = np.zeros((ny, nx), dtype=bool)
    for c in range(nx):
        if not np.any(raw_sat_2d[:, c]):
            continue
        out = col_outputs[c]
        c_down = (c - 1) if OUTPUT_REVERSE_READOUT[out] else (c + 1)
        if 0 <= c_down < nx and col_outputs[c_down] == out:
            downstream[:, c_down] |= raw_sat_2d[:, c]

    return dilated | downstream


def _estimate_alternate_smear_stamp(
    *,
    sci_pre_smear: np.ndarray,
    sci_sigma_pre: np.ndarray,
    flat_2d: np.ndarray,
    collat_smear_stamp: np.ndarray | None,
    collat_smear_err_stamp: np.ndarray | None,
    collat_sm_sat_cols: np.ndarray | None,
    off_target_smear: np.ndarray | None,
    off_target_smear_err: np.ndarray | None,
    raw_sat_2d: np.ndarray,
    sat_bleed_2d: np.ndarray,
    strap_col_mask: np.ndarray,
    source_mask_2d: np.ndarray,
    bleed_dilation_rows: int = 4,
    min_clean_rows: int = 5,
    flank_width: int = 3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate concavity-guarded collateral-anchored smear ``M_hat`` across ``(N_t, nx)``.

    Works identically in ``e-/frame`` (for :func:`recalibrate_tpf`) and in
    calibrated ``e-/s`` with ``flat_2d = 1`` and ``collat_smear_stamp = 0`` (for
    :func:`desmear_calibrated_cutout`).

    Returns
    -------
    m_hat : np.ndarray
        Estimated per-cadence, per-column smear ``(N_t, nx)`` in same units as ``sci_pre_smear``.
    m_err : np.ndarray
        Column-coherent smear uncertainty ``(N_t, nx)`` in same units as ``sci_pre_smear``.
    unresolved_col_mask : np.ndarray
        Boolean ``(nx,)`` mask of bleed columns with ``< min_clean_rows`` clean rows.
    """
    n_t, ny, nx = sci_pre_smear.shape
    safe_flat = np.where(np.isfinite(flat_2d) & (flat_2d > 0.0), flat_2d, 1.0)

    # Clean unsaturated rows after >= 4-row vertical bleed-tip dilation + SOURCE_MASK
    dilated_bleed = _dilate_vertical_bool(raw_sat_2d, radius=bleed_dilation_rows)
    clean_rows_2d = (~dilated_bleed) & (~source_mask_2d)  # (ny, nx)

    col_has_sci_sat = np.any(raw_sat_2d, axis=0)
    col_has_sm_sat = (
        np.asarray(collat_sm_sat_cols, dtype=bool)
        if collat_sm_sat_cols is not None
        else np.zeros(nx, dtype=bool)
    )

    # Detect columns where collateral smrow bloomed (drives clean science rows negative)
    col_has_bloom = np.zeros(nx, dtype=bool)
    if collat_smear_stamp is not None:
        for c in range(nx):
            r_mask = clean_rows_2d[:, c]
            if not np.any(r_mask):
                col_has_bloom[c] = True
                continue
            sub = sci_pre_smear[:, r_mask, c] - collat_smear_stamp[:, None, c]
            finite_t = np.any(np.isfinite(sub), axis=1)
            if not np.any(finite_t):
                continue
            min_per_t = np.nanmin(sub[finite_t], axis=1)
            med_min = float(np.nanmedian(min_per_t))
            sig_med = float(np.nanmedian(sci_sigma_pre[finite_t][:, r_mask, c]))
            if med_min < -2.5 * max(sig_med, 1e-3):
                col_has_bloom[c] = True
    else:
        # FFI cutout mode (no collateral smrow): bloomed columns have negative clean rows
        for c in range(nx):
            r_mask = clean_rows_2d[:, c]
            if not np.any(r_mask):
                col_has_bloom[c] = True
                continue
            sub = sci_pre_smear[:, r_mask, c]
            finite_t = np.any(np.isfinite(sub), axis=1)
            if not np.any(finite_t):
                continue
            min_per_t = np.nanmin(sub[finite_t], axis=1)
            med_min = float(np.nanmedian(min_per_t))
            sig_med = float(np.nanmedian(sci_sigma_pre[finite_t][:, r_mask, c]))
            if med_min < -2.5 * max(sig_med, 1e-3):
                col_has_bloom[c] = True

    c_bleed_mask = col_has_sci_sat | col_has_sm_sat | col_has_bloom  # (nx,)

    # Initialize m_hat / m_err with collateral smear (or 0.0 in FFI desmear mode)
    if collat_smear_stamp is not None:
        m_hat = np.array(collat_smear_stamp, dtype=np.float64, copy=True)
        m_err = (
            np.array(collat_smear_err_stamp, dtype=np.float64, copy=True)
            if collat_smear_err_stamp is not None
            else np.full((n_t, nx), 1.0, dtype=np.float64)
        )
    else:
        m_hat = np.zeros((n_t, nx), dtype=np.float64)
        m_err = np.zeros((n_t, nx), dtype=np.float64)

    if not np.any(c_bleed_mask):
        return m_hat, m_err, np.zeros(nx, dtype=bool)

    # Check if ANY non-bleeding anchor column exists in the stamp (excluding STRAP_COL)
    non_bleed_non_strap = (~c_bleed_mask) & (~strap_col_mask)

    # Narrow-stamp upper-register bloom fallback (e.g. Sirius S6 11-col stamp where
    # smrow blooms across all 11 columns, while outer science columns are unsaturated):
    # anchor unsaturated outer columns with off-target clean collateral smrow baseline!
    if (
        not np.any(non_bleed_non_strap)
        and off_target_smear is not None
        and np.any(~col_has_sci_sat)
    ):
        edge_candidates = (~col_has_sci_sat) & (~strap_col_mask)
        if not np.any(edge_candidates):
            edge_candidates = ~col_has_sci_sat
        for c in np.flatnonzero(edge_candidates):
            m_hat[:, c] = off_target_smear
            if off_target_smear_err is not None:
                m_err[:, c] = off_target_smear_err
            c_bleed_mask[c] = False

    if not np.any(c_bleed_mask):
        return m_hat, m_err, np.zeros(nx, dtype=bool)

    # Precompute flat-corrected anchor signal U_{t,r,c'} = (I_{t,r,c'} - M_{t,c'}) / F_{r,c'}
    u_anchor_all = (sci_pre_smear - m_hat[:, None, :]) / safe_flat[None, :, :]
    anchor_pixel_ok = (
        (~c_bleed_mask)[None, :]
        & (~strap_col_mask)[None, :]
        & (~sat_bleed_2d)
        & (~source_mask_2d)
    )
    # If every non-bleed column happened to fall on a strap, relax strap exclusion rather than fail
    if not np.any(anchor_pixel_ok) and np.any(~c_bleed_mask):
        anchor_pixel_ok = (~c_bleed_mask)[None, :] & (~sat_bleed_2d) & (~source_mask_2d)

    bleed_cols = np.flatnonzero(c_bleed_mask)
    unresolved_col_mask = np.zeros(nx, dtype=bool)
    resolvable_cols: list[int] = []

    for c in bleed_cols:
        n_clean = int(np.sum(clean_rows_2d[:, c]))
        if n_clean < min_clean_rows:
            unresolved_col_mask[int(c)] = True
        else:
            resolvable_cols.append(int(c))

    if resolvable_cols:
        n_res = len(resolvable_cols)
        j_hat = np.zeros((n_t, ny, n_res), dtype=np.float64)
        flank_smear_var = np.zeros((n_t, n_res), dtype=np.float64)
        clamp_grad_var = np.zeros((n_t, n_res), dtype=np.float64)
        clean_sub = clean_rows_2d[:, resolvable_cols]  # (ny, n_res)

        for idx_b, c in enumerate(resolvable_cols):
            # Vectorized nearest-flank_width anchor selection across all ny rows:
            # Left flank: columns c-1, c-2, ..., 0 (nearest first)
            if c > 0:
                ok_l_rev = anchor_pixel_ok[:, c - 1 :: -1]  # (ny, c)
                sel_l = ok_l_rev & (np.cumsum(ok_l_rev, axis=1) <= flank_width)
                n_l = np.sum(sel_l, axis=1)  # (ny,)
                has_l = n_l > 0
                xs_l_rev = np.arange(c - 1, -1, -1, dtype=np.float64)
                c_l = np.sum(np.where(sel_l, xs_l_rev[None, :], 0.0), axis=1) / np.maximum(
                    n_l, 1
                )
                u_l_rev = u_anchor_all[:, :, c - 1 :: -1]  # (n_t, ny, c)
                u_l = np.sum(
                    np.where(sel_l[None, :, :], u_l_rev, 0.0), axis=2
                ) / np.maximum(n_l[None, :], 1)
                m_var_l_rev = m_err[:, c - 1 :: -1] ** 2  # (n_t, c)
                v_l = np.sum(
                    np.where(sel_l[None, :, :], m_var_l_rev[:, None, :], 0.0), axis=2
                ) / np.maximum(n_l[None, :] * n_l[None, :], 1)
            else:
                ok_l_rev = np.zeros((ny, 0), dtype=bool)
                xs_l_rev = np.zeros(0, dtype=np.float64)
                u_l_rev = np.zeros((n_t, ny, 0), dtype=np.float64)
                has_l = np.zeros(ny, dtype=bool)
                c_l = np.zeros(ny, dtype=np.float64)
                u_l = np.zeros((n_t, ny), dtype=np.float64)
                v_l = np.zeros((n_t, ny), dtype=np.float64)

            # Right flank: columns c+1, c+2, ..., nx-1 (nearest first)
            if c + 1 < nx:
                ok_r = anchor_pixel_ok[:, c + 1 :]  # (ny, nx - c - 1)
                sel_r = ok_r & (np.cumsum(ok_r, axis=1) <= flank_width)
                n_r = np.sum(sel_r, axis=1)  # (ny,)
                has_r = n_r > 0
                xs_r = np.arange(c + 1, nx, dtype=np.float64)
                c_r = np.sum(np.where(sel_r, xs_r[None, :], 0.0), axis=1) / np.maximum(
                    n_r, 1
                )
                u_r_sub = u_anchor_all[:, :, c + 1 :]
                u_r = np.sum(
                    np.where(sel_r[None, :, :], u_r_sub, 0.0), axis=2
                ) / np.maximum(n_r[None, :], 1)
                m_var_r = m_err[:, c + 1 :] ** 2
                v_r = np.sum(
                    np.where(sel_r[None, :, :], m_var_r[:, None, :], 0.0), axis=2
                ) / np.maximum(n_r[None, :] * n_r[None, :], 1)
            else:
                ok_r = np.zeros((ny, 0), dtype=bool)
                xs_r = np.zeros(0, dtype=np.float64)
                u_r_sub = np.zeros((n_t, ny, 0), dtype=np.float64)
                has_r = np.zeros(ny, dtype=bool)
                c_r = np.zeros(ny, dtype=np.float64)
                u_r = np.zeros((n_t, ny), dtype=np.float64)
                v_r = np.zeros((n_t, ny), dtype=np.float64)

            both = has_l & has_r
            only_l = has_l & (~has_r)
            only_r = has_r & (~has_l)
            row_has_anchor = clean_sub[:, idx_b] & (has_l | has_r)

            if not np.any(row_has_anchor):
                # No anchor columns available anywhere in stamp: fallback to off-target smear
                clean_sub[:, idx_b] = False
                if off_target_smear is not None:
                    m_hat[:, c] = off_target_smear
                    m_err[:, c] = (
                        off_target_smear_err
                        if off_target_smear_err is not None
                        else np.full(n_t, 5.0, dtype=np.float64)
                    )
                else:
                    unresolved_col_mask[c] = True
                continue

            clean_sub[:, idx_b] &= row_has_anchor
            if int(np.sum(clean_sub[:, idx_b])) < min_clean_rows:
                clean_sub[:, idx_b] = False
                unresolved_col_mask[c] = True
                continue

            denom_lr = np.where(both, c_r - c_l, 1.0)
            alpha = (float(c) - c_l) / denom_lr
            u_both = (1.0 - alpha[None, :]) * u_l + alpha[None, :] * u_r
            v_both = ((1.0 - alpha[None, :]) ** 2) * v_l + (alpha[None, :] ** 2) * v_r

            u_comb = np.where(
                both[None, :],
                u_both,
                np.where(only_l[None, :], u_l, np.where(only_r[None, :], u_r, 0.0)),
            )
            j_hat[:, :, idx_b] = safe_flat[:, c][None, :] * u_comb
            row_flank_var = np.where(
                both[None, :],
                v_both,
                np.where(only_l[None, :], v_l, np.where(only_r[None, :], v_r, 0.0)),
            )
            row_clamp_dist = np.where(
                only_l,
                np.abs(float(c) - c_l),
                np.where(only_r, np.abs(float(c) - c_r), 0.0),
            )

            # Vectorized 4-column horizontal gradient on one-sidedly clamped rows
            row_grad = np.zeros((n_t, ny), dtype=np.float64)
            r_ok = clean_sub[:, idx_b]
            clamped_l_rows = np.flatnonzero(r_ok & only_l)
            if clamped_l_rows.size > 0 and c > 0:
                ok4 = ok_l_rev[clamped_l_rows]
                sel4 = ok4 & (np.cumsum(ok4, axis=1) <= 4)
                cnt4 = np.sum(sel4, axis=1)
                xm = np.sum(np.where(sel4, xs_l_rev[None, :], 0.0), axis=1) / np.maximum(
                    cnt4, 1
                )
                dx = np.where(sel4, xs_l_rev[None, :] - xm[:, None], 0.0)
                den = np.sum(dx * dx, axis=1)
                slp = np.sum(
                    np.where(sel4[None, :, :], u_l_rev[:, clamped_l_rows, :], 0.0)
                    * dx[None, :, :],
                    axis=2,
                ) / np.where(den[None, :] > 0.0, den[None, :], 1.0)
                has2 = (cnt4 >= 2) & (den > 0.0)
                row_grad[:, clamped_l_rows] = np.where(
                    has2[None, :],
                    np.abs(slp) * safe_flat[clamped_l_rows, c][None, :],
                    0.05 * np.abs(j_hat[:, clamped_l_rows, idx_b]),
                )

            clamped_r_rows = np.flatnonzero(r_ok & only_r)
            if clamped_r_rows.size > 0 and (c + 1 < nx):
                ok4 = ok_r[clamped_r_rows]
                sel4 = ok4 & (np.cumsum(ok4, axis=1) <= 4)
                cnt4 = np.sum(sel4, axis=1)
                xm = np.sum(np.where(sel4, xs_r[None, :], 0.0), axis=1) / np.maximum(
                    cnt4, 1
                )
                dx = np.where(sel4, xs_r[None, :] - xm[:, None], 0.0)
                den = np.sum(dx * dx, axis=1)
                slp = np.sum(
                    np.where(sel4[None, :, :], u_r_sub[:, clamped_r_rows, :], 0.0)
                    * dx[None, :, :],
                    axis=2,
                ) / np.where(den[None, :] > 0.0, den[None, :], 1.0)
                has2 = (cnt4 >= 2) & (den > 0.0)
                row_grad[:, clamped_r_rows] = np.where(
                    has2[None, :],
                    np.abs(slp) * safe_flat[clamped_r_rows, c][None, :],
                    0.05 * np.abs(j_hat[:, clamped_r_rows, idx_b]),
                )

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                flank_smear_var[:, idx_b] = np.nanmedian(row_flank_var[:, r_ok], axis=1)
                clamped_rows = r_ok & (row_clamp_dist > 0.0)
                if np.any(clamped_rows):
                    # Evaluate gradient inflation over the furthest low-flux clamped rows
                    j_med_r = np.nanmedian(j_hat[:, clamped_rows, idx_b], axis=0)
                    if np.any(np.isfinite(j_med_r)):
                        q_thresh = float(np.nanpercentile(j_med_r, 35.0))
                        far_clamped = np.flatnonzero(clamped_rows)[j_med_r <= q_thresh]
                    else:
                        far_clamped = np.flatnonzero(clamped_rows)
                    if far_clamped.size == 0:
                        far_clamped = np.flatnonzero(clamped_rows)
                    grad_term = (
                        row_clamp_dist[None, far_clamped] * row_grad[:, far_clamped]
                    ) ** 2
                    clamp_grad_var[:, idx_b] = np.nanmedian(grad_term, axis=1)

        active_rel = [i for i in range(n_res) if np.any(clean_sub[:, i])]
        if active_rel:
            act_cols = [resolvable_cols[i] for i in active_rel]
            j_act = j_hat[:, :, active_rel]
            clean_act = clean_sub[:, active_rel]
            i_act = sci_pre_smear[:, :, act_cols]
            sig_act = np.maximum(sci_sigma_pre[:, :, act_cols], 1e-6)

            # Inverse-wing weights: omega_{t,r,c} = 1 / ((J_hat - min_{r'} J_hat)^2 + sigma_I^2)
            j_masked = np.where(clean_act[None, :, :], j_act, np.inf)
            j_min = np.min(j_masked, axis=1, keepdims=True)
            j_min = np.where(np.isfinite(j_min), j_min, 0.0)
            wing_excess = np.where(clean_act[None, :, :], j_act - j_min, 0.0)
            inv_wing_w = np.where(
                clean_act[None, :, :],
                1.0 / (wing_excess * wing_excess + sig_act * sig_act),
                0.0,
            )

            residuals = i_act - j_act
            loc_est, loc_err = _solve_huber_smear_location(
                residuals,
                inv_wing_w,
                sig_act,
                clean_act,
            )

            for k_idx, col_idx in enumerate(act_cols):
                rel_i = active_rel[k_idx]
                m_hat[:, col_idx] = loc_est[:, k_idx]
                tot_var = (
                    loc_err[:, k_idx] ** 2
                    + flank_smear_var[:, rel_i]
                    + clamp_grad_var[:, rel_i]
                )
                m_err[:, col_idx] = np.sqrt(np.maximum(tot_var, 1e-12))

    # Fill any SMEAR_UNRESOLVED columns conservatively so FLUX/SMEAR_ERR remain finite
    if np.any(unresolved_col_mask):
        resolved_cols = np.flatnonzero(~unresolved_col_mask)
        for c_unres in np.flatnonzero(unresolved_col_mask):
            if off_target_smear is not None:
                m_hat[:, c_unres] = off_target_smear
            elif resolved_cols.size > 0:
                nearest = int(resolved_cols[np.argmin(np.abs(resolved_cols - c_unres))])
                m_hat[:, c_unres] = m_hat[:, nearest]
            finite_cad = np.any(np.isfinite(sci_pre_smear[:, :, c_unres]), axis=1)
            fallback_err = (
                float(np.nanmedian(m_err[finite_cad])) * 10.0 + 50.0
                if np.any(finite_cad) and np.any(np.isfinite(m_err[finite_cad]))
                else 100.0
            )
            m_err[:, c_unres] = np.where(finite_cad, fallback_err, np.nan)

    return m_hat, m_err, unresolved_col_mask


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
    """Stream ``tvcol`` in ``chunk_rows`` slices and publish ``tvcol-...-blf.fits``."""
    tv_hdr0, tv_hdr1, tv_time, tv_cad = read_collateral_headers_and_index(tv_col_path)
    n_total = len(tv_cad)
    n_eff_tv = effective_frame_count(
        tv_hdr0, tv_hdr1, fallback_hdr0=sci_hdr0, fallback_hdr1=sci_hdr1
    )

    fxdoff = float(
        tv_hdr1["FXDOFF"]
        if "FXDOFF" in tv_hdr1
        else (sci_hdr1["FXDOFF"] if sci_hdr1 is not None and "FXDOFF" in sci_hdr1 else 209700.0)
    )
    meanblc_key = f"MEANBLC{output}"
    if meanblc_key in tv_hdr1:
        meanblc = float(tv_hdr1[meanblc_key])
    elif sci_hdr1 is not None and meanblc_key in sci_hdr1:
        meanblc = float(sci_hdr1[meanblc_key])
    else:
        raise KeyError(f"Missing {meanblc_key} in collateral and science FITS headers.")

    _, _, tx1, tx2 = OUTPUT_BOUNDS_0BASED[output]
    b2d_tv = np.asarray(black2d_img[:, tx1:tx2], dtype=np.float64)  # (2078, 11)
    bias_shift = -fxdoff + n_eff_tv * (meanblc - b2d_tv)[np.newaxis, :, :]

    b1d_all = np.full((n_total, 2078), np.nan, dtype=np.float64)
    b1d_err_all = np.full((n_total, 2078), np.nan, dtype=np.float32)

    for t0, t1, _t_chunk, _c_chunk, raw_chunk in iter_collateral_chunks(
        tv_col_path, "TVCOL_RAW", chunk_rows=chunk_rows
    ):
        debiased = raw_chunk + bias_shift
        b1d_chunk, err_chunk = fit_onedblack_legendre(debiased)
        b1d_all[t0:t1] = np.asarray(b1d_chunk, dtype=np.float64)
        b1d_err_all[t0:t1] = np.asarray(err_chunk, dtype=np.float32)

    time_unit = str(tv_hdr1.get("TUNIT1", "BJD - 2457000, days"))
    return write_onedblack_blf(
        blf_dst,
        time=tv_time,
        cadenceno=tv_cad,
        black1d=b1d_all,
        black1d_err=b1d_err_all,
        n_eff=n_eff_tv,
        time_unit=time_unit,
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calibrate a ``(N_c, 10, 512)`` ``SMROW_RAW`` chunk to ``(smear, smear_err, sm_sat_col)``."""
    scale_1d = n_eff_sm / n_eff_tv
    debiased = (
        sm_raw_chunk
        - fxdoff
        + n_eff_sm * (meanblc - b2d_sm)[np.newaxis, :, :]
        - scale_1d * b1d_sm[:, :, np.newaxis]
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
        smear_med = np.nanmedian(sm_e_np, axis=1)  # (N_c, 512)
        mad_rows = 1.4826 * np.nanmedian(
            np.abs(sm_e_np - smear_med[:, np.newaxis, :]), axis=1
        )

    n_valid = np.sum(np.isfinite(sm_e_np), axis=1)  # (N_c, 512)
    emp_var = (np.pi / (2.0 * np.maximum(n_valid, 1))) * (mad_rows * mad_rows)
    row_shot_read = (np.maximum(smear_med, 0.0) + read_noise_e * read_noise_e) / n_eff_sm
    theo_var = (np.pi / (2.0 * np.maximum(n_valid, 1))) * row_shot_read
    if b1d_err_sm is not None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean_g = np.nanmedian(sm_gain_np, axis=1)
            b1d_sm_med = np.nanmean(b1d_err_sm, axis=1)[:, np.newaxis]
        mean_g = np.where(np.isfinite(mean_g), mean_g, gain)
        b1d_var = ((mean_g * b1d_sm_med) / n_eff_tv) ** 2
        b1d_var = np.where(np.isfinite(b1d_var), b1d_var, 0.0)
    else:
        b1d_var = 0.0

    smear_err = np.sqrt(np.maximum(emp_var, theo_var) + b1d_var)
    smear_err = np.where(n_valid > 0, smear_err, np.nan)
    sm_sat_cols = np.any(sm_sat_np, axis=(0, 1))  # (512,)
    return smear_med, smear_err, sm_sat_cols


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
    """Stream ``smrow`` in ``chunk_rows`` slices and publish standard ``smrow-...-blf.fits``."""
    sm_hdr0, sm_hdr1, sm_time, sm_cad = read_collateral_headers_and_index(sm_col_path)
    n_total = len(sm_cad)
    n_eff_sm = effective_frame_count(
        sm_hdr0, sm_hdr1, fallback_hdr0=sci_hdr0, fallback_hdr1=sci_hdr1
    )

    with fits.open(tv_blf_path, memmap=True) as tv_hdul:
        tv_hdu1 = tv_hdul[1]
        tv_cad = np.asarray(tv_hdu1.data["CADENCENO"], dtype=np.int64)
        tv_b1d_sm = np.asarray(tv_hdu1.data["1DBLACK"][:, 2058:2068], dtype=np.float64)
        has_err = "1DBLACK_ERR" in (tv_hdu1.data.names or ())
        tv_err_sm = (
            np.asarray(tv_hdu1.data["1DBLACK_ERR"][:, 2058:2068], dtype=np.float64)
            if has_err
            else None
        )
        n_eff_tv = float(tv_hdu1.header.get("NEFF", n_eff_sm))

    matched_tv, tv_idx = match_cadenceno(
        sm_cad, tv_cad, cadence_match=cadence_match, label=f"tvcol-blf ({output})"
    )

    def _hdr_float(key: str, default: float | None = None) -> float:
        for hdr in (sm_hdr1, sci_hdr1):
            if hdr is not None and key in hdr:
                return float(hdr[key])
        if default is not None:
            return default
        raise KeyError(f"Missing {key} in collateral and science FITS headers.")

    fxdoff = _hdr_float("FXDOFF", 209700.0)
    meanblc = _hdr_float(f"MEANBLC{output}")
    gain = _hdr_float(f"GAIN{output}")
    read_noise_e = _hdr_float(f"READNOI{output}", 8.5)
    int_time = _hdr_float("INT_TIME", 1.98)

    sx1, sx2, _, _ = OUTPUT_BOUNDS_0BASED[output]
    b2d_sm = np.asarray(black2d_img[2058:2068, sx1:sx2], dtype=np.float64)

    smear_all = np.full((n_total, 512), np.nan, dtype=np.float64)
    smear_err_all = np.full((n_total, 512), np.nan, dtype=np.float32)

    for t0, t1, _t_chunk, _c_chunk, raw_chunk in iter_collateral_chunks(
        sm_col_path, "SMROW_RAW", chunk_rows=chunk_rows
    ):
        chunk_matched = matched_tv[t0:t1]
        if not np.any(chunk_matched):
            continue
        safe_tv_idx = np.where(chunk_matched, tv_idx[t0:t1], 0)
        b1d_chunk = np.where(chunk_matched[:, None], tv_b1d_sm[safe_tv_idx], np.nan)
        err_chunk = (
            np.where(chunk_matched[:, None], tv_err_sm[safe_tv_idx], np.nan)
            if tv_err_sm is not None
            else None
        )
        sm_med, sm_err, _ = _calibrate_smrow_chunk_e_per_frame(
            raw_chunk,
            output=output,
            fxdoff=fxdoff,
            meanblc=meanblc,
            gain=gain,
            read_noise_e=read_noise_e,
            n_eff_sm=n_eff_sm,
            n_eff_tv=n_eff_tv,
            b2d_sm=b2d_sm,
            b1d_sm=b1d_chunk,
            b1d_err_sm=err_chunk,
            u_coeff=u_coeff,
            lin_channel=lin_channel,
        )
        smear_all[t0:t1] = np.where(chunk_matched[:, None], sm_med, np.nan)
        smear_err_all[t0:t1] = np.where(
            chunk_matched[:, None], sm_err.astype(np.float32), np.float32(np.nan)
        )

    time_unit = str(sm_hdr1.get("TUNIT1", "BJD - 2457000, days"))
    return write_smear_blf(
        blf_dst,
        time=sm_time,
        cadenceno=sm_cad,
        smear=smear_all,
        smear_err=smear_err_all,
        smear_mode="standard",
        int_time=int_time,
        time_unit=time_unit,
        stage_dir=stage_dir,
    )


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
    cleanup_staged_raw: bool = True,
    s3_fs: Any = None,
    verbose: bool = False,
) -> tuple[Path, Path]:
    """Build or reuse ``tvcol`` and ``smrow`` ``-blf.fits`` caches for one CCD output."""
    out_upper = output.upper()
    if out_upper not in OUTPUT_ORDER:
        raise ValueError(f"Invalid CCD output {output!r}; expected one of {OUTPUT_ORDER}.")

    resolved_cache = (
        Path(cache_dir).expanduser()
        if cache_dir is not None
        else resolve_cache_dir(create=True)
    )
    sec_dir = resolved_cache / f"s{sector:04d}"

    tv_cached = find_cached_blf(
        sector=sector,
        camera=camera,
        ccd=ccd,
        output=out_upper,
        kind="tvcol",
        cadence=cadence,
        smear_mode="standard",
        ddir=ddir,
        cache_dir=resolved_cache,
    )
    sm_cached = find_cached_blf(
        sector=sector,
        camera=camera,
        ccd=ccd,
        output=out_upper,
        kind="smrow",
        cadence=cadence,
        smear_mode="standard",
        ddir=ddir,
        cache_dir=resolved_cache,
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
        tv_raw_path = resolve_collateral_file(
            sector=sector,
            camera=camera,
            ccd=ccd,
            output=out_upper,
            kind="tvcol",
            cadence=cadence,
            ddir=ddir,
            cache_dir=resolved_cache,
            offline=offline,
            stage_dir=stage_dir,
            keep_raw_in_stage_only=cleanup_staged_raw,
            s3_fs=s3_fs,
            verbose=verbose,
        )
        tv_dst = collateral_blf_cache_path(
            sec_dir,
            sector=sector,
            camera=camera,
            ccd=ccd,
            output=out_upper,
            kind="tvcol",
            cadence=cadence,
            col_stem=tv_raw_path.stem,
        )
        tv_cached = _calibrate_tvcol_raw_to_blf(
            tv_raw_path,
            tv_dst,
            output=out_upper,
            black2d_img=black2d_img,
            chunk_rows=chunk_rows,
            stage_dir=stage_dir,
        )
        if (
            cleanup_staged_raw
            and stage_root is not None
            and tv_raw_path.resolve(strict=False).parent == stage_root
        ):
            with contextlib.suppress(OSError):
                tv_raw_path.unlink()

    if sm_cached is None:
        sm_raw_path = resolve_collateral_file(
            sector=sector,
            camera=camera,
            ccd=ccd,
            output=out_upper,
            kind="smrow",
            cadence=cadence,
            ddir=ddir,
            cache_dir=resolved_cache,
            offline=offline,
            stage_dir=stage_dir,
            keep_raw_in_stage_only=cleanup_staged_raw,
            s3_fs=s3_fs,
            verbose=verbose,
        )
        sm_dst = collateral_blf_cache_path(
            sec_dir,
            sector=sector,
            camera=camera,
            ccd=ccd,
            output=out_upper,
            kind="smrow",
            cadence=cadence,
            smear_mode="standard",
            col_stem=sm_raw_path.stem,
        )
        sm_cached = _calibrate_smrow_raw_to_standard_blf(
            sm_raw_path,
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
            and sm_raw_path.resolve(strict=False).parent == stage_root
        ):
            with contextlib.suppress(OSError):
                sm_raw_path.unlink()

    return tv_cached, sm_cached


def _smear_cache_matches_cutout(
    hdr: fits.Header,
    smear_mode: str,
    cutout_bounds: tuple[int, int, int, int],
) -> bool:
    """Verify a cached smear ``-blf.fits`` matches ``smear_mode`` and cutout bounds."""
    if "SMEARMOD" in hdr and str(hdr["SMEARMOD"]).strip() != smear_mode:
        return False
    if smear_mode == "alternate":
        crval1p, cutnx, crval2p, cutny = cutout_bounds
        expected = (
            ("CRVAL1P", crval1p),
            ("CUTNX", cutnx),
            ("CRVAL2P", crval2p),
            ("CUTNY", cutny),
        )
        for key, val in expected:
            if key in hdr and int(hdr[key]) != int(val):
                return False
    return True


def _extract_raw_cnts_float64(hdu1: fits.BinTableHDU) -> np.ndarray:
    """Extract ``RAW_CNTS`` as ``float64`` with ``-1`` / ``TNULL`` sentinels masked to ``NaN``."""
    names = hdu1.data.names or ()
    if "RAW_CNTS" not in names:
        raise ValueError(
            "TPF HDU[1] lacks a 'RAW_CNTS' column; use mundey.desmear_calibrated_cutout(tpf) "
            "for TessCut / FFI cutouts without collateral data."
        )
    raw_arr = hdu1.data["RAW_CNTS"]
    raw_f64 = np.asarray(raw_arr, dtype=np.float64)

    null_val: int = -1
    for idx, name in enumerate(hdu1.columns.names or (), start=1):
        if name == "RAW_CNTS":
            tnull_key = f"TNULL{idx}"
            if tnull_key in hdu1.header:
                with contextlib.suppress(TypeError, ValueError):
                    null_val = int(hdu1.header[tnull_key])
            break

    bad = (raw_f64 < 0.0) | (~np.isfinite(raw_f64))
    if np.issubdtype(raw_arr.dtype, np.integer):
        bad = bad | (raw_arr == null_val) | (raw_arr == -1)
    if np.any(bad):
        raw_f64 = raw_f64.copy()
        raw_f64[bad] = np.nan
    return raw_f64


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
    target_obj: Any = None,
) -> None:
    """Execute the v2.0 calibration pipeline on ``hdul`` in place."""
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

    hdr0 = hdul[0].header
    hdu1 = hdul[1]
    hdr1 = hdu1.header
    hdr2 = hdul[2].header

    sector = int(hdr0["SECTOR"])
    camera = int(hdr0["CAMERA"])
    ccd = int(hdr0["CCD"])
    ticid = int(hdr0["TICID"]) if "TICID" in hdr0 else None

    raw_f64 = _extract_raw_cnts_float64(hdu1)
    n_tpf, ny, nx = raw_f64.shape

    col0_1based = int(hdr2["CRVAL1P"])
    row0_1based = int(hdr2["CRVAL2P"])
    x1_ccd0 = col0_1based - 1
    y1_ccd0 = row0_1based - 1
    cutout_bounds = (col0_1based, nx, row0_1based, ny)

    cols_1based = np.arange(col0_1based, col0_1based + nx, dtype=np.int64)
    col_outputs = tuple(output_for_rawx_col(int(c)) for c in cols_1based)
    unique_outputs = tuple(out for out in OUTPUT_ORDER if out in col_outputs)

    if row0_1based < 1 or y1_ccd0 + ny > 2048:
        raise ValueError(
            f"TPF rows {row0_1based}..{row0_1based + ny - 1} fall outside "
            "TESS imaging rows 1..2048."
        )

    eff_smear: Literal["standard", "alternate"] = (
        ("alternate" if (y1_ccd0 + ny >= 2048) else "standard")
        if smear == "automatic"
        else smear  # type: ignore[assignment]
    )

    # Resolve cache / write directory: when ddir is explicitly supplied without cache_dir,
    # write -blf.fits into ddir for seamless drop-in compatibility; otherwise use cache_dir.
    resolved_cache: Path | None
    if cache_dir is not None:
        resolved_cache = resolve_cache_dir(cache_dir, strict_hpc=strict_hpc, create=True)
    elif ddir is not None:
        if strict_hpc:
            resolve_cache_dir(ddir, strict_hpc=True, create=False)
        resolved_cache = None
    else:
        resolved_cache = resolve_cache_dir(None, strict_hpc=strict_hpc, create=True)

    blf_write_base = (
        (resolved_cache / f"s{sector:04d}")
        if resolved_cache is not None
        else Path(ddir).expanduser()  # type: ignore[arg-type]
    )

    cadence_kind = infer_cadence_kind(
        tpf_path=source_path,
        hdu0_header=hdr0,
        hdu1_header=hdr1,
    )
    tpf_parts = source_path.name.split("-") if source_path is not None else None
    if tpf_parts is not None and len(tpf_parts) < 4:
        tpf_parts = None

    if verbose:
        print(
            f"Recalibrating Sector {sector} Camera {camera} CCD {ccd} "
            f"Outputs {list(unique_outputs)} (smear={eff_smear!r}, cadence={cadence_kind!r})"
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
            flat_hdul[1].data[y1_ccd0 : y1_ccd0 + ny, x1_ccd0 : x1_ccd0 + nx],
            dtype=np.float64,
        )

    n_eff_sci = effective_frame_count(hdr0, hdr1)
    fxdoff_sci = float(hdr1.get("FXDOFF", 209700.0))
    int_time = float(hdr1["INT_TIME"])
    tpf_time = np.asarray(hdu1.data["TIME"], dtype=np.float64)
    tpf_cad = np.asarray(hdu1.data["CADENCENO"], dtype=np.int64)
    time_unit = str(hdr1.get("TUNIT1", "BJD - 2457000, days"))

    # Preserve incoming SOURCE_MASK if caller pre-populated PIXEL_QUALITY
    existing_pq = (
        np.asarray(hdul[PIXEL_QUALITY_EXTNAME].data, dtype=np.uint16)
        if PIXEL_QUALITY_EXTNAME in hdul
        else np.zeros((ny, nx), dtype=np.uint16)
    )
    source_mask_2d = (existing_pq & int(PixelQuality.SOURCE_MASK)) != 0
    strap_col_mask = strap_mask_for_rawx_columns(cols_1based, buffer_px=2)

    # Per-output collateral loading, CADENCENO alignment, and pre-smear science calibration
    combined_matched = np.ones(n_tpf, dtype=bool)
    sci_e_frame = np.full((n_tpf, ny, nx), np.nan, dtype=np.float64)
    sci_sigma_pre = np.full((n_tpf, ny, nx), np.nan, dtype=np.float64)
    sci_sat_3d = np.zeros((n_tpf, ny, nx), dtype=bool)

    collat_smear_stamp = np.full((n_tpf, nx), np.nan, dtype=np.float64)
    collat_smear_err_stamp = np.full((n_tpf, nx), np.nan, dtype=np.float64)
    collat_sm_sat_cols = np.zeros(nx, dtype=bool)

    alt_cached_smear_stamp = np.full((n_tpf, nx), np.nan, dtype=np.float64)
    alt_cached_smear_err_stamp = np.full((n_tpf, nx), np.nan, dtype=np.float64)
    all_outputs_have_alt_cache = eff_smear == "alternate"

    # Store per-output 512-col collateral smear for alternate-mode cache writing & off-target anchor
    per_out_collat_512: dict[str, np.ndarray] = {}
    per_out_collat_err_512: dict[str, np.ndarray] = {}
    off_target_smears: list[np.ndarray] = []
    off_target_errs: list[np.ndarray] = []

    for out in unique_outputs:
        local_cols = np.flatnonzero(np.asarray([o == out for o in col_outputs], dtype=bool))
        c0_loc, c1_loc = int(local_cols[0]), int(local_cols[-1]) + 1
        sm_col_start = int(cols_1based[c0_loc] - 1 - OUTPUT_BOUNDS_0BASED[out][0])
        sm_col_end = sm_col_start + (c1_loc - c0_loc)

        lin_channel = load_linearity_channel(models.linearity_path, camera, ccd, out)
        u_coeff = load_undershoot_coeff(camera, ccd, out)
        gain_sci = float(hdr1[f"GAIN{out}"])
        meanblc_sci = float(hdr1[f"MEANBLC{out}"])
        readnoi_sci = float(hdr1.get(f"READNOI{out}", 8.5))

        # 1. Load or build tvcol 1D-black -blf.fits (always on collateral CADENCENO grid)
        tv_blf = find_cached_blf(
            sector=sector,
            camera=camera,
            ccd=ccd,
            output=out,
            kind="tvcol",
            cadence=cadence_kind,
            smear_mode="standard",
            ddir=ddir,
            cache_dir=resolved_cache,
            tpf_parts=tpf_parts,
        )
        if tv_blf is None:
            tv_raw_path = resolve_collateral_file(
                sector=sector,
                camera=camera,
                ccd=ccd,
                output=out,
                kind="tvcol",
                cadence=cadence_kind,
                ddir=ddir,
                cache_dir=resolved_cache,
                offline=offline,
                tpf_parts=tpf_parts,
                s3_fs=s3_fs,
                verbose=verbose,
            )
            tv_dst = collateral_blf_cache_path(
                blf_write_base,
                sector=sector,
                camera=camera,
                ccd=ccd,
                output=out,
                kind="tvcol",
                cadence=cadence_kind,
                tpf_parts=tpf_parts,
                col_stem=tv_raw_path.stem if tpf_parts is None else None,
            )
            tv_blf = _calibrate_tvcol_raw_to_blf(
                tv_raw_path,
                tv_dst,
                output=out,
                black2d_img=black2d_img,
                sci_hdr0=hdr0,
                sci_hdr1=hdr1,
                chunk_rows=chunk_rows,
            )

        with fits.open(tv_blf, memmap=True) as tv_hdul:
            tv_hdu1 = tv_hdul[1]
            tv_cad = np.asarray(tv_hdu1.data["CADENCENO"], dtype=np.int64)
            tv_matched, tv_idx = match_cadenceno(
                tpf_cad,
                tv_cad,
                cadence_match=cadence_match,
                label=f"tvcol ({out})",
            )
            combined_matched &= tv_matched
            safe_tv_idx = np.where(tv_matched, tv_idx, 0)
            b1d_full = np.asarray(tv_hdu1.data["1DBLACK"], dtype=np.float64)[safe_tv_idx]
            has_tv_err = "1DBLACK_ERR" in (tv_hdu1.data.names or ())
            if has_tv_err:
                b1d_err_full = np.asarray(tv_hdu1.data["1DBLACK_ERR"], dtype=np.float64)[
                    safe_tv_idx
                ]
            else:
                b1d_err_full = np.full_like(b1d_full, 0.5)
            n_eff_tv = float(tv_hdu1.header.get("NEFF", n_eff_sci))

        b1d_sci_rows = np.where(
            tv_matched[:, None], b1d_full[:, y1_ccd0 : y1_ccd0 + ny], np.nan
        )
        b1d_err_sci_rows = np.where(
            tv_matched[:, None], b1d_err_full[:, y1_ccd0 : y1_ccd0 + ny], np.nan
        )

        # 2. Pre-smear science calibration for this output slice:
        #    de-bias -> causal LDE undershoot in ADU -> per-frame DN -> clamped linearity & gain
        b2d_sci_slice = black2d_img[
            y1_ccd0 : y1_ccd0 + ny,
            x1_ccd0 + c0_loc : x1_ccd0 + c1_loc,
        ]
        scale_1d_sci = n_eff_sci / n_eff_tv
        sci_debiased = (
            raw_f64[:, :, c0_loc:c1_loc]
            - fxdoff_sci
            + n_eff_sci * (meanblc_sci - b2d_sci_slice)[np.newaxis, :, :]
            - scale_1d_sci * b1d_sci_rows[:, :, np.newaxis]
        )
        e_out_np, sat_out_np, g_out_np = calibrate_pre_smear_e_per_frame(
            sci_debiased,
            u_coeff=u_coeff,
            n_eff=n_eff_sci,
            channel=lin_channel,
            gain=gain_sci,
            reverse_readout=OUTPUT_REVERSE_READOUT[out],
        )

        # Marginal pre-smear pixel variance in (e-/frame)^2
        n_loc = c1_loc - c0_loc
        u_mult = np.full(n_loc, 1.0 / (1.0 - u_coeff * u_coeff), dtype=np.float64)
        first_read_idx = (n_loc - 1) if OUTPUT_REVERSE_READOUT[out] else 0
        u_mult[first_read_idx] = 1.0 / ((1.0 - u_coeff) ** 2)

        shot_read_var = (np.maximum(e_out_np, 0.0) + readnoi_sci * readnoi_sci) / n_eff_sci
        b1d_var_sci = ((g_out_np * b1d_err_sci_rows[:, :, np.newaxis]) / n_eff_tv) ** 2
        b1d_var_sci = np.where(np.isfinite(b1d_var_sci), b1d_var_sci, 0.0)
        sig_out = np.sqrt(u_mult[None, None, :] * (shot_read_var + b1d_var_sci))

        sci_e_frame[:, :, c0_loc:c1_loc] = e_out_np
        sci_sigma_pre[:, :, c0_loc:c1_loc] = sig_out
        sci_sat_3d[:, :, c0_loc:c1_loc] = sat_out_np

        # 3. Smear cache resolution:
        # Check target-scoped -alt-*-blf.fits first when eff_smear == 'alternate'
        used_alt_cache_this_out = False
        if eff_smear == "alternate":
            alt_blf = find_cached_blf(
                sector=sector,
                camera=camera,
                ccd=ccd,
                output=out,
                kind="smrow",
                cadence=cadence_kind,
                smear_mode="alternate",
                ddir=ddir,
                cache_dir=resolved_cache,
                tpf_parts=tpf_parts,
                ticid=ticid,
            )
            if alt_blf is not None:
                with fits.open(alt_blf, memmap=True) as alt_hdul:
                    alt_hdu1 = alt_hdul[1]
                    if _smear_cache_matches_cutout(
                        alt_hdu1.header, "alternate", cutout_bounds
                    ):
                        alt_cad = np.asarray(alt_hdu1.data["CADENCENO"], dtype=np.int64)
                        alt_matched, alt_idx = match_cadenceno(
                            tpf_cad,
                            alt_cad,
                            cadence_match=cadence_match,
                            label=f"smrow-alt-blf ({out})",
                        )
                        combined_matched &= alt_matched
                        safe_alt_idx = np.where(alt_matched, alt_idx, 0)
                        alt_sm_512 = np.asarray(alt_hdu1.data["SMEAR"], dtype=np.float64)[
                            safe_alt_idx
                        ]
                        if "SMEAR_ERR" in (alt_hdu1.data.names or ()):
                            alt_err_512 = np.asarray(
                                alt_hdu1.data["SMEAR_ERR"], dtype=np.float64
                            )[safe_alt_idx]
                        else:
                            alt_err_512 = np.sqrt(
                                (np.pi / 20.0)
                                * (np.maximum(alt_sm_512, 0.0) + readnoi_sci**2)
                                / n_eff_sci
                            )
                        alt_cached_smear_stamp[:, c0_loc:c1_loc] = np.where(
                            alt_matched[:, None],
                            alt_sm_512[:, sm_col_start:sm_col_end],
                            np.nan,
                        )
                        alt_cached_smear_err_stamp[:, c0_loc:c1_loc] = np.where(
                            alt_matched[:, None],
                            alt_err_512[:, sm_col_start:sm_col_end],
                            np.nan,
                        )
                        used_alt_cache_this_out = True

        if not used_alt_cache_this_out:
            all_outputs_have_alt_cache = False

            # Check standard collateral smrow -blf.fits cache
            std_sm_blf = find_cached_blf(
                sector=sector,
                camera=camera,
                ccd=ccd,
                output=out,
                kind="smrow",
                cadence=cadence_kind,
                smear_mode="standard",
                ddir=ddir,
                cache_dir=resolved_cache,
                tpf_parts=tpf_parts,
            )
            if std_sm_blf is not None:
                with fits.open(std_sm_blf, memmap=True) as sm_hdul:
                    if not _smear_cache_matches_cutout(
                        sm_hdul[1].header, "standard", cutout_bounds
                    ):
                        std_sm_blf = None

            if std_sm_blf is None and eff_smear == "standard":
                sm_raw_path = resolve_collateral_file(
                    sector=sector,
                    camera=camera,
                    ccd=ccd,
                    output=out,
                    kind="smrow",
                    cadence=cadence_kind,
                    ddir=ddir,
                    cache_dir=resolved_cache,
                    offline=offline,
                    tpf_parts=tpf_parts,
                    s3_fs=s3_fs,
                    verbose=verbose,
                )
                sm_dst = collateral_blf_cache_path(
                    blf_write_base,
                    sector=sector,
                    camera=camera,
                    ccd=ccd,
                    output=out,
                    kind="smrow",
                    cadence=cadence_kind,
                    smear_mode="standard",
                    tpf_parts=tpf_parts,
                    col_stem=sm_raw_path.stem if tpf_parts is None else None,
                )
                std_sm_blf = _calibrate_smrow_raw_to_standard_blf(
                    sm_raw_path,
                    tv_blf,
                    sm_dst,
                    output=out,
                    black2d_img=black2d_img,
                    lin_channel=lin_channel,
                    u_coeff=u_coeff,
                    sci_hdr0=hdr0,
                    sci_hdr1=hdr1,
                    cadence_match=cadence_match,
                    chunk_rows=chunk_rows,
                )

            if std_sm_blf is not None:
                with fits.open(std_sm_blf, memmap=True) as sm_hdul:
                    sm_hdu1 = sm_hdul[1]
                    sm_cad = np.asarray(sm_hdu1.data["CADENCENO"], dtype=np.int64)
                    sm_matched, sm_idx = match_cadenceno(
                        tpf_cad,
                        sm_cad,
                        cadence_match=cadence_match,
                        label=f"smrow ({out})",
                    )
                    combined_matched &= sm_matched
                    safe_sm_idx = np.where(sm_matched, sm_idx, 0)
                    sm_512 = np.asarray(sm_hdu1.data["SMEAR"], dtype=np.float64)[safe_sm_idx]
                    if "SMEAR_ERR" in (sm_hdu1.data.names or ()):
                        sm_err_512 = np.asarray(
                            sm_hdu1.data["SMEAR_ERR"], dtype=np.float64
                        )[safe_sm_idx]
                    else:
                        sm_err_512 = np.sqrt(
                            (np.pi / 20.0)
                            * (np.maximum(sm_512, 0.0) + readnoi_sci**2)
                            / n_eff_sci
                        )
                    sm_512 = np.where(sm_matched[:, None], sm_512, np.nan)
                    sm_err_512 = np.where(sm_matched[:, None], sm_err_512, np.nan)
                    max_sm_e = lin_channel.max_domain * gain_sci
                    sm_sat_512 = np.any(
                        np.isfinite(sm_512) & (sm_512 > 0.95 * max_sm_e), axis=0
                    )
            else:
                # Alternate mode without pre-existing standard smrow-blf: stream only matched
                # chunks from raw smrow without polluting standard smrow-blf path.
                sm_raw_path = resolve_collateral_file(
                    sector=sector,
                    camera=camera,
                    ccd=ccd,
                    output=out,
                    kind="smrow",
                    cadence=cadence_kind,
                    ddir=ddir,
                    cache_dir=resolved_cache,
                    offline=offline,
                    tpf_parts=tpf_parts,
                    s3_fs=s3_fs,
                    verbose=verbose,
                )
                sm_hdr0, sm_hdr1, _sm_time, sm_cad = read_collateral_headers_and_index(
                    sm_raw_path
                )
                sm_matched, sm_idx = match_cadenceno(
                    tpf_cad,
                    sm_cad,
                    cadence_match=cadence_match,
                    label=f"smrow ({out})",
                )
                combined_matched &= sm_matched
                n_eff_sm = effective_frame_count(
                    sm_hdr0, sm_hdr1, fallback_hdr0=hdr0, fallback_hdr1=hdr1
                )
                fxdoff_sm = float(sm_hdr1.get("FXDOFF", fxdoff_sci))
                meanblc_sm = float(sm_hdr1.get(f"MEANBLC{out}", meanblc_sci))
                gain_sm = float(sm_hdr1.get(f"GAIN{out}", gain_sci))
                readnoi_sm = float(sm_hdr1.get(f"READNOI{out}", readnoi_sci))
                sx1, sx2, _, _ = OUTPUT_BOUNDS_0BASED[out]
                b2d_sm = black2d_img[2058:2068, sx1:sx2]

                sm_512 = np.full((n_tpf, 512), np.nan, dtype=np.float64)
                sm_err_512 = np.full((n_tpf, 512), np.nan, dtype=np.float64)
                sm_sat_512 = np.zeros(512, dtype=bool)
                matched_tpf_rows = np.flatnonzero(sm_matched)
                needed_col_rows = sm_idx[matched_tpf_rows]

                for t0, t1, _tc, _cc, raw_chunk in iter_collateral_chunks(
                    sm_raw_path,
                    "SMROW_RAW",
                    chunk_rows=chunk_rows,
                    row_indices=needed_col_rows,
                ):
                    in_chunk = (sm_idx >= t0) & (sm_idx < t1) & sm_matched
                    if not np.any(in_chunk):
                        continue
                    tpf_sel = np.flatnonzero(in_chunk)
                    rel_sel = sm_idx[tpf_sel] - t0
                    sm_med_c, sm_err_c, sm_sat_c = _calibrate_smrow_chunk_e_per_frame(
                        raw_chunk[rel_sel],
                        output=out,
                        fxdoff=fxdoff_sm,
                        meanblc=meanblc_sm,
                        gain=gain_sm,
                        read_noise_e=readnoi_sm,
                        n_eff_sm=n_eff_sm,
                        n_eff_tv=n_eff_tv,
                        b2d_sm=b2d_sm,
                        b1d_sm=b1d_full[tpf_sel, 2058:2068],
                        b1d_err_sm=b1d_err_full[tpf_sel, 2058:2068],
                        u_coeff=u_coeff,
                        lin_channel=lin_channel,
                    )
                    sm_512[tpf_sel] = sm_med_c
                    sm_err_512[tpf_sel] = sm_err_c
                    sm_sat_512 |= sm_sat_c

            per_out_collat_512[out] = sm_512
            per_out_collat_err_512[out] = sm_err_512
            collat_smear_stamp[:, c0_loc:c1_loc] = sm_512[:, sm_col_start:sm_col_end]
            collat_smear_err_stamp[:, c0_loc:c1_loc] = sm_err_512[:, sm_col_start:sm_col_end]
            collat_sm_sat_cols[c0_loc:c1_loc] = sm_sat_512[sm_col_start:sm_col_end]

            # Off-target clean collateral smear baseline outside [sm_col_start-10, sm_col_end+10]
            cols_512 = np.arange(512)
            off_mask = (
                ((cols_512 < sm_col_start - 10) | (cols_512 >= sm_col_end + 10))
                & (~sm_sat_512)
            )
            if np.sum(off_mask) >= 8:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    off_med = np.nanmedian(sm_512[:, off_mask], axis=1)
                    off_mad = 1.4826 * np.nanmedian(
                        np.abs(sm_512[:, off_mask] - off_med[:, None]), axis=1
                    )
                off_target_smears.append(off_med)
                off_target_errs.append(np.maximum(off_mad, 0.25))

    # Enforce missing_cadence policy ("nan" warning vs "error" ValueError)
    n_missing = _enforce_missing_cadence_policy(
        tpf_cad,
        combined_matched,
        missing_cadence=missing_cadence,
    )

    if not np.all(combined_matched):
        sci_e_frame[~combined_matched] = np.nan
        sci_sigma_pre[~combined_matched] = np.nan
        sci_sat_3d[~combined_matched] = False

    # Construct spatial PIXEL_QUALITY mask
    raw_sat_2d = np.any(sci_sat_3d, axis=0)
    sat_bleed_2d = _build_sat_bleed_mask_2d(sci_sat_3d, col_outputs)
    pq_out = np.zeros((ny, nx), dtype=np.uint16)
    pq_out[:, strap_col_mask] |= np.uint16(PixelQuality.STRAP_COL)
    pq_out[sat_bleed_2d] |= np.uint16(PixelQuality.SAT_BLEED)
    pq_out[source_mask_2d] |= np.uint16(PixelQuality.SOURCE_MASK)

    # Compute final smear M_hat (e-/frame) and M_err (e-/frame) across the stamp
    if eff_smear == "standard":
        m_hat_stamp = collat_smear_stamp
        m_err_stamp = collat_smear_err_stamp
    elif all_outputs_have_alt_cache:
        m_hat_stamp = alt_cached_smear_stamp
        m_err_stamp = alt_cached_smear_err_stamp
        # Still evaluate SMEAR_UNRESOLVED flag from clean row count on bleeding columns
        dilated_bleed = _dilate_vertical_bool(raw_sat_2d, radius=4)
        clean_rows_2d = (~dilated_bleed) & (~source_mask_2d)
        for c in np.flatnonzero(np.any(raw_sat_2d, axis=0)):
            if int(np.sum(clean_rows_2d[:, c])) < 5:
                pq_out[:, c] |= np.uint16(PixelQuality.SMEAR_UNRESOLVED)
    else:
        if off_target_smears:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                off_sm_vec = np.nanmean(np.vstack(off_target_smears), axis=0)
                off_err_vec = np.nanmean(np.vstack(off_target_errs), axis=0)
        else:
            off_sm_vec = None
            off_err_vec = None

        m_hat_stamp, m_err_stamp, unresolved_cols = _estimate_alternate_smear_stamp(
            sci_pre_smear=sci_e_frame,
            sci_sigma_pre=sci_sigma_pre,
            flat_2d=flat_cutout,
            collat_smear_stamp=collat_smear_stamp,
            collat_smear_err_stamp=collat_smear_err_stamp,
            collat_sm_sat_cols=collat_sm_sat_cols,
            off_target_smear=off_sm_vec,
            off_target_smear_err=off_err_vec,
            raw_sat_2d=raw_sat_2d,
            sat_bleed_2d=sat_bleed_2d,
            strap_col_mask=strap_col_mask,
            source_mask_2d=source_mask_2d,
        )
        pq_out[:, unresolved_cols] |= np.uint16(PixelQuality.SMEAR_UNRESOLVED)

        # Persist target-scoped -alt-*-blf.fits caches for each output without cached alt
        for out in unique_outputs:
            if out not in per_out_collat_512:
                continue
            local_cols = np.flatnonzero(np.asarray([o == out for o in col_outputs], dtype=bool))
            c0_loc, c1_loc = int(local_cols[0]), int(local_cols[-1]) + 1
            sm_col_start = int(cols_1based[c0_loc] - 1 - OUTPUT_BOUNDS_0BASED[out][0])
            sm_col_end = sm_col_start + (c1_loc - c0_loc)

            alt_sm_512 = per_out_collat_512[out].copy()
            alt_err_512 = per_out_collat_err_512[out].copy()
            alt_sm_512[:, sm_col_start:sm_col_end] = m_hat_stamp[:, c0_loc:c1_loc]
            alt_err_512[:, sm_col_start:sm_col_end] = m_err_stamp[:, c0_loc:c1_loc]

            alt_dst = collateral_blf_cache_path(
                blf_write_base,
                sector=sector,
                camera=camera,
                ccd=ccd,
                output=out,
                kind="smrow",
                cadence=cadence_kind,
                smear_mode="alternate",
                tpf_parts=tpf_parts,
                col_stem=(
                    f"tess-s{sector:04d}-smrow-{camera}-{ccd}-{out.lower()}-"
                    f"{int(ticid):016d}"
                    if ticid is not None and tpf_parts is None
                    else None
                ),
            )
            write_smear_blf(
                alt_dst,
                time=tpf_time,
                cadenceno=tpf_cad,
                smear=alt_sm_512,
                smear_err=alt_err_512.astype(np.float32),
                smear_mode="alternate",
                cutout_bounds=cutout_bounds,
                int_time=int_time,
                time_unit=time_unit,
            )

    # Subtract smear, divide by flat field and single-frame exposure time INT_TIME -> e-/s
    safe_flat = np.where(np.isfinite(flat_cutout) & (flat_cutout > 0.0), flat_cutout, np.nan)
    denom_rate = safe_flat[np.newaxis, :, :] * int_time

    cal_flux = (sci_e_frame - m_hat_stamp[:, np.newaxis, :]) / denom_rate
    cal_flux_err = (
        np.sqrt(sci_sigma_pre * sci_sigma_pre + (m_err_stamp * m_err_stamp)[:, np.newaxis, :])
        / denom_rate
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        col_mean_flat = np.nanmean(safe_flat, axis=0)
    col_mean_flat = np.where(np.isfinite(col_mean_flat) & (col_mean_flat > 0.0), col_mean_flat, 1.0)
    cal_smear_err = m_err_stamp / (col_mean_flat[np.newaxis, :] * int_time)

    if not np.all(combined_matched):
        cal_flux[~combined_matched] = np.nan
        cal_flux_err[~combined_matched] = np.nan
        cal_smear_err[~combined_matched] = np.nan

    # Write calibrated columns & HALOCVER = 1 provenance across all N_raw rows of HDU 1
    col_names = hdul[1].data.names or ()
    hdul[1].data["FLUX"] = cal_flux.astype(np.float32)
    if "FLUX_ERR" in col_names:
        hdul[1].data["FLUX_ERR"] = cal_flux_err.astype(np.float32)
    if "FLUX_BKG" in col_names:
        hdul[1].data["FLUX_BKG"] = np.float32(0.0)
    if "FLUX_BKG_ERR" in col_names:
        hdul[1].data["FLUX_BKG_ERR"] = np.float32(0.0)

    append_smear_err_column(hdul, cal_smear_err.astype(np.float32))
    write_pixel_quality(hdul, pq_out)
    write_mundey_provenance(target_obj if target_obj is not None else hdul, bkgsub=False)

    hdul[0].header["CADMATCH"] = (cadence_match, "Collateral cadence pairing mode")
    hdul[0].header["CADMISS"] = (int(n_missing), "Unmatched TPF cadences masked to NaN")
    if target_obj is not None:
        meta = getattr(target_obj, "meta", None)
        meta_dict = getattr(meta, "_dict", None)
        if isinstance(meta_dict, dict):
            meta_dict["CADMATCH"] = cadence_match
            meta_dict["CADMISS"] = int(n_missing)


def _clone_input_tpf(
    tpf: str | os.PathLike[str] | lk.TargetPixelFile,
) -> tuple[fits.HDUList, Path | None, Any]:
    """Return a deep-copied ``HDUList``, optional source path, and ``quality_bitmask``."""
    if isinstance(tpf, (str, os.PathLike)):
        src_path = Path(tpf).expanduser()
        with fits.open(src_path, mode="readonly", memmap=False) as src_hdul:
            hdul = fits.HDUList([h.copy() for h in src_hdul])
        return hdul, src_path, "default"

    if hasattr(tpf, "hdu") and isinstance(tpf.hdu, fits.HDUList):
        hdul = fits.HDUList([h.copy() for h in tpf.hdu])
        raw_path = getattr(tpf, "path", None)
        src_path = Path(raw_path).expanduser() if raw_path else None
        qmask = getattr(tpf, "quality_bitmask", "default") or "default"
        return hdul, src_path, qmask

    raise TypeError(
        f"Expected file path or lightkurve.TargetPixelFile, got {type(tpf)!r}."
    )


def recalibrate_tpf(
    tpf: str | os.PathLike[str] | lk.TargetPixelFile,
    *,
    cache_dir: str | os.PathLike[str] | None = None,
    ddir: str | os.PathLike[str] | None = None,
    smear: SmearMode = "automatic",
    cadence_match: CadenceMatchMode = "cadenceno",
    missing_cadence: MissingCadenceMode = "nan",
    offline: bool = False,
    strict_hpc: bool = False,
    verbose: bool = False,
) -> lk.TargetPixelFile:
    """Recalibrate a TESS Target Pixel File from ``RAW_CNTS`` under ``HALOCVER = 1``.

    Pure functional transformation: never mutates an in-memory ``TargetPixelFile``
    argument. Returns a calibrated ``TessTargetPixelFile`` with un-background-subtracted
    ``FLUX`` (``e-/s``), total marginal ``FLUX_ERR``, column-coherent ``SMEAR_ERR`` in
    ``hdu[1]``, zeroed ``FLUX_BKG`` / ``FLUX_BKG_ERR``, spatial ``PIXEL_QUALITY``
    ``ImageHDU``, and ``HALOCVER = 1``, ``BKGSUB = False`` provenance cards.
    """
    hdul, src_path, quality_bitmask = _clone_input_tpf(tpf)
    _recalibrate_hdulist_inplace(
        hdul,
        source_path=src_path,
        cache_dir=cache_dir,
        ddir=ddir,
        smear=smear,
        cadence_match=cadence_match,
        missing_cadence=missing_cadence,
        offline=offline,
        strict_hpc=strict_hpc,
        verbose=verbose,
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
    """Repair collateral-bleed smear corruption on calibrated FFI cutouts (e.g. ``TessCut``).

    Designed for full-frame-image cutouts where raw ``RAW_CNTS`` and collateral
    ``_col.fits`` tables are unavailable. Restores any previously subtracted
    ``FLUX_BKG`` when ``BKGSUB`` / ``BACKAPP`` is ``True``, bridges clean flank
    columns across bloomed/saturated columns with inverse-wing-weighted Huber
    row-location estimation, writes ``SMEAR_ERR`` and ``PIXEL_QUALITY``, zeroes
    ``FLUX_BKG`` / ``FLUX_BKG_ERR``, and sets ``HALOCVER = 1``, ``BKGSUB = False``.
    """
    hdul, src_path, quality_bitmask = _clone_input_tpf(tpf)
    hdr0 = hdul[0].header
    hdu1 = hdul[1]
    hdr1 = hdu1.header
    hdr2 = hdul[2].header

    names = hdu1.data.names or ()
    if "FLUX" not in names:
        raise ValueError("Cutout HDU[1] lacks a 'FLUX' column.")

    flux_in = np.asarray(hdu1.data["FLUX"], dtype=np.float64)
    n_tpf, ny, nx = flux_in.shape

    # Restore un-background-subtracted flux if cutout had background subtracted
    bkgsub_hdr = hdr0.get("BKGSUB", hdr1.get("BKGSUB", None))
    backapp_hdr = bool(hdr1.get("BACKAPP", hdr0.get("BACKAPP", False)))
    has_bkg_col = "FLUX_BKG" in names
    if has_bkg_col:
        bkg_arr = np.asarray(hdu1.data["FLUX_BKG"], dtype=np.float64)
        bkg_nonzero = bool(np.any(np.isfinite(bkg_arr) & (bkg_arr != 0.0)))
    else:
        bkg_arr = np.zeros_like(flux_in)
        bkg_nonzero = False

    should_restore_bkg = (
        (bkgsub_hdr is True) or (bkgsub_hdr is None and backapp_hdr and bkg_nonzero)
    )
    if should_restore_bkg:
        u_in = flux_in + np.where(np.isfinite(bkg_arr), bkg_arr, 0.0)
    else:
        u_in = flux_in.copy()

    if "FLUX_ERR" in names:
        err_in = np.asarray(hdu1.data["FLUX_ERR"], dtype=np.float64)
        if should_restore_bkg and "FLUX_BKG_ERR" in names:
            berr = np.asarray(hdu1.data["FLUX_BKG_ERR"], dtype=np.float64)
            berr = np.where(np.isfinite(berr), berr, 0.0)
            sigma_in = np.sqrt(np.maximum(err_in * err_in - berr * berr, 1e-6))
        else:
            sigma_in = np.where(np.isfinite(err_in) & (err_in > 0.0), err_in, 1.0)
    else:
        sigma_in = np.sqrt(np.maximum(np.abs(u_in), 1.0) / 360.0)

    col0_1based = int(hdr2.get("CRVAL1P", 45))
    cols_1based = np.arange(col0_1based, col0_1based + nx, dtype=np.int64)
    in_sci_range = (cols_1based >= 45) & (cols_1based <= 2092)
    if np.all(in_sci_range):
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

    sci_sat_3d = np.isfinite(u_in) & (u_in > float(saturation_flux_threshold))
    raw_sat_2d = np.any(sci_sat_3d, axis=0) | ((existing_pq & int(PixelQuality.SAT_BLEED)) != 0)
    sat_bleed_2d = _build_sat_bleed_mask_2d(sci_sat_3d, col_outputs) | raw_sat_2d

    if verbose:
        print(f"Desmearing calibrated cutout of shape {(n_tpf, ny, nx)}")

    delta_m, delta_m_err, unresolved_cols = _estimate_alternate_smear_stamp(
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
    )

    pq_out = np.zeros((ny, nx), dtype=np.uint16)
    pq_out[:, strap_col_mask] |= np.uint16(PixelQuality.STRAP_COL)
    pq_out[sat_bleed_2d] |= np.uint16(PixelQuality.SAT_BLEED)
    pq_out[:, unresolved_cols] |= np.uint16(PixelQuality.SMEAR_UNRESOLVED)
    pq_out[source_mask_2d] |= np.uint16(PixelQuality.SOURCE_MASK)

    flux_out = u_in - delta_m[:, np.newaxis, :]
    flux_err_out = np.sqrt(sigma_in * sigma_in + (delta_m_err * delta_m_err)[:, np.newaxis, :])

    finite_cad = np.any(np.isfinite(u_in), axis=(1, 2))
    if not np.all(finite_cad):
        flux_out[~finite_cad] = np.nan
        flux_err_out[~finite_cad] = np.nan
        delta_m_err[~finite_cad] = np.nan

    hdul[1].data["FLUX"] = flux_out.astype(np.float32)
    if "FLUX_ERR" in names:
        hdul[1].data["FLUX_ERR"] = flux_err_out.astype(np.float32)
    if "FLUX_BKG" in names:
        hdul[1].data["FLUX_BKG"] = np.float32(0.0)
    if "FLUX_BKG_ERR" in names:
        hdul[1].data["FLUX_BKG_ERR"] = np.float32(0.0)

    append_smear_err_column(hdul, delta_m_err.astype(np.float32))
    write_pixel_quality(hdul, pq_out)
    write_mundey_provenance(hdul, bkgsub=False)

    out_tpf = lk.TessTargetPixelFile(hdul, quality_bitmask=quality_bitmask)
    if src_path is not None and hasattr(out_tpf, "path"):
        out_tpf.path = str(src_path)
    return out_tpf
