"""Concavity-guarded alternate smear estimation, Huber location solver, and bloom detection."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
from scipy.special import ndtri

from ._precision import in_float64

TESS_PARALLEL_SHIFTS: int = 2058  # 2048 imaging rows + 10 lower buffer rows into frame-store
_HUBER_ROW_BUCKETS: tuple[int, ...] = (64, 256, 2048)

# Exact variance of the sample median of n i.i.d. N(0, 1) variates for 1 <= n <= 20.
# See https://en.wikipedia.org/wiki/Median#Sampling_distribution
_UNIT_NORMAL_MEDIAN_VAR_1_TO_20: np.ndarray = np.array(
    [
        1.0,
        1.000000000000,
        0.500000000000,
        0.448671104578,
        0.298199618435,
        0.286833661606,
        0.214742666582,
        0.210446861535,
        0.168180855159,
        0.166101281359,
        0.138326435836,
        0.137162433548,
        0.117516187119,
        0.116798995001,
        0.102167836796,
        0.101694652082,
        0.090375211225,
        0.090046581423,
        0.081028510053,
        0.080790975073,
        0.073437027705,
    ],
    dtype=np.float64,
)


@dataclass(frozen=True)
class AlternateSmearConfig:
    """Parameters controlling concavity-guarded alternate smear estimation."""

    row0_1based: int = 1
    row_smear_weight: float = 0.02 / (1.98 * float(TESS_PARALLEL_SHIFTS))
    saturation_limit: float | None = None
    bleed_dilation_rows: int = 4
    min_clean_rows: int = 5
    flank_width: int = 3


def _unit_normal_median_variance(n: np.ndarray | int) -> np.ndarray:
    """Return exact unit-normal sample-median variance for ``1 <= n <= 20`` and ``pi / (2n)``."""
    n_arr = np.maximum(np.asarray(n, dtype=np.int64), 1)
    exact = _UNIT_NORMAL_MEDIAN_VAR_1_TO_20[np.minimum(n_arr, 20)]
    asymptotic = np.pi / (2.0 * n_arr.astype(np.float64))
    return np.where(n_arr <= 20, exact, asymptotic)


def _select_bucket(size: int, buckets: tuple[int, ...]) -> int:
    for b in buckets:
        if size <= b:
            return b
    return size


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

    w_max = jnp.max(w_base, axis=1, keepdims=True)
    far_mask = valid & (w_base >= 0.05 * w_max)
    scale_mask = jnp.where(jnp.sum(far_mask, axis=1, keepdims=True) >= 5, far_mask, valid)

    # Exact weighted-median initialization over scale_mask (immune to extreme outliers)
    init_mask = scale_mask & (w_base > 0.0)
    order = jnp.argsort(jnp.where(init_mask, residuals, jnp.inf), axis=1)
    r_sorted = jnp.take_along_axis(residuals, order, axis=1)
    w_sorted = jnp.take_along_axis(jnp.where(init_mask, w_base, 0.0), order, axis=1)
    w_cum = jnp.cumsum(w_sorted, axis=1)
    w_tot = w_cum[:, -1:]
    half_w = 0.5 * w_tot
    idx_lo = jnp.argmax(w_cum >= half_w, axis=1)
    idx_hi = jnp.where(jnp.any(w_cum > half_w, axis=1), jnp.argmax(w_cum > half_w, axis=1), idx_lo)
    loc_lo = jnp.take_along_axis(r_sorted, idx_lo[:, None], axis=1)[:, 0]
    loc_hi = jnp.take_along_axis(r_sorted, idx_hi[:, None], axis=1)[:, 0]
    loc = jnp.where(w_tot[:, 0] > 0.0, 0.5 * (loc_lo + loc_hi), 0.0)

    sigma_eff = jnp.where(w_base > 0.0, jax.lax.rsqrt(jnp.maximum(w_base, 1e-30)), 1.0)
    c_huber = jnp.asarray(1.345, dtype=residuals.dtype)

    r0 = jnp.where(valid, residuals - loc[:, None], 0.0)
    std_r0 = jnp.where(scale_mask, jnp.abs(r0) / jnp.maximum(sigma_eff, 1e-12), jnp.nan)
    tau0 = jnp.maximum(1.0, 1.4826 * jnp.nanmedian(std_r0, axis=1))
    tau0 = jnp.where(jnp.isfinite(tau0), tau0, 1.0)
    inv_denom_scale = 1.0 / jnp.maximum(tau0[:, None] * sigma_eff, 1e-12)
    res_valid = jnp.where(valid, residuals, 0.0)

    def _irls_step(_idx: int, m_curr: jax.Array) -> jax.Array:
        r = jnp.where(valid, residuals - m_curr[:, None], 0.0)
        u = jnp.abs(r) * inv_denom_scale
        w_rob = jnp.where(u <= c_huber, 1.0, c_huber / jnp.maximum(u, 1e-12))
        w_step = w_base * w_rob
        denom = jnp.sum(w_step, axis=1)
        m_next = jnp.sum(w_step * res_valid, axis=1) / jnp.maximum(denom, 1e-30)
        return jnp.where(denom > 0.0, m_next, m_curr)

    loc_final = jax.lax.fori_loop(0, 12, _irls_step, loc)

    # M-estimator sandwich variance: Hessian bread sum(w_base * I(|u| <= c))
    r_fin = jnp.where(valid, residuals - loc_final[:, None], 0.0)
    u_fin = jnp.abs(r_fin) * inv_denom_scale
    in_quad = valid & (u_fin <= c_huber)
    w_rob_fin = jnp.where(u_fin <= c_huber, 1.0, c_huber / jnp.maximum(u_fin, 1e-12))
    w_tot_fin = w_base * w_rob_fin

    hess_fin = jnp.sum(jnp.where(in_quad, w_base, 0.0), axis=1)
    denom_tot = jnp.sum(w_tot_fin, axis=1)
    denom_fin = jnp.where(
        jnp.sum(in_quad, axis=1) >= 2,
        hess_fin,
        jnp.maximum(hess_fin, 0.25 * denom_tot),
    )

    emp_meat = jnp.sum((w_tot_fin * r_fin) ** 2, axis=1)
    theo_meat = jnp.sum((w_tot_fin * jnp.where(valid, sigma_pix, 0.0)) ** 2, axis=1)
    loc_err = jnp.sqrt(jnp.maximum(jnp.maximum(emp_meat, theo_meat), 0.0)) / jnp.maximum(
        denom_fin, 1e-30
    )

    has_any = jnp.any(valid, axis=1)
    return jnp.where(has_any, loc_final, jnp.nan), jnp.where(has_any, loc_err, jnp.nan)


@in_float64
def _solve_huber_static_chunk(
    res_buf: Any,
    w_buf: Any,
    sig_buf: Any,
    mask_buf: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    loc, err = _huber_weighted_smear_location_2d_jit(
        res_buf, w_buf, sig_buf, jnp.asarray(mask_buf, dtype=bool)
    )
    return np.asarray(loc, dtype=np.float64), np.asarray(err, dtype=np.float64)


def _solve_huber_smear_location(
    residuals: Any,
    inv_wing_weights: Any,
    sigma_pix: Any,
    clean_row_mask: Any,
) -> tuple[np.ndarray, np.ndarray]:
    """Run shape-bucketed ``@in_float64`` MAD-scaled Huber smear solver."""
    res_np = np.asarray(residuals, dtype=np.float64)
    w_np = np.asarray(inv_wing_weights, dtype=np.float64)
    sig_np = np.asarray(sigma_pix, dtype=np.float64)
    mask_np = np.asarray(clean_row_mask, dtype=bool)

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
            out |= np.broadcast_to(np.any(mask_2d, axis=0, keepdims=True), (ny, nx))
            break
    return out


def _build_clean_rows_2d(
    raw_sat_2d: np.ndarray,
    sat_bleed_2d: np.ndarray,
    source_mask_2d: np.ndarray,
    *,
    bleed_dilation_rows: int = 4,
    min_clean_rows: int = 5,
) -> np.ndarray:
    """Exclude dilated bleed, ``sat_bleed_2d``, and ``±1``-col lateral bleed neighbours.

    Falls back per column to single-column dilated bleed exclusion only if stricter
    lateral/downstream bleed masking leaves fewer than ``min_clean_rows`` clean rows.
    """
    _, nx = raw_sat_2d.shape
    dilated = _dilate_vertical_bool(raw_sat_2d, radius=bleed_dilation_rows)
    base_clean = (~dilated) & (~source_mask_2d)

    lat_bleed = np.zeros_like(raw_sat_2d, dtype=bool)
    if nx > 1:
        lat_bleed[:, 1:] |= raw_sat_2d[:, :-1]
        lat_bleed[:, :-1] |= raw_sat_2d[:, 1:]

    strict_clean = base_clean & (~sat_bleed_2d) & (~lat_bleed)
    clean_rows = strict_clean.copy()
    for c in range(nx):
        if int(np.sum(clean_rows[:, c])) < min_clean_rows:
            clean_rows[:, c] = base_clean[:, c]
    return clean_rows


def _eval_bloom_from_stats(
    w_num: np.ndarray,
    w_sum: np.ndarray,
    n_pix_2d: np.ndarray,
    min_z_2d: np.ndarray,
    clean_rows_2d: np.ndarray,
    *,
    fwer_alpha: float = 1.0e-4,
    coherent_spatial_err: float = 0.0,
) -> np.ndarray:
    """Evaluate exact Bonferroni persistent + transient bloom tests from 2D sufficient stats."""
    _, nx = clean_rows_2d.shape
    col_has_bloom = np.zeros(nx, dtype=bool)
    coh_var = max(float(coherent_spatial_err), 0.0) ** 2

    for c in range(nx):
        r_mask = clean_rows_2d[:, c]
        if not np.any(r_mask):
            col_has_bloom[c] = True
            continue

        ws = w_sum[r_mask, c]
        valid_r = ws > 0.0
        n_r_tests = int(np.sum(valid_r))
        if n_r_tests > 0:
            mean_r = w_num[r_mask, c][valid_r] / ws[valid_r]
            se_r = np.sqrt(1.0 / ws[valid_r] + coh_var)
            p_row = min(max(0.5 * fwer_alpha / float(n_r_tests), 1.0e-300), 0.49)
            if float(np.min(mean_r / se_r)) < float(ndtri(p_row)):
                col_has_bloom[c] = True
                continue

        n_pix_tests = int(np.sum(n_pix_2d[r_mask, c]))
        if n_pix_tests > 0:
            p_pix = min(max(0.5 * fwer_alpha / float(n_pix_tests), 1.0e-300), 0.49)
            if float(np.min(min_z_2d[r_mask, c])) < float(ndtri(p_pix)):
                col_has_bloom[c] = True

    return col_has_bloom


def _detect_bloom_columns(
    sub_cube: np.ndarray,
    sigma_cube: np.ndarray,
    clean_rows_2d: np.ndarray,
    *,
    fwer_alpha: float = 1.0e-4,
    coherent_spatial_err: float = 0.0,
) -> np.ndarray:
    """Detect columns whose clean rows are driven negative by upper-register smear bloom."""
    _, ny, nx = sub_cube.shape
    coh_err = max(float(coherent_spatial_err), 0.0)
    valid = np.isfinite(sub_cube) & np.isfinite(sigma_cube)
    sig_safe = np.maximum(np.where(valid, sigma_cube, 1.0), 1e-6)
    w = np.where(valid, 1.0 / (sig_safe * sig_safe), 0.0)
    w_sum = np.sum(w, axis=0)
    w_num = np.sum(np.where(valid, w * sub_cube, 0.0), axis=0)
    n_pix_2d = np.sum(valid, axis=0)
    sig_pix = np.hypot(sig_safe, coh_err) if coh_err > 0.0 else sig_safe
    z_cube = np.where(valid, sub_cube / sig_pix, np.inf)
    min_z_2d = np.min(z_cube, axis=0) if z_cube.shape[0] > 0 else np.full((ny, nx), np.inf)
    return _eval_bloom_from_stats(
        w_num,
        w_sum,
        n_pix_2d,
        min_z_2d,
        clean_rows_2d,
        fwer_alpha=fwer_alpha,
        coherent_spatial_err=coh_err,
    )


def _fit_upper_column_floor_line_vec(
    x_1d: np.ndarray,
    col_floors: np.ndarray,
) -> np.ndarray:
    """Vectorized Theil-Sen upper-envelope linear trend across per-column floors ``(m, nx)``."""
    m, nx = col_floors.shape
    pred_x = np.zeros((m, nx), dtype=np.float64)
    if nx == 0 or m == 0:
        return pred_x

    sort_vals = np.where(np.isfinite(col_floors), col_floors, -np.inf)
    order = np.argsort(sort_vals, axis=1)
    n_fin = np.sum(np.isfinite(col_floors), axis=1)

    for t in np.flatnonzero(n_fin > 0):
        nf = int(n_fin[t])
        if nf == 1:
            pred_x[t] = max(float(col_floors[t, order[t, -1]]), 0.0)
            continue
        k = min(nf, max(3, int(np.ceil(0.45 * nf))))
        top_idx = order[t, nx - k :]
        xs, ys = x_1d[top_idx], col_floors[t, top_idx]
        iu, ju = np.triu_indices(k, k=1)
        dx = xs[ju] - xs[iu]
        ok = np.abs(dx) > 1e-12
        slope = float(np.median((ys[ju[ok]] - ys[iu[ok]]) / dx[ok])) if np.any(ok) else 0.0
        intercept = float(np.median(ys - slope * xs))
        pred_x[t] = np.maximum(intercept + slope * x_1d, 0.0)
    return pred_x


def _estimate_cutout_sky_plane(
    u_in: np.ndarray,
    sigma_in: np.ndarray,
    clean_rows_2d: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Estimate a non-negative 2D sky baseline ``(N_t, ny, nx)`` and coherent spatial error."""
    n_t, ny, nx = u_in.shape
    npix = ny * nx
    if npix == 0 or n_t == 0 or not np.any(clean_rows_2d):
        return np.zeros_like(u_in, dtype=np.float64), 0.0

    x_1d = np.linspace(-1.0, 1.0, nx, dtype=np.float64) if nx > 1 else np.zeros(1)
    y_1d = np.linspace(-1.0, 1.0, ny, dtype=np.float64) if ny > 1 else np.zeros(1)
    yy, xx = np.meshgrid(y_1d, x_1d, indexing="ij")
    xr, yr = xx.ravel(), yy.ravel()

    use_quad = (nx >= 4) and (ny >= 4)
    design_lin = np.column_stack([np.ones(npix), xr, yr])
    design_quad = np.column_stack([np.ones(npix), xr, yr, xr * xr, xr * yr, yr * yr])
    ridge_lin = np.diag([1e-12, 1e-6, 1e-6])
    ridge_quad = np.diag([1e-12, 1e-6, 1e-6, 1e-5, 1e-5, 1e-5])
    clean_flat = clean_rows_2d.ravel()

    out = np.zeros((n_t, npix), dtype=np.float64)
    inlier_count = np.zeros(npix, dtype=np.int64)
    for t0 in range(0, n_t, 256):
        t1 = min(t0 + 256, n_t)
        m = t1 - t0
        u_chunk = u_in[t0:t1]
        u_flat = u_chunk.reshape(m, npix)
        sig_c = sigma_in[t0:t1]
        sig_flat = np.maximum(np.where(np.isfinite(sig_c), sig_c, 1.0), 1e-6).reshape(m, npix)
        valid = clean_flat[None, :] & np.isfinite(u_flat)
        u_safe = np.where(valid, u_flat, 0.0)

        col_floors = np.full((m, nx), np.nan, dtype=np.float64)
        vert_mads = np.full((m, nx), np.nan, dtype=np.float64)
        for c in range(nx):
            r_ok = clean_rows_2d[:, c]
            if not np.any(r_ok):
                continue
            uc = u_chunk[:, r_ok, c]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                p15 = np.nanpercentile(uc, 15.0, axis=1)
                p30 = np.nanpercentile(uc, 30.0, axis=1)
                low_m = np.isfinite(uc) & (uc <= p30[:, None])
                low_med = np.nanmedian(np.where(low_m, uc, np.nan), axis=1)
                vmad = 1.4826 * np.nanmedian(
                    np.where(low_m, np.abs(uc - low_med[:, None]), np.nan), axis=1
                )
            col_floors[:, c] = p15
            vert_mads[:, c] = vmad

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            cad_vert_mad = np.nanmedian(vert_mads, axis=1)
        cad_vert_mad = np.where(np.isfinite(cad_vert_mad), cad_vert_mad, 0.0)

        pred_x = _fit_upper_column_floor_line_vec(x_1d, col_floors)
        has_valid = np.any(valid, axis=1) & np.any(np.isfinite(col_floors), axis=1)
        pred = np.broadcast_to(pred_x[:, None, :], (m, ny, nx)).reshape(m, npix).copy()

        final_w_rob = np.zeros((m, npix), dtype=np.float64)
        for it in range(5):
            res = u_safe - pred
            if it == 0:
                scale = np.maximum(np.maximum(sig_flat, cad_vert_mad[:, None]), 1e-6)
            else:
                sky_cand = valid & (res >= -4.0 * sig_flat)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    mad = 1.4826 * np.nanmedian(np.where(sky_cand, np.abs(res), np.nan), axis=1)
                mad = np.where(np.isfinite(mad), mad, 0.0)
                mad_cap = np.minimum(mad[:, None], 3.0 * sig_flat)
                scale = np.maximum(np.maximum(sig_flat, mad_cap), 1e-6)

            z = res / scale
            if it < 2:
                w_pos = np.where(z <= 2.5, 1.0, 2.5 / np.maximum(z, 1e-6))
            else:
                w_pos = np.where(z <= 2.0, 1.0, np.where(z <= 3.5, 2.0 / np.maximum(z, 1e-6), 0.0))
            w_neg = np.where(z >= -2.0, 1.0, np.where(z >= -4.0, 2.0 / np.maximum(-z, 1e-6), 0.0))
            final_w_rob = np.where(z >= 0.0, w_pos, w_neg)
            w = np.where(valid, final_w_rob / (sig_flat * sig_flat), 0.0)

            too_few = has_valid & (np.sum(w, axis=1) <= 0.0)
            if np.any(too_few):
                w[too_few] = np.where(
                    valid[too_few],
                    1.0 / (np.maximum(np.abs(z[too_few]), 1.0) * (sig_flat[too_few] ** 2)),
                    0.0,
                )

            inlier_cols = np.any((w > 0.0).reshape(m, ny, nx), axis=1)
            n_inlier_cols = (
                int(np.min(np.sum(inlier_cols[has_valid], axis=1))) if np.any(has_valid) else 0
            )
            if use_quad and n_inlier_cols >= 4:
                des, rdg = design_quad, ridge_quad
            else:
                des, rdg = design_lin, ridge_lin

            atw = des.T[None, :, :] * w[:, None, :]
            ata = np.einsum("tip,pj->tij", atw, des) + rdg[None, :, :]
            atu = np.einsum("tip,tp->ti", atw, u_safe)
            beta = np.linalg.solve(ata, atu[..., None])[..., 0]
            pred = np.where(has_valid[:, None], beta @ des.T, 0.0)

        out[t0:t1] = np.maximum(pred, 0.0)
        inlier_count += np.sum(valid & (final_w_rob > 0.5), axis=0)

    sky_cube = out.reshape(n_t, ny, nx)
    inlier_2d = (inlier_count.reshape(ny, nx) >= max(1, n_t // 2)) & clean_rows_2d
    coherent_spatial_err = 0.0
    if int(np.sum(inlier_2d)) >= 6:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            vals = np.nanmean(u_in - sky_cube, axis=0)[inlier_2d]
        vals = vals[np.isfinite(vals)]
        if vals.size >= 6:
            med_v = float(np.median(vals))
            coherent_spatial_err = float(1.4826 * np.median(np.abs(vals - med_v)))
    return sky_cube, coherent_spatial_err


def _flank_gradient_4col(
    ok_side: np.ndarray,
    xs_side: np.ndarray,
    u_side: np.ndarray,
    rows: np.ndarray,
    flat_col: np.ndarray,
    j_col: np.ndarray,
) -> np.ndarray:
    """Vectorized 4-column horizontal OLS gradient on one-sidedly clamped anchor rows."""
    ok4 = ok_side[rows]
    sel4 = ok4 & (np.cumsum(ok4, axis=1) <= 4)
    cnt4 = np.sum(sel4, axis=1)
    xm = np.sum(np.where(sel4, xs_side[None, :], 0.0), axis=1) / np.maximum(cnt4, 1)
    dx = np.where(sel4, xs_side[None, :] - xm[:, None], 0.0)
    den = np.sum(dx * dx, axis=1)
    slp = np.sum(
        np.where(sel4[None, :, :], u_side[:, rows, :], 0.0) * dx[None, :, :],
        axis=2,
    ) / np.where(den[None, :] > 0.0, den[None, :], 1.0)
    has2 = (cnt4 >= 2) & (den > 0.0)
    return np.where(
        has2[None, :],
        np.abs(slp) * flat_col[rows][None, :],
        0.05 * np.abs(j_col[:, rows]),
    )


def _select_flank_side(
    ok_side: np.ndarray,
    xs_side: np.ndarray,
    u_side: np.ndarray,
    indep_var_side: np.ndarray,
    shared_side: np.ndarray,
    flank_width: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Select up to ``flank_width`` nearest valid anchors on one side of a bleed column."""
    n_t, ny, n_side = u_side.shape
    if n_side == 0:
        return (
            np.zeros(ny, dtype=bool),
            np.zeros(ny, dtype=np.float64),
            np.zeros((n_t, ny), dtype=np.float64),
            np.zeros((n_t, ny), dtype=np.float64),
            np.zeros(ny, dtype=np.float64),
        )
    sel = ok_side & (np.cumsum(ok_side, axis=1) <= flank_width)
    n_s = np.sum(sel, axis=1)
    n_safe = np.maximum(n_s, 1)
    c_pos = np.sum(np.where(sel, xs_side[None, :], 0.0), axis=1) / n_safe
    u_mean = np.sum(np.where(sel[None, :, :], u_side, 0.0), axis=2) / n_safe[None, :]
    v_indep = np.sum(
        np.where(sel[None, :, :], indep_var_side[:, None, :], 0.0), axis=2
    ) / (n_safe[None, :] ** 2)
    s_frac = np.sum(np.where(sel, shared_side[None, :], 0.0), axis=1) / n_safe
    return n_s > 0, c_pos, u_mean, v_indep, s_frac


def _anchor_all_bloomed_columns(
    *,
    sci_pre_smear: np.ndarray,
    sci_sigma_pre: np.ndarray,
    safe_flat: np.ndarray,
    off_target_smear: np.ndarray | None,
    off_target_smear_err: np.ndarray | None,
    raw_sat_2d: np.ndarray,
    sat_bleed_2d: np.ndarray,
    strap_col_mask: np.ndarray,
    source_mask_2d: np.ndarray,
    clean_rows_2d: np.ndarray,
    col_has_sci_sat: np.ndarray,
    c_bleed_mask: np.ndarray,
    m_hat: np.ndarray,
    m_err: np.ndarray,
    indep_m_var: np.ndarray,
    shared_off_var: np.ndarray,
    col_uses_shared: np.ndarray,
    config: AlternateSmearConfig,
    bloom_detect_cube: np.ndarray | None,
) -> bool:
    """Anchor unsaturated edge columns when every stamp column bloomed."""
    n_t, ny, nx = sci_pre_smear.shape
    if np.any((~c_bleed_mask) & (~strap_col_mask)):
        return False

    if off_target_smear is not None and np.any(~col_has_sci_sat):
        edge_candidates = (~col_has_sci_sat) & (~strap_col_mask)
        if not np.any(edge_candidates):
            edge_candidates = ~col_has_sci_sat
        edge_cols = np.flatnonzero(edge_candidates)

        u0_edge = (
            sci_pre_smear[:, :, edge_cols] - off_target_smear[:, None, None]
        ) / safe_flat[None, :, edge_cols]
        clean_edge = (~sat_bleed_2d[:, edge_cols]) & (~source_mask_2d[:, edge_cols])
        if not np.any(clean_edge):
            clean_edge = ~raw_sat_2d[:, edge_cols]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            sky_floor = np.nanpercentile(
                np.where(clean_edge[None, :, :], u0_edge, np.nan).reshape(n_t, -1), 10.0, axis=1
            )
        sky_floor = np.where(np.isfinite(sky_floor), np.maximum(sky_floor, 0.0), 0.0)

        n_below = max(int(config.row0_1based) - 1, 0)
        n_above = max(2048 - (int(config.row0_1based) - 1 + ny), 0)
        n_edge_r = min(3, ny)
        var_off = (
            off_target_smear_err * off_target_smear_err
            if off_target_smear_err is not None
            else np.zeros(n_t, dtype=np.float64)
        )
        shared_off_var[:] = var_off

        for idx_e, c in enumerate(edge_cols):
            u0_c = u0_edge[:, :, idx_e]
            star_ex = np.where(np.isfinite(u0_c), np.maximum(u0_c - sky_floor[:, None], 0.0), 0.0)
            delta_m_in = float(config.row_smear_weight) * np.sum(star_ex, axis=1)
            sig_u_c = np.where(
                np.isfinite(sci_sigma_pre[:, :, c]),
                sci_sigma_pre[:, :, c] / safe_flat[None, :, c],
                0.0,
            )
            var_in = (float(config.row_smear_weight) ** 2) * np.sum(sig_u_c * sig_u_c, axis=1)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                edge_below = np.nanmedian(star_ex[:, :n_edge_r], axis=1)
                edge_above = np.nanmedian(star_ex[:, ny - n_edge_r :], axis=1)
            sigma_out = float(config.row_smear_weight) * (
                float(n_below) * np.where(np.isfinite(edge_below), edge_below, 0.0)
                + float(n_above) * np.where(np.isfinite(edge_above), edge_above, 0.0)
            )
            col_indep_var = var_in + sigma_out * sigma_out
            m_hat[:, c] = off_target_smear + delta_m_in
            indep_m_var[:, c] = col_indep_var
            m_err[:, c] = np.sqrt(np.maximum(var_off + col_indep_var, 1e-12))
            col_uses_shared[c] = True
            c_bleed_mask[c] = False
        return False

    if off_target_smear is not None:
        return False

    resolvable_any = np.sum(clean_rows_2d, axis=0) >= config.min_clean_rows
    edge_candidates = (~col_has_sci_sat) & (~strap_col_mask) & resolvable_any
    if not np.any(edge_candidates):
        edge_candidates = (~col_has_sci_sat) & resolvable_any
    if not np.any(edge_candidates) and np.any(resolvable_any):
        cand_idx = np.flatnonzero((~strap_col_mask) & resolvable_any)
        if cand_idx.size == 0:
            cand_idx = np.flatnonzero(resolvable_any)
        edge_candidates = np.zeros(nx, dtype=bool)
        edge_candidates[cand_idx[0]] = True
        edge_candidates[cand_idx[-1]] = True

    edge_cols = np.flatnonzero(edge_candidates)
    if edge_cols.size == 0:
        return False

    floor_src = bloom_detect_cube if bloom_detect_cube is not None else sci_pre_smear
    z_q10 = -float(ndtri(0.10))
    clamped_without_bkg = False
    for c in edge_cols:
        r_ok = clean_rows_2d[:, c]
        u_c_clean = (floor_src[:, :, c] / safe_flat[None, :, c])[:, r_ok]
        sig_dim = np.where(
            np.isfinite(sci_sigma_pre[:, :, c]),
            sci_sigma_pre[:, :, c] / safe_flat[None, :, c],
            0.0,
        )[:, r_ok]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            p20 = np.nanpercentile(u_c_clean, 20.0, axis=1)
            p25 = np.nanpercentile(u_c_clean, 25.0, axis=1)
            p10 = np.nanpercentile(u_c_clean, 10.0, axis=1)
            dim_mask = np.isfinite(u_c_clean) & (u_c_clean <= p20[:, None])
            col_floor = np.nanmedian(np.where(dim_mask, u_c_clean, np.nan), axis=1)
            sig_dim_med = np.nanmedian(np.where(dim_mask, sig_dim, np.nan), axis=1)
        sig_dim_med = np.where(np.isfinite(sig_dim_med), sig_dim_med, 0.0)
        n_dim = np.maximum(np.sum(dim_mask, axis=1), 1)

        col_floor_debiased = col_floor + z_q10 * sig_dim_med
        floor_est = np.where(
            np.isfinite(col_floor_debiased),
            np.minimum(col_floor_debiased, 0.0),
            0.0,
        )
        if bloom_detect_cube is None and np.any(floor_est < 0.0):
            clamped_without_bkg = True

        med_var_dim = float(_unit_normal_median_variance(int(np.median(n_dim))))
        theo_floor_var = (
            med_var_dim * np.sum(np.where(dim_mask, sig_dim * sig_dim, 0.0), axis=1) / n_dim
        )
        sigma_ped = np.where(
            np.isfinite(p25) & np.isfinite(p10), np.maximum(p25 - p10, 0.0), 0.0
        )
        tot_var = np.maximum(
            theo_floor_var + sig_dim_med * sig_dim_med + sigma_ped * sigma_ped, 1e-12
        )
        m_hat[:, c] = floor_est
        indep_m_var[:, c] = tot_var
        m_err[:, c] = np.sqrt(tot_var)
        c_bleed_mask[c] = False

    return clamped_without_bkg


def _bridge_bleed_columns(
    *,
    resolvable_cols: list[int],
    sci_pre_smear: np.ndarray,
    sci_sigma_pre: np.ndarray,
    safe_flat: np.ndarray,
    u_anchor_all: np.ndarray,
    anchor_pixel_ok: np.ndarray,
    clean_rows_2d: np.ndarray,
    indep_m_var: np.ndarray,
    shared_off_var: np.ndarray,
    col_uses_shared: np.ndarray,
    off_target_smear: np.ndarray | None,
    off_target_smear_err: np.ndarray | None,
    unres_smear_scale: float,
    m_hat: np.ndarray,
    m_err: np.ndarray,
    unresolved_col_mask: np.ndarray,
    use_good: np.ndarray,
    cfg: AlternateSmearConfig,
    precomputed_wing_excess: np.ndarray | None,
) -> None:
    """Bridge clean flank columns across bleeding columns and solve Huber smear locations."""
    n_t, ny, nx = sci_pre_smear.shape
    n_res = len(resolvable_cols)
    j_hat = np.zeros((n_t, ny, n_res), dtype=np.float64)
    flank_smear_var = np.zeros((n_t, n_res), dtype=np.float64)
    clamp_grad_var = np.zeros((n_t, n_res), dtype=np.float64)
    clean_sub = clean_rows_2d[:, resolvable_cols].copy()

    empty_u = np.zeros((n_t, ny, 0), dtype=np.float64)
    empty_iv = np.zeros((n_t, 0), dtype=np.float64)
    empty_ok = np.zeros((ny, 0), dtype=bool)
    empty_1d = np.zeros(0, dtype=np.float64)

    for idx_b, c in enumerate(resolvable_cols):
        if c > 0:
            ok_l_rev = anchor_pixel_ok[:, c - 1 :: -1]
            xs_l_rev = np.arange(c - 1, -1, -1, dtype=np.float64)
            u_l_rev = u_anchor_all[:, :, c - 1 :: -1]
            iv_l_rev = indep_m_var[:, c - 1 :: -1]
            sh_l_rev = col_uses_shared[c - 1 :: -1].astype(np.float64)
        else:
            ok_l_rev, xs_l_rev, u_l_rev, iv_l_rev, sh_l_rev = (
                empty_ok,
                empty_1d,
                empty_u,
                empty_iv,
                empty_1d,
            )
        has_l, c_l, u_l, v_l, s_l = _select_flank_side(
            ok_l_rev, xs_l_rev, u_l_rev, iv_l_rev, sh_l_rev, cfg.flank_width
        )

        if c + 1 < nx:
            ok_r = anchor_pixel_ok[:, c + 1 :]
            xs_r = np.arange(c + 1, nx, dtype=np.float64)
            u_r_sub = u_anchor_all[:, :, c + 1 :]
            iv_r = indep_m_var[:, c + 1 :]
            sh_r = col_uses_shared[c + 1 :].astype(np.float64)
        else:
            ok_r, xs_r, u_r_sub, iv_r, sh_r = (
                empty_ok,
                empty_1d,
                empty_u,
                empty_iv,
                empty_1d,
            )
        has_r, c_r, u_r, v_r, s_r = _select_flank_side(
            ok_r, xs_r, u_r_sub, iv_r, sh_r, cfg.flank_width
        )

        both = has_l & has_r
        only_l = has_l & (~has_r)
        only_r = has_r & (~has_l)
        row_has_anchor = clean_sub[:, idx_b] & (has_l | has_r)

        if not np.any(row_has_anchor):
            clean_sub[:, idx_b] = False
            if off_target_smear is not None:
                m_hat[:, c] = off_target_smear
                base_off = (
                    np.where(np.isfinite(off_target_smear_err), off_target_smear_err, 0.0)
                    if off_target_smear_err is not None
                    else np.zeros(n_t, dtype=np.float64)
                )
                m_err[:, c] = np.hypot(base_off, unres_smear_scale)
            else:
                unresolved_col_mask[c] = True
            continue

        clean_sub[:, idx_b] &= row_has_anchor
        if int(np.sum(clean_sub[:, idx_b])) < cfg.min_clean_rows:
            clean_sub[:, idx_b] = False
            unresolved_col_mask[c] = True
            continue

        alpha = (float(c) - c_l) / np.where(both, c_r - c_l, 1.0)
        u_comb = np.where(
            both[None, :],
            (1.0 - alpha[None, :]) * u_l + alpha[None, :] * u_r,
            np.where(only_l[None, :], u_l, np.where(only_r[None, :], u_r, 0.0)),
        )
        j_hat[:, :, idx_b] = safe_flat[:, c][None, :] * u_comb

        v_indep_comb = np.where(
            both[None, :],
            ((1.0 - alpha[None, :]) ** 2) * v_l + (alpha[None, :] ** 2) * v_r,
            np.where(only_l[None, :], v_l, np.where(only_r[None, :], v_r, 0.0)),
        )
        s_comb = np.where(
            both,
            (1.0 - alpha) * s_l + alpha * s_r,
            np.where(only_l, s_l, np.where(only_r, s_r, 0.0)),
        )
        row_flank_var = v_indep_comb + (s_comb[None, :] ** 2) * shared_off_var[:, None]

        row_clamp_dist = np.where(
            only_l, np.abs(float(c) - c_l), np.where(only_r, np.abs(float(c) - c_r), 0.0)
        )
        row_grad = np.zeros((n_t, ny), dtype=np.float64)
        r_ok = clean_sub[:, idx_b]
        clamped_l = np.flatnonzero(r_ok & only_l)
        if clamped_l.size > 0 and c > 0:
            row_grad[:, clamped_l] = _flank_gradient_4col(
                ok_l_rev, xs_l_rev, u_l_rev, clamped_l, safe_flat[:, c], j_hat[:, :, idx_b]
            )
        clamped_r = np.flatnonzero(r_ok & only_r)
        if clamped_r.size > 0 and (c + 1 < nx):
            row_grad[:, clamped_r] = _flank_gradient_4col(
                ok_r, xs_r, u_r_sub, clamped_r, safe_flat[:, c], j_hat[:, :, idx_b]
            )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            flank_smear_var[:, idx_b] = np.nanmedian(row_flank_var[:, r_ok], axis=1)
            clamped_rows = r_ok & (row_clamp_dist > 0.0)
            if np.any(clamped_rows):
                j_med_r = np.nanmedian(j_hat[:, clamped_rows, idx_b], axis=0)
                if np.any(np.isfinite(j_med_r)):
                    far_clamped = np.flatnonzero(clamped_rows)[
                        j_med_r <= float(np.nanpercentile(j_med_r, 35.0))
                    ]
                else:
                    far_clamped = np.flatnonzero(clamped_rows)
                if far_clamped.size == 0:
                    far_clamped = np.flatnonzero(clamped_rows)
                grad_term = (row_clamp_dist[None, far_clamped] * row_grad[:, far_clamped]) ** 2
                clamp_grad_var[:, idx_b] = np.nanmedian(grad_term, axis=1)

    active_rel = [i for i in range(n_res) if np.any(clean_sub[:, i])]
    if not active_rel:
        return

    act_cols = [resolvable_cols[i] for i in active_rel]
    j_act = j_hat[:, :, active_rel]
    clean_act = clean_sub[:, active_rel]
    i_act = sci_pre_smear[:, :, act_cols]
    sig_act = np.maximum(sci_sigma_pre[:, :, act_cols], 1e-6)

    if precomputed_wing_excess is not None:
        wing_excess_2d = precomputed_wing_excess[:, act_cols]
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            j_med_2d = np.nanmedian(j_act[use_good], axis=0)
            j_floor = np.nanpercentile(
                np.where(clean_act, j_med_2d, np.nan), 10.0, axis=0, keepdims=True
            )
        j_floor = np.where(np.isfinite(j_floor), j_floor, 0.0)
        wing_excess_2d = np.where(
            clean_act & np.isfinite(j_med_2d), np.maximum(j_med_2d - j_floor, 0.0), 0.0
        )

    inv_wing_w = np.where(
        clean_act[None, :, :], 1.0 / (wing_excess_2d[None, :, :] ** 2 + sig_act * sig_act), 0.0
    )
    loc_est, loc_err = _solve_huber_smear_location(i_act - j_act, inv_wing_w, sig_act, clean_act)
    for k_idx, col_idx in enumerate(act_cols):
        rel_i = active_rel[k_idx]
        m_hat[:, col_idx] = loc_est[:, k_idx]
        tot_var = (
            loc_err[:, k_idx] ** 2 + flank_smear_var[:, rel_i] + clamp_grad_var[:, rel_i]
        )
        m_err[:, col_idx] = np.sqrt(np.maximum(tot_var, 1e-12))


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
    config: AlternateSmearConfig | None = None,
    bloom_detect_cube: np.ndarray | None = None,
    good_cad_mask: np.ndarray | None = None,
    precomputed_bleed_mask: np.ndarray | None = None,
    precomputed_wing_excess: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    """Estimate concavity-guarded collateral-anchored smear ``M_hat`` across ``(N_t, nx)``."""
    cfg = config if config is not None else AlternateSmearConfig()
    n_t, ny, nx = sci_pre_smear.shape
    safe_flat = np.where(np.isfinite(flat_2d) & (flat_2d > 0.0), flat_2d, 1.0)

    sat_lim = cfg.saturation_limit
    if sat_lim is not None and np.isfinite(sat_lim) and sat_lim > 0.0:
        sat_scale = float(sat_lim)
    else:
        fin_sci = sci_pre_smear[np.isfinite(sci_pre_smear)]
        sat_scale = max(float(np.max(np.abs(fin_sci))) if fin_sci.size > 0 else 1.0, 1.0)
    unres_smear_scale = float(TESS_PARALLEL_SHIFTS) * float(cfg.row_smear_weight) * sat_scale

    clean_rows_2d = _build_clean_rows_2d(
        raw_sat_2d,
        sat_bleed_2d,
        source_mask_2d,
        bleed_dilation_rows=cfg.bleed_dilation_rows,
        min_clean_rows=cfg.min_clean_rows,
    )
    use_good = (
        good_cad_mask
        if (good_cad_mask is not None and np.any(good_cad_mask))
        else np.ones(n_t, dtype=bool)
    )
    col_has_sci_sat = np.any(raw_sat_2d, axis=0)
    col_has_sm_sat = (
        np.asarray(collat_sm_sat_cols, dtype=bool)
        if collat_sm_sat_cols is not None
        else np.zeros(nx, dtype=bool)
    )

    if precomputed_bleed_mask is not None:
        c_bleed_mask = np.array(precomputed_bleed_mask, dtype=bool, copy=True)
    else:
        if collat_smear_stamp is not None:
            sub_for_bloom = sci_pre_smear[use_good] - collat_smear_stamp[use_good, None, :]
            col_has_bloom = _detect_bloom_columns(
                sub_for_bloom, sci_sigma_pre[use_good], clean_rows_2d
            )
        else:
            bloom_src = (
                bloom_detect_cube[use_good]
                if bloom_detect_cube is not None
                else sci_pre_smear[use_good]
            )
            col_has_bloom = _detect_bloom_columns(
                bloom_src, sci_sigma_pre[use_good], clean_rows_2d
            )
            if bloom_detect_cube is None and nx >= 3 and not np.all(col_has_bloom):
                sky_clean_2d = clean_rows_2d & (~col_has_bloom)[None, :]
                if np.any(sky_clean_2d):
                    sky_base, coh_err = _estimate_cutout_sky_plane(
                        sci_pre_smear[use_good], sci_sigma_pre[use_good], sky_clean_2d
                    )
                    if np.any(sky_base > 0.0):
                        col_has_bloom |= _detect_bloom_columns(
                            sci_pre_smear[use_good] - sky_base,
                            sci_sigma_pre[use_good],
                            clean_rows_2d,
                            coherent_spatial_err=coh_err,
                        )
        c_bleed_mask = col_has_sci_sat | col_has_sm_sat | col_has_bloom

    if collat_smear_stamp is not None:
        m_hat = np.array(collat_smear_stamp, dtype=np.float64, copy=True)
        m_err = (
            np.array(collat_smear_err_stamp, dtype=np.float64, copy=True)
            if collat_smear_err_stamp is not None
            else np.sqrt(
                float(_unit_normal_median_variance(10)) * np.maximum(np.abs(m_hat), 1.0)
            )
        )
    else:
        m_hat = np.zeros((n_t, nx), dtype=np.float64)
        m_err = np.zeros((n_t, nx), dtype=np.float64)

    if not np.any(c_bleed_mask):
        return m_hat, m_err, np.zeros(nx, dtype=bool), False

    indep_m_var = m_err * m_err
    shared_off_var = np.zeros(n_t, dtype=np.float64)
    col_uses_shared = np.zeros(nx, dtype=bool)
    zero_floor_clamped = _anchor_all_bloomed_columns(
        sci_pre_smear=sci_pre_smear,
        sci_sigma_pre=sci_sigma_pre,
        safe_flat=safe_flat,
        off_target_smear=off_target_smear,
        off_target_smear_err=off_target_smear_err,
        raw_sat_2d=raw_sat_2d,
        sat_bleed_2d=sat_bleed_2d,
        strap_col_mask=strap_col_mask,
        source_mask_2d=source_mask_2d,
        clean_rows_2d=clean_rows_2d,
        col_has_sci_sat=col_has_sci_sat,
        c_bleed_mask=c_bleed_mask,
        m_hat=m_hat,
        m_err=m_err,
        indep_m_var=indep_m_var,
        shared_off_var=shared_off_var,
        col_uses_shared=col_uses_shared,
        config=cfg,
        bloom_detect_cube=bloom_detect_cube,
    )
    if not np.any(c_bleed_mask):
        return m_hat, m_err, np.zeros(nx, dtype=bool), zero_floor_clamped

    u_anchor_all = (sci_pre_smear - m_hat[:, None, :]) / safe_flat[None, :, :]
    anchor_pixel_ok = (
        (~c_bleed_mask)[None, :] & (~strap_col_mask)[None, :] & (~sat_bleed_2d) & (~source_mask_2d)
    )
    if not np.any(anchor_pixel_ok) and np.any(~c_bleed_mask):
        anchor_pixel_ok = (~c_bleed_mask)[None, :] & (~sat_bleed_2d) & (~source_mask_2d)

    unresolved_col_mask = np.zeros(nx, dtype=bool)
    resolvable_cols: list[int] = []
    for c in np.flatnonzero(c_bleed_mask):
        if int(np.sum(clean_rows_2d[:, c])) < cfg.min_clean_rows:
            unresolved_col_mask[int(c)] = True
        else:
            resolvable_cols.append(int(c))

    if resolvable_cols:
        _bridge_bleed_columns(
            resolvable_cols=resolvable_cols,
            sci_pre_smear=sci_pre_smear,
            sci_sigma_pre=sci_sigma_pre,
            safe_flat=safe_flat,
            u_anchor_all=u_anchor_all,
            anchor_pixel_ok=anchor_pixel_ok,
            clean_rows_2d=clean_rows_2d,
            indep_m_var=indep_m_var,
            shared_off_var=shared_off_var,
            col_uses_shared=col_uses_shared,
            off_target_smear=off_target_smear,
            off_target_smear_err=off_target_smear_err,
            unres_smear_scale=unres_smear_scale,
            m_hat=m_hat,
            m_err=m_err,
            unresolved_col_mask=unresolved_col_mask,
            use_good=use_good,
            cfg=cfg,
            precomputed_wing_excess=precomputed_wing_excess,
        )

    if np.any(unresolved_col_mask):
        resolved_cols = np.flatnonzero(~unresolved_col_mask)
        for c_unres in np.flatnonzero(unresolved_col_mask):
            if off_target_smear is not None:
                m_hat[:, c_unres] = off_target_smear
                base_err = (
                    np.asarray(off_target_smear_err, dtype=np.float64)
                    if off_target_smear_err is not None
                    else np.zeros(n_t, dtype=np.float64)
                )
            elif resolved_cols.size > 0:
                nearest = int(resolved_cols[np.argmin(np.abs(resolved_cols - c_unres))])
                m_hat[:, c_unres] = m_hat[:, nearest]
                base_err = m_err[:, nearest].copy()
            else:
                base_err = np.zeros(n_t, dtype=np.float64)
            finite_cad = np.any(np.isfinite(sci_pre_smear[:, :, c_unres]), axis=1)
            safe_base = np.where(np.isfinite(base_err), base_err, 0.0)
            m_err[:, c_unres] = np.where(
                finite_cad, np.hypot(safe_base, unres_smear_scale), np.nan
            )

    return m_hat, m_err, unresolved_col_mask, zero_floor_clamped
