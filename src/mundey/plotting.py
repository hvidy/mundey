"""Lightweight visual diagnostics for ``mundey`` recalibration results."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np

from ._contract import PixelQuality, read_pixel_quality, read_smear_err


def _extract_median_stamp(tpf: Any) -> np.ndarray:
    flux = getattr(tpf, "flux", None)
    arr = np.asarray(getattr(flux, "value", flux), dtype=np.float64)
    if arr.ndim == 2:
        return arr
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(arr, axis=0)


def _signed_asinh_image(img: np.ndarray) -> tuple[np.ndarray, float]:
    """Compute a soft signed ``arcsinh`` display transform and characteristic scale."""
    fin = img[np.isfinite(img)]
    if fin.size == 0:
        return np.zeros_like(img, dtype=np.float64), 1.0
    scale = float(np.nanpercentile(np.abs(fin), 80.0))
    scale = max(scale, 1.0)
    return np.arcsinh(img / scale), scale


def plot_diagnostics(
    cal_tpf: Any,
    raw_tpf: Any | None = None,
    ax: Any | None = None,
) -> tuple[Any, np.ndarray]:
    """Plot median calibrated stamp, ``PIXEL_QUALITY`` overlay, and per-column smear profile.

    Parameters
    ----------
    cal_tpf : TargetPixelFile
        Recalibrated or desmeared TPF returned by :func:`~mundey.recalibrate_tpf` or
        :func:`~mundey.desmear_calibrated_cutout`.
    raw_tpf : TargetPixelFile or None, optional
        Original uncalibrated/SPOC TPF for side-by-side before/after comparison.
    ax : matplotlib Axes or array of Axes, optional
        Existing axes to draw into (2 panels without ``raw_tpf``, 3 panels with ``raw_tpf``).

    Returns
    -------
    fig, axes
        Matplotlib ``Figure`` and 1D array of ``Axes`` objects.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "mundey.plot_diagnostics requires matplotlib (`pip install 'mundey[plot]'`)."
        ) from exc

    cal_med = _extract_median_stamp(cal_tpf)
    ny, nx = cal_med.shape
    pq = read_pixel_quality(cal_tpf)
    smear_err = read_smear_err(cal_tpf)

    n_panels = 3 if raw_tpf is not None else 2
    if ax is None:
        fig, axes = plt.subplots(
            1, n_panels, figsize=(4.5 * n_panels, 4.5), constrained_layout=True
        )
        axes_arr = np.atleast_1d(axes).ravel()
    else:
        axes_arr = np.atleast_1d(ax).ravel()
        if axes_arr.size < n_panels:
            raise ValueError(
                f"Expected at least {n_panels} matplotlib Axes, got {axes_arr.size}."
            )
        fig = axes_arr[0].figure

    aspect = "auto" if (max(ny, nx) / max(min(ny, nx), 1) > 4.0) else "equal"
    cmap = plt.get_cmap("viridis").with_extremes(bad="lightgray")

    hdr0 = cal_tpf.hdu[0].header if hasattr(cal_tpf, "hdu") else {}
    smear_mode = str(hdr0.get("SMEARMOD", "unknown"))
    cad_miss = int(hdr0.get("CADMISS", 0))

    panel_idx = 0
    if raw_tpf is not None:
        raw_med = _extract_median_stamp(raw_tpf)
        raw_disp, raw_scale = _signed_asinh_image(raw_med)
        im_raw = axes_arr[panel_idx].imshow(
            raw_disp, origin="lower", aspect=aspect, cmap=cmap, interpolation="nearest"
        )
        axes_arr[panel_idx].set_title(f"Before (asinh, scale={raw_scale:.1f} e-/s)")
        axes_arr[panel_idx].set_xlabel("Stamp column (px)")
        axes_arr[panel_idx].set_ylabel("Stamp row (px)")
        cbar_raw = fig.colorbar(im_raw, ax=axes_arr[panel_idx], fraction=0.046, pad=0.04)
        cbar_raw.set_label("asinh(flux / scale)")
        panel_idx += 1

    cal_disp, cal_scale = _signed_asinh_image(cal_med)
    ax_cal = axes_arr[panel_idx]
    im_cal = ax_cal.imshow(
        cal_disp, origin="lower", aspect=aspect, cmap=cmap, interpolation="nearest"
    )
    ax_cal.set_title(f"After ({smear_mode}, CADMISS={cad_miss})")
    ax_cal.set_xlabel("Stamp column (px)")
    ax_cal.set_ylabel("Stamp row (px)")
    cbar_cal = fig.colorbar(im_cal, ax=ax_cal, fraction=0.046, pad=0.04)
    cbar_cal.set_label(f"asinh(flux / {cal_scale:.1f} e-/s)")

    sat_mask = (pq & int(PixelQuality.SAT_BLEED)) != 0
    unres_cols = np.any((pq & int(PixelQuality.SMEAR_UNRESOLVED)) != 0, axis=0)
    strap_cols = np.any((pq & int(PixelQuality.STRAP_COL)) != 0, axis=0)
    if np.any(sat_mask) and ny >= 2 and nx >= 2:
        ax_cal.contour(
            sat_mask.astype(float), levels=[0.5], colors=["red"], linewidths=0.8
        )
    for col_i in np.flatnonzero(unres_cols):
        ax_cal.axvline(float(col_i), color="orange", ls="--", lw=0.8, alpha=0.8)
    for col_i in np.flatnonzero(strap_cols):
        ax_cal.axvline(float(col_i), color="cyan", ls=":", lw=0.8, alpha=0.8)
    panel_idx += 1

    # Panel: per-column outer-row median flux & SMEAR_ERR
    ax_col = axes_arr[panel_idx]
    cols = np.arange(nx)
    n_edge = max(1, min(20, ny // 4))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        edge_flux = np.nanmedian(cal_med[:n_edge, :], axis=0)
        se_med = (
            np.nanmedian(smear_err, axis=0)
            if smear_err is not None
            else np.full(nx, np.nan, dtype=np.float64)
        )

    ax_col.plot(cols, edge_flux, "o-", ms=3, lw=1.2, label=f"Bottom {n_edge}-row median flux")
    if np.any(np.isfinite(se_med)):
        ax_col.fill_between(
            cols,
            edge_flux - se_med,
            edge_flux + se_med,
            alpha=0.25,
            label="± median SMEAR_ERR",
        )
    if np.any(unres_cols):
        ax_col.scatter(
            cols[unres_cols],
            edge_flux[unres_cols],
            color="orange",
            marker="x",
            s=36,
            zorder=5,
            label="SMEAR_UNRESOLVED",
        )
    ax_col.axhline(0.0, color="0.4", ls="--", lw=0.8)
    ax_col.set_xlabel("Stamp column (px)")
    ax_col.set_ylabel("Flux (e-/s)")
    ax_col.set_title("Outer-row flux & SMEAR_ERR")
    ax_col.legend(fontsize="x-small", loc="best")

    return fig, axes_arr
