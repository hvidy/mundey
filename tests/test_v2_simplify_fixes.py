"""Regression tests for mundey v2 simplification, scientific bugfixes, chunking, and plotting."""

from __future__ import annotations

from pathlib import Path

import jax.numpy as jnp
import lightkurve as lk
import numpy as np
import pytest
from astropy.io import fits

import mundey
from mundey import (
    CalibrationConfig,
    PixelQuality,
    desmear_calibrated_cutout,
    plot_diagnostics,
    read_pixel_quality,
    read_smear_err,
    recalibrate_tpf,
)
from mundey._precision import run_in
from mundey.calibration import (
    AlternateSmearConfig,
    _estimate_alternate_smear_stamp,
    _invert_lde_undershoot_forward_impl,
    _solve_huber_smear_location,
)


def test_sirius_s06_bleed_neighbour_trench_fix(tmp_path: Path) -> None:
    """[MUND-ADV-1] Masking sat_bleed_2d and ±1-col bleed neighbours eliminates Sirius trenches."""
    nb_dir = Path(__file__).resolve().parent.parent / "notebooks"
    sirius_tpf = nb_dir / "tess2018349182459-s0006-0000000322899250-0126-s_tp.fits"
    tpf_sub = lk.read(str(sirius_tpf))[:20]
    tpf_sub.path = str(sirius_tpf)

    cal = recalibrate_tpf(
        tpf_sub,
        ddir=nb_dir,
        cache_dir=tmp_path / "sirius_trench",
        smear="automatic",
        offline=True,
    )
    med = np.nanmedian(cal.flux.value, axis=0)
    col4_bot50 = float(np.nanmedian(med[:50, 4]))
    col6_bot50 = float(np.nanmedian(med[:50, 6]))
    # Previously biased down to 53.4 e-/s (col 4) and 63.9 e-/s (col 6)
    assert col4_bot50 > 90.0
    assert col6_bot50 > 110.0


def test_inv_wing_w_time_median_eliminates_selection_bias() -> None:
    """[MUND-ADV-2] Time-medianed inv_wing_w must remove +1.5 sigma/cadence selection bias."""
    rng = np.random.default_rng(42)
    n_t, ny, nx = 200, 1287, 11
    sigma = 10.0
    sci_sky = 700.0 + rng.normal(0.0, sigma, size=(n_t, ny, nx))
    rs = np.zeros((ny, nx), dtype=bool)
    rs[640:645, 5] = True
    sci_sky[:, rs] = 2.0e5
    csm = np.full((n_t, nx), 500.0)
    csm[:, 5] = 5500.0

    m_hat, m_err, _, _ = _estimate_alternate_smear_stamp(
        sci_pre_smear=sci_sky,
        sci_sigma_pre=np.full((n_t, ny, nx), sigma),
        flat_2d=np.ones((ny, nx)),
        collat_smear_stamp=csm,
        collat_smear_err_stamp=np.full((n_t, nx), 2.0),
        collat_sm_sat_cols=np.zeros(nx, dtype=bool),
        off_target_smear=np.full(n_t, 500.0),
        off_target_smear_err=np.full(n_t, 1.0),
        raw_sat_2d=rs,
        sat_bleed_2d=rs,
        strap_col_mask=np.zeros(nx, dtype=bool),
        source_mask_2d=np.zeros((ny, nx), dtype=bool),
    )
    bias = m_hat[:, 5] - 500.0
    se_mean = float(np.std(bias, ddof=1) / np.sqrt(n_t))
    # Previously +1.34 e-/frame (+60 sigma); now unbiased within 3 SE
    assert abs(float(np.mean(bias))) < max(0.25, 3.5 * se_mean)
    assert float(np.median(m_err[:, 5])) > 0.0


def test_lde_undershoot_nan_gap_and_leading_nan_reset() -> None:
    """[MUND-ADV-3] Causal LDE inversion must reset to steady state after any NaN pixel."""
    u_c = 0.00231
    x_t = np.array([[40000.0, 50.0, 50.0, 50.0, 50.0, 50.0]])
    w_o = x_t.copy()
    w_o[:, 1:] = x_t[:, 1:] - u_c * x_t[:, :-1]
    w_o[:, 0] = x_t[:, 0] * (1.0 - u_c)
    w_gap = w_o.copy()
    w_gap[:, 1:5] = np.nan

    w_n0 = np.array([[np.nan, 10000.0 * (1.0 - u_c), 10000.0 * (1.0 - u_c)]])

    with run_in("float64"):
        x_gap = np.asarray(
            _invert_lde_undershoot_forward_impl(jnp.asarray(w_gap), jnp.asarray(u_c))
        )
        x_n0 = np.asarray(
            _invert_lde_undershoot_forward_impl(jnp.asarray(w_n0), jnp.asarray(u_c))
        )

    assert np.all(np.isnan(x_gap[:, 1:5]))
    # Previously +92.28 DN teleportation error; now exact steady-state restart
    assert abs(float(x_gap[0, 5] - 50.0)) < 1e-10
    assert np.isnan(x_n0[0, 0])
    # Previously -23.10 DN col-0 NaN bias; now exact 1/(1-u) restart
    np.testing.assert_allclose(x_n0[0, 1:], 10000.0, rtol=1e-12)


def test_smear_automatic_row500_saturated_star_stays_standard(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """[MUND-ADV-4] Saturated star far below row 2048 with clean smrow must stay 'standard'."""
    bundle = tpf_builder(
        tmp_path / "row500_sat",
        col0_1based=675,
        row0_1based=500,
        nx=15,
        ny=15,
        n_cadences=12,
        include_bleed=True,
        bleed_in_smrow=False,
        flag_bad_cadence=False,
        noise_sigma=2.0,
        sky_dark_level=50.0,
        seed=4004,
    )
    cal_auto = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=tmp_path / "c_auto",
        smear="automatic",
        offline=True,
    )
    cal_std = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=tmp_path / "c_std",
        smear="standard",
        offline=True,
    )
    assert cal_auto.hdu[0].header["SMEARMOD"] == "standard"
    pq = read_pixel_quality(cal_auto)
    assert np.any((pq & int(PixelQuality.SAT_BLEED)) != 0)
    assert not np.any((pq & int(PixelQuality.SMEAR_UNRESOLVED)) != 0)
    np.testing.assert_allclose(read_smear_err(cal_auto), read_smear_err(cal_std), rtol=1e-6)


def test_quality_flagged_saturating_cadence_does_not_poison_sector(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """[MUND-ADV-5] A single QUALITY=32 saturating cadence must not flag SAT_BLEED or bloom."""
    bundle = tpf_builder(
        tmp_path / "qual_poison",
        col0_1based=675,
        row0_1based=500,
        nx=11,
        ny=11,
        n_cadences=16,
        include_bleed=False,
        bleed_in_smrow=False,
        flag_bad_cadence=False,
        noise_sigma=1.0,
        sky_dark_level=50.0,
        seed=5005,
    )
    # Corrupt cadence 0 with QUALITY=32 (coarse point) and full-stamp saturation
    with fits.open(bundle.tpf_path, mode="update", memmap=False) as hdul:
        hdul[1].data["QUALITY"][0] = 32
        hdul[1].data["RAW_CNTS"][0] = 15_000_000
        hdul.flush()

    tpf_lk = lk.TessTargetPixelFile(str(bundle.tpf_path), quality_bitmask="default")
    cal = recalibrate_tpf(
        tpf_lk,
        ddir=bundle.ddir,
        cache_dir=tmp_path / "c_qual",
        smear="automatic",
        offline=True,
    )
    pq = read_pixel_quality(cal)
    assert cal.hdu[0].header["SMEARMOD"] == "standard"
    assert not np.any((pq & int(PixelQuality.SAT_BLEED)) != 0)
    assert not np.any((pq & int(PixelQuality.SMEAR_UNRESOLVED)) != 0)
    # All N_raw rows remain calibrated in HDU[1]
    assert len(cal.hdu[1].data["FLUX"]) == 16


def test_huber_weighted_median_and_shared_off_target_smear_variance() -> None:
    """[MUND-ADV-6] Huber resists +1e5 outliers and preserves shared off_target_smear variance."""
    rng = np.random.default_rng(123)
    res_o = rng.normal(0.0, 1.0, size=(60, 200, 1))
    res_o[:, 0, 0] = 1.0e5
    lh, eh = _solve_huber_smear_location(
        res_o, np.ones_like(res_o), np.ones_like(res_o), np.ones((200, 1), dtype=bool)
    )
    # Previously loc=+5.02, err=5.02; now unbiased around 0 with SE ~0.075
    assert abs(float(np.mean(lh))) < 0.1
    assert 0.05 < float(np.median(eh)) < 0.15

    # Shared off_target_smear_err = 50.0 must not shrink by sqrt(2 * flank_width)
    rs6 = np.zeros((60, 11), dtype=bool)
    rs6[25:35, 5] = True
    sci6 = np.full((20, 60, 11), 200.0)
    sci6[:, rs6] = 2.0e5
    _, me6, _, _ = _estimate_alternate_smear_stamp(
        sci_pre_smear=sci6,
        sci_sigma_pre=np.full((20, 60, 11), 1.0),
        flat_2d=np.ones((60, 11)),
        collat_smear_stamp=np.full((20, 11), 5000.0),
        collat_smear_err_stamp=np.full((20, 11), 2.0),
        collat_sm_sat_cols=np.zeros(11, dtype=bool),
        off_target_smear=np.full(20, 100.0),
        off_target_smear_err=np.full(20, 50.0),
        raw_sat_2d=rs6,
        sat_bleed_2d=rs6,
        strap_col_mask=np.zeros(11, dtype=bool),
        source_mask_2d=np.zeros((60, 11), dtype=bool),
        config=AlternateSmearConfig(),
    )
    assert float(me6[0, 5]) >= 49.5


def test_desmear_all_bloomed_positive_sky_warns_and_marks_bkgsub() -> None:
    """[MUND-ADV-7] All-bloomed BACKAPP=False cutout without FLUX_BKG warns & sets BKGSUB=True."""
    rng = np.random.default_rng(77)
    n_t, ny, nx = 16, 40, 11
    yy, xx = np.meshgrid(np.arange(ny) - ny // 2, np.arange(nx) - 5, indexing="ij")
    star = 1.5e5 / ((xx * xx + yy * yy + 0.5) ** 1.5)
    fl = 100.0 - 150.0 + star[None, :, :] + rng.normal(0.0, 2.0, (n_t, ny, nx))

    h0 = fits.PrimaryHDU()
    h0.header.update({"SECTOR": 6, "CAMERA": 2, "CCD": 3, "BACKAPP": False})
    h1 = fits.BinTableHDU.from_columns(
        fits.ColDefs(
            [
                fits.Column(name="TIME", format="D", array=1500.0 + np.arange(n_t) * 0.01),
                fits.Column(name="CADENCENO", format="J", array=np.arange(n_t, dtype=np.int32)),
                fits.Column(
                    name="FLUX",
                    format=f"{ny * nx}E",
                    dim=f"({nx},{ny})",
                    array=fl.astype(np.float32),
                ),
                fits.Column(
                    name="FLUX_ERR",
                    format=f"{ny * nx}E",
                    dim=f"({nx},{ny})",
                    array=np.full((n_t, ny, nx), 2.0, dtype=np.float32),
                ),
                fits.Column(name="QUALITY", format="J", array=np.zeros(n_t, dtype=np.int32)),
            ]
        ),
        name="PIXELS",
    )
    h1.header.update(
        {"1CRV5P": 675, "2CRV5P": 1000, "INT_TIME": 1.98, "READTIME": 0.02, "BACKAPP": False}
    )
    h2 = fits.ImageHDU(data=np.ones((ny, nx), dtype=np.int32), name="APERTURE")
    h2.header.update({"CRVAL1P": 675, "CRVAL2P": 1000})
    tpf = lk.TessTargetPixelFile(fits.HDUList([h0, h1, h2]), quality_bitmask="none")

    with pytest.warns(UserWarning, match="degenerate with uniform smear bloom"):
        out = desmear_calibrated_cutout(tpf)
    assert out.meta["BKGSUB"] is True
    assert out.hdu[1].header["BACKAPP"] is True


def test_chunked_calibration_config_and_plot_diagnostics(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """Multi-chunk science calibration, CalibrationConfig, and plot_diagnostics must work."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bundle = tpf_builder(
        tmp_path / "chunk_test",
        col0_1based=210,
        row0_1based=1993,
        nx=10,
        ny=56,
        n_cadences=10,
        bleed_in_smrow=True,
        include_bleed=True,
        flag_bad_cadence=False,
        noise_sigma=1.0,
        sky_dark_level=120.0,
        seed=808,
    )
    raw_tpf = lk.read(bundle.tpf_path)

    cfg_full = CalibrationConfig(
        smear="automatic",
        ddir=bundle.ddir,
        cache_dir=tmp_path / "c_full",
        offline=True,
        chunk_rows=4096,
    )
    cfg_chunked = CalibrationConfig(
        smear="automatic",
        ddir=bundle.ddir,
        cache_dir=tmp_path / "c_chunk",
        offline=True,
        chunk_rows=4,  # forces 3 science chunks (4 + 4 + 2)
    )
    cal_full = recalibrate_tpf(bundle.tpf_path, config=cfg_full)
    cal_chunked = recalibrate_tpf(bundle.tpf_path, config=cfg_chunked)
    np.testing.assert_allclose(cal_chunked.flux.value, cal_full.flux.value, rtol=1e-2, atol=2.0)

    # Standard mode chunking is bit-for-bit identical
    cal_std_full = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=tmp_path / "cs1",
        smear="standard",
        offline=True,
    )
    cal_std_ch = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=tmp_path / "cs2",
        smear="standard",
        offline=True,
        chunk_rows=3,
    )
    np.testing.assert_allclose(cal_std_ch.flux.value, cal_std_full.flux.value, rtol=1e-6)

    # Diagnostic plotting (2-panel and 3-panel)
    fig2, axes2 = plot_diagnostics(cal_full)
    assert len(axes2) == 2
    plt.close(fig2)

    fig3, axes3 = plot_diagnostics(cal_full, raw_tpf=raw_tpf)
    assert len(axes3) == 3
    out_png = tmp_path / "diag.png"
    fig3.savefig(out_png)
    assert out_png.is_file() and out_png.stat().st_size > 1000
    plt.close(fig3)
    assert "plot_diagnostics" in mundey.__all__
