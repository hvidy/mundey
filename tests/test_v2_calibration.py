"""Tests for mundey v2.0 physics, @in_float64 precision, alternate smear, and CADENCENO pairing."""

from __future__ import annotations

from pathlib import Path

import jax
import lightkurve as lk
import numpy as np
import pytest
from astropy.io import fits

from mundey import (
    PixelQuality,
    desmear_calibrated_cutout,
    read_pixel_quality,
    read_smear_err_column,
    recalibrate_tpf,
)
from mundey.calibration import _detect_bloom_columns, fit_onedblack_legendre


def test_unbleeding_roundtrip_1_2_and_4_outputs_within_1e4(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """Unbleeding 1-, 2-, and 4-output TPFs must recover S + B within < 1e-4 relative error."""
    # 1. 1-output (A) & 2-output (A/B) high-flux unbleeding round-trips
    for label, col0, nx, expected_outs in (
        ("1out_a", 210, 8, ("A",)),
        ("2out_ab", 553, 8, ("A", "B")),
    ):
        bundle = tpf_builder(
            tmp_path / label,
            col0_1based=col0,
            row0_1based=1400,
            nx=nx,
            ny=10,
            n_cadences=6,
            include_bleed=False,
            flag_bad_cadence=False,
            noise_sigma=0.0,
            sky_dark_level=450.0,
            seed=11,
        )
        assert bundle.outputs == expected_outs
        orig_tpf = lk.read(bundle.tpf_path)
        orig_flux_snapshot = orig_tpf.hdu[1].data["FLUX"].copy()

        cal = recalibrate_tpf(
            orig_tpf,
            ddir=bundle.ddir,
            smear="standard",
            offline=True,
        )
        # Pure functional guarantee: caller's TPF is never mutated
        np.testing.assert_array_equal(orig_tpf.hdu[1].data["FLUX"], orig_flux_snapshot)
        assert "SMEAR_ERR" not in (orig_tpf.hdu[1].data.names or ())
        assert jax.config.jax_enable_x64 is False

        non_cr = bundle.cr_cube == 0.0
        truth = bundle.true_star_cube + bundle.true_sky_cube
        rel_err = np.abs(cal.flux.value[non_cr] - truth[non_cr]) / truth[non_cr]
        assert float(np.max(rel_err)) < 1e-4

    # 2. Wide 4-output stamp straddling A (<=556), B (557..1068), C (1069..1580), D (>=1581)
    bundle_4 = tpf_builder(
        tmp_path / "4out_abcd",
        col0_1based=555,
        row0_1based=1200,
        nx=1030,  # cols 555..1584 -> spans A, B, C, and D!
        ny=4,
        n_cadences=2,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        sky_dark_level=850.0,
        seed=19,
    )
    assert bundle_4.outputs == ("A", "B", "C", "D")
    cal_4 = recalibrate_tpf(bundle_4.tpf_path, ddir=bundle_4.ddir, smear="standard", offline=True)
    non_cr_4 = bundle_4.cr_cube == 0.0
    truth_4 = bundle_4.true_star_cube + bundle_4.true_sky_cube
    rel_err_4 = np.abs(cal_4.flux.value[non_cr_4] - truth_4[non_cr_4]) / truth_4[non_cr_4]
    assert float(np.max(rel_err_4)) < 1e-4


def test_degree6_maxdomain_clamping_and_undershoot_first_preserves_halo(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """Degree-6 >maxDomain bleed must not corrupt downstream readout halo columns."""
    # Output A reads left-to-right (0 -> nx-1); Output B reads right-to-left (nx-1 -> 0)
    for out_tag, col0, downstream_offset in (("outA", 200, +1), ("outB", 700, -1)):
        bundle = tpf_builder(
            tmp_path / out_tag,
            col0_1based=col0,
            row0_1based=1500,
            nx=10,
            ny=14,
            n_cadences=4,
            include_bleed=True,
            degree6_linearity=True,
            flag_bad_cadence=False,
            noise_sigma=0.0,
            sky_dark_level=300.0,
            seed=55,
        )
        cal = recalibrate_tpf(bundle.tpf_path, ddir=bundle.ddir, smear="standard", offline=True)

        bleed_col = 5
        downstream_col = bleed_col + downstream_offset
        bleed_rows = bundle.bleed_mask[:, bleed_col]
        assert np.any(bleed_rows)

        # Adjacent downstream halo column at the exact rows of saturated bleed is uncorrupted!
        got_halo = cal.flux.value[:, bleed_rows, downstream_col]
        exp_halo = (bundle.true_star_cube + bundle.true_sky_cube)[:, bleed_rows, downstream_col]
        np.testing.assert_allclose(got_halo, exp_halo, rtol=5e-4, atol=0.08)

        # Saturated bleed pixels stay positive (not flipped negative by unclamped degree-6 poly)
        assert np.all(cal.flux.value[:, bleed_rows, bleed_col] > 1.0e5)

        # PIXEL_QUALITY flags both saturated pixels AND 1-px downstream undershoot neighbour
        pq = read_pixel_quality(cal)
        assert np.all((pq[bleed_rows, bleed_col] & int(PixelQuality.SAT_BLEED)) != 0)
        assert np.all((pq[bleed_rows, downstream_col] & int(PixelQuality.SAT_BLEED)) != 0)


def test_raw_cnts_minus_one_sentinel_and_20s_fast_collateral(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """RAW_CNTS == -1 sentinels must become NaN and 20-s CRMITEN=False must scale 10x B_2D."""
    bundle = tpf_builder(
        tmp_path / "fast20s",
        col0_1based=220,
        row0_1based=1450,
        nx=8,
        ny=10,
        n_cadences=6,
        fast_20s=True,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        sky_dark_level=400.0,
        seed=73,
    )
    assert bundle.tpf_path.name.endswith("-a_fast-tp.fits")

    # Inject RAW_CNTS == -1 missing cadence (row 2) and isolated pixel sentinel (row 4, 3, 3)
    with fits.open(bundle.tpf_path, mode="update") as hdul:
        hdul[1].data["RAW_CNTS"][2, :, :] = -1
        hdul[1].data["RAW_CNTS"][4, 3, 3] = -1

    cal = recalibrate_tpf(bundle.tpf_path, ddir=bundle.ddir, smear="standard", offline=True)
    flux = cal.hdu[1].data["FLUX"]
    assert np.all(np.isnan(flux[2]))
    assert np.isnan(flux[4, 3, 3])

    valid = np.isfinite(flux) & (bundle.cr_cube == 0.0)
    truth = bundle.true_star_cube + bundle.true_sky_cube
    np.testing.assert_allclose(flux[valid], truth[valid], rtol=5e-4, atol=0.1)


def test_concavity_guarded_alternate_smear_pulsation_and_edge_cases(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """smear='alternate' must preserve pulsation amplitude (<0.2% bias) & handle edge cases."""
    # 1. Strong 2.5% pulsation on tall stamp reaching row 2048 with strap columns in stamp
    # RAWX 125..148 includes CCID-80 straps at 125, 126, 129, 130 (cols 0..7) while leaving
    # clean non-strap flanks on both sides of bleed_col = 12 (RAWX 137)!
    bundle = tpf_builder(
        tmp_path / "alt_pulsation",
        col0_1based=125,
        row0_1based=1985,
        nx=24,
        ny=64,
        n_cadences=24,
        bleed_in_smrow=True,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        pulsation_amp=0.025,
        sky_dark_level=50.0,
        star_scale=250.0,
        seed=88,
    )
    cal = recalibrate_tpf(bundle.tpf_path, ddir=bundle.ddir, smear="automatic", offline=True)
    pq = read_pixel_quality(cal)
    assert np.any((pq & int(PixelQuality.STRAP_COL)) != 0)
    assert cal.meta["BKGSUB"] is False
    assert np.all(cal.hdu[1].data["FLUX_BKG"] == 0.0)

    # Measure fractional pulsation amplitude recovery on mid-halo rows of the bleed column!
    bleed_col = 12
    halo_row = 28  # clean unsaturated wing row below bleed (bleed starts at row 38)
    assert not bundle.bleed_mask[halo_row, bleed_col]
    dark_slice = slice(0, 12)  # constant sky pedestal -> pure pulsation + transit modulation
    rec_curve = cal.flux.value[dark_slice, halo_row, bleed_col]
    true_curve = (bundle.true_star_cube + bundle.true_sky_cube)[dark_slice, halo_row, bleed_col]

    # Fit linear slope of recovered vs true modulation to check fractional self-subtraction
    tc = true_curve - np.mean(true_curve)
    rc = rec_curve - np.mean(rec_curve)
    amp_ratio = float(np.dot(rc, tc) / np.dot(tc, tc))
    assert abs(amp_ratio - 1.0) < 2e-3  # < 0.2% self-subtraction!

    # 2. Boundary-clamped bleed at col 0 under a horizontal gradient must inflate SMEAR_ERR[:, 0]
    bundle_edge = tpf_builder(
        tmp_path / "alt_edge",
        col0_1based=240,
        row0_1based=2001,
        nx=10,
        ny=48,
        n_cadences=4,
        bleed_in_smrow=False,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        sky_dark_level=800.0,
        seed=91,
    )
    # Inject saturated bleed at boundary col 0 (rows 30..47) AND corrupt smrow col 0
    with fits.open(bundle_edge.tpf_path, mode="update") as hdul:
        hdul[1].data["RAW_CNTS"][:, 30:, 0] = 12_000_000
    sm_path = next(bundle_edge.ddir.glob("*-smrow-*-s_col.fits"))
    sm_c0 = 240 - 1 - 44
    with fits.open(sm_path, mode="update") as sh:
        sh[1].data["SMROW_RAW"][:, :, sm_c0] = 12_000_000

    cal_edge = recalibrate_tpf(
        bundle_edge.tpf_path, ddir=bundle_edge.ddir, smear="alternate", offline=True
    )
    serr_edge = read_smear_err_column(cal_edge, quality_masked=False)
    assert float(np.median(serr_edge[:, 0])) > 2.0 * float(np.median(serr_edge[:, 5]))

    # 3. Bleed leaving < 5 clean rows after 4-row dilation must set PixelQuality.SMEAR_UNRESOLVED
    bundle_unres = tpf_builder(
        tmp_path / "alt_unres",
        col0_1based=260,
        row0_1based=2033,
        nx=10,
        ny=16,
        n_cadences=4,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        seed=92,
    )
    with fits.open(bundle_unres.tpf_path, mode="update") as hdul:
        # Saturate rows 6..15 -> after 4-row dilation, rows 2..15 are masked -> only 2 clean rows!
        hdul[1].data["RAW_CNTS"][:, 6:, 4] = 12_000_000
    cal_unres = recalibrate_tpf(
        bundle_unres.tpf_path, ddir=bundle_unres.ddir, smear="alternate", offline=True
    )
    pq_unres = read_pixel_quality(cal_unres)
    assert np.all((pq_unres[:, 4] & int(PixelQuality.SMEAR_UNRESOLVED)) != 0)
    assert np.all((pq_unres[:, 0] & int(PixelQuality.SMEAR_UNRESOLVED)) == 0)


def test_desmear_calibrated_cutout_tesscut_and_bkg_restoration(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """desmear_calibrated_cutout must repair FFI cutouts without RAW_CNTS or _col.fits."""
    bundle = tpf_builder(
        tmp_path / "ffi_cutout",
        col0_1based=300,
        row0_1based=1993,
        nx=12,
        ny=56,
        n_cadences=8,
        write_collateral=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        sky_dark_level=120.0,
        seed=105,
    )

    truth = bundle.true_star_cube + bundle.true_sky_cube
    bleed_col = 6
    # Build TessCut/HLSP-like FITS without RAW_CNTS where col 6 has -1,850 e-/s smear bloom
    with fits.open(bundle.tpf_path) as src:
        hdu0 = src[0].copy()
        hdu2 = src[2].copy()
        kept = [c for c in src[1].columns if c.name != "RAW_CNTS"]
        hdu1 = fits.BinTableHDU.from_columns(
            fits.ColDefs(kept), header=src[1].header, name="PIXELS"
        )

        corrupted = truth.copy()
        corrupted[:, :, bleed_col] -= 1_850.0 * bundle.true_signal[:, None]
        fake_bkg = bundle.true_sky_cube.astype(np.float32)
        hdu1.data["FLUX"] = (corrupted - fake_bkg).astype(np.float32)
        hdu1.data["FLUX_BKG"] = fake_bkg
        hdu1.data["FLUX_BKG_ERR"] = np.full_like(fake_bkg, 0.5)
        hdu1.header["BACKAPP"] = True

        ffi_path = tmp_path / "synthetic_tesscut.fits"
        fits.HDUList([hdu0, hdu1, hdu2]).writeto(ffi_path, overwrite=True)

    cutout_in = lk.read(ffi_path)
    assert "RAW_CNTS" not in (cutout_in.hdu[1].data.names or ())

    repaired = desmear_calibrated_cutout(cutout_in)
    assert repaired.meta["HALOCVER"] == 1
    assert repaired.meta["BKGSUB"] is False
    assert repaired.hdu[1].header["BACKAPP"] is False
    assert np.all(repaired.hdu[1].data["FLUX_BKG"] == 0.0)

    unsat = np.broadcast_to(~bundle.bleed_mask, repaired.flux.shape)
    np.testing.assert_allclose(
        repaired.flux.value[unsat],
        truth[unsat],
        rtol=5e-3,
        atol=2.0,
    )
    serr = read_smear_err_column(repaired)
    assert serr.shape == (8, 12)
    assert np.all(np.isfinite(serr))


def test_cadenceno_subset_thermal_step_and_blf_cache_sharing(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """QUALITY-filtered TPFs (N_tpf < N_col) across a thermal step must pair by CADENCENO."""
    base_dir = tmp_path / "thermal_sector"
    # Inject a steep post-gap thermal step (+400 ADU 1D black, +90 e-/frame smear) at cads >= 16
    bundle = tpf_builder(
        base_dir,
        col0_1based=210,
        row0_1based=1500,
        nx=8,
        ny=10,
        n_cadences=32,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        sky_dark_level=350.0,
        thermal_step_cadence=16,
        thermal_step_b1d_adu=400.0,
        thermal_step_smear_e=90.0,
        seed=404,
    )
    sm_path = next(base_dir.glob("*-smrow-*-s_col.fits"))
    tv_path = next(base_dir.glob("*-tvcol-*-s_col.fits"))

    # Write two differently QUALITY-filtered subset TPFs (N_tpf < N_col = 32):
    # TPF A drops cadences 0..1 and 12..19 (22 rows); TPF B drops cadences 5..12 (24 rows)
    keep_a = np.r_[2:12, 20:32]
    keep_b = np.r_[0:5, 13:32]
    tpf_a_path = tmp_path / "tess2020001000000-s0020-0000000100000001-0165-s_tp.fits"
    tpf_b_path = tmp_path / "tess2020001000000-s0020-0000000100000002-0165-s_tp.fits"

    for subset_idx, dst_p in ((keep_a, tpf_a_path), (keep_b, tpf_b_path)):
        with fits.open(bundle.tpf_path) as src:
            h0 = src[0].copy()
            h1 = fits.BinTableHDU(data=src[1].data[subset_idx], header=src[1].header, name="PIXELS")
            h2 = src[2].copy()
            fits.HDUList([h0, h1, h2]).writeto(dst_p, overwrite=True)

    shared_cache = tmp_path / "shared_cache"
    cal_a = recalibrate_tpf(
        tpf_a_path,
        ddir=base_dir,
        cache_dir=shared_cache,
        smear="standard",
        cadence_match="cadenceno",
        offline=True,
    )
    # Remove raw _col.fits before TPF B to prove TPF B reuses the 32-row -blf.fits cache!
    sm_path.unlink()
    tv_path.unlink()

    cal_b = recalibrate_tpf(
        tpf_b_path,
        ddir=base_dir,
        cache_dir=shared_cache,
        smear="standard",
        cadence_match="cadenceno",
        offline=True,
    )

    truth_full = bundle.true_star_cube + bundle.true_sky_cube
    non_cr_a = bundle.cr_cube[keep_a] == 0.0
    non_cr_b = bundle.cr_cube[keep_b] == 0.0

    # Zero thermal step artefact across the dropped cadence gap for BOTH targets!
    np.testing.assert_allclose(
        cal_a.flux.value[non_cr_a],
        truth_full[keep_a][non_cr_a],
        rtol=5e-4,
        atol=0.08,
    )
    np.testing.assert_allclose(
        cal_b.flux.value[non_cr_b],
        truth_full[keep_b][non_cr_b],
        rtol=5e-4,
        atol=0.08,
    )

    # Verify -blf.fits stored all 32 collateral CADENCENO rows (not 22)
    blfs = sorted((shared_cache / "s0020").glob("*-blf.fits"))
    assert len(blfs) == 2
    for bp in blfs:
        with fits.open(bp) as bh:
            assert len(bh[1].data["CADENCENO"]) == 32

    # Mismatched sector CADENCENO (n_matched == 0) must raise clear ValueError
    bad_tpf = tmp_path / "mismatched_cad.fits"
    with fits.open(tpf_a_path) as src:
        h0, h1, h2 = src[0].copy(), src[1].copy(), src[2].copy()
        h1.data["CADENCENO"] = np.arange(900_000, 900_000 + len(keep_a), dtype=np.int32)
        fits.HDUList([h0, h1, h2]).writeto(bad_tpf, overwrite=True)

    with pytest.raises(ValueError, match="Zero matching CADENCENO"):
        recalibrate_tpf(bad_tpf, ddir=base_dir, cache_dir=shared_cache, offline=True)


def test_masked_onedblack_nan_gap_avoids_gibbs_ringing_and_uses_n_finite() -> None:
    """fit_onedblack_legendre must solve masked normal equations across interior NaN row gaps."""
    rows = np.arange(2078, dtype=np.float64)
    x_norm = (rows[:2048] - 1023.5) / 1023.5
    true_sci = (
        180.0
        + 35.0 * x_norm
        - 18.0 * (x_norm**2)
        + 9.5 * (x_norm**5)
        - 4.2 * (x_norm**8)
    )
    true_ov = np.full(30, 192.5, dtype=np.float64)
    true_profile = np.concatenate([true_sci, true_ov])

    # 1. Noiseless 120-row interior NaN dropout (rows 800..919) + 5 overclock NaNs
    tv_noiseless = np.broadcast_to(true_profile[None, :, None], (2, 2078, 11)).copy()
    tv_noiseless[1, 800:920, :] = np.nan
    tv_noiseless[1, 2050:2055, :] = np.nan

    b1d_clean, err_clean = fit_onedblack_legendre(tv_noiseless)
    # Zero Gibbs ringing across ALL 2,048 rows (including inside the 120-row NaN gap!)
    np.testing.assert_allclose(b1d_clean[0], true_profile, atol=1e-7)
    np.testing.assert_allclose(b1d_clean[1], true_profile, atol=1e-7)
    assert np.all(np.isfinite(err_clean))

    # 2. With Gaussian noise, uncertainty inside the unobserved gap exceeds observed rows
    rng = np.random.default_rng(2026)
    tv_noisy = tv_noiseless + rng.normal(0.0, 0.8, size=tv_noiseless.shape)
    tv_noisy[1, 800:920, :] = np.nan
    tv_noisy[1, 2050:2055, :] = np.nan
    _, err_noisy = fit_onedblack_legendre(tv_noisy)
    gap_err = float(np.median(err_noisy[1, 830:890]))
    obs_err = float(np.median(err_noisy[1, 200:600]))
    assert gap_err > 1.2 * obs_err


def test_tall_zero_sky_column_not_falsely_flagged_as_bloom(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """Tall (ny=1287) zero-sky noise columns must not trigger extreme-value bloom false alarms."""
    rng = np.random.default_rng(2026)
    sigma_cube = np.full((40, 1287, 6), 5.0, dtype=np.float64)
    sub_cube = rng.normal(0.0, 5.0, size=sigma_cube.shape)
    clean_rows = np.ones((1287, 6), dtype=bool)

    flagged = _detect_bloom_columns(sub_cube, sigma_cube, clean_rows)
    assert not np.any(flagged)

    # End-to-end desmear_calibrated_cutout on a tall zero-sky unbloomed FFI cutout
    bundle = tpf_builder(
        tmp_path / "tall_zerosky",
        col0_1based=300,
        row0_1based=762,
        nx=6,
        ny=1287,
        n_cadences=16,
        write_collateral=False,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=2.0,
        sky_dark_level=0.0,
        star_scale=0.0,
        seed=42,
    )
    with fits.open(bundle.tpf_path) as src:
        h0, h2 = src[0].copy(), src[2].copy()
        kept = [c for c in src[1].columns if c.name != "RAW_CNTS"]
        h1 = fits.BinTableHDU.from_columns(
            fits.ColDefs(kept), header=src[1].header, name="PIXELS"
        )
        zero_sky_flux = rng.normal(0.0, 2.0, size=h1.data["FLUX"].shape).astype(np.float32)
        h1.data["FLUX"] = zero_sky_flux
        h1.data["FLUX_ERR"] = np.full_like(zero_sky_flux, 2.0)
        h1.data["FLUX_BKG"] = np.zeros_like(zero_sky_flux)
        h1.data["FLUX_BKG_ERR"] = np.zeros_like(zero_sky_flux)
        h1.header["BACKAPP"] = False
        cut_path = tmp_path / "tall_zerosky_cutout.fits"
        fits.HDUList([h0, h1, h2]).writeto(cut_path, overwrite=True)

    repaired = desmear_calibrated_cutout(cut_path)
    np.testing.assert_allclose(repaired.flux.value, zero_sky_flux, atol=1e-5)


def test_partial_multioutput_alt_cache_no_nans_and_warm_pixel_quality_parity(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """Warm -alt- reload must preserve SMEAR_UNRESOLVED & partial multi-output cache has no NaNs."""
    bundle = tpf_builder(
        tmp_path / "partial_2out",
        col0_1based=553,
        row0_1based=2033,
        nx=8,  # cols 553..556 (Output A) & 557..560 (Output B)
        ny=16,
        n_cadences=4,
        bleed_in_smrow=True,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        seed=612,
    )
    assert bundle.outputs == ("A", "B")
    # Saturate Output A col 1 rows 6..15 -> < 5 clean rows after 4-row dilation -> SMEAR_UNRESOLVED
    with fits.open(bundle.tpf_path, mode="update") as hdul:
        hdul[1].data["RAW_CNTS"][:, 6:, 1] = 12_000_000

    cache_dir = tmp_path / "alt_cache"
    cal_cold = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=cache_dir,
        smear="alternate",
        offline=True,
    )
    pq_cold = read_pixel_quality(cal_cold)
    assert np.all((pq_cold[:, 1] & int(PixelQuality.SMEAR_UNRESOLVED)) != 0)

    # 1. Full warm reload from both A & B -alt- caches preserves PIXEL_QUALITY bit-for-bit!
    cal_warm = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=cache_dir,
        smear="alternate",
        offline=True,
    )
    np.testing.assert_array_equal(read_pixel_quality(cal_warm), pq_cold)
    np.testing.assert_allclose(cal_warm.flux.value, cal_cold.flux.value, rtol=1e-5, atol=1e-4)

    # 2. Delete ONLY Output B's -alt- cache -> partial hit must still load standard smrow for all
    # outputs and calibrate with zero NaNs and identical flux!
    out_b_alt = next((cache_dir / "s0020").glob("*-smrow-2-3-b-*-alt-s_col-blf.fits"))
    out_b_alt.unlink()

    cal_partial = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=cache_dir,
        smear="alternate",
        offline=True,
    )
    assert np.all(np.isfinite(cal_partial.flux.value))
    assert np.all(np.isfinite(read_smear_err_column(cal_partial)))
    np.testing.assert_array_equal(read_pixel_quality(cal_partial), pq_cold)
    np.testing.assert_allclose(cal_partial.flux.value, cal_cold.flux.value, rtol=1e-5, atol=1e-4)


def test_upper_register_bloom_stellar_wing_smear_and_desmear_exposure_scaling(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """All-column smrow bloom subtracts in-cutout stellar smear; desmear scales with exposure."""
    # Part A: Narrow 4-column stamp with bright unsaturated stellar wing & all-column smrow bloom
    bundle = tpf_builder(
        tmp_path / "upper_bloom",
        col0_1based=220,
        row0_1based=1993,
        nx=4,
        ny=56,
        n_cadences=4,
        include_bleed=False,
        bleed_in_smrow=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        sky_dark_level=120.0,
        star_scale=2500.0,
        seed=718,
    )
    sm_path = next(bundle.ddir.glob("*-smrow-*-s_col.fits"))
    sm_c0 = 220 - 1 - 44
    with fits.open(sm_path, mode="update") as sh:
        sh[1].data["SMROW_RAW"][:, :, sm_c0 : sm_c0 + 4] = 12_000_000

    cache_dir = tmp_path / "ub_cache"
    cal_cold = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=cache_dir,
        smear="alternate",
        offline=True,
    )
    # Delete -alt- cache so second run reloads standard smrow -blf.fits + SMROW_SAT extension!
    for alt_p in (cache_dir / "s0020").glob("*-alt-s_col-blf.fits"):
        alt_p.unlink()
    cal_warm_sat = recalibrate_tpf(
        bundle.tpf_path,
        ddir=bundle.ddir,
        cache_dir=cache_dir,
        smear="alternate",
        offline=True,
    )
    np.testing.assert_allclose(
        cal_warm_sat.flux.value, cal_cold.flux.value, rtol=1e-5, atol=1e-4
    )
    np.testing.assert_allclose(
        read_smear_err_column(cal_warm_sat),
        read_smear_err_column(cal_cold),
        rtol=1e-5,
        atol=1e-4,
    )

    truth = bundle.true_star_cube + bundle.true_sky_cube
    non_cr = bundle.cr_cube == 0.0
    np.testing.assert_allclose(cal_cold.flux.value[non_cr], truth[non_cr], rtol=8e-3, atol=2.5)

    # Part B: desmear_calibrated_cutout without FLUX_ERR scales Poisson fallback with exposure time
    # and scales SMEAR_UNRESOLVED uncertainty covariantly with flux units!
    bundle_b = tpf_builder(
        tmp_path / "desmear_exp_src",
        col0_1based=220,
        row0_1based=1993,
        nx=8,
        ny=56,
        n_cadences=4,
        write_collateral=False,
        include_bleed=False,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        seed=719,
    )

    def _make_errless_cutout(path_out: Path, num_frm: int, flux_scale: float) -> Path:
        with fits.open(bundle_b.tpf_path) as src:
            h0, h2 = src[0].copy(), src[2].copy()
            kept = [c for c in src[1].columns if c.name not in ("RAW_CNTS", "FLUX_ERR")]
            h1 = fits.BinTableHDU.from_columns(
                fits.ColDefs(kept), header=src[1].header, name="PIXELS"
            )
            h0.header["NUM_FRM"] = num_frm
            h0.header["CRMITEN"] = False
            h1.header["NUM_FRM"] = num_frm
            h1.header["CRMITEN"] = False
            base_flux = np.full(h1.data["FLUX"].shape, 200.0 * flux_scale, dtype=np.float32)
            # Col 1 has symmetric clean flanks (col 0 & cols 2..4) -> pure shot-noise SMEAR_ERR
            base_flux[:, :, 1] = np.float32(-500.0 * flux_scale)
            # Col 5 is an unresolved saturated column (all 56 rows > sat limit)
            base_flux[:, :, 5] = np.float32(2.0e6 * flux_scale)
            h1.data["FLUX"] = base_flux
            h1.data["FLUX_BKG"] = np.zeros_like(base_flux)
            h1.data["FLUX_BKG_ERR"] = np.zeros_like(base_flux)
            h1.header["BACKAPP"] = False
            fits.HDUList([h0, h1, h2]).writeto(path_out, overwrite=True)
        return path_out

    cut_20s = _make_errless_cutout(tmp_path / "cut_20s.fits", num_frm=10, flux_scale=1.0)
    cut_30m = _make_errless_cutout(tmp_path / "cut_30m.fits", num_frm=900, flux_scale=1.0)
    cut_20s_10x = _make_errless_cutout(tmp_path / "cut_20s_10x.fits", num_frm=10, flux_scale=10.0)

    rep_20s = desmear_calibrated_cutout(cut_20s)
    rep_30m = desmear_calibrated_cutout(cut_30m)
    rep_20s_10x = desmear_calibrated_cutout(cut_20s_10x, saturation_flux_threshold=1.0e6)

    serr_20s = read_smear_err_column(rep_20s, quality_masked=False)
    serr_30m = read_smear_err_column(rep_30m, quality_masked=False)
    serr_20s_10x = read_smear_err_column(rep_20s_10x, quality_masked=False)

    # Active-repair bloom column 1: 20-s exposure has ~sqrt(900/10) = 9.5x larger shot-noise error!
    assert np.all(serr_20s[:, 0] == 0.0)
    ratio_exp = float(np.median(serr_20s[:, 1]) / np.median(serr_30m[:, 1]))
    assert 8.0 < ratio_exp < 11.0

    # Unresolved column 5 scales linearly (10x) with flux scale (zero magic additive floor!)
    pq_20s = read_pixel_quality(rep_20s)
    assert np.all((pq_20s[:, 5] & int(PixelQuality.SMEAR_UNRESOLVED)) != 0)
    ratio_scale = float(np.median(serr_20s_10x[:, 5]) / np.median(serr_20s[:, 5]))
    np.testing.assert_allclose(ratio_scale, 10.0, rtol=0.05)
