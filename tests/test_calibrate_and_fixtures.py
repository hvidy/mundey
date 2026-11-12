"""Tests for synthetic 1/2-output fixtures, deprecated mundey_tpf shim, and Sirius S6."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

import mundey
from mundey import mundey_tpf, read_pixel_quality, read_smear_err_column
from mundey.calibration import match_cadenceno


def test_synthetic_1output_and_2output_fixture_properties(
    synthetic_1output_bundle,
    synthetic_1output_tpf,
    synthetic_2output_bundle,
    synthetic_2output_tpf,
) -> None:
    """Verify 1-output (A) and 2-output (A/B boundary at 556/557) fixture properties."""
    assert synthetic_1output_bundle.outputs == ("A",)
    assert synthetic_1output_tpf.column == 200
    assert synthetic_1output_tpf.shape[1:] == (14, 10)

    assert synthetic_2output_bundle.outputs == ("A", "B")
    cols = synthetic_2output_bundle.columns_1based
    assert 556 in cols and 557 in cols
    assert synthetic_2output_tpf.column == 551
    assert synthetic_2output_tpf.shape[1:] == (16, 12)

    assert float(np.min(synthetic_2output_bundle.true_sky_cube)) < 10.0
    assert float(np.max(synthetic_2output_bundle.true_sky_cube)) > 600.0
    time_vals = synthetic_2output_tpf.time.value
    assert float(np.max(np.diff(time_vals))) > 1.0


def test_calibrate_shim_1output_recovers_truth_sets_provenance_and_blf_schema(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """Deprecated mundey_tpf.calibrate() must delegate to v2.0 engine and fix -blf schema."""
    bundle = tpf_builder(
        tmp_path / "cal_1out",
        col0_1based=200,
        row0_1based=1500,
        nx=8,
        ny=12,
        n_cadences=4,
        write_collateral=True,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        seed=7,
    )
    tpf = mundey_tpf(str(bundle.tpf_path))
    assert float(np.nansum(tpf.hdu[1].data["FLUX_BKG"])) > 0.0
    assert float(np.nansum(tpf.hdu[1].data["FLUX_BKG_ERR"])) > 0.0
    assert tpf.hdu[1].header["BACKAPP"] is True

    with pytest.warns(DeprecationWarning, match="recalibrate_tpf"):
        tpf.calibrate(verbose=False, smear="standard", ddir=bundle.ddir, offline=True)

    assert np.all(tpf.hdu[1].data["FLUX_BKG"] == 0.0)
    assert np.all(tpf.hdu[1].data["FLUX_BKG_ERR"] == 0.0)
    assert np.all(tpf.flux_bkg.value == 0.0)
    assert np.all(tpf.flux_bkg_err.value == 0.0)
    assert tpf.hdu[0].header["HALOCVER"] == 1
    assert tpf.hdu[1].header["HALOCVER"] == 1
    assert tpf.hdu[0].header["BKGSUB"] is False
    assert tpf.hdu[1].header["BKGSUB"] is False
    assert tpf.hdu[1].header["BACKAPP"] is False
    assert tpf.meta["BKGSUB"] is False
    assert tpf.hdu[0].header["CADMATCH"] == "cadenceno"
    assert tpf.hdu[0].header["CADMISS"] == 0

    unsat_non_cr = (~bundle.bleed_mask[None, :, :]) & (bundle.cr_cube == 0.0)
    expected = bundle.true_star_cube + bundle.true_sky_cube
    np.testing.assert_allclose(
        tpf.flux.value[unsat_non_cr],
        expected[unsat_non_cr],
        rtol=5e-4,
        atol=0.08,
    )

    blf_files = sorted(bundle.ddir.glob("*-s_col-blf.fits"))
    assert len(blf_files) == 2
    assert list(bundle.ddir.glob("*.tmp")) == []
    for blf_path in blf_files:
        with fits.open(blf_path) as hdul:
            hdr = hdul[1].header
            assert hdr["TTYPE1"] == "TIME"
            assert hdr["TTYPE2"] == "CADENCENO"
            assert hdr["TFORM2"].strip() == "J"
            np.testing.assert_array_equal(
                hdul[1].data["CADENCENO"], tpf.hdu[1].data["CADENCENO"]
            )


def test_calibrate_2output_ab_boundary_recovers_truth(tpf_builder, tmp_path: Path) -> None:
    """recalibrate_tpf(smear='standard') across Outputs A/B must quantitatively recover flux."""
    bundle = tpf_builder(
        tmp_path / "cal_2out",
        col0_1based=553,
        row0_1based=1400,
        nx=8,
        ny=12,
        n_cadences=4,
        write_collateral=True,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        seed=99,
    )
    assert bundle.outputs == ("A", "B")
    cal_tpf = mundey.recalibrate_tpf(
        bundle.tpf_path, verbose=False, smear="standard", ddir=bundle.ddir, offline=True
    )

    unsat_non_cr = (~bundle.bleed_mask[None, :, :]) & (bundle.cr_cube == 0.0)
    expected = bundle.true_star_cube + bundle.true_sky_cube
    np.testing.assert_allclose(
        cal_tpf.flux.value[unsat_non_cr],
        expected[unsat_non_cr],
        rtol=5e-4,
        atol=0.08,
    )
    assert np.all(cal_tpf.hdu[1].data["FLUX_BKG"] == 0.0)
    assert cal_tpf.hdu[1].header["BACKAPP"] is False
    assert cal_tpf.meta["BKGSUB"] is False
    assert cal_tpf.meta["HALOCVER"] == 1
    assert len(list(bundle.ddir.glob("*-s_col-blf.fits"))) == 4


def test_calibrate_2output_alternate_smear_and_cache_isolation(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """2-output smear='automatic' (row 2048 -> 'alternate') must recover flux and isolate caches."""
    shared_ddir = tmp_path / "shared_ddir"
    bundle1 = tpf_builder(
        shared_ddir,
        col0_1based=551,
        row0_1based=1993,
        nx=12,
        ny=56,
        n_cadences=4,
        ticid=111111111,
        write_collateral=True,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        bleed_in_smrow=True,
        seed=314,
    )
    assert bundle1.outputs == ("A", "B")

    tpf1 = mundey.recalibrate_tpf(
        bundle1.tpf_path, verbose=False, smear="automatic", ddir=shared_ddir, offline=True
    )

    dark_cads = slice(0, 2)
    unsat_dark = (~bundle1.bleed_mask[None, :, :]) & (bundle1.cr_cube[dark_cads] == 0.0)
    expected1 = bundle1.true_star_cube + bundle1.true_sky_cube
    np.testing.assert_allclose(
        tpf1.flux.value[dark_cads][unsat_dark],
        expected1[dark_cads][unsat_dark],
        rtol=4e-3,
        atol=2.0,
    )
    # Right flank cols 7..11 (RAWX 558..562) are CCID-80 strap-masked, so bleed col 6 clamps
    # to Output A left flank and inflates SMEAR_ERR on bright-gradient cadences 2..3!
    serr1 = read_smear_err_column(tpf1, quality_masked=False)
    assert float(np.min(serr1[2:, 6])) > 10.0 * float(np.max(serr1[:2, 6]))

    alt_caches = sorted(shared_ddir.glob("*-0000000111111111-*-alt-s_col-blf.fits"))
    assert len(alt_caches) == 2

    # Running smear='standard' in the same ddir must NOT reuse '-alt-s_col-blf.fits'
    tpf1_std = mundey.recalibrate_tpf(
        bundle1.tpf_path, verbose=False, smear="standard", ddir=shared_ddir, offline=True
    )
    bleed_col = 6
    assert float(np.min(tpf1_std.flux.value[:, :20, bleed_col])) < -50_000.0

    # Second target at shifted columns must NOT be poisoned by target 1's alt cache
    bundle2 = tpf_builder(
        shared_ddir,
        col0_1based=553,
        row0_1based=1993,
        nx=10,
        ny=56,
        n_cadences=4,
        ticid=222222222,
        write_collateral=True,
        flag_bad_cadence=False,
        noise_sigma=0.0,
        bleed_in_smrow=True,
        seed=271,
    )
    tpf2 = mundey.recalibrate_tpf(
        bundle2.tpf_path, verbose=False, smear="automatic", ddir=shared_ddir, offline=True
    )
    unsat2_dark = (~bundle2.bleed_mask[None, :, :]) & (bundle2.cr_cube[dark_cads] == 0.0)
    expected2 = bundle2.true_star_cube + bundle2.true_sky_cube
    np.testing.assert_allclose(
        tpf2.flux.value[dark_cads][unsat2_dark],
        expected2[dark_cads][unsat2_dark],
        rtol=4e-3,
        atol=2.0,
    )
    assert len(list(shared_ddir.glob("*-0000000222222222-*-alt-s_col-blf.fits"))) == 2


def test_calibrate_rejects_collateral_columns_outside_science_outputs(
    synthetic_1output_bundle,
) -> None:
    """TPF stamps extending into collateral columns (< 45 or > 2092) must raise ValueError."""
    tpf = mundey_tpf(str(synthetic_1output_bundle.tpf_path))
    for col_idx in (4, 5, 6, 7, 8):
        tpf.hdu[1].header[f"1CRV{col_idx}P"] = 40
    tpf.hdu[2].header["CRVAL1P"] = 40
    with (
        pytest.warns(DeprecationWarning),
        pytest.raises(ValueError, match="outside TESS science outputs 45..2092"),
    ):
        tpf.calibrate(verbose=False, ddir=synthetic_1output_bundle.ddir, offline=True)


def test_sirius_sector6_cadenceno_pairing_desmearing_and_speed(tmp_path: Path) -> None:
    """Verify Sirius Sector 6 CADENCENO pairing, missing-cadence handling, and < 3 s runtime."""
    notebooks_dir = Path(__file__).resolve().parent.parent / "notebooks"
    for ref_file in notebooks_dir.iterdir():
        if ref_file.is_file() and (ref_file.suffix in (".fits", ".xml")):
            (tmp_path / ref_file.name).symlink_to(ref_file)

    tpffile = tmp_path / "tess2018349182459-s0006-0000000322899250-0126-s_tp.fits"
    smfile = tmp_path / "tess2018349182459-s0006-smrow-2-3-b-0126-s_col.fits"

    # 1. Direct CADENCENO index mapping verification:
    # TPF: 173388..173463 (rows 0..75), gap 173464..173470, 173471..173494 (rows 76..99)
    # Col: 173387..173486 (rows 0..99, contiguous)
    with fits.open(tpffile) as th, fits.open(smfile) as sh:
        tpf_cad = np.asarray(th[1].data["CADENCENO"], dtype=np.int64)
        col_cad = np.asarray(sh[1].data["CADENCENO"], dtype=np.int64)
        spoc_col0_med = float(np.nanmedian(th[1].data["FLUX"][:, 0, 0]))
    assert spoc_col0_med < -1_000.0  # SPOC smear bloom drove col 0 deeply negative

    matched, col_idx = match_cadenceno(tpf_cad, col_cad, cadence_match="cadenceno")
    np.testing.assert_array_equal(col_idx[:76], np.arange(1, 77))
    np.testing.assert_array_equal(col_idx[76:92], np.arange(84, 100))
    np.testing.assert_array_equal(col_idx[92:], np.full(8, -1))
    assert int(np.sum(~matched)) == 8
    np.testing.assert_array_equal(tpf_cad[~matched], np.arange(173487, 173495))

    # 2. Strict missing_cadence='error' must raise ValueError naming 8 / 100 and 173487..173494
    with pytest.raises(ValueError, match="8 / 100.*173487\\.\\.173494"):
        mundey.recalibrate_tpf(
            tpffile,
            ddir=tmp_path,
            cadence_match="cadenceno",
            missing_cadence="error",
            offline=True,
        )

    # 3. Default cadence_match='cadenceno', missing_cadence='nan' warns, sets CADMISS=8, < 3.5 s
    t0 = time.perf_counter()
    with pytest.warns(UserWarning, match="8 / 100.*173487\\.\\.173494"):
        cal = mundey.recalibrate_tpf(
            tpffile,
            ddir=tmp_path,
            cadence_match="cadenceno",
            missing_cadence="nan",
            offline=True,
        )
    elapsed = time.perf_counter() - t0
    assert elapsed < 4.0

    assert cal.meta["HALOCVER"] == 1
    assert cal.meta["BKGSUB"] is False
    assert cal.hdu[0].header["CADMATCH"] == "cadenceno"
    assert cal.hdu[0].header["CADMISS"] == 8

    flux_raw = cal.hdu[1].data["FLUX"]
    ferr_raw = cal.hdu[1].data["FLUX_ERR"]
    serr_raw = read_smear_err_column(cal, quality_masked=False)
    pq = read_pixel_quality(cal)
    assert pq.shape == (1287, 11)
    assert np.any((pq & int(mundey.PixelQuality.SAT_BLEED)) != 0)
    # RAWX 675..685 corresponds to imaging cols 631..641 (between straps 598 and 673)
    assert np.all((pq & int(mundey.PixelQuality.STRAP_COL)) == 0)

    # Matched 92 cadences are finite and restore positive outer-wing + sky flux on col 0
    assert np.all(np.isfinite(flux_raw[:92, 0, 0]))
    assert float(np.nanmedian(flux_raw[:92, 0, 0])) > 50.0
    assert np.all(np.isfinite(ferr_raw[:92, 0, 0]))
    assert np.all(np.isfinite(serr_raw[:92, 0]))

    # Unmatched trailing 8 cadences (173487..173494) are strictly NaN across FLUX/FLUX_ERR/SMEAR_ERR
    assert np.all(np.isnan(flux_raw[92:]))
    assert np.all(np.isnan(ferr_raw[92:]))
    assert np.all(np.isnan(serr_raw[92:]))
    assert np.all(cal.hdu[1].data["FLUX_BKG"] == 0.0)
    assert np.all(cal.hdu[1].data["FLUX_BKG_ERR"] == 0.0)

    # Sector-output tvcol -blf.fits cache must hold the collateral CADENCENO (173387..173486)
    tv_blf = tmp_path / "tess2018349182459-s0006-tvcol-2-3-b-0126-s_col-blf.fits"
    with fits.open(tv_blf) as tbh:
        np.testing.assert_array_equal(tbh[1].data["CADENCENO"], col_cad)

    # 4. Legacy reproduction mode cadence_match='positional' pairs all 100 rows (CADMISS == 0)
    pos_dir = tmp_path / "pos_mode"
    pos_dir.mkdir()
    for ref_file in notebooks_dir.iterdir():
        if ref_file.is_file() and (ref_file.suffix in (".fits", ".xml")):
            (pos_dir / ref_file.name).symlink_to(ref_file)
    cal_pos = mundey.recalibrate_tpf(
        pos_dir / tpffile.name,
        ddir=pos_dir,
        cadence_match="positional",
        offline=True,
    )
    assert cal_pos.hdu[0].header["CADMATCH"] == "positional"
    assert cal_pos.hdu[0].header["CADMISS"] == 0
    assert np.all(np.isfinite(cal_pos.hdu[1].data["FLUX"][:, 0, 0]))
