"""Tests for synthetic 1/2-output fixtures, calibrate provenance, and -blf.fits schema."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from mundey import mundey_tpf


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

    # 2-output fixture must straddle Outputs A (<=556) and B (>=557)
    assert synthetic_2output_bundle.outputs == ("A", "B")
    cols = synthetic_2output_bundle.columns_1based
    assert 556 in cols and 557 in cols
    assert synthetic_2output_tpf.column == 551
    assert synthetic_2output_tpf.shape[1:] == (16, 12)

    # Sky sweeps dark (<10 e-/s) and bright (>600 e-/s) regimes
    assert float(np.min(synthetic_2output_bundle.true_sky_cube)) < 10.0
    assert float(np.max(synthetic_2output_bundle.true_sky_cube)) > 600.0
    # Mid-sector gap > 1 day present
    time_vals = synthetic_2output_tpf.time.value
    assert float(np.max(np.diff(time_vals))) > 1.0


def test_calibrate_1output_recovers_truth_sets_provenance_and_fixes_blf_schema(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """calibrate() on 1-output TPF must quantitatively recover truth and fix -blf schema."""
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

    unsat_non_cr = (~bundle.bleed_mask[None, :, :]) & (bundle.cr_cube == 0.0)
    expected = bundle.true_star_cube + bundle.true_sky_cube
    np.testing.assert_allclose(
        tpf.flux.value[unsat_non_cr],
        expected[unsat_non_cr],
        rtol=2e-3,
        atol=0.35,
    )

    # Inspect cached -blf.fits files: column 2 must be CADENCENO ('J'), NOT TIMECORR ('E')
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
    """calibrate(smear='standard') across Outputs A/B must quantitatively recover flux."""
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
    tpf = mundey_tpf(str(bundle.tpf_path))
    tpf.calibrate(verbose=False, smear="standard", ddir=bundle.ddir, offline=True)

    unsat_non_cr = (~bundle.bleed_mask[None, :, :]) & (bundle.cr_cube == 0.0)
    expected = bundle.true_star_cube + bundle.true_sky_cube
    np.testing.assert_allclose(
        tpf.flux.value[unsat_non_cr],
        expected[unsat_non_cr],
        rtol=2e-3,
        atol=0.35,
    )
    assert np.all(tpf.hdu[1].data["FLUX_BKG"] == 0.0)
    assert tpf.hdu[1].header["BACKAPP"] is False
    assert tpf.meta["BKGSUB"] is False
    assert tpf.meta["HALOCVER"] == 1
    assert len(list(bundle.ddir.glob("*-s_col-blf.fits"))) == 4


def test_calibrate_2output_alternate_smear_and_cache_isolation(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """2-output smear='automatic' (row 2048 -> 'alternate') must recover flux and isolate caches."""
    shared_ddir = tmp_path / "shared_ddir"
    # row0_1based=1993, ny=56 -> top row is 2048 -> smear='automatic' selects 'alternate'
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

    tpf1 = mundey_tpf(str(bundle1.tpf_path))
    tpf1.calibrate(verbose=False, smear="automatic", ddir=shared_ddir, offline=True)

    # Dark-sky cadences (first half) have negligible 2D sky spread across columns (< 1.5 e-/s),
    # so alternate smear recovers unsaturated flux across BOTH Output A and Output B
    dark_cads = slice(0, 2)
    unsat_dark = (~bundle1.bleed_mask[None, :, :]) & (bundle1.cr_cube[dark_cads] == 0.0)
    expected_dark = (bundle1.true_star_cube + bundle1.true_sky_cube)[dark_cads]
    np.testing.assert_allclose(
        tpf1.flux.value[dark_cads][unsat_dark],
        expected_dark[unsat_dark],
        rtol=5e-3,
        atol=2.5,
    )

    alt_caches = sorted(shared_ddir.glob("*-0000000111111111-*-alt-s_col-blf.fits"))
    assert len(alt_caches) == 2

    # Running smear='standard' in the same ddir must NOT reuse '-alt-s_col-blf.fits':
    # collateral bleed corruption in SMROW_RAW drives the bleed column deeply negative.
    tpf1_std = mundey_tpf(str(bundle1.tpf_path))
    tpf1_std.calibrate(verbose=False, smear="standard", ddir=shared_ddir, offline=True)
    bleed_col = 6
    assert float(np.min(tpf1_std.flux.value[:, :20, bleed_col])) < -50_000.0

    # Now both 'standard' (-s_col-blf.fits) and target 1's 'alternate' caches exist in shared_ddir.
    # Calibrating a second target at different columns must NOT be poisoned by either cache.
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
    tpf2 = mundey_tpf(str(bundle2.tpf_path))
    tpf2.calibrate(verbose=False, smear="automatic", ddir=shared_ddir, offline=True)

    unsat_dark2 = (~bundle2.bleed_mask[None, :, :]) & (bundle2.cr_cube[dark_cads] == 0.0)
    expected_dark2 = (bundle2.true_star_cube + bundle2.true_sky_cube)[dark_cads]
    np.testing.assert_allclose(
        tpf2.flux.value[dark_cads][unsat_dark2],
        expected_dark2[unsat_dark2],
        rtol=5e-3,
        atol=2.5,
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
    with pytest.raises(ValueError, match="outside TESS science outputs 45..2092"):
        tpf.calibrate(verbose=False, ddir=synthetic_1output_bundle.ddir, offline=True)


def test_sirius_sector6_golden_checksum_and_contract(tmp_path: Path) -> None:
    """Verify Sirius Sector 6 100-cadence golden checksum offline in an isolated temp dir."""
    notebooks_dir = Path(__file__).resolve().parent.parent / "notebooks"
    for ref_file in notebooks_dir.iterdir():
        if ref_file.is_file() and (ref_file.suffix in (".fits", ".xml")):
            (tmp_path / ref_file.name).symlink_to(ref_file)

    tpffile = tmp_path / "tess2018349182459-s0006-0000000322899250-0126-s_tp.fits"
    tpf = mundey_tpf(str(tpffile))
    # OBJECT = 'TIC 322899250' is in hdu[0].header, so use_simbad defaults to False automatically
    tpf.calibrate(verbose=False, ddir=tmp_path)

    flux_sum = np.nansum(tpf.flux.value)
    assert np.abs(flux_sum - np.float32(41293316000.0)) < 1.0
    assert np.all(tpf.hdu[1].data["FLUX_BKG"] == 0.0)
    assert np.all(tpf.hdu[1].data["FLUX_BKG_ERR"] == 0.0)
    assert tpf.meta["HALOCVER"] == 1
    assert tpf.meta["BKGSUB"] is False
