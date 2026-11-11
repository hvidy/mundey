"""Tests for synthetic 1/2-output fixtures, calibrate provenance, and -blf.fits schema."""

from __future__ import annotations

from pathlib import Path

import numpy as np
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


def test_calibrate_1output_zeroes_bkg_sets_provenance_and_fixes_blf_schema(
    tpf_builder,
    tmp_path: Path,
) -> None:
    """calibrate() must zero FLUX_BKG/ERR, set BKGSUB=False, and write CADENCENO ('J') in -blf."""
    bundle = tpf_builder(
        tmp_path / "cal_1out",
        col0_1based=200,
        row0_1based=1500,
        nx=8,
        ny=12,
        n_cadences=4,
        write_collateral=True,
        flag_bad_cadence=False,
        seed=7,
    )
    tpf = mundey_tpf(str(bundle.tpf_path))
    assert float(np.nansum(tpf.hdu[1].data["FLUX_BKG"])) > 0.0
    assert float(np.nansum(tpf.hdu[1].data["FLUX_BKG_ERR"])) > 0.0

    tpf.calibrate(verbose=False, smear="standard", ddir=bundle.ddir, offline=True)

    assert np.all(tpf.hdu[1].data["FLUX_BKG"] == 0.0)
    assert np.all(tpf.hdu[1].data["FLUX_BKG_ERR"] == 0.0)
    assert np.all(tpf.flux_bkg.value == 0.0)
    assert np.all(tpf.flux_bkg_err.value == 0.0)
    assert tpf.hdu[0].header["HALOCVER"] == 1
    assert tpf.hdu[1].header["HALOCVER"] == 1
    assert tpf.hdu[0].header["BKGSUB"] is False
    assert tpf.hdu[1].header["BKGSUB"] is False
    assert tpf.meta["BKGSUB"] is False
    assert np.all(np.isfinite(tpf.flux.value))

    # Inspect cached -blf.fits files: column 2 must be CADENCENO ('J'), NOT TIMECORR ('E')
    blf_files = sorted(bundle.ddir.glob("*-s_col-blf.fits"))
    assert len(blf_files) == 2
    for blf_path in blf_files:
        with fits.open(blf_path) as hdul:
            hdr = hdul[1].header
            assert hdr["TTYPE1"] == "TIME"
            assert hdr["TTYPE2"] == "CADENCENO"
            assert hdr["TFORM2"].strip() == "J"
            np.testing.assert_array_equal(
                hdul[1].data["CADENCENO"], tpf.hdu[1].data["CADENCENO"]
            )


def test_calibrate_2output_ab_boundary_runs_cleanly(tpf_builder, tmp_path: Path) -> None:
    """calibrate(smear='standard') on a 2-output A/B cutout must succeed without crashing."""
    bundle = tpf_builder(
        tmp_path / "cal_2out",
        col0_1based=553,
        row0_1based=1400,
        nx=8,
        ny=10,
        n_cadences=4,
        write_collateral=True,
        flag_bad_cadence=False,
        seed=99,
    )
    assert bundle.outputs == ("A", "B")
    tpf = mundey_tpf(str(bundle.tpf_path))
    tpf.calibrate(verbose=False, smear="standard", ddir=bundle.ddir, offline=True)

    assert np.all(np.isfinite(tpf.flux.value))
    assert np.all(tpf.hdu[1].data["FLUX_BKG"] == 0.0)
    assert tpf.meta["BKGSUB"] is False
    assert tpf.meta["HALOCVER"] == 1
    assert len(list(bundle.ddir.glob("*-s_col-blf.fits"))) == 4


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
