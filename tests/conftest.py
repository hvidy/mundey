"""Deterministic synthetic 1-output and 2-output TESS TPF fixture builders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import lightkurve as lk
import numpy as np
import pytest
from astropy.io import fits

from mundey.mundey import load_undershoot_xml


@dataclass(frozen=True)
class SyntheticTPFBundle:
    """Container for a synthetic TESS TPF file, collateral directory, and truth arrays."""

    tpf_path: Path
    ddir: Path
    true_signal: np.ndarray
    true_star_cube: np.ndarray
    true_sky_cube: np.ndarray
    cr_cube: np.ndarray
    injected_flux_cube: np.ndarray
    bleed_mask: np.ndarray
    columns_1based: np.ndarray
    outputs: tuple[str, ...]


def _output_for_col_1based(col: int) -> str:
    if 45 <= col <= 556:
        return "A"
    if 557 <= col <= 1068:
        return "B"
    if 1069 <= col <= 1580:
        return "C"
    if 1581 <= col <= 2092:
        return "D"
    raise ValueError(f"Column {col} outside TESS science outputs 45..2092.")


def _write_linearity_xml(path: Path, camera: int, ccd: int) -> None:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<ns1:linearityModel dataSetId="41006" linearityModelId="100" '
        'xmlns:ns1="http://tess.mit.edu/linearity-model">',
    ]
    for output in ("A", "B", "C", "D"):
        lines.extend(
            [
                f'  <channel cameraNumber="{camera}" ccdNumber="{ccd}" ccdOutput="{output}">',
                '    <linearityPoly maxDomain="42000.0" offsetx="0.0" order="2" '
                'originx="850.0" scalex="2.5e-05" type="standard" xindex="0">',
                "      <coeffs>",
                '        <coeff value="1.0"/>',
                '        <coeff value="-0.02"/>',
                '        <coeff value="0.03"/>',
                "      </coeffs>",
                "    </linearityPoly>",
                "  </channel>",
            ]
        )
    lines.append("</ns1:linearityModel>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _invert_linearity_poly(
    dn_linear: np.ndarray,
    *,
    originx: float = 850.0,
    scalex: float = 2.5e-5,
    coeffs: tuple[float, float, float] = (1.0, -0.02, 0.03),
) -> np.ndarray:
    """Invert ``d * p((d - originx) * scalex) = dn_linear`` via Newton-Raphson."""
    d = np.array(dn_linear, dtype=np.float64, copy=True)
    c0, c1, c2 = coeffs
    for _ in range(10):
        s = (d - originx) * scalex
        ps = c0 + c1 * s + c2 * (s * s)
        dps_dd = (c1 + 2.0 * c2 * s) * scalex
        f_val = d * ps - dn_linear
        df_val = ps + d * dps_dd
        d = d - f_val / df_val
    return d


def _apply_forward_undershoot(
    e_per_frame: np.ndarray,
    output: str,
    u_coeff: float,
) -> np.ndarray:
    """Apply exact forward serial LDE undershoot inverse to ``mundey_tpf.undershoot``."""
    w = np.asarray(e_per_frame, dtype=np.float64)
    z = np.empty_like(w)
    m = w.shape[-1]
    if output in ("A", "C"):
        z[..., 0] = w[..., 0] / (1.0 + u_coeff)
        for j in range(1, m):
            z[..., j] = w[..., j] - u_coeff * z[..., j - 1]
    else:
        z[..., -1] = w[..., -1] / (1.0 + u_coeff)
        for j in range(m - 2, -1, -1):
            z[..., j] = w[..., j] - u_coeff * z[..., j + 1]
    return z


def _true_1dblack_profile(output: str, rows_0based: np.ndarray) -> np.ndarray:
    """Evaluate output-specific quadratic 1D black profile (ADU summed over cadence)."""
    offsets = {"A": 4.0, "B": 5.5, "C": 3.5, "D": 4.5}
    x_norm = np.asarray(rows_0based, dtype=np.float64) / 2078.0 - 0.5
    return offsets[output] + 3.0 * x_norm + 0.8 * (x_norm * x_norm)


def _flat_field_subarray(rows_0based: np.ndarray, cols_0based: np.ndarray) -> np.ndarray:
    """Evaluate deterministic smooth flat-field subarray matching ``flat.fits``."""
    row_grid = np.linspace(0.0, 3.0, 2048, dtype=np.float64)[rows_0based][:, None]
    col_grid = np.linspace(0.0, 4.0, 2136, dtype=np.float64)[cols_0based][None, :]
    return (1.0 + 0.015 * np.sin(row_grid) * np.cos(col_grid)).astype(np.float32)


def build_synthetic_tess_tpf(
    out_dir: Path,
    *,
    col0_1based: int = 200,
    row0_1based: int = 1800,
    nx: int = 10,
    ny: int = 16,
    n_cadences: int = 44,
    camera: int = 2,
    ccd: int = 3,
    sector: int = 20,
    ticid: int = 123456789,
    include_object_header: bool = True,
    flag_bad_cadence: bool = True,
    write_collateral: bool = True,
    noise_sigma: float = 2.0,
    bleed_in_smrow: bool = False,
    seed: int = 42,
) -> SyntheticTPFBundle:
    """Build a deterministic SPOC-schema TESS TPF (1-output or 2-output A/B).

    Synthesizes ``RAW_CNTS`` using the exact forward CCD readout physics matching
    the generated static ``2dblack``, ``flat``, ``linearity.xml``, bundled
    ``tess-undershoot.xml``, and collateral ``SMROW_RAW`` / ``TVCOL_RAW`` tables.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    cols_1based = np.arange(col0_1based, col0_1based + nx, dtype=int)
    col_outputs = tuple(_output_for_col_1based(int(c)) for c in cols_1based)
    unique_outputs = tuple(o for o in ("A", "B", "C", "D") if o in col_outputs)

    # Time axis with a 1.2-day mid-sector downlink gap
    dt_days = 2.0 / (24.0 * 60.0)
    half = n_cadences // 2
    t_first = 1900.0 + np.arange(half, dtype=np.float64) * dt_days
    t_second = t_first[-1] + 1.2 + np.arange(n_cadences - half, dtype=np.float64) * dt_days
    time = np.concatenate([t_first, t_second])
    cadenceno = np.arange(450_000, 450_000 + n_cadences, dtype=np.int32)
    timecorr = np.full(n_cadences, 1.5e-4, dtype=np.float32)

    # Coherent asteroseismic pulsation + box transit dip
    pulsation = 2.5e-3 * np.sin(2.0 * np.pi * 6.0 * (time - time[0]))
    transit = np.where((np.arange(n_cadences) >= 8) & (np.arange(n_cadences) < 13), -6.0e-3, 0.0)
    true_signal = 1.0 + pulsation + transit

    # Spatial Moffat halo + central saturated bleed column (> maxDomain = 42,000 DN/frame)
    yy, xx = np.meshgrid(
        np.arange(ny, dtype=np.float64),
        np.arange(nx, dtype=np.float64),
        indexing="ij",
    )
    cy = 0.72 * (ny - 1)
    cx = 0.5 * (nx - 1)
    r2 = (yy - cy) ** 2 + (xx - cx) ** 2
    psf_halo = 2_400.0 / (1.0 + r2 / 3.5) ** 2.0
    bleed_col = nx // 2
    bleed_mask = np.zeros((ny, nx), dtype=bool)
    bleed_mask[max(1, int(0.55 * ny)) :, bleed_col] = True

    true_star_cube = true_signal[:, None, None] * psf_halo[None, :, :]
    # Inject saturated core/bleed exceeding maxDomain in electrons/s
    true_star_cube[:, bleed_mask] = 1.35e6 * true_signal[:, None]

    # Sky background sweeping from dark (< 10 e-/s) to bright Earthshine (> 600 e-/s)
    xi = np.linspace(-1.0, 1.0, nx)[None, :]
    eta = np.linspace(-1.0, 1.0, ny)[:, None]
    plane = 1.0 + 0.12 * xi - 0.08 * eta
    sky_pedestal = np.where(
        np.arange(n_cadences) < half,
        6.0,  # dark high-latitude sky (< 10 e-/s)
        680.0 + 60.0 * np.linspace(0.0, 1.0, n_cadences),  # bright Earthshine (> 600 e-/s)
    )
    true_sky_cube = sky_pedestal[:, None, None] * plane[None, :, :]

    # Inject deterministic single-pixel cosmic-ray impulses
    cr_cube = np.zeros((n_cadences, ny, nx), dtype=np.float64)
    if n_cadences >= 6 and ny >= 4 and nx >= 4:
        cr_cube[3, 1, 1] = 450.0
        cr_cube[min(n_cadences - 2, half + 2), ny - 2, nx - 2] = 620.0

    int_time = 1.98
    num_frm = 60
    crmiten = True
    crblksz = 10
    n_eff = num_frm * (crblksz - 2) / crblksz  # 48 effective frames
    fxdoff = 209_700
    black2d_val = 1.5
    smear_base_e = 24.0
    gains = {"A": 5.28, "B": 5.26, "C": 5.25, "D": 5.21}
    meanblc = {"A": 6812, "B": 6724, "C": 6753, "D": 6433}
    sx1_map = {"A": 44, "B": 556, "C": 1068, "D": 1580}

    noise_cube = (
        rng.normal(0.0, noise_sigma, size=(n_cadences, ny, nx))
        if noise_sigma > 0.0
        else np.zeros((n_cadences, ny, nx), dtype=np.float64)
    )
    injected_flux_cube = true_star_cube + true_sky_cube + cr_cube + noise_cube

    rows_0based = np.arange(row0_1based - 1, row0_1based - 1 + ny, dtype=int)
    cols_0based = cols_1based - 1
    flat_cutout = _flat_field_subarray(rows_0based, cols_0based).astype(np.float64)

    true_smear_cutout = np.full((n_cadences, nx), smear_base_e, dtype=np.float64)
    if bleed_in_smrow:
        true_smear_cutout[:, bleed_col] += 180.0 * true_signal

    # Forward readout chain per output (flat -> smear -> undershoot -> nonlinearity -> blacks)
    underroot = load_undershoot_xml()
    u_coeffs: dict[str, float] = {}
    raw_cnts = np.zeros((n_cadences, ny, nx), dtype=np.int32)
    observed_e_per_s = np.zeros((n_cadences, ny, nx), dtype=np.float64)

    for out in unique_outputs:
        idx = [i for i, o in enumerate(col_outputs) if o == out]
        u_node = underroot.find(
            f".//*[@cameraNumber='{camera}'][@ccdNumber='{ccd}'][@ccdOutput='{out}']"
        )
        assert u_node is not None
        u_coeff = float(u_node.attrib["undershoot"])
        u_coeffs[out] = u_coeff

        sci_e_per_frame = injected_flux_cube[:, :, idx] * int_time * flat_cutout[None, :, idx]
        w_e_per_frame = sci_e_per_frame + true_smear_cutout[:, None, idx]
        z_e_per_frame = _apply_forward_undershoot(w_e_per_frame, out, u_coeff)
        observed_e_per_s[:, :, idx] = z_e_per_frame / int_time

        dn_per_frame = _invert_linearity_poly(z_e_per_frame / gains[out])
        b1d_cutout = _true_1dblack_profile(out, rows_0based)
        counts = np.round(
            dn_per_frame * n_eff
            + b1d_cutout[None, :, None]
            + black2d_val
            + fxdoff
            - meanblc[out]
        ).astype(np.int32)
        raw_cnts[:, :, idx] = np.clip(counts, 0, 2_000_000_000)

    # Stale SPOC background (non-zero, to verify mundey zeroes FLUX_BKG & FLUX_BKG_ERR)
    spoc_flux_bkg = (true_sky_cube + 12.5).astype(np.float32)
    spoc_flux_bkg_err = np.full((n_cadences, ny, nx), 1.75, dtype=np.float32)
    spoc_flux = (observed_e_per_s - spoc_flux_bkg).astype(np.float32)
    spoc_flux_err = np.sqrt(np.maximum(observed_e_per_s, 1.0) / (int_time * n_eff)).astype(
        np.float32
    )

    quality = np.zeros(n_cadences, dtype=np.int32)
    if flag_bad_cadence and n_cadences > 5:
        quality[4] = 128  # manual exclude / attitude tweak flag

    pos_corr1 = np.zeros(n_cadences, dtype=np.float32)
    pos_corr2 = np.zeros(n_cadences, dtype=np.float32)

    # Primary HDU [0]
    hdu0 = fits.PrimaryHDU()
    hdr0 = hdu0.header
    hdr0["ORIGIN"] = "NASA/Ames"
    hdr0["CREATOR"] = "TargetPixelExporterPipelineModule"
    hdr0["TELESCOP"] = "TESS"
    hdr0["INSTRUME"] = "TESS Photometer"
    if include_object_header:
        hdr0["OBJECT"] = f"TIC {ticid}"
    hdr0["TICID"] = ticid
    hdr0["SECTOR"] = sector
    hdr0["CAMERA"] = camera
    hdr0["CCD"] = ccd
    hdr0["RADESYS"] = "ICRS"
    hdr0["RA_OBJ"] = 101.287
    hdr0["DEC_OBJ"] = -16.716
    hdr0["EQUINOX"] = 2000.0
    hdr0["CRMITEN"] = crmiten
    hdr0["CRBLKSZ"] = crblksz

    # PIXELS BinTableHDU [1]
    npix = ny * nx
    dim_str = f"({nx},{ny})"
    cols = fits.ColDefs(
        [
            fits.Column(
                name="TIME",
                format="D",
                unit="BJD - 2457000, days",
                disp="D14.7",
                array=time,
            ),
            fits.Column(
                name="TIMECORR",
                format="E",
                unit="d",
                disp="E14.7",
                array=timecorr,
            ),
            fits.Column(
                name="CADENCENO",
                format="J",
                disp="I10",
                array=cadenceno,
            ),
            fits.Column(
                name="RAW_CNTS",
                format=f"{npix}J",
                unit="count",
                null=-1,
                disp="I8",
                dim=dim_str,
                array=raw_cnts,
            ),
            fits.Column(
                name="FLUX",
                format=f"{npix}E",
                unit="e-/s",
                disp="E14.7",
                dim=dim_str,
                array=spoc_flux,
            ),
            fits.Column(
                name="FLUX_ERR",
                format=f"{npix}E",
                unit="e-/s",
                disp="E14.7",
                dim=dim_str,
                array=spoc_flux_err,
            ),
            fits.Column(
                name="FLUX_BKG",
                format=f"{npix}E",
                unit="e-/s",
                disp="E14.7",
                dim=dim_str,
                array=spoc_flux_bkg,
            ),
            fits.Column(
                name="FLUX_BKG_ERR",
                format=f"{npix}E",
                unit="e-/s",
                disp="E14.7",
                dim=dim_str,
                array=spoc_flux_bkg_err,
            ),
            fits.Column(name="QUALITY", format="J", disp="B16.16", array=quality),
            fits.Column(
                name="POS_CORR1",
                format="E",
                unit="pixel",
                disp="E14.7",
                array=pos_corr1,
            ),
            fits.Column(
                name="POS_CORR2",
                format="E",
                unit="pixel",
                disp="E14.7",
                array=pos_corr2,
            ),
        ]
    )
    hdu1 = fits.BinTableHDU.from_columns(cols, name="PIXELS")
    hdr1 = hdu1.header
    hdr1["INHERIT"] = True
    hdr1["TELESCOP"] = "TESS"
    hdr1["INSTRUME"] = "TESS Photometer"
    if include_object_header:
        hdr1["OBJECT"] = f"TIC {ticid}"
    hdr1["TICID"] = ticid
    hdr1["RADESYS"] = "ICRS"
    hdr1["RA_OBJ"] = 101.287
    hdr1["DEC_OBJ"] = -16.716
    hdr1["EQUINOX"] = 2000.0
    hdr1["TIMEREF"] = "SOLARSYSTEM"
    hdr1["TASSIGN"] = "SPACECRAFT"
    hdr1["TIMESYS"] = "TDB"
    hdr1["BJDREFI"] = 2457000
    hdr1["BJDREFF"] = 0.0
    hdr1["TIMEUNIT"] = "d"
    hdr1["INT_TIME"] = int_time
    hdr1["READTIME"] = 0.02
    hdr1["FRAMETIM"] = 2.0
    hdr1["NUM_FRM"] = num_frm
    hdr1["TIMEDEL"] = dt_days
    hdr1["BACKAPP"] = True
    hdr1["NREADOUT"] = int(n_eff)
    hdr1["FXDOFF"] = fxdoff
    for out in ("A", "B", "C", "D"):
        hdr1[f"GAIN{out}"] = gains[out]
        hdr1[f"READNOI{out}"] = 8.5
        hdr1[f"MEANBLC{out}"] = meanblc[out]
    for col_idx in (4, 5, 6, 7, 8):
        hdr1[f"WCSN{col_idx}P"] = "PHYSICAL"
        hdr1[f"WCAX{col_idx}P"] = 2
        hdr1[f"1CTY{col_idx}P"] = "RAWX"
        hdr1[f"2CTY{col_idx}P"] = "RAWY"
        hdr1[f"1CUN{col_idx}P"] = "PIXEL"
        hdr1[f"2CUN{col_idx}P"] = "PIXEL"
        hdr1[f"1CRV{col_idx}P"] = col0_1based
        hdr1[f"2CRV{col_idx}P"] = row0_1based
        hdr1[f"1CDL{col_idx}P"] = 1.0
        hdr1[f"2CDL{col_idx}P"] = 1.0
        hdr1[f"1CRP{col_idx}P"] = 1
        hdr1[f"2CRP{col_idx}P"] = 1

    # APERTURE ImageHDU [2]
    aperture = np.ones((ny, nx), dtype=np.int32)
    aperture[1:-1, 1:-1] = 3
    hdu2 = fits.ImageHDU(data=aperture, name="APERTURE")
    hdr2 = hdu2.header
    hdr2["INHERIT"] = True
    hdr2["TELESCOP"] = "TESS"
    hdr2["INSTRUME"] = "TESS Photometer"
    hdr2["RADESYS"] = "ICRS"
    hdr2["EQUINOX"] = 2000.0
    hdr2["WCSAXES"] = 2
    hdr2["CTYPE1"] = "RA---TAN"
    hdr2["CTYPE2"] = "DEC--TAN"
    hdr2["CRPIX1"] = nx / 2.0
    hdr2["CRPIX2"] = ny / 2.0
    hdr2["CRVAL1"] = 101.287
    hdr2["CRVAL2"] = -16.716
    hdr2["CUNIT1"] = "deg"
    hdr2["CUNIT2"] = "deg"
    hdr2["CDELT1"] = -0.0057
    hdr2["CDELT2"] = 0.0057
    hdr2["PC1_1"] = 1.0
    hdr2["PC1_2"] = 0.0
    hdr2["PC2_1"] = 0.0
    hdr2["PC2_2"] = 1.0
    hdr2["WCSNAMEP"] = "PHYSICAL"
    hdr2["WCSAXESP"] = 2
    hdr2["CTYPE1P"] = "RAWX"
    hdr2["CUNIT1P"] = "PIXEL"
    hdr2["CRPIX1P"] = 1
    hdr2["CRVAL1P"] = col0_1based
    hdr2["CDELT1P"] = 1.0
    hdr2["CTYPE2P"] = "RAWY"
    hdr2["CUNIT2P"] = "PIXEL"
    hdr2["CRPIX2P"] = 1
    hdr2["CRVAL2P"] = row0_1based
    hdr2["CDELT2P"] = 1.0

    prefix = f"tess2020001000000-s{sector:04d}"
    run_id = "0165"
    tpf_path = out_dir / f"{prefix}-{ticid:016d}-{run_id}-s_tp.fits"
    fits.HDUList([hdu0, hdu1, hdu2]).writeto(tpf_path, overwrite=True)

    if write_collateral:
        # Static 2D black (2078 x 2136) and flat field (2048 x 2136)
        black2d_path = out_dir / f"tess2018323-{camera}-{ccd}-2dblack.fits"
        if not black2d_path.exists():
            black2d = np.full((2078, 2136), black2d_val, dtype=np.float32)
            fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=black2d)]).writeto(
                black2d_path, overwrite=True
            )
        flat_path = out_dir / f"tess2018323-{camera}-{ccd}-flat.fits"
        if not flat_path.exists():
            flat = _flat_field_subarray(np.arange(2048), np.arange(2136))
            fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=flat)]).writeto(
                flat_path, overwrite=True
            )
        lin_path = out_dir / "tess2018143203310-41006_100-linearity.xml"
        if not lin_path.exists():
            _write_linearity_xml(lin_path, camera=camera, ccd=ccd)

        all_rows = np.arange(2078, dtype=int)
        sm_rows = np.arange(2058, 2068, dtype=int)
        bleed_out = col_outputs[bleed_col]
        bleed_sm_col = int(cols_1based[bleed_col] - 1 - sx1_map[bleed_out])

        # Collateral smrow (N, 10, 512) and tvcol (N, 2078, 11) per output
        for out in unique_outputs:
            base_black = fxdoff - meanblc[out]
            sm_w = np.full((n_cadences, 10, 512), smear_base_e, dtype=np.float64)
            if bleed_in_smrow and out == bleed_out and 0 <= bleed_sm_col < 512:
                sm_w[:, :, bleed_sm_col] += 250_000.0
            sm_z = _apply_forward_undershoot(sm_w, out, u_coeffs[out])
            sm_dn = _invert_linearity_poly(sm_z / gains[out])
            b1d_sm = _true_1dblack_profile(out, sm_rows)
            smrow_raw = np.round(
                sm_dn * n_eff + b1d_sm[None, :, None] + black2d_val + base_black
            ).astype(np.int32)
            smrow_col = fits.Column(
                name="SMROW_RAW",
                format=f"{10 * 512}J",
                dim="(512,10)",
                array=smrow_raw,
            )
            smrow_path = (
                out_dir / f"{prefix}-smrow-{camera}-{ccd}-{out.lower()}-{run_id}-s_col.fits"
            )
            fits.HDUList(
                [fits.PrimaryHDU(), fits.BinTableHDU.from_columns([smrow_col])]
            ).writeto(smrow_path, overwrite=True)

            b1d_all = _true_1dblack_profile(out, all_rows)
            tv_sigma = 0.25 if noise_sigma > 0.0 else 0.0
            tv_noise = (
                rng.normal(0.0, tv_sigma, size=(n_cadences, 2078, 11))
                if tv_sigma > 0.0
                else np.zeros((n_cadences, 2078, 11), dtype=np.float64)
            )
            tvcol_raw = np.round(
                base_black + black2d_val + b1d_all[None, :, None] + tv_noise
            ).astype(np.int32)
            tvcol_col = fits.Column(
                name="TVCOL_RAW",
                format=f"{2078 * 11}J",
                dim="(11,2078)",
                array=tvcol_raw,
            )
            tvcol_path = (
                out_dir / f"{prefix}-tvcol-{camera}-{ccd}-{out.lower()}-{run_id}-s_col.fits"
            )
            fits.HDUList(
                [fits.PrimaryHDU(), fits.BinTableHDU.from_columns([tvcol_col])]
            ).writeto(tvcol_path, overwrite=True)

    return SyntheticTPFBundle(
        tpf_path=tpf_path,
        ddir=out_dir,
        true_signal=true_signal,
        true_star_cube=true_star_cube,
        true_sky_cube=true_sky_cube,
        cr_cube=cr_cube,
        injected_flux_cube=injected_flux_cube,
        bleed_mask=bleed_mask,
        columns_1based=cols_1based,
        outputs=unique_outputs,
    )


@pytest.fixture
def tpf_builder():
    """Return the deterministic synthetic TESS TPF builder callable."""
    return build_synthetic_tess_tpf


@pytest.fixture
def synthetic_1output_bundle(tmp_path: Path) -> SyntheticTPFBundle:
    """Synthetic 1-output TESS TPF bundle purely on Output A (cols 200..209)."""
    return build_synthetic_tess_tpf(
        tmp_path / "one_output",
        col0_1based=200,
        row0_1based=1800,
        nx=10,
        ny=14,
        n_cadences=40,
        write_collateral=False,
        seed=101,
    )


@pytest.fixture
def synthetic_1output_tpf(synthetic_1output_bundle: SyntheticTPFBundle) -> lk.TessTargetPixelFile:
    """Loaded lightkurve TessTargetPixelFile for the 1-output fixture."""
    return lk.read(synthetic_1output_bundle.tpf_path)


@pytest.fixture
def synthetic_2output_bundle(tmp_path: Path) -> SyntheticTPFBundle:
    """Synthetic 2-output TESS TPF bundle straddling Outputs A/B at cols 556/557."""
    return build_synthetic_tess_tpf(
        tmp_path / "two_output",
        col0_1based=551,
        row0_1based=1800,
        nx=12,
        ny=16,
        n_cadences=44,
        write_collateral=False,
        seed=202,
    )


@pytest.fixture
def synthetic_2output_tpf(synthetic_2output_bundle: SyntheticTPFBundle) -> lk.TessTargetPixelFile:
    """Loaded lightkurve TessTargetPixelFile straddling Outputs A and B."""
    return lk.read(synthetic_2output_bundle.tpf_path)
