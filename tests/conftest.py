"""Deterministic synthetic 1/2/3/4-output TESS TPF fixture builders (physical v2.0 chain)."""

from __future__ import annotations

import io
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


def _write_linearity_xml(
    path: Path,
    camera: int,
    ccd: int,
    *,
    degree6: bool = False,
) -> tuple[float, float, tuple[float, ...]]:
    if degree6:
        # Monotone on [0, maxDomain=42000], diverging strongly for unclamped d > 42000:
        originx = 10_000.0
        scalex = 2.0e-5
        coeffs = (1.0, -0.01, 0.015, -0.02, 0.03, -0.05, -0.80)
        order = 6
    else:
        originx = 850.0
        scalex = 2.5e-5
        coeffs = (1.0, -0.02, 0.03)
        order = 2

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<ns1:linearityModel dataSetId="41006" linearityModelId="100" '
        'xmlns:ns1="http://tess.mit.edu/linearity-model">',
    ]
    for output in ("A", "B", "C", "D"):
        lines.extend(
            [
                f'  <channel cameraNumber="{camera}" ccdNumber="{ccd}" ccdOutput="{output}">',
                f'    <linearityPoly maxDomain="42000.0" offsetx="0.0" order="{order}" '
                f'originx="{originx}" scalex="{scalex}" type="standard" xindex="0">',
                "      <coeffs>",
            ]
        )
        for c_val in coeffs:
            lines.append(f'        <coeff value="{c_val}"/>')
        lines.extend(
            [
                "      </coeffs>",
                "    </linearityPoly>",
                "  </channel>",
            ]
        )
    lines.append("</ns1:linearityModel>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return originx, scalex, coeffs


def _invert_clamped_linearity_poly(
    dn_linear: np.ndarray,
    *,
    originx: float = 850.0,
    scalex: float = 2.5e-5,
    offsetx: float = 0.0,
    max_domain: float = 42_000.0,
    coeffs: tuple[float, ...] = (1.0, -0.02, 0.03),
) -> np.ndarray:
    """Invert ``d * P(offsetx + (clip(d, 0, maxDomain) - originx) * scalex) = dn_linear``."""
    coeffs_arr = np.asarray(coeffs, dtype=np.float64)
    coeffs_desc = coeffs_arr[::-1]
    deriv_desc = (
        (coeffs_arr[1:] * np.arange(1, len(coeffs_arr), dtype=np.float64))[::-1]
        if len(coeffs_arr) > 1
        else np.zeros(1, dtype=np.float64)
    )

    s_max = offsetx + (max_domain - originx) * scalex
    p_max = float(np.polyval(coeffs_desc, s_max))
    linear_at_max = max_domain * p_max

    target = np.asarray(dn_linear, dtype=np.float64)
    d = np.clip(target, 0.0, max_domain)
    for _ in range(14):
        s = offsetx + (d - originx) * scalex
        ps = np.polyval(coeffs_desc, s)
        dps_dd = np.polyval(deriv_desc, s) * scalex
        f_val = d * ps - target
        df_val = ps + d * dps_dd
        d = np.clip(d - f_val / df_val, 0.0, max_domain)

    # Above maxDomain, forward charge exceeds maxDomain linearly at saturated gain factor p_max
    d_sat = target / p_max
    return np.where(target > linear_at_max, d_sat, d)


def _apply_forward_causal_undershoot_adu(
    x_adu: np.ndarray,
    output: str,
    u_coeff: float,
) -> np.ndarray:
    """Apply physical forward causal serial LDE undershoot ``y_m = x_m - u * x_{m-1}`` in ADU."""
    x = np.asarray(x_adu, dtype=np.float64)
    y = np.empty_like(x)
    m = x.shape[-1]
    if m == 0:
        return y
    if output in ("A", "C"):
        y[..., 0] = (1.0 - u_coeff) * x[..., 0]
        if m > 1:
            y[..., 1:] = x[..., 1:] - u_coeff * x[..., :-1]
    else:
        y[..., -1] = (1.0 - u_coeff) * x[..., -1]
        if m > 1:
            y[..., :-1] = x[..., :-1] - u_coeff * x[..., 1:]
    return y


def _true_1dblack_profile(output: str, rows_0based: np.ndarray) -> np.ndarray:
    """Evaluate output-specific quadratic 1D black profile (ADU summed over cadence)."""
    offsets = {"A": 4.0, "B": 5.5, "C": 3.5, "D": 4.5}
    x_norm = np.asarray(rows_0based, dtype=np.float64) / 2078.0 - 0.5
    return offsets[output] + 3.0 * x_norm + 0.8 * (x_norm * x_norm)


def _spatial_2dblack_full_ccd(meanblc: dict[str, int]) -> np.ndarray:
    """Build realistic per-frame ``2dblack.fits`` ``(2078, 2136)`` around ``MEANBLC{out}``."""
    r_grid = np.linspace(-1.0, 1.0, 2078, dtype=np.float64)[:, None]
    c_grid = np.linspace(-1.0, 1.0, 2136, dtype=np.float64)[None, :]
    delta = 2.4 + 0.9 * np.sin(2.5 * r_grid) * np.cos(3.0 * c_grid)

    img = np.full((2078, 2136), float(meanblc["A"]), dtype=np.float64)
    bounds = {
        "A": ((44, 556), (2092, 2103)),
        "B": ((556, 1068), (2103, 2114)),
        "C": ((1068, 1580), (2114, 2125)),
        "D": ((1580, 2092), (2125, 2136)),
    }
    for out, ((sx1, sx2), (tx1, tx2)) in bounds.items():
        img[:, sx1:sx2] = float(meanblc[out]) + delta[:, sx1:sx2]
        img[:, tx1:tx2] = float(meanblc[out]) + delta[:, tx1:tx2]
    return img.astype(np.float32)


def _flat_field_subarray(rows_0based: np.ndarray, cols_0based: np.ndarray) -> np.ndarray:
    """Evaluate deterministic smooth flat-field subarray matching ``flat.fits``."""
    row_grid = np.linspace(0.0, 3.0, 2048, dtype=np.float64)[rows_0based][:, None]
    col_grid = np.linspace(0.0, 4.0, 2136, dtype=np.float64)[cols_0based][None, :]
    return (1.0 + 0.015 * np.sin(row_grid) * np.cos(col_grid)).astype(np.float32)


_STATIC_ARRAY_CACHE: dict[str, object] = {}


def _get_cached_black2d(meanblc: dict[str, int]) -> tuple[np.ndarray, bytes]:
    key = f"black2d:{tuple(sorted(meanblc.items()))}"
    cached = _STATIC_ARRAY_CACHE.get(key)
    if cached is not None:
        return cached  # type: ignore[return-value]
    img_f32 = _spatial_2dblack_full_ccd(meanblc)
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=img_f32)]).writeto(buf)
    result = (img_f32.astype(np.float64), buf.getvalue())
    _STATIC_ARRAY_CACHE[key] = result
    return result


def _get_cached_flat_bytes() -> bytes:
    cached = _STATIC_ARRAY_CACHE.get("flat_bytes")
    if cached is not None:
        return cached  # type: ignore[return-value]
    flat = _flat_field_subarray(np.arange(2048), np.arange(2136))
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=flat)]).writeto(buf)
    data = buf.getvalue()
    _STATIC_ARRAY_CACHE["flat_bytes"] = data
    return data


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
    include_bleed: bool = True,
    degree6_linearity: bool = False,
    fast_20s: bool = False,
    pulsation_amp: float = 2.5e-3,
    sky_dark_level: float = 6.0,
    star_scale: float = 1.0,
    thermal_step_cadence: int | None = None,
    thermal_step_b1d_adu: float = 0.0,
    thermal_step_smear_e: float = 0.0,
    seed: int = 42,
) -> SyntheticTPFBundle:
    """Build a deterministic SPOC-schema TESS TPF (1/2/3/4-output, 2-min or 20-s).

    Synthesizes ``RAW_CNTS`` using the physical v2.0 spacecraft readout chain:
    flat field -> shutterless smear -> inverse nonlinearity -> ADU causal LDE
    undershoot -> 1D black + frame-scaled 2D black + ``FXDOFF - N_eff * MEANBLC``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    cols_1based = np.arange(col0_1based, col0_1based + nx, dtype=int)
    col_outputs = tuple(_output_for_col_1based(int(c)) for c in cols_1based)
    unique_outputs = tuple(o for o in ("A", "B", "C", "D") if o in col_outputs)

    dt_days = (20.0 if fast_20s else 120.0) / 86400.0
    half = max(1, n_cadences // 2)
    t_first = 1900.0 + np.arange(half, dtype=np.float64) * dt_days
    t_second = t_first[-1] + 1.2 + np.arange(n_cadences - half, dtype=np.float64) * dt_days
    time = np.concatenate([t_first, t_second])
    cadenceno = np.arange(450_000, 450_000 + n_cadences, dtype=np.int32)
    timecorr = np.full(n_cadences, 1.5e-4, dtype=np.float32)

    pulsation = pulsation_amp * np.sin(2.0 * np.pi * 6.0 * (time - time[0]))
    transit = np.where((np.arange(n_cadences) >= 8) & (np.arange(n_cadences) < 13), -6.0e-3, 0.0)
    true_signal = 1.0 + pulsation + transit

    yy, xx = np.meshgrid(
        np.arange(ny, dtype=np.float64),
        np.arange(nx, dtype=np.float64),
        indexing="ij",
    )
    cy = 0.78 * (ny - 1)
    cx = 0.5 * (nx - 1)
    r2 = (yy - cy) ** 2 + (xx - cx) ** 2
    psf_halo = (3_200.0 * float(star_scale)) / (1.0 + r2 / 4.0) ** 1.5
    bleed_col = nx // 2
    bleed_mask = np.zeros((ny, nx), dtype=bool)
    if include_bleed:
        bleed_mask[max(1, int(0.60 * ny)) :, bleed_col] = True

    true_star_cube = true_signal[:, None, None] * psf_halo[None, :, :]
    if include_bleed:
        # Exceeds maxDomain = 42,000 DN/frame (~56,700 DN/frame < 65,535 16-bit ADC limit)
        true_star_cube[:, bleed_mask] = 1.42e5 * true_signal[:, None]

    xi = np.linspace(-1.0, 1.0, nx)[None, :]
    eta = np.linspace(-1.0, 1.0, ny)[:, None]
    plane = 1.0 + 0.12 * xi - 0.08 * eta
    sky_pedestal = np.where(
        np.arange(n_cadences) < half,
        sky_dark_level,
        680.0 + 60.0 * np.linspace(0.0, 1.0, n_cadences),
    )
    true_sky_cube = sky_pedestal[:, None, None] * plane[None, :, :]

    cr_cube = np.zeros((n_cadences, ny, nx), dtype=np.float64)
    if n_cadences >= 6 and ny >= 4 and nx >= 4:
        cr_cube[3, 1, 1] = 450.0
        cr_cube[min(n_cadences - 2, half + 2), ny - 2, nx - 2] = 620.0

    int_time = 1.98
    num_frm = 10 if fast_20s else 60
    crmiten = not fast_20s
    crblksz = 10
    n_eff = float(num_frm) * ((crblksz - 2.0) / crblksz if crmiten else 1.0)
    fxdoff = 209_700
    smear_base_e = 24.0
    gains = {"A": 5.28, "B": 5.26, "C": 5.25, "D": 5.21}
    meanblc = {"A": 6812, "B": 6724, "C": 6753, "D": 6433}
    sx1_map = {"A": 44, "B": 556, "C": 1068, "D": 1580}
    tx1_map = {"A": 2092, "B": 2103, "C": 2114, "D": 2125}

    noise_cube = (
        rng.normal(0.0, noise_sigma, size=(n_cadences, ny, nx))
        if noise_sigma > 0.0
        else np.zeros((n_cadences, ny, nx), dtype=np.float64)
    )
    injected_flux_cube = true_star_cube + true_sky_cube + cr_cube + noise_cube

    rows_0based = np.arange(row0_1based - 1, row0_1based - 1 + ny, dtype=int)
    cols_0based = cols_1based - 1
    flat_cutout = _flat_field_subarray(rows_0based, cols_0based).astype(np.float64)
    black2d_full, black2d_bytes = _get_cached_black2d(meanblc)

    if degree6_linearity:
        lin_originx = 10_000.0
        lin_scalex = 2.0e-5
        lin_coeffs = (1.0, -0.01, 0.015, -0.02, 0.03, -0.05, -0.80)
    else:
        lin_originx = 850.0
        lin_scalex = 2.5e-5
        lin_coeffs = (1.0, -0.02, 0.03)

    step_mask = (
        (np.arange(n_cadences) >= int(thermal_step_cadence))
        if thermal_step_cadence is not None
        else np.zeros(n_cadences, dtype=bool)
    )
    b1d_step_vec = np.where(step_mask, float(thermal_step_b1d_adu), 0.0)
    smear_step_vec = np.where(step_mask, float(thermal_step_smear_e), 0.0)

    true_smear_cutout = (smear_base_e + smear_step_vec[:, None]) * np.ones(
        (n_cadences, nx), dtype=np.float64
    )
    if bleed_in_smrow and include_bleed:
        true_smear_cutout[:, bleed_col] += 180.0 * true_signal

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
        observed_e_per_s[:, :, idx] = w_e_per_frame / int_time

        dn_per_frame = _invert_clamped_linearity_poly(
            w_e_per_frame / gains[out],
            originx=lin_originx,
            scalex=lin_scalex,
            coeffs=lin_coeffs,
        )
        x_adu = dn_per_frame * n_eff
        y_adu = _apply_forward_causal_undershoot_adu(x_adu, out, u_coeff)

        b1d_cutout = _true_1dblack_profile(out, rows_0based)[None, :, None] + b1d_step_vec[
            :, None, None
        ]
        b2d_slice = black2d_full[np.ix_(rows_0based, cols_0based[idx])]
        counts = np.round(
            y_adu
            + b1d_cutout
            + n_eff * (b2d_slice[None, :, :] - float(meanblc[out]))
            + fxdoff
        ).astype(np.int32)
        raw_cnts[:, :, idx] = np.clip(counts, 0, 2_000_000_000)

    spoc_flux_bkg = (true_sky_cube + 12.5).astype(np.float32)
    spoc_flux_bkg_err = np.full((n_cadences, ny, nx), 1.75, dtype=np.float32)
    spoc_flux = (observed_e_per_s - spoc_flux_bkg).astype(np.float32)
    spoc_flux_err = np.sqrt(np.maximum(observed_e_per_s, 1.0) / (int_time * n_eff)).astype(
        np.float32
    )

    quality = np.zeros(n_cadences, dtype=np.int32)
    if flag_bad_cadence and n_cadences > 5:
        quality[4] = 128

    pos_corr1 = np.zeros(n_cadences, dtype=np.float32)
    pos_corr2 = np.zeros(n_cadences, dtype=np.float32)

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
    tp_tag = "a_fast-tp" if fast_20s else "s_tp"
    col_tag = "a_fast-col" if fast_20s else "s_col"
    tpf_path = out_dir / f"{prefix}-{ticid:016d}-{run_id}-{tp_tag}.fits"
    fits.HDUList([hdu0, hdu1, hdu2]).writeto(tpf_path, overwrite=True)

    if write_collateral:
        black2d_path = out_dir / f"tess2018323-{camera}-{ccd}-2dblack.fits"
        if not black2d_path.exists():
            black2d_path.write_bytes(black2d_bytes)

        flat_path = out_dir / f"tess2018323-{camera}-{ccd}-flat.fits"
        if not flat_path.exists():
            flat_path.write_bytes(_get_cached_flat_bytes())

        lin_path = out_dir / "tess2018143203310-41006_100-linearity.xml"
        _write_linearity_xml(lin_path, camera=camera, ccd=ccd, degree6=degree6_linearity)

        all_rows = np.arange(2078, dtype=int)
        sm_rows = np.arange(2058, 2068, dtype=int)
        bleed_out = col_outputs[bleed_col]
        bleed_sm_col = int(cols_1based[bleed_col] - 1 - sx1_map[bleed_out])

        for out in unique_outputs:
            sx1 = sx1_map[out]
            tx1 = tx1_map[out]
            b2d_sm = black2d_full[2058:2068, sx1 : sx1 + 512]
            b2d_tv = black2d_full[:, tx1 : tx1 + 11]

            sm_w = (smear_base_e + smear_step_vec[:, None, None]) * np.ones(
                (n_cadences, 10, 512), dtype=np.float64
            )
            if bleed_in_smrow and include_bleed and out == bleed_out and 0 <= bleed_sm_col < 512:
                sm_w[:, :, bleed_sm_col] += 250_000.0
            sm_dn = _invert_clamped_linearity_poly(
                sm_w / gains[out],
                originx=lin_originx,
                scalex=lin_scalex,
                coeffs=lin_coeffs,
            )
            sm_y_adu = _apply_forward_causal_undershoot_adu(
                sm_dn * n_eff, out, u_coeffs[out]
            )
            b1d_sm = _true_1dblack_profile(out, sm_rows)[None, :, None] + b1d_step_vec[
                :, None, None
            ]
            smrow_raw = np.round(
                sm_y_adu
                + b1d_sm
                + n_eff * (b2d_sm[None, :, :] - float(meanblc[out]))
                + fxdoff
            ).astype(np.int32)

            col_hdr0 = fits.Header()
            col_hdr0["SECTOR"] = sector
            col_hdr0["CAMERA"] = camera
            col_hdr0["CCD"] = ccd
            col_hdr0["CRMITEN"] = crmiten
            col_hdr0["CRBLKSZ"] = crblksz

            sm_cols = fits.ColDefs(
                [
                    fits.Column(
                        name="TIME",
                        format="D",
                        unit="BJD - 2457000, days",
                        array=time,
                    ),
                    fits.Column(name="CADENCENO", format="J", array=cadenceno),
                    fits.Column(
                        name="SMROW_RAW",
                        format=f"{10 * 512}J",
                        null=-1,
                        dim="(512,10)",
                        array=smrow_raw,
                    ),
                ]
            )
            sm_hdu1 = fits.BinTableHDU.from_columns(sm_cols, name="COLLATERAL")
            sm_hdu1.header["NUM_FRM"] = num_frm
            sm_hdu1.header["INT_TIME"] = int_time
            sm_hdu1.header["FXDOFF"] = fxdoff
            sm_hdu1.header[f"MEANBLC{out}"] = meanblc[out]
            sm_hdu1.header[f"GAIN{out}"] = gains[out]
            sm_hdu1.header[f"READNOI{out}"] = 8.5

            smrow_path = (
                out_dir / f"{prefix}-smrow-{camera}-{ccd}-{out.lower()}-{run_id}-{col_tag}.fits"
            )
            fits.HDUList([fits.PrimaryHDU(header=col_hdr0), sm_hdu1]).writeto(
                smrow_path, overwrite=True
            )

            b1d_all = _true_1dblack_profile(out, all_rows)[None, :, None] + b1d_step_vec[
                :, None, None
            ]
            tv_sigma = 0.25 if noise_sigma > 0.0 else 0.0
            tv_noise = (
                rng.normal(0.0, tv_sigma, size=(n_cadences, 2078, 11))
                if tv_sigma > 0.0
                else np.zeros((n_cadences, 2078, 11), dtype=np.float64)
            )
            tv_adu = (
                b1d_all
                + n_eff * (b2d_tv[None, :, :] - float(meanblc[out]))
                + fxdoff
                + tv_noise
            )
            # Inject settling ramp in cols 0..2, step in col 10, and 5-sigma CR in col 5
            tv_adu[:, :, 0:3] += 380.0
            tv_adu[:, :, 10] -= 260.0
            tv_adu[:, 150, 5] += 1_200.0
            tvcol_raw = np.round(tv_adu).astype(np.int32)

            tv_cols = fits.ColDefs(
                [
                    fits.Column(
                        name="TIME",
                        format="D",
                        unit="BJD - 2457000, days",
                        array=time,
                    ),
                    fits.Column(name="CADENCENO", format="J", array=cadenceno),
                    fits.Column(
                        name="TVCOL_RAW",
                        format=f"{2078 * 11}J",
                        null=-1,
                        dim="(11,2078)",
                        array=tvcol_raw,
                    ),
                ]
            )
            tv_hdu1 = fits.BinTableHDU.from_columns(tv_cols, name="COLLATERAL")
            tv_hdu1.header["NUM_FRM"] = num_frm
            tv_hdu1.header["INT_TIME"] = int_time
            tv_hdu1.header["FXDOFF"] = fxdoff
            tv_hdu1.header[f"MEANBLC{out}"] = meanblc[out]
            tv_hdu1.header[f"GAIN{out}"] = gains[out]
            tv_hdu1.header[f"READNOI{out}"] = 8.5

            tvcol_path = (
                out_dir / f"{prefix}-tvcol-{camera}-{ccd}-{out.lower()}-{run_id}-{col_tag}.fits"
            )
            fits.HDUList([fits.PrimaryHDU(header=col_hdr0), tv_hdu1]).writeto(
                tvcol_path, overwrite=True
            )

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
