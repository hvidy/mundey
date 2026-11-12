# Changelog

All notable changes to `mundey` are documented in this file.

## [2.0.0a0] - 2026-11-12

### Breaking & Default Behaviour Changes
- **Functional `lightkurve` API (`mundey.recalibrate_tpf`, `mundey.desmear_calibrated_cutout`)**:
  `recalibrate_tpf(tpf, ...)` returns a new calibrated `TessTargetPixelFile` without mutating caller inputs. Legacy `mundey_tpf` and `mundey_tpf.calibrate` remain available as `@deprecated` compatibility shims delegating to the v2.0 engine (`use_simbad` is retired; `mundey` v2.0 leaves sky background intact and leaves background modelling to `tessbkgd` v2.0).
- **Exact `CADENCENO` collateral pairing (`cadence_match="cadenceno"`, `missing_cadence="nan"`)**:
  TPF cadences are now paired against `smrow`, `tvcol`, and `-blf.fits` tables by exact integer `CADENCENO` equality rather than row-index position, preventing cadence-shift thermal step artefacts and `ValueError` broadcast crashes on `QUALITY`-filtered TPFs (`N_tpf < N_col`). Unmatched TPF cadences emit a `UserWarning`, set `FLUX` / `FLUX_ERR` / `SMEAR_ERR` to `NaN`, and record `CADMATCH` and `CADMISS` in `HDU[0]` (`missing_cadence="error"` raises `ValueError`; `cadence_match="positional"` reproduces legacy row-by-row pairing).
- **Physical 64-bit JAX calibration order & non-truncated black subtraction**:
  - Masked `RAW_CNTS == -1` and integer `TNULLn` sentinels to `NaN` in `float64` prior to arithmetic.
  - Scaled 2D black by effective co-added frame count `N_eff = factor * NUM_FRM` read independently from science and collateral headers (supporting both 2-minute `CRMITEN=True` and 20-second `CRMITEN=False` `_fast-col.fits`).
  - Fitted 1D black on clean virtual columns (`cols 3..9`) in locally scoped `@in_float64` JAX using masked discrete orthonormal Legendre least squares with cancellation-free orthogonal-complement AIC degree selection (`k in [20, 49]`), without mutating global `jax_enable_x64` or truncating calibrations to `int64`.
  - Inverted causal serial LDE undershoot in ADU via `jax.lax.associative_scan` along physical readout direction (`Outputs A/C` vs `Outputs B/D`) *before* evaluating `[0, maxDomain]`-clamped `linearity.xml` polynomials and channel gain.
- **`HALOCVER = 1` FITS contract (`BKGSUB = False`, `SMEAR_ERR`, `PIXEL_QUALITY`)**:
  - Calibrated outputs preserve un-background-subtracted `FLUX` (`e-/s`), zero `FLUX_BKG` and `FLUX_BKG_ERR`, stamp `HALOCVER = 1`, `BKGSUB = False`, and `BACKAPP = False`, append column-coherent `(N_t, N_x)` `SMEAR_ERR` (`e-/s`) to `HDU[1]`, and write spatial `uint16` `PIXEL_QUALITY` (`STRAP_COL`, `SAT_BLEED`, `SMEAR_UNRESOLVED`, `SOURCE_MASK`, `BKG_USED`).
  - `-blf.fits` caches store integer `CADENCENO` (`TFORM2 = 'J'`), calibrated `SMEAR` and `SMEAR_ERR` in `e-/frame`, per-cadence `SMROW_SAT` saturation masks, and target-scoped `ALT_UNRES` masks.

### Added
- Unified 2,048-column multi-output (`'A'`, `'B'`, `'C'`, `'D'`) TPF calibration across arbitrary 1-, 2-, 3-, and 4-output cutouts.
- Concavity-guarded collateral-anchored `smear="alternate"` (`smear="automatic"` selects `"alternate"` when stamps reach CCD row 2,048) with Bonferroni-calibrated bloom detection, inverse-wing-weighted Huber location estimation, and physical stellar wing smear accounting.
- `mundey.desmear_calibrated_cutout` for calibrated FFI / `TessCut` cutouts lacking `RAW_CNTS` and collateral `_col.fits` tables.
- `mundey.archive` streaming collateral reader, anonymous MAST S3 resolver (`s3://stpubdata/tess/public/`), `mundey prefetch` / `mundey calibrate` CLI, and Swinburne Ngarrgu Tindebeek Slurm + Apptainer templates (`hpc/`).
