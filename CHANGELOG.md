# Changelog

All notable changes to `mundey` are documented in this file.

## [2.0.0] - 2026-11-14

### Breaking & Default Behaviour Changes
- **Functional `lightkurve` API (`mundey.recalibrate_tpf`, `mundey.desmear_calibrated_cutout`)**:
  `recalibrate_tpf(tpf, ...)` returns a new calibrated `TessTargetPixelFile` without mutating caller inputs, preserving sliced in-memory TPFs (`tpf[10:38]`, `tpf[mask]`, where `tpf.path` is an `astropy.io.fits.HDUList`) and falsy `quality_bitmask` values (`0` and `None`) verbatim. Legacy `mundey_tpf` and `mundey_tpf.calibrate` remain available as `@deprecated` compatibility shims delegating to the v2.0 engine (`use_simbad` is retired; `mundey` v2.0 leaves sky background intact and leaves background modelling to `tessbkgd` v2.0).
- **Exact `CADENCENO` collateral pairing (`cadence_match="cadenceno"`, `missing_cadence="nan"`)**:
  TPF cadences are now paired against `smrow`, `tvcol`, and `-blf.fits` tables by exact integer `CADENCENO` equality rather than row-index position, preventing cadence-shift thermal step artefacts and `ValueError` broadcast crashes on `QUALITY`-filtered TPFs (`N_tpf < N_col`). Unmatched TPF cadences emit a `UserWarning`, set `FLUX` / `FLUX_ERR` / `SMEAR_ERR` to `NaN`, and record `CADMATCH` and `CADMISS` in `HDU[0]` (`missing_cadence="error"` raises `ValueError`; `cadence_match="positional"` reproduces legacy row-by-row pairing).
- **Physical 64-bit JAX calibration order & order-statistic uncertainty propagation**:
  - Masked `RAW_CNTS == -1` and integer `TNULLn` sentinels to `NaN` in `float64` prior to arithmetic.
  - Scaled 2D black by effective co-added frame count `N_eff = factor * NUM_FRM` read independently from science and collateral headers (supporting both 2-minute `CRMITEN=True` and 20-second `CRMITEN=False` `_fast-col.fits`).
  - Fitted 1D black on clean virtual columns (`cols 3..9`) in locally scoped `@in_float64` JAX using masked discrete orthonormal Legendre least squares with cancellation-free orthogonal-complement AIC degree selection (`k in [20, 49]`), without mutating global `jax_enable_x64` or truncating calibrations to `int64`.
  - Inverted causal serial LDE undershoot in ADU via `jax.lax.associative_scan` along physical readout direction (`Outputs A/C` vs `Outputs B/D`) *before* evaluating `[0, maxDomain]`-clamped `linearity.xml` polynomials and channel gain.
  - Calibrated `SMEAR_ERR` with exact finite-sample unit-normal sample-median variances for `1 <= n <= 20` (`0.138326` at `n = 10` smear rows vs asymptotic `pi / 20 = 0.15708`) and gated empirical 10-row MAD variance (`emp_var > 2.5 * theo_var`) to prevent upward rectification bias on clean cadences while retaining unmodelled bleed/strap scatter.
  - Propagated the exact Sarhan & Greenberg (1956) symmetric 8-of-10 truncated-mean variance ratio (`8 * 0.105356 = 0.84285`) into 2-minute `CRMITEN=True` (`CRBLKSZ=10`) science `FLUX_ERR` and `smrow` `SMEAR_ERR` while preserving exact `N_eff = 0.8 * NUM_FRM` flux scaling.
- **`HALOCVER = 1` FITS contract (`BKGSUB = False`, `SMEAR_ERR`, `PIXEL_QUALITY`) & HPC `/fred` guard**:
  - Calibrated outputs preserve un-background-subtracted `FLUX` (`e-/s`), zero `FLUX_BKG` and `FLUX_BKG_ERR`, stamp `HALOCVER = 1`, `BKGSUB = False`, `BACKAPP = False`, and `SMEARMOD` (`"standard"`, `"alternate"`, or `"cutout"`) in `HDU[0]`, append column-coherent `(N_t, N_x)` `SMEAR_ERR` (`e-/s`) to `HDU[1]` (readable via `read_smear_err` / `read_smear_err_column`), and write spatial `uint16` `PIXEL_QUALITY` (`STRAP_COL`, `SAT_BLEED`, `SMEAR_UNRESOLVED`, `SOURCE_MASK`, `BKG_USED`).
  - Standardized `resolve_cache_dir(cache_dir=None, *, strict_hpc=True, create=False)` across `mundey`, `tessbkgd`, and `halophot` with `$MUNDEY_CACHE_DIR` -> `$TESS_BRIGHT_CACHE_DIR` -> `$SCRATCH/.cache/tess_bright` -> `$XDG_CACHE_HOME/tess_bright` -> `~/.cache/tess_bright` precedence and symlink-resolved `/fred` protection enabled by default (`strict_hpc=True` raises `ValueError`; `strict_hpc=False` emits `UserWarning`).
  - `-blf.fits` caches store integer `CADENCENO` (`TFORM2 = 'J'`), calibrated `SMEAR` and `SMEAR_ERR` in `e-/frame`, per-cadence `SMROW_SAT` saturation masks, and target-scoped `ALT_UNRES` masks.

### Added
- Unified 2,048-column multi-output (`'A'`, `'B'`, `'C'`, `'D'`) TPF calibration across arbitrary 1-, 2-, 3-, and 4-output cutouts.
- Concavity-guarded collateral-anchored `smear="alternate"` (`smear="automatic"` selects `"alternate"` when stamps reach CCD row 2,048) with Bonferroni-calibrated bloom detection, inverse-wing-weighted Huber location estimation, and physical stellar wing smear accounting.
- `mundey.desmear_calibrated_cutout` for calibrated FFI / `TessCut` cutouts lacking `RAW_CNTS` and collateral `_col.fits` tables, stamping `CADMATCH = "none"`, `CADMISS = 0`, and `SMEARMOD = "cutout"` in `HDU[0]`.
- `mundey.archive` streaming collateral reader, anonymous MAST S3 resolver (`s3://stpubdata/tess/public/`), `mundey prefetch` / `mundey calibrate` CLI, and Swinburne Ngarrgu Tindebeek Slurm + Apptainer Stage-1 prefetch template (`hpc/mundey_prefetch.sbatch` with `halo.sif`, `module -q load apptainer`, `set -eo pipefail` for Lmod batch compatibility, `/aphid/scratch-3month/%u/logs/` Slurm logs, and `trevor` outbound-internet submission instructions).
