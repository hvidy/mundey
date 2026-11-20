# Changelog

All notable changes to `mundey` are documented in this file.

## [2.0.1] - 2026-11-20

### Fixed
- **Physical shutterless row-smear weight (`_frame_transfer_row_weight`, `[MUND-01]`)**:
  Compute `READTIME / (INT_TIME * 2058)` directly whenever positive `READTIME` and `INT_TIME` cards exist in `HDU[1]`/`HDU[0]`, instead of plugging SPOC's `CRMITEN` duty-cycle `DEADC = 0.792` into `(1 - DEADC) / (DEADC * 2058)`. Only fall back to `DEADC` when `READTIME`/`INT_TIME` are absent, removing the `(CRBLKSZ - 2) / CRBLKSZ` cosmic-ray trimming factor when `CRMITEN` is `True`. Eliminates the `26.0x` row-smear weight inflation (`-56.7 e-/s` sky over-subtraction on Sirius Sector 6 and `26x` inflated `SMEAR_UNRESOLVED` fallback error).
- **Deferred bloom/saturation detection in `smear="automatic"` (`[MUND-02]`)**:
  Defer `eff_smear` resolution for stamps ending below CCD row 2048 until after pre-smear calibration and collateral `smrow` loading, automatically switching to `"alternate"` whenever `y1_ccd0 + ny >= 2048`, `collat_sm_sat_cols` is saturated, `raw_sat_2d` is saturated, or `_detect_bloom_columns` flags negative smear bloom.
- **`-blf.fits` cache provenance validation (`[MUND-03]`, `[MUND-P2]`)**:
  Stamp `CADMATCH` (`cadenceno`/`positional`) and a deterministic 16-hex-character SHA-256 `SRCMASK` digest (`"NONE"` when no `SOURCE_MASK` bits are set) in `write_smear_blf`, require exact provenance match in `_find_valid_smrow_blf` and `_alt_blf_has_more_cadences`, persist only matched cadences into `-alt-*-blf.fits`, and tolerate read-only cache directories with a `UserWarning` when persisting optional `-alt-` caches. Added `--build-alt` to `mundey prefetch --tpf` to pre-stage target-scoped alternate BLFs for read-only compute nodes.
- **Calibrated cutout desmearing on positive-sky and all-bloomed stamps (`desmear_calibrated_cutout`, `[MUND-04]`)**:
  Detect negative bloom streaks before restoring `FLUX_BKG` and against a robust non-negative 2D linear sky plane when `sky > |bloom|`, anchor all-bloomed collateral-free cutouts to their unsaturated clean-row lower envelope instead of leaving 100% of columns `SMEAR_UNRESOLVED`, and preserve `BACKAPP`/`BKGSUB` when `FLUX_BKG` is all-`NaN`.
- **`Lightkurve.cutout()` CCD coordinate handling (`_stamp_ccd_origin_1based`, `[MUND-05]`)**:
  Prefer `HDU[1]` `1CRV5P`/`2CRV5P` (and `1CRV4P`/`2CRV4P`) over stale `HDU[2]` `CRVAL1P`/`CRVAL2P` across `recalibrate_tpf`, `desmear_calibrated_cutout`, and `mundey prefetch --tpf`, synchronizing `HDU[2]` origin cards on output.
- **`HALOCVER = 1` contract validation (`[TBKG-03a]`)**:
  Reject booleans and `HALOCVER != 1` in `write_mundey_provenance`/`read_mundey_provenance`, and reject negative `SMEAR_ERR` values in `append_smear_err_column` and `read_smear_err_column`.
- **MAST S3 collateral and model selection (`[MUND-P1]`)**:
  Filter S3 collateral matches by `tpf_parts` run ID when available and select the latest reprocessing timestamp (`sorted(...)[-1]`) instead of the oldest (`[0]`).
- **Sirius example notebook (`notebooks/sirius_example.ipynb`, `[MUND-06]`)**:
  Updated Cell 8 to the `HALOCVER = 1` v2 functional workflow (`mundey.recalibrate_tpf` -> `tessbkgd.subtract_background` -> `halophot.extract_halo_lightcurve`) and executed all cells cleanly.

## [2.0.0] - 2026-11-19

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
  - Calibrated outputs preserve un-background-subtracted `FLUX` (`e-/s`), zero `FLUX_BKG` and `FLUX_BKG_ERR`, stamp `HALOCVER = 1`, `BKGSUB = False`, `BACKAPP = False`, and `SMEARMOD` (`"standard"`, `"alternate"`, or `"cutout"`) in `HDU[0]`, append column-coherent `(N_t, N_x)` `SMEAR_ERR` (`e-/s`) to `HDU[1]` (readable from `TargetPixelFile`, `HDUList`, or FITS paths via `read_smear_err` / `read_smear_err_column`), and write spatial `uint16` `PIXEL_QUALITY` (`STRAP_COL`, `SAT_BLEED`, `SMEAR_UNRESOLVED`, `SOURCE_MASK`, `BKG_USED`).
  - Standardized `resolve_cache_dir(cache_dir=None, *, strict_hpc=True, create=False)` across `mundey`, `tessbkgd`, and `halophot` with `$MUNDEY_CACHE_DIR` -> `$TESS_BRIGHT_CACHE_DIR` -> `$SCRATCH/.cache/tess_bright` -> `$XDG_CACHE_HOME/tess_bright` -> `~/.cache/tess_bright` precedence (including automatic rerouting when `ddir` resides under `/fred`) and symlink-resolved `/fred` protection enabled by default (`strict_hpc=True` raises `ValueError`; `strict_hpc=False` emits `UserWarning`).
  - `-blf.fits` caches store integer `CADENCENO` (`TFORM2 = 'J'`), calibrated `SMEAR` and `SMEAR_ERR` in `e-/frame`, per-cadence `SMROW_SAT` saturation masks, and target-scoped `ALT_UNRES` masks with full `CADENCENO` coverage validation across sliced and full-sector TPFs.

### Added
- Unified 2,048-column multi-output (`'A'`, `'B'`, `'C'`, `'D'`) TPF calibration across arbitrary 1-, 2-, 3-, and 4-output cutouts.
- Concavity-guarded collateral-anchored `smear="alternate"` (`smear="automatic"` selects `"alternate"` when stamps reach CCD row 2,048) with Bonferroni-calibrated bloom detection, inverse-wing-weighted Huber location estimation, and physical stellar wing smear accounting.
- `mundey.desmear_calibrated_cutout` for calibrated FFI / `TessCut` cutouts lacking `RAW_CNTS` and collateral `_col.fits` tables, stamping `CADMATCH = "none"`, `CADMISS = 0`, and `SMEARMOD = "cutout"` in `HDU[0]`.
- `mundey.archive` streaming collateral reader, anonymous MAST S3 resolver (`s3://stpubdata/tess/public/`), `mundey prefetch` / `mundey calibrate` CLI, and Swinburne Ngarrgu Tindebeek Slurm + Apptainer Stage-1 prefetch template (`hpc/mundey_prefetch.sbatch` with `halo.sif`, `module -q load apptainer`, `set -eo pipefail` for Lmod batch compatibility, `/aphid/scratch-3month/%u/logs/` Slurm logs, and `trevor` outbound-internet submission instructions).
