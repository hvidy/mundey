#!/bin/bash
#SBATCH --job-name=mundey-prefetch
#SBATCH --account=oz440
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --tmp=25G
#SBATCH --time=01:00:00
#SBATCH --array=0-15%4
#SBATCH --output=/aphid/scratch-3month/%u/logs/%x-%A_%a.out
#SBATCH --error=/aphid/scratch-3month/%u/logs/%x-%A_%a.err
#
# Swinburne Ngarrgu Tindebeek (NT) / OzSTAR Stage-1 collateral prefetcher.
# Deduplicates the 16 (camera, ccd) collateral sets for a TESS sector via
# anonymous MAST S3 into local NVMe $JOBFS, streams 4,096-cadence chunks (< 2 GB
# RSS) through @in_float64 JAX kernels, and publishes lockless -blf.fits caches
# to $SCRATCH/.cache/tess_bright (never /fred).
#
# IMPORTANT: Swinburne NT compute nodes have NO outbound internet access.
# Stage-1 MAST S3 prefetching (mundey_prefetch.sbatch) MUST be submitted from a
# `trevor` node (creating the scratch Slurm log directory first so array logs
# never consume /fred/oz440 inodes or /home quota; cd to your mundey checkout if
# outside $HOME):
#   ssh trevor 'mkdir -p /aphid/scratch-3month/$USER/logs && sbatch --export=ALL,SECTOR=6 hpc/mundey_prefetch.sbatch'
#
# Or directly from a shell on trevor for Sector 6 (2-min) or Sector 27 (both 2-min & 20-s):
#   mkdir -p /aphid/scratch-3month/$USER/logs
#   sbatch --export=ALL,SECTOR=6 hpc/mundey_prefetch.sbatch
#   sbatch --export=ALL,SECTOR=27,CADENCE=both hpc/mundey_prefetch.sbatch

# Omit -u: Lmod's `module load` reads unset PS1 in non-interactive Slurm shells on NT.
set -eo pipefail

: "${SECTOR:?Set SECTOR=<int> via sbatch --export=ALL,SECTOR=...}"
: "${USER:=$(id -un)}"
CADENCE="${CADENCE:-2min}"

export SCRATCH="${SCRATCH:-/aphid/scratch-3month/$USER}"
export MYFRED="${MYFRED:-/fred/oz440/$USER}"
export SIF="${SIF:-$MYFRED/containers/halo.sif}"

# Zero-/fred cache hierarchy on /aphid scratch (including Apptainer cache & session tmpdir)
export TESS_BRIGHT_CACHE_DIR="${TESS_BRIGHT_CACHE_DIR:-$SCRATCH/.cache/tess_bright}"
export MUNDEY_CACHE_DIR="${MUNDEY_CACHE_DIR:-$TESS_BRIGHT_CACHE_DIR}"
export JAX_COMPILATION_CACHE_DIR="${JAX_COMPILATION_CACHE_DIR:-$SCRATCH/jax_cache}"
export LIGHTKURVE_CACHE_DIR="${SCRATCH}/.cache/lightkurve"
export ASTROPY_CACHE_DIR="${SCRATCH}/.cache/astropy"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$SCRATCH/.apptainer/cache}"
mkdir -p "$TESS_BRIGHT_CACHE_DIR" "$MUNDEY_CACHE_DIR" "$JAX_COMPILATION_CACHE_DIR" "$LIGHTKURVE_CACHE_DIR" "$ASTROPY_CACHE_DIR" "$APPTAINER_CACHEDIR"

# CPU-only JAX configuration & strict thread pinning (omit --nv)
CPUS="${SLURM_CPUS_PER_TASK:-4}"
export JAX_PLATFORMS=cpu
export JAX_DEFAULT_MATMUL_PRECISION=highest
export OMP_NUM_THREADS="$CPUS"
export OPENBLAS_NUM_THREADS="$CPUS"
export MKL_NUM_THREADS="$CPUS"
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads=${CPUS}"

STAGE_DIR="${JOBFS:-$SCRATCH/tmp/${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}}"
export APPTAINER_TMPDIR="${JOBFS:-$STAGE_DIR}"
export APPTAINER_BIND="/fred,/aphid/scratch-3month${JOBFS:+,$JOBFS}"
export APPTAINER_BINDPATH="/fred,/aphid/scratch-3month${JOBFS:+,$JOBFS}"

OUT="${STAGE_DIR}/task_manifest"
ARCHIVE_DIR="${MYFRED}/outputs/${SLURM_JOB_NAME:-mundey-prefetch}-${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
TASK_TAG="${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$STAGE_DIR" "$APPTAINER_TMPDIR" "$OUT"

# Mandatory EXIT-trapped single-inode .tar archival to /fred
trap 'mkdir -p "$ARCHIVE_DIR" && tar -cf "${ARCHIVE_DIR}/task_${TASK_TAG}.tar" -C "$OUT" .' EXIT

echo "sector=${SECTOR} task_id=${TASK_TAG} cadence=${CADENCE} host=$(hostname)" > "${OUT}/prefetch_meta.txt"

if [[ -f "$SIF" ]]; then
    module -q load apptainer
    apptainer exec "$SIF" mundey prefetch \
        --sector "$SECTOR" \
        --task-id "$TASK_TAG" \
        --cadence "$CADENCE" \
        --cache-dir "$MUNDEY_CACHE_DIR" \
        --stage-dir "$STAGE_DIR" \
        --chunk-rows 4096 \
        --strict-hpc \
        --verbose | tee "${OUT}/prefetch.log"
else
    mundey prefetch \
        --sector "$SECTOR" \
        --task-id "$TASK_TAG" \
        --cadence "$CADENCE" \
        --cache-dir "$MUNDEY_CACHE_DIR" \
        --stage-dir "$STAGE_DIR" \
        --chunk-rows 4096 \
        --strict-hpc \
        --verbose | tee "${OUT}/prefetch.log"
fi
