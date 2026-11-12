#!/usr/bin/env bash
#SBATCH --job-name=mundey-prefetch
#SBATCH --account=oz440
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --tmp=25G
#SBATCH --time=01:00:00
#SBATCH --array=0-15%4
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#
# Swinburne Ngarrgu Tindebeek (NT) / OzSTAR Stage-1 collateral prefetcher.
# Deduplicates the 16 (camera, ccd) collateral sets for a TESS sector via
# anonymous MAST S3 into local NVMe $JOBFS, streams 4,096-cadence chunks (< 2 GB
# RSS) through @in_float64 JAX kernels, and publishes lockless -blf.fits caches
# to $SCRATCH/.cache/tess_bright (never /fred).
#
# Submit for Sector 6 (2-min) or Sector 27 (both 2-min & 20-s):
#   mkdir -p logs
#   sbatch --export=ALL,SECTOR=6 hpc/mundey_prefetch.sbatch
#   sbatch --export=ALL,SECTOR=27,CADENCE=both hpc/mundey_prefetch.sbatch

set -euo pipefail

: "${SECTOR:?Set SECTOR=<int> via sbatch --export=ALL,SECTOR=...}"
CADENCE="${CADENCE:-2min}"

export SCRATCH="${SCRATCH:-/aphid/scratch-3month/$USER}"
export MYFRED="${MYFRED:-/fred/oz440/$USER}"
export SIF="${SIF:-$MYFRED/containers/virgil.sif}"

# Zero-/fred cache hierarchy on /aphid scratch
export MUNDEY_CACHE_DIR="${MUNDEY_CACHE_DIR:-$SCRATCH/.cache/tess_bright}"
export JAX_COMPILATION_CACHE_DIR="${JAX_COMPILATION_CACHE_DIR:-$SCRATCH/jax_cache}"
export LIGHTKURVE_CACHE_DIR="${SCRATCH}/.cache/lightkurve"
export ASTROPY_CACHE_DIR="${SCRATCH}/.cache/astropy"
mkdir -p "$MUNDEY_CACHE_DIR" "$JAX_COMPILATION_CACHE_DIR" "$LIGHTKURVE_CACHE_DIR" "$ASTROPY_CACHE_DIR"

# EPYC Milan CPU-only JAX configuration & strict thread pinning (omit --nv)
CPUS="${SLURM_CPUS_PER_TASK:-4}"
export JAX_PLATFORMS=cpu
export JAX_DEFAULT_MATMUL_PRECISION=highest
export OMP_NUM_THREADS="$CPUS"
export OPENBLAS_NUM_THREADS="$CPUS"
export MKL_NUM_THREADS="$CPUS"
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads=${CPUS}"

export APPTAINER_BINDPATH="/fred,/aphid${JOBFS:+,$JOBFS}"

STAGE_DIR="${JOBFS:-$SCRATCH/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
OUT="${STAGE_DIR}/task_manifest"
ARCHIVE_DIR="${MYFRED}/outputs/${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"
TASK_TAG="${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$STAGE_DIR" "$OUT"

# Mandatory EXIT-trapped single-inode .tar archival to /fred
trap 'mkdir -p "$ARCHIVE_DIR" && tar -cf "${ARCHIVE_DIR}/task_${TASK_TAG}.tar" -C "$OUT" .' EXIT

echo "sector=${SECTOR} task_id=${TASK_TAG} cadence=${CADENCE} host=$(hostname)" > "${OUT}/prefetch_meta.txt"

if [[ -f "$SIF" ]]; then
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
