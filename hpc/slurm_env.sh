#!/usr/bin/env bash
# Shared SLURM environment bootstrap for SomnoDynamics jobs.
#
# Override on the cluster before sbatch if needed, for example:
#   export SOMNO_REPO=/home/bmilinkovic/somnodynamics
#   export SOMNODYNAMICS_DATA_ROOT=/scratch/bmilinkovic/data_somnoemergence
#   export MATLAB_BIN=/path/to/matlab

set -euo pipefail

export SOMNO_REPO="${SOMNO_REPO:-/home/bmilinkovic/somnodynamics}"
export SOMNO_CONDA_ENV="${SOMNO_CONDA_ENV:-somnodynamics}"
export SOMNO_CONDA_PREFIX="${SOMNO_CONDA_PREFIX:-${HOME}/.conda/envs/${SOMNO_CONDA_ENV}}"
export SOMNODYNAMICS_DATA_ROOT="${SOMNODYNAMICS_DATA_ROOT:-${HOME}/data_somnoemergence}"
export MATLAB_BIN="${MATLAB_BIN:-matlab}"

# Keep accidental BLAS/OpenMP oversubscription under control. Submit scripts can
# raise these deliberately when MATLAB parpool workers are not used.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export MATLAB_NUM_THREADS="${MATLAB_NUM_THREADS:-1}"

export MPLCONFIGDIR="${MPLCONFIGDIR:-${SOMNO_REPO}/.cache/matplotlib}"
mkdir -p "${MPLCONFIGDIR}" "${SOMNO_REPO}/hpc/logs"
cd "${SOMNO_REPO}"

if command -v module >/dev/null 2>&1; then
  module load matlab 2>/dev/null || true
  module load anaconda/3.11 2>/dev/null || true
fi

if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck source=/dev/null
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck source=/dev/null
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
fi

if command -v conda >/dev/null 2>&1; then
  if [ -x "${SOMNO_CONDA_PREFIX}/bin/python" ]; then
    conda activate "${SOMNO_CONDA_PREFIX}"
  else
    conda activate "${SOMNO_CONDA_ENV}" 2>/dev/null || true
  fi
fi

echo "SOMNO_REPO=${SOMNO_REPO}"
echo "SOMNODYNAMICS_DATA_ROOT=${SOMNODYNAMICS_DATA_ROOT}"
echo "MATLAB_BIN=${MATLAB_BIN}"
echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
echo "MKL_NUM_THREADS=${MKL_NUM_THREADS}"
echo "OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS}"
echo "MATLAB_NUM_THREADS=${MATLAB_NUM_THREADS}"
if command -v python >/dev/null 2>&1; then
  python --version
  which python
fi
command -v "${MATLAB_BIN}" >/dev/null 2>&1 || {
  echo "ERROR: MATLAB_BIN not found: ${MATLAB_BIN}" >&2
  exit 2
}

require_somno_data_root() {
  if [ ! -d "${SOMNODYNAMICS_DATA_ROOT}/results" ]; then
    echo "ERROR: expected SomnoDynamics data root with results/:" >&2
    echo "       ${SOMNODYNAMICS_DATA_ROOT}" >&2
    echo "Set SOMNODYNAMICS_DATA_ROOT to the transferred data_somnoemergence root." >&2
    return 7
  fi
}

run_matlab_batch() {
  local code="$1"
  "${MATLAB_BIN}" -batch "maxNumCompThreads(str2double(getenv('MATLAB_NUM_THREADS'))); ${code}"
}
