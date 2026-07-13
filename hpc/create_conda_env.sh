#!/usr/bin/env bash
# Optional one-time setup for Python plotting/statistics on the cluster.
set -euo pipefail

ENV_NAME="${SOMNO_CONDA_ENV:-somnodynamics}"
ENV_PREFIX="${SOMNO_CONDA_PREFIX:-${HOME}/.conda/envs/${ENV_NAME}}"

if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck source=/dev/null
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck source=/dev/null
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
else
  if command -v module >/dev/null 2>&1; then
    module load anaconda/3.11 2>/dev/null || true
  fi
  CONDA_BASE="$(conda info --base)"
  # shellcheck source=/dev/null
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
fi

if [ -d "${ENV_PREFIX}" ]; then
  echo "Using existing Conda environment: ${ENV_PREFIX}"
else
  conda create -y -p "${ENV_PREFIX}" python=3.11 pip
fi

conda activate "${ENV_PREFIX}"
python -m pip install --upgrade pip setuptools wheel
python -m pip install numpy scipy pandas matplotlib seaborn statsmodels mne h5py pymatreader
python - <<'PY'
import sys
print(sys.version)
for mod in ("numpy", "scipy", "pandas", "matplotlib", "mne", "h5py"):
    m = __import__(mod)
    print(mod, getattr(m, "__version__", "unknown"))
PY
