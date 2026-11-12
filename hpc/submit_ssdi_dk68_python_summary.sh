#!/usr/bin/env bash
#SBATCH --job-name=somno-dd-dk68-summary
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=hpc/logs/%x-%j.out
#SBATCH --error=hpc/logs/%x-%j.err

set -euo pipefail
source hpc/slurm_env.sh
require_somno_data_root

SSDI_DIR="${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz"
SSDI_COUNT="$(find "${SSDI_DIR}" -maxdepth 2 -name '*_ssdi_spectral_domain_results.mat' -type f | wc -l | tr -d ' ')"
if [ "${SSDI_COUNT}" -ne 143 ]; then
  echo "ERROR: Python DK68 DD results incomplete: ${SSDI_COUNT}/143." >&2
  exit 9
fi
echo "Python DK68 DD results complete: ${SSDI_COUNT}/143."

python code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/summarize_ssdi_results_dk.py \
  --spectral-dir "${SSDI_DIR}" --out-dir "${SSDI_DIR}"
python code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/plot_ssdi_stage_figures_dk.py \
  --spectral-dir "${SSDI_DIR}" \
  --long-csv "${SSDI_DIR}/ssdi_subject_long.csv" \
  --out-dir "${SSDI_DIR}/figures" \
  --pairwise-root "${SSDI_DIR}/downstream_stage_pairwise"
