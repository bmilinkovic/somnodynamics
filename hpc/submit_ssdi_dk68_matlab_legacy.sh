#!/usr/bin/env bash
#SBATCH --job-name=somno-ssdi-dk68
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --exclusive
#SBATCH --time=14-00:00:00
#SBATCH --output=hpc/logs/%x-%j.out
#SBATCH --error=hpc/logs/%x-%j.err

set -euo pipefail
source hpc/slurm_env.sh
require_somno_data_root

# One thread per MATLAB worker avoids nested oversubscription.
export MATLAB_NUM_THREADS="${MATLAB_NUM_THREADS:-1}"
export SSDI_MAX_WORKERS="${SSDI_MAX_WORKERS:-${SLURM_CPUS_PER_TASK:-64}}"

MODEL_COUNT="$(find "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz/model_inputs/subject_measures" -name '*_subject_measures.mat' 2>/dev/null | wc -l | tr -d ' ')"
if [ "${MODEL_COUNT}" -lt 143 ]; then
  echo "ERROR: DK68 SSDI model inputs incomplete: ${MODEL_COUNT}/143." >&2
  echo "Run hpc/submit_ssdi_dk68_model_inputs_lrt30.sh first, then resubmit this job." >&2
  exit 8
fi

run_matlab_batch "addpath('code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi'); run_ssdi_sleep_stage_pipeline_dk('mode','stage_ssdi','atlas','dk68','recompute',false,'parallel_enable',true,'parallel_workers',${SSDI_MAX_WORKERS},'nrunsp',100,'nitersp',10000,'niterso',10000);"

SSDI_DIR="${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz"
SSDI_COUNT="$(find "${SSDI_DIR}" -name '*_ssdi_spectral_domain_results.mat' -type f | wc -l | tr -d ' ')"
if [ "${SSDI_COUNT}" -ne 143 ]; then
  echo "ERROR: DK68 SSDI results incomplete: ${SSDI_COUNT}/143." >&2
  exit 9
fi
echo "DK68 SSDI results complete: ${SSDI_COUNT}/143."

python code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/summarize_ssdi_results_dk.py \
  --spectral-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz" \
  --out-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz"
python code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/plot_ssdi_stage_figures_dk.py \
  --spectral-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz" \
  --long-csv "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz/ssdi_subject_long.csv" \
  --out-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz/figures" \
  --pairwise-root "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz/downstream_stage_pairwise"
