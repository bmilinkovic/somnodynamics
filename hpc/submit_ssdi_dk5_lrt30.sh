#!/usr/bin/env bash
#SBATCH --job-name=somno-ssdi-dk5
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0
#SBATCH --time=14-00:00:00
#SBATCH --output=hpc/logs/%x-%j.out
#SBATCH --error=hpc/logs/%x-%j.err

set -euo pipefail
source hpc/slurm_env.sh
require_somno_data_root

export MATLAB_NUM_THREADS="${MATLAB_NUM_THREADS:-1}"
export SSDI_MAX_WORKERS="${SSDI_MAX_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"

run_matlab_batch "addpath('code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi'); run_ssdi_sleep_stage_pipeline_dk('mode','stage_ssdi','atlas','dk5','recompute',false,'parallel_enable',true,'parallel_workers',${SSDI_MAX_WORKERS},'nrunsp',100,'nitersp',10000,'niterso',10000);"

python code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/summarize_ssdi_results_dk.py \
  --spectral-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz" \
  --out-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz"
python code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/plot_ssdi_stage_figures_dk.py \
  --spectral-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz" \
  --long-csv "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/ssdi_subject_long.csv" \
  --out-dir "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/figures" \
  --pairwise-root "${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/downstream_stage_pairwise"
