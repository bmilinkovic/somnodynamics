#!/usr/bin/env bash
#SBATCH --job-name=somno-dk68-var
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

# DK68 VAR/LRT fitting is memory-heavy. On HPC we use the high-memory node
# aggressively, but keep one computational thread per worker to avoid nested
# BLAS/MATLAB oversubscription.
export MATLAB_NUM_THREADS="${MATLAB_NUM_THREADS:-1}"
WORKERS="${DK68_MODEL_INPUT_WORKERS:-${SLURM_CPUS_PER_TASK:-64}}"

if [ "${WORKERS}" -le 1 ]; then
  PARALLEL_ENABLE="false"
else
  PARALLEL_ENABLE="true"
fi

run_matlab_batch "addpath('code/run_all'); addpath('code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/stage_ssdi_dk68'); prepare_dk68_ssdi_model_inputs('parallel_enable',${PARALLEL_ENABLE},'parallel_workers',${WORKERS},'recompute',false,'momax',30,'mosel','LRT');"
