#!/usr/bin/env bash
#SBATCH --job-name=somno-gc-dk5
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
WORKERS="${GC_WORKERS:-${SLURM_CPUS_PER_TASK:-48}}"

run_matlab_batch "addpath('code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc'); run_gc_sleep_stage_pipeline_dk5('mode','stage_mvgc2','atlas','dk5','recompute_subject_measures',false,'recompute_group_stats',true,'parallel_enable',true,'parallel_workers',${WORKERS},'momax',30,'mosel','LRT');"

python code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc/plot_gc_stage_figures_dk5.py
