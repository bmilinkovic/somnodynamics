#!/usr/bin/env bash
#SBATCH --job-name=somno-dd-dk68
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --exclusive
#SBATCH --time=14-00:00:00
#SBATCH --array=0-142%4
#SBATCH --output=hpc/logs/%x-%A_%a.out
#SBATCH --error=hpc/logs/%x-%A_%a.err

set -euo pipefail
source hpc/slurm_env.sh
require_somno_data_root

# The vendored package travels with the SomnoDynamics repository to the HPC.
export PYTHONPATH="${SOMNO_REPO}/code/external/emergence-multiscale/src:${PYTHONPATH:-}"
python -c 'from multiscale import fit_var_auto, run_grid_parallel'

WORKERS="${SSDI_MAX_WORKERS:-${SLURM_CPUS_PER_TASK:-64}}"
if [ "${WORKERS}" -gt 50 ]; then
  WORKERS=50
fi
TASK_INDEX="${SLURM_ARRAY_TASK_ID:?submit this script as a Slurm array}"

echo "DK68 Python DD: task ${TASK_INDEX}/142, ${WORKERS} block workers"
python -u code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/run_python_dd_dk68.py \
  --task-index "${TASK_INDEX}" \
  --expected-count 143 \
  --fs 256 --fres 256 --pmax 30 --criterion LRT --regmode LWR \
  --dimensions 1 2 3 4 5 6 7 8 9 10 \
  --runs 100 --pre-iters 10000 --final-iters 10000 \
  --workers "${WORKERS}"
