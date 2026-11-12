#!/usr/bin/env bash
# Run on the login node after transferring repo/data. Does not launch heavy jobs.
set -euo pipefail

source hpc/slurm_env.sh
require_somno_data_root

python -m py_compile \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc/plot_gc_stage_figures_dk5.py \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/plot_ssdi_stage_figures_dk.py \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/summarize_ssdi_results_dk.py \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/run_python_dd_dk68.py \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/fit_python_var_dk68.py

export PYTHONPATH="${SOMNO_REPO}/code/external/emergence-multiscale/src:${PYTHONPATH:-}"
python -c 'from multiscale import fit_var_auto, run_grid_parallel; print("Python DD engine imports OK")'

MODEL_ROOT="${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz/model_inputs"
if [ ! -d "${MODEL_ROOT}" ]; then
  echo "ERROR: compact DK68 model directory is missing: ${MODEL_ROOT}" >&2
  exit 7
fi
MODEL_COUNT="$(find "${MODEL_ROOT}" -name '*_var_model.npz' -type f | wc -l | tr -d ' ')"
echo "DK68 fitted VAR models: ${MODEL_COUNT}"
if [ "${MODEL_COUNT}" -ne 143 ]; then
  echo "ERROR: expected 143 compact DK68 models before submitting the DD array." >&2
  exit 7
fi

if [ "${1:-}" = "--matlab" ]; then
  run_matlab_batch "addpath('code/run_all'); setup_matlab_toolboxes(pwd); disp(sd_external_results_root(pwd));"
fi

echo "Smoke test complete. The 143 compact models are ready for DD."
