#!/usr/bin/env bash
# Run on the login node after transferring repo/data. Does not launch heavy jobs.
set -euo pipefail

source hpc/slurm_env.sh
require_somno_data_root

python -m py_compile \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc/plot_gc_stage_figures_dk5.py \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/plot_ssdi_stage_figures_dk.py \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/summarize_ssdi_results_dk.py \
  code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/run_python_dd_dk68.py

export PYTHONPATH="${SOMNO_REPO}/code/external/emergence-multiscale/src:${PYTHONPATH:-}"
python -c 'from multiscale import fit_var_auto, run_grid_parallel; print("Python DD engine imports OK")'

find "${SOMNODYNAMICS_DATA_ROOT}/results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk5" \
  -name '*_source_time_series.mat' | wc -l | awk '{print "DK5 source files:", $1}'
find "${SOMNODYNAMICS_DATA_ROOT}/results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk68" \
  -name '*_source_time_series_68.mat' | wc -l | awk '{print "DK68 source files:", $1}'

if [ "${1:-}" = "--matlab" ]; then
  run_matlab_batch "addpath('code/run_all'); setup_matlab_toolboxes(pwd); disp(sd_external_results_root(pwd));"
fi

echo "Smoke test complete. Check that DK68 source files total 143 before submitting."
