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

DK5_DIR="${SOMNODYNAMICS_DATA_ROOT}/results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk5"
DK68_DIR="${SOMNODYNAMICS_DATA_ROOT}/results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk68"
if [ -d "${DK5_DIR}" ]; then
  find "${DK5_DIR}" -name '*_source_time_series.mat' -type f |
    wc -l | awk '{print "DK5 source files:", $1}'
fi
if [ ! -d "${DK68_DIR}" ]; then
  echo "ERROR: DK68 source directory is missing: ${DK68_DIR}" >&2
  exit 7
fi
DK68_COUNT="$(find "${DK68_DIR}" -name '*_source_time_series_68.mat' -type f | wc -l | tr -d ' ')"
echo "DK68 source files: ${DK68_COUNT}"
if [ "${DK68_COUNT}" -ne 143 ]; then
  echo "ERROR: expected 143 DK68 source files before submitting the DD array." >&2
  exit 7
fi

if [ "${1:-}" = "--matlab" ]; then
  run_matlab_batch "addpath('code/run_all'); setup_matlab_toolboxes(pwd); disp(sd_external_results_root(pwd));"
fi

echo "Smoke test complete. Check that DK68 source files total 143 before submitting."
