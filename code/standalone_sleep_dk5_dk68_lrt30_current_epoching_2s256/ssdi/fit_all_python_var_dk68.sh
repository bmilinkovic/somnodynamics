#!/usr/bin/env bash
# Fit all DK68 source files serially, one Python process per recording.
# Separate processes release each subject-stage recording before loading the
# next one. Completed compact models are skipped on rerun.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${REPO_ROOT}"
export SOMNODYNAMICS_DATA_ROOT="${SOMNODYNAMICS_DATA_ROOT:-/Volumes/ex_data/data_somnoemergence}"
export OPENBLAS_NUM_THREADS="${SOMNO_FIT_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${SOMNO_FIT_THREADS:-1}"
export OMP_NUM_THREADS="${SOMNO_FIT_THREADS:-1}"
export MKL_NUM_THREADS="${SOMNO_FIT_THREADS:-1}"

SCRIPT="code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/fit_python_var_dk68.py"
failed=0
for index in $(seq 0 142); do
  echo "===== DK68 model $(($index + 1))/143; task index $index ====="
  if ! python -u "${SCRIPT}" --task-index "${index}" --expected-count 143; then
    echo "FAILED task index ${index}; continuing with the next recording" >&2
    failed=$((failed + 1))
  fi
done

MODEL_ROOT="${SOMNODYNAMICS_DATA_ROOT}/results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz/model_inputs"
MODEL_COUNT="$(find "${MODEL_ROOT}" -name '*_var_model.npz' -type f | wc -l | tr -d ' ')"
echo "Finished DK68 model files: ${MODEL_COUNT}/143"
echo "Failed task attempts: ${failed}"
[ "${MODEL_COUNT}" -eq 143 ]
