# SomnoDynamics HPC Submitters

These scripts mirror the style of `/Users/borjan/CNRS/projects/TVBToolkit/hpc`:
shared `slurm_env.sh`, logs in `hpc/logs`, long `workq` jobs, explicit
thread caps, and resumable analysis commands.

## Required Cluster Variables

Set these before `sbatch` if the defaults do not match the institutional
machine:

```bash
export SOMNO_REPO=/home/bmilinkovic/somnodynamics
export SOMNODYNAMICS_DATA_ROOT=/path/to/data_somnoemergence
export MATLAB_BIN=matlab
export SOMNO_CONDA_ENV=somnodynamics
```

`SOMNODYNAMICS_DATA_ROOT` must contain the transferred external-drive-style
layout, including:

```text
results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk5
results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk68
```

## Smoke Test

Run on the login node after transfer:

```bash
bash hpc/smoke_test_login.sh
```

## Submit Order

For the current DK5/DK68 LRT30 branch:

```bash
sbatch hpc/submit_gc_dk5_lrt30.sh
sbatch hpc/submit_ssdi_dk5_lrt30.sh
```

DK68 DD now runs directly from source reconstructed epochs using the vendored
Emergence Engine Python package. No MATLAB model-input or SSDI job is needed
for this branch. After transferring the repository and source data to the
cluster, run `bash hpc/create_conda_env.sh` once, then submit from the
repository root:

```bash
mkdir -p hpc/logs
sbatch --array=0-0 hpc/submit_ssdi_dk68_lrt30.sh
```

This first task is a pilot: inspect its log, selected VAR order, memory use
and wall time. Once it finishes satisfactorily, submit the full array and
summary:

```bash
mkdir -p hpc/logs
dd_job=$(sbatch --parsable hpc/submit_ssdi_dk68_lrt30.sh)
summary_job=$(sbatch --parsable --dependency=afterok:${dd_job} hpc/submit_ssdi_dk68_python_summary.sh)
echo "DK68 DD array: ${dd_job}; summary: ${summary_job}"
squeue -j "${dd_job},${summary_job}"
```

The array has 143 zero-based tasks and runs up to four tasks concurrently by
default. Each task fits one subject-stage VAR model and optimises its 50
dimension × band blocks. The runner checks that all 143 source files are
present before starting and skips existing Python result files on resubmission.
If a task stops part-way through, it reuses its fitted VAR model and finished
dimension–band checkpoints on the next submission.
The summary job starts only after the array succeeds and checks for 143
results. Logs in `hpc/logs/` show the VAR order and progress through the DD
blocks. To use one node at a time, submit with `sbatch --array=0-142%1 ...`.

Python results are separate from the MATLAB archive under
`$SOMNODYNAMICS_DATA_ROOT/results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz`.
Each stage has the same `*_ssdi_spectral_domain_results.mat` layout consumed
by the existing summary and figure scripts. The old MATLAB submitter is
`submit_ssdi_dk68_matlab_legacy.sh`; the MATLAB model-input submitter remains
available for historical comparison.

DK68 GC is provided but is not part of the currently preferred local branch:

```bash
sbatch hpc/submit_gc_dk68_lrt30.sh
```

## Resource Logic

- Each Python DD array task requests a high-memory exclusive node. The 50
  independent dimension × band blocks can use up to 50 Python worker processes.
- `SSDI_MAX_WORKERS` limits processes per array task; BLAS/OpenMP threads are
  capped at one to avoid oversubscription.
- LRT/LWR VAR fitting occurs once per subject-stage file. The proxy and
  canonical spectral optimisers then run for every requested block.
- DK68 with 100 restarts and 10,000 iterations per stage is computationally
  substantial. Start with one array task and inspect its runtime before
  releasing the full array if your allocation is limited.

## Current Analysis Parameters

- Epoching/source data: current 2 s, 256 Hz source-reconstructed branch.
- GC/VAR model selection: `momax=30`, `mosel=LRT`.
- SSDI: Emergence Engine proxy followed by canonical spectral DD.
- DK68 macro dimensions: `m=1:10`.
- Frequency bands: delta, theta, alpha, sigma, beta; no gamma.
