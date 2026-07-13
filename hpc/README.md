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
sbatch hpc/submit_ssdi_dk68_model_inputs_lrt30.sh
sbatch hpc/submit_ssdi_dk68_lrt30.sh
```

DK68 GC is provided but is not part of the currently preferred local branch:

```bash
sbatch hpc/submit_gc_dk68_lrt30.sh
```

## Resource Logic

- DK68 model input fitting is the memory-heavy VAR/LRT stage.
- The DK68 jobs request a high-memory exclusive node, matching the aggressive
  style used in the TVBToolkit scripts.
- DK68 model input fitting defaults to `DK68_MODEL_INPUT_WORKERS=$SLURM_CPUS_PER_TASK`.
- DK68 SSDI defaults to `SSDI_MAX_WORKERS=$SLURM_CPUS_PER_TASK`.
- Each MATLAB worker is kept to one computational thread by default. That is
  not conservative; it is the realistic way to use a full node without nested
  BLAS oversubscription turning 64 workers into hundreds of fighting threads.
- All heavy scripts use `recompute=false`, so completed subject-stage files are
  skipped when jobs are resumed.

## Current Analysis Parameters

- Epoching/source data: current 2 s, 256 Hz source-reconstructed branch.
- GC/VAR model selection: `momax=30`, `mosel=LRT`.
- SSDI: canonical SomnoDynamics/DMT-faithful spectral SSDI.
- DK68 SSDI macro dimensions: configured in MATLAB as `m=1:10`.
- Frequency bands: delta, theta, alpha, sigma, beta; no gamma.
