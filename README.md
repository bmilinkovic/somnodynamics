# SomnoDynamics

SomnoDynamics is a reproducible multiscale EEG analysis toolbox for healthy sleep data.

## Scope

SomnoDynamics provides:
- preprocessing
- sensor-space spectral power analysis (absolute + relative)
- source reconstruction
- canonical coarse-graining
- stage-based Granger causality (GC)
- stage-based dynamical dependence / SSDI
- GC-SSDI cross-analysis
- publication-oriented figures and machine-readable statistics

## Raw Data Policy

Raw data stays outside the repository at:

`/Volumes/ex_data/data_somnoemergence/raw`

SomnoDynamics reads from that location and writes all derived outputs under this repository.

## Stage Labels

Internal stage labels:
- `W_`
- `N1`
- `N2`
- `N3`
- `R_`

Standardized display labels used in outputs:
- `W_ -> Wake`
- `N1 -> N1`
- `N2 -> N2`
- `N3 -> N3`
- `R_ -> REM`

## Repository Layout

- `code/legacy/` preserved legacy healthy-sleep scripts
- `code/preprocessing/` sensor preprocessing
- `code/sensor_power/` sensor spectral power analysis
- `code/src_recon/` source reconstruction wrappers
  - primary pipeline: anatomical46 reconstruction then anatomical7 coarse-graining
- `code/coarse_graining/` canonical region mapping assets
- `code/gc/` GC pipeline wrappers
  - primary GC module: `code/gc/stage_mvgc2_anatomical7/`
- `code/external/` vendored external toolboxes (includes `MVGC2`)
- `code/ssdi/` SSDI pipeline and summaries
  - primary SSDI module: `code/ssdi/stage_ssdi_anatomical7/`
- `code/stats/` repeated-measures statistics
- `code/cross_analysis/` GC-SSDI cross-analysis
- `code/run_all/` orchestration + shared config
- `results/` stage-specific outputs
- `docs/` migration and methodology docs

## Run

Stage runner:

```bash
python3 code/run_all/run_pipeline.py all
```

Examples:

```bash
python3 code/run_all/run_pipeline.py preprocessing
python3 code/run_all/run_pipeline.py sensor_power
python3 code/run_all/run_pipeline.py source_recon
python3 code/run_all/run_pipeline.py all --skip-matlab
```

The DK68 stage-wise DD run uses the vendored Emergence Engine Python package
from source reconstructed epochs. Its HPC array and dependent summary job are
documented in [hpc/README.md](hpc/README.md). The DK5 and legacy GC/SSDI
paths still use MATLAB and initialise the existing toolboxes.

## Key Docs

- `MIGRATION_PLAN_sleep_to_dmt_style.md`
- `docs/PREPROCESSING.md`
- `docs/SOURCE_RECONSTRUCTION.md`
- `docs/CANONICAL_REGION_MAPPING.md`
- `docs/FIGURES.md`
- `docs/VALIDATION.md`
- `docs/IMPLEMENTATION_SUMMARY.md`
- `docs/gc_pipeline_dmt_audit.md`
- `docs/ssdi_pipeline_dmt_audit.md`
