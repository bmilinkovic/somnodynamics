# SomnoDynamics Stage SSDI (Anatomical7)

This module ports DMT `ssdi_anatomical7` methodology to SomnoDynamics sleep stages.

## Design
- Input models: `A`, `V` from stage-wise GC subject measures
- Stages: `W_`, `N1`, `N2`, `N3`, `R_`
- Source space: anatomical7
- Optimizer: `opt_gd_dds_mruns` (SSDI toolbox)

## Files
- `sd_ssdi_config_anatomical7.m`
- `run_stage_ssdi_anatomical7_pipeline.m`
- `sd_ssdi01_compute_stage_spectral.m`

## Outputs
- `results/ssdi/stage_ssdi_anatomical7/spectral/*.mat`
- `results/ssdi/stage_ssdi_anatomical7/logs/*.csv`
- `results/ssdi/stage_ssdi_anatomical7/reports/*`
- `results/ssdi/stage_ssdi_anatomical7/diagnostics/*`

## Run
From MATLAB:

```matlab
addpath('code/ssdi');
run_ssdi_sleep_stage_pipeline('mode','stage_ssdi');
```

Or through stage runner:

```bash
python3 code/run_all/run_pipeline.py ssdi
```
