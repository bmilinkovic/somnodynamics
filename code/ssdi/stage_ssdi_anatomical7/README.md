# SomnoDynamics Stage SSDI (Anatomical7)

This module is a DMT-faithful canonical SSDI port (`Eq. 25/26`) adapted to SomnoDynamics sleep stages.

## Design
- Input models: `A`, `V` from stage-wise GC subject measures
- Stages: `W_`, `N1`, `N2`, `N3`, `R_`
- Source space: anatomical7
- Optimizer (canonical, two-stage):
  - pre-stage proxy DD: `canonical_opt_ddx_mruns` + clustering
  - final canonical spectral DD: `canonical_opt_dds_mruns` (Eq. 25/26)

## Key parity with DMT SSDI
- Same optimisation structure (`preopt -> Lcluster -> second stage`) in canonical form.
- Same output family (`*_preopt_m*`, `*_Loptp_m*`, `*_dopto_m*`, `*_diag_m*`, node/beta tables).
- Same DD diagnostics (canonical band objective, full-spectrum DD check, negativity debug summaries).
- Same default optimiser scales (`nrunsp=100`, `nitersp=10000`, `niterso=10000`).

Sleep-specific adaptation:
- Conditions become stages (`W_`, `N1`, `N2`, `N3`, `R_`).
- Non-overlapping sleep bands: delta/theta/alpha/sigma/beta/gamma.

## Files
- `sd_ssdi_config_anatomical7.m`
- `run_stage_ssdi_anatomical7_pipeline.m`
- `sd_ssdi01_compute_stage_spectral.m`

## Outputs
- `results/ssdi/stage_ssdi_anatomical7/<stage>/*_ssdi_spectral_domain_results.mat`
- `results/ssdi/stage_ssdi_anatomical7/logs/*.csv`
- `results/ssdi/stage_ssdi_anatomical7/summaries/*`
- `results/ssdi/stage_ssdi_anatomical7/figures/*`
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

Legacy fast/simplified implementation is preserved under:
- `code/ssdi/stage_ssdi_anatomical7_legacy_fast`
