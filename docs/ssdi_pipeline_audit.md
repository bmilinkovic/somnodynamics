# SSDI Pipeline Audit: SomnoDynamics vs DMT (Canonical Requirement)

Date: 2026-03-17

## Scope
This audit compares:
- DMT SSDI reference (`/Users/borjan/code/dmt_project/code/ssdi_anatomical7`)
- Prior SomnoDynamics fast SSDI implementation (now preserved as legacy)
- Current SomnoDynamics corrected SSDI implementation

Primary question: why the previous run was too fast, and whether current SSDI now follows canonical spectral DD computation.

## Files Audited

### DMT reference (non-canonical + canonical)
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/run_spectral_domain_optimisation_dmt_anatomical7.m`
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/run_spectral_domain_optimisation_dmt_anatomical7_canonical.m`
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/canonical_opt_ddx_mruns.m`
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/canonical_opt_dds_mruns.m`
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/canonical_band_objective_and_grad.m`
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/ssdi_anatomical7_config.m`

### SomnoDynamics previous fast implementation (preserved)
- `code/ssdi/stage_ssdi_anatomical7_legacy_fast/sd_ssdi_config_anatomical7.m`
- `code/ssdi/stage_ssdi_anatomical7_legacy_fast/sd_ssdi01_compute_stage_spectral.m`

### SomnoDynamics current corrected canonical implementation
- `code/ssdi/stage_ssdi_anatomical7/sd_ssdi_config_anatomical7.m`
- `code/ssdi/stage_ssdi_anatomical7/sd_ssdi01_compute_stage_spectral.m`
- `code/ssdi/stage_ssdi_anatomical7/canonical_opt_ddx_mruns.m`
- `code/ssdi/stage_ssdi_anatomical7/canonical_opt_dds_mruns.m`
- `code/ssdi/stage_ssdi_anatomical7/canonical_band_objective_and_grad.m`
- `code/ssdi/stage_ssdi_anatomical7/run_stage_ssdi_anatomical7_pipeline.m`
- `code/ssdi/run_ssdi_sleep_stage_pipeline.m`

## Why the Prior SomnoDynamics SSDI Was Too Fast
The prior implementation was materially simplified versus DMT.

### Evidence
1. Reduced optimizer scale:
- `nrunsp=30`, `nitersp=3000` (no second-stage `niterso` in config)
  in legacy config.

2. Single-stage optimization path:
- One `opt_gd_dds_mruns` call per band/m and direct save, with no DMT-style preopt clustering + second stage.

3. Reduced outputs:
- Missing DMT families such as `*_preopt_m*`, `*_Loptp_m*`, `*_diag_m*`.

### Audit conclusion
The previous implementation could not be considered faithful to DMT SSDI optimization depth; fast runtime was expected from this simplification.

## Canonical Requirement: Current State
Per your instruction, pipeline computation is now canonical-only.

### What changed
Current stage SSDI now uses canonical optimization functions:
- Pre-stage: `canonical_opt_ddx_mruns` (proxy DD on `CAK`)
- Final-stage: `canonical_opt_dds_mruns` (Eq. 25/26 canonical band objective)

in `sd_ssdi01_compute_stage_spectral.m`.

### Canonical objective implementation
- `canonical_band_objective_and_grad.m` computes per-frequency canonical spectral DD and analytic gradient.
- `canonical_opt_dds_mruns.m` runs multi-restart gradient-descent optimization on the canonical band objective.

### Output parity maintained
Current outputs keep DMT-style result organization and key families (`*_preopt_m*`, `*_Loptp_m*`, `*_dopto_m*`, `*_diag_m*`, node/beta outputs), with sleep-stage adaptation (`W_`, `N1`, `N2`, `N3`, `R_`).

## Sleep-Specific Adaptations (explicit and intended)
- Condition axis replaced with stage axis.
- Bands are non-overlapping sleep bands requested for SomnoDynamics:
  - delta 1-4, theta 5-8, alpha 9-12, sigma 13-16, beta 17-30, gamma 31-44.
- Output path is stage-structured:
  - `results/ssdi/stage_ssdi_anatomical7/<stage>/<subject>_ssdi_spectral_domain_results.mat`

## Legacy Preservation
Fast/simplified prior SSDI was preserved and clearly marked non-authoritative:
- `code/ssdi/stage_ssdi_anatomical7_legacy_fast`
- `code/ssdi/stage_ssdi_anatomical7_legacy_fast/README_LEGACY_FAST.md`

## Validation Snapshot (Canonical Run)
Representative canonical run executed:
- Command:
  - `run_ssdi_sleep_stage_pipeline('mode','stage_ssdi','recompute',true,'max_files',1)`
- Case:
  - `EPCTL01 | N1`
- Runtime:
  - ~112.1 s
- Status table captured:
  - `results/ssdi/stage_ssdi_anatomical7/logs/stage1_ssdi_status.csv`

This runtime increase versus the earlier ~14 s run is consistent with switching from non-canonical simplified flow to canonical optimization.

Validation summary:
- `docs/ssdi_validation_summary.md`
