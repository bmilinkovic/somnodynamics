# Standalone Audit And Amendments (DK5/DK68 + LRT)

## Epoching decision (locked)
- Kept current established epoching unchanged:
  - sensor preprocessing epochs: 10 s (existing)
  - source windows: 2 s non-overlapping (existing)
- No 3 s epoching introduced in this standalone branch.

## What was already correct
1. Beamformer family and core settings were already aligned:
   - LCMV with NAI, `pick_ori='max-power'`, `reg=0.05`
   - fsaverage + template BEM workflow
2. Pooled normalisation for multi-trial VAR inputs was already implemented and retained.
3. Stage naming and repeated-measures downstream logic were already in place.

## What needed amending
1. Atlas mismatch:
   - active pipeline was anatomical46->anatomical7, not DK5/DK68.
   - Added standalone source reconstruction producing DK5 and DK68 outputs.
2. Model-order selection mismatch:
   - active GC branch selected AIC.
   - Standalone branch now computes AIC/BIC/HQC/LRT and selects LRT by default, with HQC stored as a check.
3. Sampling defaults mismatch:
   - several branches defaulted to 128 Hz.
   - Standalone DK branches set to 256 Hz.
4. Band definition mismatch:
   - removed sigma from this branch to match requested delta/theta/alpha/beta/gamma set.

## Files added/modified in this standalone package
- `src_recon/run_source_reconstruction_dk5_dk68.py`
- `gc/run_gc_sleep_stage_pipeline_dk5.m`
- `gc/plot_gc_stage_figures_dk5.py`
- `gc/stage_mvgc2_dk5/*`
- `gc/stage_mvgc2_dk68/*`
- `ssdi/run_ssdi_sleep_stage_pipeline_dk.m`
- `ssdi/summarize_ssdi_results_dk.py`
- `ssdi/plot_ssdi_stage_figures_dk.py`
- `ssdi/stage_ssdi_dk5/*`
- `ssdi/stage_ssdi_dk68/*`

## Output roots used by this standalone branch
- Source recon:
  - `results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/`
- GC DK5:
  - `results/gc/stage_mvgc2_dk5_lrt30_current_epoching_2s_256hz/`
- GC DK68:
  - `results/gc/stage_mvgc2_dk68_lrt30_current_epoching_2s_256hz/`
- SSDI DK5:
  - `results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/`
- SSDI DK68:
  - `results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz/`
