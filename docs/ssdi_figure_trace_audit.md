# SSDI Figure Trace Audit (DMT -> SomnoDynamics Stage Pipeline)

Date: 2026-03-19

## Scope
This audit traces the exact DMT scripts that generate the requested SSDI reference figures, maps their SomnoDynamics equivalents, and records the patches made so SomnoDynamics can generate sleep-stage analogues with DMT-parity logic.

## Phase 1: DMT Figure Provenance

### 1) `sf_index_pooled_main.pdf`
- Figure path:
  - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_canonical/figures/similarity_sf_lopto/sf_index_pooled_main.pdf`
- Generating script:
  - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/ssdi_anatomical7_canonical_similarity_sf_plots.py`
- Key save lines:
  - around lines 876-877 (`sf_index_pooled_main.png/.pdf`)
- Upstream inputs:
  - canonical SSDI subject outputs under:
    - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_canonical/{rest,dmt}/*_ssdi_spectral_domain_canonical.mat`
  - optional subject-measure covariance fallback source:
    - `/Users/borjan/code/dmt_project/results/source_reconstruction_anatomical7/gc_results/data_intermediate/subject_measures`
- Selection/ordering logic:
  - `BANDS = [delta, theta, alpha, beta, gamma]`
  - `MVALS = [2,3,4,5,6]`
  - condition order rendered as `dmt` then `rest` in matrix supplement; pooled line figure shows both conditions.
  - SF computed as proportion of similarity distances below threshold (`SSDI_SF_THRESHOLD`, default `1e-4`).
- Style settings:
  - palette: rest `#2EC4B6`, dmt `#5E60CE`
  - high-contrast typography, 90/95% CI shading, BH-FDR stars.
  - right-side “Structured/Fragmented” intuition arrow.

### 2) `canonical_dd_min_grand_average.png` and `canonical_dd_min_violin_by_macro.png`
- Figure paths:
  - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_canonical/figures/canonical_dd_min_grand_average.png`
  - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_canonical/figures/canonical_dd_min_violin_by_macro.png`
- Generating script:
  - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/ssdi_anatomical7_canonical_stats_plots.py`
- Key save lines:
  - around lines 286-287 (violin-by-macro)
  - around lines 434-435 (grand-average)
- Upstream inputs:
  - canonical per-subject outputs in:
    - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_canonical/{rest,dmt}`
- Selection/ordering logic:
  - same `BANDS` and `MVALS` as above
  - `dd_min = min(canonical_band_dd_m*)` per subject/condition/band/m
  - paired Wilcoxon; BH-FDR within macro across bands.
- Style settings:
  - paired violin + strip overlays
  - rest/dmt palette as above
  - significance stars in panel top rows
  - grand figure includes “Scale separated/integrated” arrow.

### 3) Timewindow canonical composite/overview figures
- Figure paths:
  - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_timewindows_canonical/figures/SF_index_4x5_timewindow_band_composite.png`
  - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_timewindows_canonical/figures/grand_average_DD_composite_timewindows.png`
  - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_timewindows_canonical/downstream_timewindow_pairwise/figures/timewindow_pairwise_overview_left_column.png`
  - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_timewindows_canonical/downstream_timewindow_pairwise/figures/timewindow_pairwise_overview_right_column.png`
- Generating script:
  - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7_timewindows/ssdi_anatomical7_timewindows_canonical_stats_plots.py`
- Key save lines:
  - around lines 905-906 (`grand_average_DD_composite_timewindows`)
  - around line 1110 (`SF_index_4x5_timewindow_band_composite`)
  - around lines 1367 and 1373 (left/right overview columns)
- Upstream inputs:
  - subject-level canonical outputs:
    - `/Users/borjan/code/dmt_project/results/ssdi_anatomical7_timewindows_canonical/subject_level/...`
  - rest-baseline canonical outputs + per-window DMT outputs.
- Selection/ordering logic:
  - states: `Rest, Onset, Peak, Plateau, Recovery`
  - pairwise families are predefined state comparisons.
  - DD/SF pairwise grand tests: paired Wilcoxon, BH within comparison across bands.
- Style settings:
  - consistent rest/window palette, band-wise panel composites, starred significance overlays.
  - pairwise heatmap matrices with printed deltas and stars.
  - final left/right overview column composition from intermediate panels.

## Phase 2: SomnoDynamics Parity Audit

## Existing SomnoDynamics SSDI scripts inspected
- SSDI run entrypoint:
  - `/Users/borjan/code/somnodynamics/code/ssdi/run_ssdi_sleep_stage_pipeline.m`
- Canonical stage SSDI runner:
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/run_stage_ssdi_anatomical7_pipeline.m`
- Existing summary and plotting:
  - `/Users/borjan/code/somnodynamics/code/ssdi/summarize_ssdi_results.py`
  - `/Users/borjan/code/somnodynamics/code/ssdi/plot_ssdi_stage_figures.py` (pre-patch)

## Parity gaps found (pre-patch)
- SomnoDynamics had basic stage DD summaries but did not implement DMT-style:
  - pooled SF main figure from similarity matrices
  - stage-composite analogue of timewindow SF composite
  - stage-composite DD panel analogous to timewindow DD composite
  - downstream pairwise left/right overview columns with trajectory+heatmap composition
- Existing plotting path did not read canonical similarity matrices from SSDI `.mat` outputs to compute SF-index tables.

## Phase 3/4: Patches Implemented in SomnoDynamics

## Files changed
- `/Users/borjan/code/somnodynamics/code/ssdi/plot_ssdi_stage_figures.py`
- `/Users/borjan/code/somnodynamics/code/ssdi/summarize_ssdi_results.py`

## What changed
- Replaced the plotting script with DMT-parity stage logic:
  - canonical DD min violin-by-macro output
  - canonical DD grand-average output
  - pooled SF index main output computed from saved similarity matrices in SSDI `.mat` files
  - stage-composite SF figure (sleep-stage analogue of timewindow composite)
  - stage-composite DD figure (sleep-stage analogue of timewindow DD composite)
  - stage pairwise trajectories + band heatmaps + composed left/right overview panels
- Added stage-pairwise stats/table outputs for DD and SF with BH-FDR correction families mirroring DMT pairwise style.
- Updated summarizer invocation so figure generation receives `--spectral-dir` and can compute SF from actual SSDI outputs.

## Adaptation note (DMT vs sleep stages)
- DMT timewindow states (`Rest, Onset, Peak, Plateau, Recovery`) are replaced in SomnoDynamics by sleep-stage states:
  - `Wake, N1, N2, N3, REM`
- Wake serves the baseline role analogous to Rest for wake-vs-stage contrasts and star annotations.
- Sleep-specific bands include sigma (`delta, theta, alpha, sigma, beta, gamma`), so stage composites are rendered across six bands.

## Phase 5: Validation checklist (to be finalized after rerun)
- DMT provenance traced: complete.
- SomnoDynamics script parity patch: complete.
- Full SSDI rerun for `momax4 + 400epochs`: pending execution/verification in this run cycle.
- Final output paths and rerun confirmation: reported after job completion.
