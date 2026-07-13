# SomnoDynamics SSDI Pipeline Audit and Port (DMT -> Sleep Stages)

Date: 2026-03-17

## Scope
This report audits the DMT SSDI implementation and documents the SomnoDynamics port for stage-based sleep analysis.

Reference implementation audited:
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7`

SomnoDynamics target:
- within-participant stage analysis across `W_`, `N1`, `N2`, `N3`, `R_`
- 7-region anatomical coarse-grained source models
- SSDI optimisation from subject-level VAR/state-space models

## DMT files audited
Primary DMT SSDI files inspected:
- `ssdi_anatomical7_config.m`
- `run_spectral_domain_optimisation_dmt_anatomical7.m`
- `run_spectral_domain_optimisation_dmt_anatomical7_canonical.m`
- `run_ssdi_group_summary_anatomical7.m`
- `run_ssdi_figures_anatomical7.m`
- `run_ssdi_canonical_stats_and_figures.m`
- `README.md`

## Upstream dependencies identified
DMT SSDI depends on:
- SSDI toolbox (`ssdi-1`) functions, primarily `opt_gd_dds_mruns`
- MVGC2/state-space utilities (`transform_var`, `var_to_ss`, `ss2trfun`, `sfreqs`)
- Subject-level models (`A`, `V`) from prior GC-style subject-measure stage

## DMT execution flow (audited)
1. Load subject-condition model (`A`, `V`) from subject measures.
2. Convert VAR to state-space/transfer function.
3. Run band-limited SSDI optimisation for each macro size `m` and band.
4. Save per subject-condition optimisation outputs (`dopto`, `Lopto`, node weights, beta metrics).
5. Aggregate paired Rest vs DMT summaries and inferential stats.
6. Generate figures and reports.

## Key assumptions in original DMT SSDI code
- Two-condition paired design (`rest`, `dmt`).
- Condition-specific directories and naming.
- Group stats are paired condition comparisons, not multi-stage repeated measures.

## Required adaptation for SomnoDynamics
SomnoDynamics differs scientifically:
- comparison is across sleep stages within participant (5 levels) not rest-vs-dmt.

Adaptations implemented:
- Replace condition logic with stage logic (`W_`, `N1`, `N2`, `N3`, `R_`).
- Use SomnoDynamics GC subject-measure outputs as model inputs:
  - `results/gc/stage_mvgc2_anatomical7/data_intermediate/subject_measures/<stage>/<subject>_subject_measures.mat`
- Preserve DMT-style optimisation core (`opt_gd_dds_mruns` over spectral bands).
- Keep downstream repeated-measures inference in Python summarizer (Friedman + pairwise Wilcoxon BH-FDR).

## SomnoDynamics implementation added
New module:
- `code/ssdi/stage_ssdi_anatomical7/`

Files:
- `sd_ssdi_config_anatomical7.m`
- `run_stage_ssdi_anatomical7_pipeline.m`
- `sd_ssdi01_compute_stage_spectral.m`

Entrypoint wiring:
- `code/ssdi/run_ssdi_sleep_stage_pipeline.m`
  - default mode: `stage_ssdi`
  - optional `mode='legacy'` to run previous pwcgc-based wrapper
- legacy wrapper preserved in:
  - `code/ssdi/run_ssdi_sleep_stage_pipeline_legacy_from_pwcgc.m`

## SSDI toolbox local vendoring
SSDI toolbox is now package-local:
- `code/external/ssdi-1`

Toolbox setup updated:
- `code/run_all/setup_matlab_toolboxes.m`
  - now checks `code/external/ssdi-1` before external machine paths.

## Frequency bands
SomnoDynamics stage SSDI uses non-overlapping bands:
- delta: 1-4 Hz
- theta: 5-8 Hz
- alpha: 9-12 Hz
- sigma: 13-16 Hz
- beta: 17-30 Hz
- gamma: 31-44 Hz

## Input/output contract (SomnoDynamics SSDI)
Input:
- subject-measure MAT files containing `subject_measure.model.A`, `subject_measure.model.V`

Primary outputs:
- `results/ssdi/stage_ssdi_anatomical7/spectral/EPCTLxx_<stage>_ssdi_spectral_domain_results.mat`
- `results/ssdi/stage_ssdi_anatomical7/logs/stage1_ssdi_status.csv`
- diagnostics and reports under `results/ssdi/stage_ssdi_anatomical7/`

Downstream summary outputs (Python):
- `results/ssdi/ssdi_subject_long.csv`
- `results/ssdi/ssdi_dd_*_friedman.csv`
- `results/ssdi/ssdi_dd_*_pairwise_wilcoxon.csv`
- `results/ssdi/figures/*`

## Notes / risks
- SSDI optimisation verbosity can be high from underlying toolbox functions.
- Runtime scales strongly with number of files, restarts, and iterations.
- Current implementation ports DMT optimisation core; canonical Eq. (25/26) audit branch from DMT is not yet ported into SomnoDynamics.

## Recommendation
Use the new stage-based SSDI module as the default SomnoDynamics SSDI path. It preserves DMT optimisation behavior while correctly adapting from two-condition designs to sleep-stage within-participant organization.
