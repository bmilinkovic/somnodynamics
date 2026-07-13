# MIGRATION_PLAN_sleep_to_dmt_style

## Project Rename

The repository has been migrated to the project name **SomnoDynamics**. New code and documentation use this name. Legacy references to `healthy-sleep-project` are preserved only inside compatibility code under `code/legacy`.

## Phase 1 Audit

### SomnoDynamics audit (pre-refactor state)

Identified legacy components:
- preprocessing scripts: `code/legacy/source/01_preprocessing.py`, `code/legacy/source/load-edf-data.py`
- source reconstruction (legacy branch): `code/legacy/source/source_reconstruction_yeo7.py`
- GC pipeline: `code/legacy/run_mvgc_atlas.m`, `code/legacy/fun_conn_analysis/run_gc_mi_stats_atlas.m`
- SSDI downstream analysis: `code/legacy/ssdi_cross_frequency_tests.py`, `code/legacy/extract_spectral_dopto.py`, `code/legacy/ssdi_analysis/*`
- visualization/statistics notebooks and scripts in `code/legacy/fun_conn_analysis/*`, `code/legacy/test/*`

Observed issues before migration:
- many hard-coded paths to `/Users/borjan/code/healthy-sleep-project/...`
- mixed output locations under legacy `code/results/*`
- non-modular stage boundaries and limited shared configuration
- inconsistent sleep stage display naming in downstream outputs

### DMT reference audit

Reference GC pipeline:
- `/Users/borjan/code/dmt_project/code/paper_style_mvgc_anatomical7`
- key methods inspected:
  - `ps_config_anatomical7.m`
  - `ps01_compute_subject_measures_anatomical7.m`
  - `ps02_group_stats_paired_anatomical7.m`
  - `ps03_make_figures_anatomical7.m`

Reference SSDI pipeline:
- `/Users/borjan/code/dmt_project/code/ssdi_anatomical7`
- key methods inspected:
  - `ssdi_anatomical7_config.m`
  - `run_spectral_domain_optimisation_dmt_anatomical7.m`
  - `run_spectral_domain_optimisation_dmt_anatomical7_canonical.m`
  - `run_ssdi_group_summary_anatomical7.m`

Canonical anatomical mapping source:
- `/Users/borjan/code/dmt_project/code/paper_style_mvgc_anatomical7/export_source_data_anatomical7.py`
- concrete mapping artifact copied from DMT results:
  - `/Users/borjan/code/dmt_project/results/source_reconstruction_anatomical7/coarse_graining/hcp_reduced_to_anatomical7_mapping.csv`

## Phase 2 Restructure (implemented)

Implemented modular structure:
- `code/legacy/` (legacy code preserved)
- `code/preprocessing/`
- `code/sensor_power/`
- `code/src_recon/`
- `code/coarse_graining/`
- `code/gc/`
- `code/ssdi/`
- `code/stats/`
- `code/figures/`
- `code/cross_analysis/`
- `code/run_all/`
- `data/interim`, `data/processed`
- `results/preprocessing`, `results/sensor_power`, `results/source_recon`, `results/coarse_grained`, `results/gc`, `results/ssdi`, `results/stats`, `results/figures`, `results/cross_analysis`

Legacy scripts were moved into `code/legacy`.

## Phase 3 Integration (implemented baseline)

Added integrated stage scripts:
- preprocessing: `code/preprocessing/preprocess_sensor_space.py`
- sensor power analysis: `code/sensor_power/run_sensor_power_analysis.py`
- source reconstruction primary pipeline: `code/src_recon/run_source_reconstruction_anatomical46.py`
  - reconstruct 46-region HCP atlas first
  - then coarse-grain to anatomical7 via DMT mapping
- source reconstruction compatibility wrapper: `code/src_recon/run_source_reconstruction_yeo7.py`
- canonical mapping validation: `code/coarse_graining/validate_mapping.py`
- GC pipeline wrapper: `code/gc/run_gc_sleep_stage_pipeline.m`
- SSDI pipeline wrapper: `code/ssdi/run_ssdi_sleep_stage_pipeline.m`
- SSDI summary/statistics: `code/ssdi/summarize_ssdi_results.py`
- GC repeated-measures stats: `code/stats/repeated_measures_gc.py`
- GC-SSDI cross-analysis: `code/cross_analysis/run_gc_ssdi_cross_analysis.py`
- orchestrator: `code/run_all/run_pipeline.py`

## Phase 4 Statistics (implemented baseline)

Repeated-measures statistics across stages are implemented for:
- sensor power (`Friedman`, pairwise `Wilcoxon`, BH-FDR)
- GC node-level metrics (`Friedman`, pairwise `Wilcoxon`, BH-FDR)
- SSDI summary metrics (`Friedman`, pairwise `Wilcoxon`, BH-FDR)

## Phase 5 Figures (implemented baseline)

Figure outputs in PNG/SVG include:
- sensor power topomaps by stage
- GC stage-node heatmaps
- SSDI curves across macro sizes
- GC-SSDI scatter diagnostics

## Phase 6 Validation (implemented baseline)

Validation entrypoints:
- per-stage execution via `code/run_all/run_pipeline.py <stage>`
- full pipeline via `code/run_all/run_pipeline.py all`
- figure consolidation via `code/figures/build_publication_figures.py`
- artifact validation via `code/run_all/validate_pipeline_outputs.py`

Additional validation docs:
- `docs/PREPROCESSING.md`
- `docs/CANONICAL_REGION_MAPPING.md`
- `docs/IMPLEMENTATION_SUMMARY.md`

## Remaining Work for Full DMT-Level Equivalence

Items still requiring dataset-scale runs and numeric parity checks:
- full numeric parity benchmarking against legacy outputs
- full SSDI optimization runtime validation on all subjects/stages
- visual style fine-tuning to match DMT publication figure aesthetics exactly
- final manuscript-grade report templates with frozen versions
