# SomnoDynamics Stage MVGC2 (Anatomical7)

This module ports the DMT `paper_style_mvgc_anatomical7` methodology to SomnoDynamics sleep stages.

## Design
- Input: anatomical7 source-space time series from SomnoDynamics source reconstruction
- Comparison design: within-participant across sleep stages (`W_`, `N1`, `N2`, `N3`, `R_`)
- Not used: rest-vs-dmt and time-windowed-within-state contrasts

## Pipeline files
- `sd_gc_config_anatomical7.m`
- `run_stage_mvgc2_anatomical7_pipeline.m`
- `sd_gc01_compute_subject_measures_by_stage.m`
- `sd_gc02_group_stats_repeated_measures.m`

## Metrics
Subject-level (MVGC2):
- directed GC (time + spectral): source in/out and ROI matrices
- undirected MI-family (time + spectral): source and ROI matrices

Group-level:
- Friedman omnibus tests across stages (per feature)
- pairwise signed-rank between stage pairs
- BH primary FDR + BY sensitivity FDR within correction families

## Frequency bands
Current non-overlapping band definitions:
- delta: 1-4 Hz
- theta: 5-8 Hz
- alpha: 9-12 Hz
- sigma: 13-16 Hz
- beta: 17-30 Hz
- gamma: 31-44 Hz

## Outputs
`results/gc/stage_mvgc2_anatomical7/`
- `data_intermediate/subject_measures/<stage>/...`
- `stats/repeated_measures_stats.mat`
- `stats/omnibus_long.csv`
- `stats/pairwise_long.csv`
- `reports/repeated_measures_summary.txt`
- `logs/*.log`

## Run
From MATLAB:

```matlab
addpath('code/gc');
run_gc_sleep_stage_pipeline('mode','stage_mvgc2');
```

Through stage runner:

```bash
python3 code/run_all/run_pipeline.py gc
```
