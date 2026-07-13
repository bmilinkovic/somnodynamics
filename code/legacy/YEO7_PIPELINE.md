# Yeo7 + Legacy Dual Pipeline

This project now supports two side-by-side pipelines:

1. Legacy Desikan pipeline (existing, unchanged):
- Source/inputs as currently used by `healthy_sleep_mvgc.m`
- Outputs in existing legacy locations (for example `code/results/pwcgc*`, `code/results/new_results_GC`, `code/results/ssdi_results`)

2. New Yeo7 pipeline (parallel):
- Source reconstruction output: `code/results/source_recon/yeo7/`
- GC/MI output: `code/results/yeo7/gc/`
- SSDI output (expected): `code/results/yeo7/ssdi_results/`

## New files

- `code/source/source_reconstruction_yeo7.py`
- `code/source/validate_yeo7_reconstruction.py`
- `code/run_mvgc_atlas.m`
- `code/fun_conn_analysis/run_gc_mi_stats_atlas.m`
- `code/run_yeo7_pipeline.m`
- `code/run_dual_atlas_pipelines.m`
- `code/ssdi_analysis/run_ssdi_crossfreq_for_atlas.py`
- `code/ssdi_analysis/run_extract_spectral_dopto_for_atlas.py`

## Yeo7 source reconstruction

Run:

```bash
# Recommended: use a fully local copy (not OneDrive CloudStorage placeholders).
# Example local layout:
#   /Users/borjan/data/APHNY_data/raw/EPCTLxx/EPCTLxx.edf
#   /Users/borjan/data/APHNY_data/Co-registered_average_positions.pos
export APHNY_DATA_DIR="/Users/borjan/data/APHNY_data"
export APHNY_MONTAGE_FILE="/Users/borjan/data/APHNY_data/Co-registered_average_positions.pos"

python3 /Users/borjan/code/healthy-sleep-project/code/source/source_reconstruction_yeo7.py \
  --raw-dir "$APHNY_DATA_DIR/raw" \
  --data-dir "$APHNY_DATA_DIR" \
  --montage-file "$APHNY_MONTAGE_FILE"
```

This writes per-subject/per-stage files like:

- `code/results/source_recon/yeo7/EPCTL01/mat/EPCTL01_N1_source_time_series.mat`
- `code/results/source_recon/yeo7/EPCTL01/csv/EPCTL01_N1_source_time_series.csv`

Each `.mat` includes:
- `source_ts` (shape: regions x time x epochs)
- `labels` (deterministic 7-label order)
- `atlas = 'yeo7'`
- `metadata` (sampling rate, epoch length, inverse settings)

## Yeo7 GC + MI downstream analysis

In MATLAB:

```matlab
startup
run_yeo7_pipeline
```

Or step-by-step:

```matlab
startup
run_mvgc_atlas('atlas','yeo7')
run_gc_mi_stats_atlas('atlas','yeo7')
```

Outputs include:
- `code/results/yeo7/gc/pwcgc/pwcgc_params_EPCTLxx_<stage>.mat`
- `code/results/yeo7/gc/gc_in_freq_mw_effect_sizes.mat`
- `code/results/yeo7/gc/gc_out_freq_mw_effect_sizes.mat`
- `code/results/yeo7/gc/gc_in_time_mw_effect_sizes.mat`
- `code/results/yeo7/gc/gc_out_time_mw_effect_sizes.mat`
- `code/results/yeo7/gc/mi/mi_data.mat`
- `code/results/yeo7/gc/mi/freq_cmi_data.mat`
- z-score matrix exports for GC and MI (atlas-specific directory)

## Yeo7 validation and parity checks

Sanity checks:

```bash
python3 /Users/borjan/code/healthy-sleep-project/code/source/validate_yeo7_reconstruction.py
```

Cross-project parity check (sleep vs dmt):

```bash
python3 /Users/borjan/code/healthy-sleep-project/code/source/validate_yeo7_reconstruction.py \
  --sleep-file /path/to/sleep_subject_stage.mat \
  --dmt-file /Users/borjan/code/dmt_project/results/source_reconstruction_7regions/rest/data/SUBJ_7regions.mat
```

Printed metrics:
- label set/count consistency (expected 7)
- ROI-correlation Frobenius difference
- upper-triangle correlation between ROI correlation matrices
- pass/warn summary

## Yeo7 SSDI wrappers

Cross-frequency SSDI figures:

```bash
python3 /Users/borjan/code/healthy-sleep-project/code/ssdi_analysis/run_ssdi_crossfreq_for_atlas.py --atlas yeo7
```

Extract spectral dopto vectors:

```bash
python3 /Users/borjan/code/healthy-sleep-project/code/ssdi_analysis/run_extract_spectral_dopto_for_atlas.py --atlas yeo7
```

## Run both pipelines side-by-side

In MATLAB:

```matlab
startup
run_dual_atlas_pipelines('run_legacy',true,'run_yeo7',true)
```

To run only one:

```matlab
run_dual_atlas_pipelines('run_legacy',true,'run_yeo7',false)  % legacy only
run_dual_atlas_pipelines('run_legacy',false,'run_yeo7',true)  % yeo7 only
```

Legacy scripts are not replaced. Yeo7 outputs are isolated under `code/results/yeo7` and `code/results/source_recon/yeo7`.

## Legacy preprocessing/source scripts now default to APHNY path

- `code/source/01_preprocessing.py`
- `code/source/load-edf-data.py`

They now read:
- raw EDF/TXT from `${APHNY_DATA_DIR}/raw`
- montage from `${APHNY_DATA_DIR}/Co-registered_average_positions.pos`

Defaults:
- `APHNY_DATA_DIR=/Users/borjan/Library/CloudStorage/OneDrive-ICM/APHNY data`
- `APHNY_MONTAGE_FILE` can be set to force a local montage file path

Outputs default to project results:
- `code/results/preprocessed_sensor_space`
- `code/results/preprocessed_data`
