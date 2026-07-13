# Source Reconstruction + Preprocessing Parity Audit (Sleep vs DMT)

Date: 2026-03-18

## Scope
- Verify preprocessing steps match legacy sleep pipeline.
- Verify montage handling matches legacy sleep preprocessing.
- Verify anatomical7 coarse-graining matches DMT anatomical7 methodology.

## Legacy preprocessing parity
Authoritative legacy code:
- `code/legacy/source/01_preprocessing.py`

Matched operations in maintained preprocessing:
- Notch 50/100 Hz (`spectrum_fit`), bandpass 1-45 Hz, resample to 256 Hz.
- Channel rename (`-Ref` removed, uppercase), montage set with `on_missing='ignore'`.
- Exclude non-EEG + missing-location channels.
- Average reference with projection + `apply_proj()`.
- 10-second epoch tiling from annotation durations.
- Legacy-compatible outputs (`W`/`R` naming + flat file layout + cleaned CSV).

## Montage handling parity
Legacy and current both:
- read electrode names from montage table,
- uppercase channel labels,
- swap X/Y columns,
- create `make_dig_montage(..., coord_frame='head')`,
- apply after channel renaming.

## DMT anatomical7 coarse-graining parity
DMT source of truth:
- `/Users/borjan/code/dmt_project/code/paper_style_mvgc_anatomical7/export_source_data_anatomical7.py`

Checked parity points:
- fsaverage + ico4 source space.
- LCMV NAI settings (`reg=0.05`, `pick_ori='max-power'`, `weight_norm='nai'`, ad-hoc noise covariance).
- HCPMMP1_combined label extraction.
- Unknown `???` label removal before downstream label time-series use.
- Explicit reduced-HCP -> anatomical7 mapping with macro order:
  Frontal, Sensorimotor, Parietal, Temporal, Occipital, Cingulate, Limbic.

Mapping parity check:
- Somno mapping: `code/coarse_graining/anatomical7_mapping_from_dmt.csv`
- DMT mapping: `.../results/source_reconstruction_anatomical7/coarse_graining/hcp_reduced_to_anatomical7_mapping.csv`
- Exact match on `(label, macro_region)`: TRUE (44 rows each).

## Fix applied during this audit
- Updated `code/src_recon/run_source_reconstruction_anatomical46.py` to remove unknown `???` labels before extraction/coarse-graining and fail if any surviving label lacks mapping.

## Validation run (1 subject, all stages)
Input:
- `/Volumes/dataSets/healthy-sleep-project-datasets/data/preprocessed_sensor_space`

Result summary:
- `results/source_recon/anatomical7_legacy_input_256_test_dmtchecked/source_reconstruction_summary.csv`
- `n_fine_labels=44`, `n_macro_labels=7`, `source_sfreq_hz=256.0`, counts match expected stage windows.

Example output:
- `results/source_recon/anatomical7_legacy_input_256_test_dmtchecked/EPCTL01/mat/EPCTL01_W__source_time_series.mat`
- shape: `(7, 512, 690)`
