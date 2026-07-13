# Source Reconstruction (Anatomical46 -> Anatomical7)

SomnoDynamics source reconstruction uses a two-step branch:

1. Reconstruct fine source time series at HCPMMP1_combined label level (46-region branch used in the original pipeline context).
2. Coarse-grain those fine labels to canonical anatomical7 macro-regions using DMT-derived mapping.

Implementation details (current):
- Common LCMV spatial filter per subject across all available sleep stages.
- Source reconstruction is applied before 2-second segmentation.
- Non-overlapping 2-second source trials are produced at 128 Hz.
- Strict count check: reconstructed trial count must equal `floor(stage_duration_seconds / 2)` from raw TXT annotations.

## Entry point

```bash
python3 code/src_recon/run_source_reconstruction_anatomical46.py
```

Example (subset + overwrite):

```bash
NUMBA_DISABLE_JIT=1 /Users/borjan/miniconda3/envs/mne/bin/python \
  code/src_recon/run_source_reconstruction_anatomical46.py \
  --subject EPCTL01 --subject EPCTL03 --overwrite --batch-size 4
```

## Inputs

- Raw data root: `/Volumes/ex_data/data_somnoemergence/raw`
- Montage file: `/Volumes/ex_data/data_somnoemergence/raw/Co-registered_average_positions.pos`
- Sensor preprocessing input: `results/preprocessing/preprocessed_sensor_space/<subject>/<subject>_<stage>_epochs.fif`
- Raw stage annotations for strict count checks: `/Volumes/ex_data/data_somnoemergence/raw/<subject>/<subject>.txt`

## Outputs

Fine-level outputs:
- `results/source_recon/anatomical46/<subject>/mat/*_source_time_series_46.mat`
- `results/source_recon/anatomical46/<subject>/csv/*_source_time_series_46.csv`

Coarse anatomical7 outputs:
- `results/source_recon/anatomical7/<subject>/mat/*_source_time_series.mat`
- `results/source_recon/anatomical7/<subject>/csv/*_source_time_series.csv`
- `results/source_recon/anatomical7/source_reconstruction_summary.csv`

Mapping source:
- `code/coarse_graining/anatomical7_mapping_from_dmt.csv`

## Core parameters

- `--target-sfreq` default `128.0`
- `--target-epoch-seconds` default `2.0`
- `--batch-size` default `8` (sensor-epoch batch size for source extraction)
- `--no-strict-count-check` (optional; disables strict expected-window validation)
