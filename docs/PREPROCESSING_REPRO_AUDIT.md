# Preprocessing Reproduction Audit (Sleep Project)

## 1) Authoritative legacy preprocessing code
- File: `code/legacy/source/01_preprocessing.py`
- Historical output location inspected: `/Volumes/dataSets/healthy-sleep-project-datasets/data/preprocessed_sensor_space`

## 2) What the legacy code does
- Reads each subject EDF + TXT annotation file from raw data.
- Applies notch filtering at 50/100 Hz (`method='spectrum_fit'`).
- Applies bandpass filtering 1-45 Hz (`firwin`) and resamples to 256 Hz.
- Renames channels (`-Ref` removed, uppercase), applies montage, drops non-EEG and bad-location channels.
- Sets average EEG reference and applies projection.
- Writes `<subject>_cleaned.csv`.
- Builds 10-second events from annotation durations and saves stage FIF files as:
  - `<subject>_W_epochs.fif`, `<subject>_N1_epochs.fif`, `<subject>_N2_epochs.fif`, `<subject>_N3_epochs.fif`, `<subject>_R_epochs.fif`
- Output is a **flat directory** (not per-subject subfolders).

## 3) How current code differed before this fix
- Current maintained runner wrote outputs under per-subject subdirectories.
- It canonicalized filenames to `W_` and `R_` instead of legacy `W` and `R`.
- It did not write legacy `<subject>_cleaned.csv` output.
- These output-contract differences prevented faithful reproduction of existing legacy artifacts.

## 4) Changes made to restore faithful behavior
- Updated `code/preprocessing/preprocess_sensor_space.py` to support legacy-compatible output contract (now default):
  - Flat output layout.
  - Legacy stage filename labels (`W`, `N1`, `N2`, `N3`, `R`).
  - Legacy cleaned CSV generation (`<subject>_cleaned.csv`).
- Cleaned CSV writing is implemented in chunked streaming mode to reduce RAM spikes while preserving legacy artifact type.
- Added explicit output-label/output-path helpers so preprocessing and import modes both follow the same legacy-compatible naming rules.

## 5) Notes on intended logic preservation
- Core preprocessing logic (filter/resample/channel handling/referencing/stage epoching) was preserved.
- Fixes target output-contract mismatches so current runs reproduce the previously used preprocessing behavior rather than introducing a new workflow.
