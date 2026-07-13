# Dataset Structure and Contents Summary

This document provides a comprehensive overview of the datasets in the healthy-sleep-project, their structure, shapes, and contents.

## Overview

The project contains multiple types of data organized in different directories:
- **Raw data** (expected location: `/data/raw/`)
- **Preprocessed data** (expected location: `/data/preprocessed_data/`)
- **Results data** (location: `/code/results/`)
- **SSDI analysis data** (location: `/code/ssdi_analysis/data/`)

---

## 1. Raw Data Structure

**Expected Location:** `/Users/borjan/code/healthy-sleep-project/data/raw/`

### Structure:
```
data/raw/
├── EPCTL01/
│   ├── EPCTL01.edf          # Raw EEG data file
│   └── EPCTL01.txt          # Sleep stage annotations
├── EPCTL02/
│   ├── EPCTL02.edf
│   └── EPCTL02.txt
└── ... (EPCTL01 through EPCTL29)
```

### File Contents:

#### `.edf` files (EDF+ format)
- **Format:** European Data Format (EDF+)
- **Content:** Raw EEG recordings
- **Channels:** Multiple EEG channels, EMG, EOG, ECG, leg channels
- **Sampling rate:** Original (likely 256 Hz or higher, resampled to 256 Hz in preprocessing)
- **Processing:** Loaded using MNE Python library

#### `.txt` files (Sleep stage annotations)
- **Format:** Tab-separated values (TSV)
- **Columns:** 
  - `label`: Sleep stage label (W, N1, N2, N3, R, or L)
  - `start_time`: Start time in seconds
  - `duration`: Duration in seconds
- **Epoch length:** 30 seconds per epoch
- **Sleep stages:**
  - `W`: Wake (mapped to code 0 or 1)
  - `N1`: Stage 1 sleep (mapped to code 1 or 2)
  - `N2`: Stage 2 sleep (mapped to code 2 or 3)
  - `N3`: Stage 3/Deep sleep (mapped to code 3 or 4)
  - `R`: REM sleep (mapped to code 4 or 5)
  - `L`: Unknown/artifact (sometimes excluded)

### Additional Raw Data Files:
- **Co-registration file:** `Co-registered_average_positions.pos` (electrode positions)
- **Artifact files:** `{subject_id}_artndxn.mat` (artifact annotations per epoch)

---

## 2. Preprocessed Data Structure

**Expected Location:** `/Users/borjan/code/healthy-sleep-project/data/preprocessed_data/`

### Structure:
```
data/preprocessed_data/
├── EPCTL01/
│   ├── csv/
│   │   ├── EPCTL01_W_source_time_series.csv
│   │   ├── EPCTL01_N1_source_time_series.csv
│   │   ├── EPCTL01_N2_source_time_series.csv
│   │   ├── EPCTL01_N3_source_time_series.csv
│   │   └── EPCTL01_R_source_time_series.csv
│   └── mat/
│       ├── EPCTL01_W_source_time_series.mat
│       ├── EPCTL01_N1_source_time_series.mat
│       └── ... (one per sleep stage)
└── ... (EPCTL01 through EPCTL29)
```

### File Contents:

#### `.mat` files (MATLAB format)
- **Key variable:** `source_ts`
  - **Shape:** `(n_regions, n_timesteps, n_epochs)`
  - **Typical dimensions:** ~68 regions × 513 time points × variable epochs
  - **Content:** Source-reconstructed time series from LCMV beamformer
  - **Sampling:** 256 Hz, 2-second epochs
  - **Regions:** Parcellated using FreeSurfer aparc labels
- **Key variable:** `labels`
  - **Type:** Cell array or list of strings
  - **Content:** Region names (e.g., 'superiorfrontal-lh', 'precentral-rh')

#### `.csv` files
- **Format:** Comma-separated values
- **Shape:** `(epochs × time_points, n_regions)`
- **Content:** Flattened version of source time series
- **Columns:** One per brain region (region names as headers)

### Processing Pipeline:
1. Raw EDF → Filtered (1-45 Hz, notch 50/100 Hz) → Resampled to 256 Hz
2. Source reconstruction using LCMV beamformer
3. Parcellation to ~68 cortical regions
4. Epoch creation: 2-second epochs per sleep stage
5. Average referencing applied

---

## 3. Results Data Structure

**Location:** `/Users/borjan/code/healthy-sleep-project/code/results/`

### 3.1 Granger Causality (GC) Results

#### Time Domain GC:
- **Files:**
  - `gc_in_time_zscores_matrix.mat` - Incoming GC z-scores
  - `gc_out_time_zscores_matrix.mat` - Outgoing GC z-scores
  - `gc_in_time_mw_effect_sizes.mat` - Mann-Whitney effect sizes (incoming)
  - `gc_out_time_mw_effect_sizes.mat` - Mann-Whitney effect sizes (outgoing)
  - `gcin_time_subject_data_struct.mat` - Subject-level incoming GC data
  - `gcout_time_subject_data_struct.mat` - Subject-level outgoing GC data

#### Frequency Domain GC:
- **Files (per frequency band):**
  - `gc_in_freq_zscores_alpha.mat` - Alpha band (8-13 Hz)
  - `gc_in_freq_zscores_beta.mat` - Beta band (13-30 Hz)
  - `gc_in_freq_zscores_delta.mat` - Delta band (0.5-4 Hz)
  - `gc_in_freq_zscores_gamma.mat` - Gamma band (30-45 Hz)
  - `gc_in_freq_zscores_theta.mat` - Theta band (4-8 Hz)
  - Similar files for `gc_out_*` (outgoing)
  - `gc_in_freq_mw_effect_sizes.mat` - Effect sizes
  - `gcin_freq_subject_data_struct.mat` - Subject-level data
  - `gcout_freq_subject_data_struct.mat` - Subject-level data

**Structure:**
- Z-scores: Statistical significance matrices (regions × regions)
- Effect sizes: Effect size matrices for group comparisons
- Subject data: Individual participant GC values

### 3.2 Mutual Information (MI) Results

- **Files:**
  - `mi_data.mat` - Time domain MI results
  - `cmi_data.mat` - Conditional MI results
  - `freq_cmi_data.mat` - Frequency domain conditional MI
  - `mi_time_zscores_matrix.mat` - Time domain z-scores
  - `mi_freq_zscores_alpha.mat` - Alpha band z-scores
  - `mi_freq_zscores_beta.mat` - Beta band z-scores
  - `mi_freq_zscores_delta.mat` - Delta band z-scores
  - `mi_freq_zscores_gamma.mat` - Gamma band z-scores
  - `mi_freq_zscores_theta.mat` - Theta band z-scores
  - `mi_time_effect_sizes.mat` - Time domain effect sizes
  - `mi_freq_effect_sizes.mat` - Frequency domain effect sizes

**Structure:**
- Similar to GC results but for mutual information measures
- Symmetric matrices (undirected connectivity)

### 3.3 Pairwise Granger Causality (PWCGC) Results

**Location:** `/code/results/pwcgc/` and `/code/results/pwcgc_lowres/`

- **Files:** `pwcgc_params_{subject}_{stage}.mat`
- **Contents:**
  - `V`: Residual covariance matrix
  - `CAK`: VAR model parameters
  - `H`: Transfer function (frequency domain)
  - `Vss`: State-space covariance
  - `eweights`: Normalized GC weights matrix

**Structure:**
- One file per subject per sleep stage
- Used for SSDI (State-Space Directed Information) analysis

### 3.4 Coarse-Grained Results

**Location:** `/code/results/coarse-grained/`

- **Files:** `coarse_pwcgc_params_{subject}_{stage}.mat`
- **Contents:** Similar to PWCGC but for coarse-grained regions
- **Purpose:** Analysis at different spatial scales (grouped brain regions)

### 3.5 SSDI Results

**Location:** `/code/results/ssdi_results/`

#### Time Domain:
- **Location:** `/code/results/ssdi_results/time/`
- **Files:** `{subject}_{stage}_ssdi_time_domain_results.mat`
- **Structure:**
  - Contains `results` structure with:
    - `dopto_m2` through `dopto_m9`: Dynamical dependence matrices for macroscopic sizes 2-9
    - Each matrix contains optimal values for different parameter combinations
- **File size:** ~3-4 MB per file
- **Total:** 145 files (29 subjects × 5 sleep stages)

#### Spectral Domain:
- **Location:** `/code/results/ssdi_results/spectral/`
- **Files:** `{subject}_{stage}_ssdi_spectral_results.mat`
- **Structure:** Similar to time domain but for frequency-resolved analysis

### 3.6 PSD (Power Spectral Density) Results

**Location:** `/code/results/PSD-raw/` and `/code/results/PSD/`

- **Structure:**
  ```
  PSD-raw/
  ├── EPCTL01/
  │   ├── EPCTL01_W_PSD_avg.png
  │   ├── EPCTL01_N1_PSD_avg.png
  │   ├── EPCTL01_N2_PSD_avg.png
  │   ├── EPCTL01_N3_PSD_avg.png
  │   └── EPCTL01_R_PSD_avg.png
  └── ... (29 subjects)
  ```
- **Content:** Average power spectral density plots per sleep stage

### 3.7 Sleep Statistics Results

**Location:** `/code/results/sleep-results/`

- **hypno/:** Hypnogram statistics and plots
  - `average_sleep_stage_distribution_including_wake.png`
  - Individual subject statistics
- **bandpower/:** Band power calculations per sleep stage

### 3.8 Labels

- **File:** `/code/results/labels.mat`
- **Content:** Brain region labels used across analyses
- **Format:** Cell array of region names (e.g., 'superiorfrontal-lh', 'precentral-rh')

---

## 4. SSDI Analysis Data

**Location:** `/code/ssdi_analysis/data/`

### Files:

#### `beta_matrices.npz`
- **Format:** NumPy compressed archive
- **Content:** Beta matrices (likely regression coefficients or optimization parameters)
- **Keys:** Sleep stage labels (W, N1, N2, N3, R)

#### `dopto_vectors.npz`
- **Format:** NumPy compressed archive
- **Content:** Dynamical dependence (dopto) vectors
- **Keys:** Sleep stage labels (W_, N1, N2, N3, R_)
- **Note:** Contains object arrays (requires `allow_pickle=True`)

#### `gw_distances.npz` and `gw_distances_new.npz`
- **Format:** NumPy compressed archive
- **Content:** Gromov-Wasserstein distances (optimal transport distances)
- **Keys:** Sleep stage labels (W_, N1, N2, N3, R_)
- **Note:** Contains object arrays (requires `allow_pickle=True`)

#### `opt_similarity_matrices.npz`
- **Format:** NumPy compressed archive
- **Content:** Optimal similarity matrices
- **Keys:** Sleep stage labels

---

## 5. Data Dimensions Summary

### Source Time Series (Preprocessed):
- **Regions:** ~68 (varies by subject, excludes regions with all-zero data)
- **Time points per epoch:** 513 (2 seconds × 256 Hz + 1)
- **Epochs:** Variable (depends on sleep stage duration)
- **Example:** 68 × 513 × 74 epochs

### Granger Causality Matrices:
- **Shape:** `(n_regions, n_regions)`
- **Type:** Directed, asymmetric
- **Values:** GC strength (normalized or z-scores)

### Mutual Information Matrices:
- **Shape:** `(n_regions, n_regions)`
- **Type:** Undirected, symmetric
- **Values:** MI strength (normalized or z-scores)

### SSDI Results:
- **Macroscopic sizes:** 2 through 9 (8 different scales)
- **For each size:** Matrix of optimal values
- **Dimensions:** Vary by macroscopic size and optimization parameters

---

## 6. Participants

- **Total subjects:** 29 (EPCTL01 through EPCTL29)
- **Sleep stages analyzed:** W (Wake), N1, N2, N3, R (REM)
- **Epoch duration:** 2 seconds (for source analysis), 30 seconds (for sleep staging)

---

## 7. Key Processing Steps

1. **Raw EEG → Preprocessed:**
   - Notch filtering (50, 100 Hz)
   - Band-pass filtering (1-45 Hz)
   - Resampling to 256 Hz
   - Channel cleaning and average referencing

2. **Preprocessed → Source Space:**
   - LCMV beamformer for source reconstruction
   - Parcellation to ~68 cortical regions
   - Epoch creation (2-second epochs)

3. **Source Time Series → Connectivity:**
   - VAR model fitting
   - Granger Causality computation
   - Mutual Information computation
   - SSDI (State-Space Directed Information) analysis

4. **Connectivity → Statistics:**
   - Z-score computation
   - Mann-Whitney U tests
   - Effect size calculations
   - Multiple comparison corrections (FDR)

---

## Notes

- The raw data directory (`/data/raw/`) may not exist in the current workspace - it's referenced in code but may be located elsewhere or require setup
- Most analyses work with preprocessed data stored in `.mat` format
- Results are primarily stored as MATLAB `.mat` files
- SSDI analysis data uses NumPy `.npz` format
- All analyses are performed per sleep stage (W, N1, N2, N3, R)


