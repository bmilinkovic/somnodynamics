# SSDI Analysis Data Structure Summary

## Location
`/code/ssdi_analysis/data/`

## Overview
All files are NumPy compressed archives (`.npz` format) containing nested dictionary structures. The data is organized by sleep stage and subject ID.

---

## File 1: `beta_matrices.npz`

### Structure:
```
beta_matrices.npz
├── 'W_' (Wake)
│   └── '01' to '29' (subject IDs)
│       └── 2 to 9 (macroscopic sizes)
│           └── numpy array: (68, 100)
├── 'N1' (Stage 1)
│   └── (same structure)
├── 'N2' (Stage 2)
├── 'N3' (Stage 3/Deep sleep)
└── 'R_' (REM)
```

### Content:
- **Purpose:** Beta matrices (likely regression coefficients or optimization parameters)
- **Organization:** 
  - Top level: Sleep stages (W_, N1, N2, N3, R_)
  - Second level: Subject IDs ('01' through '29')
  - Third level: Macroscopic sizes (2, 3, 4, 5, 6, 7, 8, 9)
  - Data: NumPy arrays of shape `(68, 100)`
- **Array details:**
  - **Shape:** `(68, 100)` - 68 brain regions × 100 parameter values
  - **Dtype:** `float64`
  - **Value range:** Typically 0.0 to ~0.4
  - **Mean:** ~0.029 (varies by subject and condition)

### Access pattern:
```python
import numpy as np
data = np.load('beta_matrices.npz', allow_pickle=True)
beta_W_subject01_size2 = data['W_'].item()['01'][2]  # Shape: (68, 100)
```

---

## File 2: `dopto_vectors.npz`

### Structure:
```
dopto_vectors.npz
├── 'W_' (Wake)
│   └── '01' to '29' (subject IDs)
│       └── (nested dict structure - likely macroscopic sizes)
├── 'N1', 'N2', 'N3', 'R_'
│   └── (same structure)
```

### Content:
- **Purpose:** Dynamical dependence (dopto) vectors
- **Organization:**
  - Top level: Sleep stages (W_, N1, N2, N3, R_)
  - Second level: Subject IDs ('01' through '29')
  - Third level: Likely macroscopic sizes (2-9) or other parameters
- **Note:** Contains nested dictionaries, structure similar to beta_matrices

### Access pattern:
```python
data = np.load('dopto_vectors.npz', allow_pickle=True)
dopto_W_subject01 = data['W_'].item()['01']  # Returns nested dict
```

---

## File 3: `gw_distances.npz`

### Structure:
```
gw_distances.npz
├── 'W_' (Wake)
│   └── 'W_', 'N1', 'N2', 'N3', 'R_' (comparison sleep stages)
│       └── (nested structure)
├── 'N1', 'N2', 'N3', 'R_'
│   └── (same structure - each stage compared to all others)
```

### Content:
- **Purpose:** Gromov-Wasserstein distances between sleep stages
- **Organization:**
  - Top level: Reference sleep stage (W_, N1, N2, N3, R_)
  - Second level: Comparison sleep stage (W_, N1, N2, N3, R_)
  - This creates a 5×5 matrix of pairwise comparisons
- **Meaning:** GW distances measure optimal transport distances between probability distributions, used to compare connectivity patterns across sleep stages

### Access pattern:
```python
data = np.load('gw_distances.npz', allow_pickle=True)
gw_W_to_N1 = data['W_'].item()['N1']  # GW distance from Wake to N1
```

---

## File 4: `gw_distances_new.npz`

### Structure:
Same as `gw_distances.npz`

### Content:
- **Purpose:** Updated or recalculated Gromov-Wasserstein distances
- **Organization:** Identical structure to `gw_distances.npz`
- **Note:** Likely a newer version with updated calculations or different parameters

---

## File 5: `opt_similarity_matrices.npz`

### Structure:
```
opt_similarity_matrices.npz
├── 'W_' (Wake)
│   └── '01' to '29' (subject IDs)
│       └── (nested dict structure)
├── 'N1', 'N2', 'N3', 'R_'
│   └── (same structure)
```

### Content:
- **Purpose:** Optimal similarity matrices (likely from optimal transport or matching algorithms)
- **Organization:**
  - Top level: Sleep stages (W_, N1, N2, N3, R_)
  - Second level: Subject IDs ('01' through '29')
  - Third level: Likely macroscopic sizes or other parameters
- **Meaning:** These matrices likely represent optimal matching/coupling between brain regions or connectivity patterns

### Access pattern:
```python
data = np.load('opt_similarity_matrices.npz', allow_pickle=True)
similarity_W_subject01 = data['W_'].item()['01']  # Returns nested dict
```

---

## Common Patterns Across All Files

### 1. Sleep Stage Keys:
- `'W_'` - Wake
- `'N1'` - Stage 1 sleep
- `'N2'` - Stage 2 sleep
- `'N3'` - Stage 3/Deep sleep
- `'R_'` - REM sleep

### 2. Subject IDs:
- String format: `'01'` through `'29'` (29 subjects total)
- Corresponds to EPCTL01 through EPCTL29

### 3. Macroscopic Sizes:
- Integer keys: `2, 3, 4, 5, 6, 7, 8, 9`
- Represents different levels of spatial coarse-graining

### 4. Brain Regions:
- 68 regions (standard FreeSurfer aparc parcellation)
- Used in arrays of shape `(68, ...)`

---

## Data Relationships

1. **beta_matrices**: Optimization parameters for SSDI analysis
2. **dopto_vectors**: Dynamical dependence measures (the main outcome of SSDI)
3. **gw_distances**: Pairwise distances between sleep stages (stage-to-stage comparisons)
4. **opt_similarity_matrices**: Optimal matching matrices (likely used in GW distance calculations)

---

## Usage Notes

### Loading the data:
```python
import numpy as np

# Load any file
data = np.load('beta_matrices.npz', allow_pickle=True)

# Access structure
sleep_stage = 'W_'
subject_id = '01'
macro_size = 2

# For beta_matrices
beta = data[sleep_stage].item()[subject_id][macro_size]  # Shape: (68, 100)

# For gw_distances (different structure)
gw = data['W_'].item()['N1']  # GW distance from Wake to N1
```

### Important:
- All files require `allow_pickle=True` when loading due to object arrays
- Use `.item()` to extract dictionaries from object arrays
- Keys are strings (subject IDs) or integers (macroscopic sizes)

---

## Summary Statistics

### beta_matrices:
- **Dimensions per entry:** (68, 100)
- **Total entries:** 5 stages × 29 subjects × 8 sizes = 1,160 arrays
- **Value range:** ~0.0 to ~0.4
- **Typical mean:** ~0.029

### File sizes:
- Files are relatively compact due to NumPy compression
- Actual size depends on nesting depth and array sizes

---

## Analysis Context

These files are outputs from State-Space Directed Information (SSDI) analysis, which:
1. Analyzes directed connectivity at multiple spatial scales (macroscopic sizes 2-9)
2. Computes dynamical dependence measures (dopto)
3. Uses optimal transport (Gromov-Wasserstein) to compare connectivity patterns
4. Provides similarity matrices for matching brain regions/patterns across conditions

The data supports comparisons:
- Across sleep stages (W, N1, N2, N3, R)
- Across subjects (29 participants)
- Across spatial scales (macroscopic sizes 2-9)


