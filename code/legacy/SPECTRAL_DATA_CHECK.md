# Spectral SSDI Data Structure Check

## Summary

Based on my examination:

### 1. `dopto_vectors.npz` Structure
- **Location:** `/code/ssdi_analysis/data/dopto_vectors.npz`
- **Current structure:** 
  ```
  Stage (W_, N1, N2, N3, R_) 
    → Subject ID ('01'-'29') 
      → Macroscopic Size (2-9) 
        → Array (1, 100) - dopto values for candidates
  ```
- **Does NOT contain frequency-specific data** - only macroscopic sizes

### 2. Expected Structure for Cross-Frequency Analysis

Your `ssdi_cross_frequency_tests.py` expects:
- `dopto_vectors_{band}.npz` for each frequency band
- `beta_matrices_{band}.npz` for each frequency band
- Where `{band}` is one of: `broadband`, `delta`, `theta`, `alpha`, `sigma`, `beta`, `gamma`

### 3. Spectral Domain Results Files

**Location:** `/code/results/ssdi_results/spectral/`

- **Files:** `{subject}_{stage}_ssdi_spectral_domain_results.mat`
- **Total:** 145 files (29 subjects × 5 sleep stages)
- **Likely contains:** Frequency-specific dopto data organized by frequency bands

### 4. What Needs to Be Done

The frequency-specific dopto data is likely stored in the spectral `.mat` files, but needs to be extracted and reorganized into the format expected by `ssdi_cross_frequency_tests.py`.

**Recommended approach:**
1. Load each spectral `.mat` file
2. Extract frequency-specific dopto data (likely organized by frequency bands)
3. Reorganize into the expected structure:
   ```
   dopto_vectors_{band}.npz:
     Stage → Subject ID → Macro Size → Array of dopto values
   ```
4. Save as separate `.npz` files for each frequency band

### 5. Next Steps

To verify the structure of spectral files, you'll need to:
1. Load a sample spectral `.mat` file using `scipy.io.loadmat()`
2. Examine the `results` structure
3. Look for frequency band keys or nested frequency structures
4. Extract and reorganize the data accordingly

Would you like me to create a script to:
- Extract frequency-specific dopto data from spectral `.mat` files?
- Reorganize it into the expected `dopto_vectors_{band}.npz` format?
- Do the same for `beta_matrices_{band}.npz`?


