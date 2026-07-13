#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract frequency-specific dopto data from spectral SSDI results and reorganize
into dopto_vectors_{band}.npz files for cross-frequency analysis.

This script:
1. Loads all spectral .mat files from code/results/ssdi_results/spectral/
2. Extracts dopto data for each frequency band
3. Reorganizes into: Stage → Subject → Macro Size → Array
4. Saves as dopto_vectors_{band}.npz files in code/ssdi_analysis/data/
"""

import os
from pathlib import Path
import numpy as np

try:
    import scipy.io
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("ERROR: scipy is required but not installed.")
    print("Please install it with: pip install scipy")
    print("Then run this script again.")
    exit(1)

# Configuration
SPECTRAL_DIR = Path("code/results/ssdi_results/spectral")
OUTPUT_DIR = Path("code/ssdi_analysis/data")

# Frequency bands (excluding broadband - not in spectral files)
FREQUENCY_BANDS = ['delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma']

# Sleep stages
STAGES = ['W_', 'N1', 'N2', 'N3', 'R_']

# Macroscopic sizes
MACRO_SIZES = list(range(2, 10))  # 2 to 9

# Subject IDs
SUBJECT_IDS = [f"{i:02d}" for i in range(1, 30)]  # '01' to '29'


def extract_spectral_dopto():
    """
    Extract frequency-specific dopto data from spectral .mat files.
    
    Returns
    -------
    dopto_by_band : dict
        dopto_by_band[band][stage][subject][size] = array of dopto values
    """
    # Initialize structure: band -> stage -> subject -> size -> array
    dopto_by_band = {
        band: {stage: {} for stage in STAGES}
        for band in FREQUENCY_BANDS
    }
    
    # Track what we find
    files_processed = 0
    files_with_errors = []
    
    print("Extracting frequency-specific dopto data from spectral files...")
    print(f"Looking in: {SPECTRAL_DIR}")
    print(f"Frequency bands: {FREQUENCY_BANDS}")
    print()
    
    # Process all spectral .mat files
    for filename in sorted(os.listdir(SPECTRAL_DIR)):
        if not filename.endswith('_ssdi_spectral_domain_results.mat'):
            continue
        
        filepath = SPECTRAL_DIR / filename
        files_processed += 1
        
        # Parse filename: EPCTL{subject}_{stage}_ssdi_spectral_domain_results.mat
        parts = filename.split('_')
        if len(parts) < 3:
            print(f"Warning: Could not parse filename {filename}")
            files_with_errors.append(filename)
            continue
        
        subject_part = parts[0]  # EPCTL01, EPCTL02, etc.
        stage = parts[1]  # W, N1, N2, N3, R (may have trailing underscore)
        subject_id = subject_part[-2:]  # '01', '02', etc.
        
        # Normalize stage name (add underscore if needed)
        if stage not in STAGES:
            if stage + '_' in STAGES:
                stage = stage + '_'
            elif stage.rstrip('_') + '_' in STAGES:
                stage = stage.rstrip('_') + '_'
            else:
                print(f"Warning: Unknown stage '{stage}' in {filename}")
                files_with_errors.append(filename)
                continue
        
        # Load the .mat file
        try:
            mat_data = scipy.io.loadmat(str(filepath), simplify_cells=True)
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            files_with_errors.append(filename)
            continue
        
        # Access the 'results' structure
        if 'results' not in mat_data:
            print(f"Warning: No 'results' key in {filename}")
            files_with_errors.append(filename)
            continue
        
        results = mat_data['results']
        
        # Extract dopto data for each frequency band and macro size
        for band in FREQUENCY_BANDS:
            # Initialize subject dict if needed
            if subject_id not in dopto_by_band[band][stage]:
                dopto_by_band[band][stage][subject_id] = {}
            
            for size in MACRO_SIZES:
                # Key format: {band}_dopto_m{size}
                key = f"{band}_dopto_m{size}"
                
                if key in results:
                    value = results[key]
                    
                    # Convert to numpy array and ensure it's 1D
                    if isinstance(value, (list, np.ndarray)):
                        arr = np.asarray(value).ravel()
                    else:
                        # Single value - convert to array
                        arr = np.array([float(value)]).ravel()
                    
                    # Store in the structure
                    dopto_by_band[band][stage][subject_id][size] = arr
                else:
                    # Key not found - this is okay, some subjects may not have all data
                    pass
    
    print(f"\nProcessed {files_processed} files")
    if files_with_errors:
        print(f"Files with errors: {len(files_with_errors)}")
        for f in files_with_errors[:5]:  # Show first 5
            print(f"  - {f}")
        if len(files_with_errors) > 5:
            print(f"  ... and {len(files_with_errors) - 5} more")
    
    return dopto_by_band


def save_dopto_vectors(dopto_by_band):
    """
    Save dopto data organized by frequency band into .npz files.
    
    Parameters
    ----------
    dopto_by_band : dict
        Structure: band -> stage -> subject -> size -> array
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\nSaving dopto vectors to {OUTPUT_DIR}...")
    
    for band in FREQUENCY_BANDS:
        # Prepare data for this band in the expected format
        # Format: stage -> subject -> size -> array
        band_data = {}
        
        for stage in STAGES:
            band_data[stage] = {}
            for subject_id in SUBJECT_IDS:
                if subject_id in dopto_by_band[band][stage]:
                    band_data[stage][subject_id] = {}
                    for size in MACRO_SIZES:
                        if size in dopto_by_band[band][stage][subject_id]:
                            band_data[stage][subject_id][size] = dopto_by_band[band][stage][subject_id][size]
        
        # Convert to numpy object array format (similar to existing dopto_vectors.npz)
        # Each stage becomes an object array containing a dict
        output_dict = {}
        for stage in STAGES:
            # Create object array with dict inside
            stage_dict = band_data[stage]
            stage_array = np.empty((), dtype=object)
            stage_array[()] = stage_dict
            output_dict[stage] = stage_array
        
        # Save as .npz file
        output_file = OUTPUT_DIR / f"dopto_vectors_{band}.npz"
        np.savez(output_file, **output_dict)
        
        # Count how much data we saved
        total_entries = sum(
            len(subj_dict)
            for stage_dict in band_data.values()
            for subj_dict in stage_dict.values()
        )
        
        print(f"  Saved {output_file.name}: {total_entries} subject-size combinations across {len(STAGES)} stages")


def verify_output():
    """Verify that the output files were created correctly."""
    print("\nVerifying output files...")
    
    for band in FREQUENCY_BANDS:
        output_file = OUTPUT_DIR / f"dopto_vectors_{band}.npz"
        
        if not output_file.exists():
            print(f"  WARNING: {output_file.name} not found!")
            continue
        
        try:
            data = np.load(str(output_file), allow_pickle=True)
            print(f"  {output_file.name}:")
            print(f"    Keys: {list(data.keys())}")
            
            # Check structure of first stage
            if STAGES[0] in data:
                stage_data = data[STAGES[0]].item()
                n_subjects = len(stage_data)
                print(f"    {STAGES[0]}: {n_subjects} subjects")
                
                # Check first subject
                if n_subjects > 0:
                    first_subj = list(stage_data.keys())[0]
                    subj_data = stage_data[first_subj]
                    n_sizes = len(subj_data)
                    print(f"      First subject ({first_subj}): {n_sizes} macro sizes")
                    
                    # Check first size
                    if n_sizes > 0:
                        first_size = list(subj_data.keys())[0]
                        arr = subj_data[first_size]
                        print(f"        First size ({first_size}): array shape {arr.shape}, dtype {arr.dtype}")
        except Exception as e:
            print(f"  ERROR loading {output_file.name}: {e}")


def main():
    """Main function."""
    print("="*60)
    print("Extracting Frequency-Specific Dopto Data from Spectral Files")
    print("="*60)
    print()
    
    # Check if spectral directory exists
    if not SPECTRAL_DIR.exists():
        print(f"ERROR: Spectral directory not found: {SPECTRAL_DIR}")
        print("Please check the path and try again.")
        return
    
    # Extract data
    dopto_by_band = extract_spectral_dopto()
    
    # Save to .npz files
    save_dopto_vectors(dopto_by_band)
    
    # Verify output
    verify_output()
    
    print("\n" + "="*60)
    print("Extraction complete!")
    print(f"Output files saved to: {OUTPUT_DIR}")
    print("="*60)


if __name__ == "__main__":
    main()

