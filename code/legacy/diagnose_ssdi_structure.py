#!/usr/bin/env python3
"""Diagnose the actual structure of SSDI .mat files"""
from pathlib import Path
from scipy.io import loadmat
import numpy as np

# Check time domain file
time_file = Path("code/results/ssdi_results/time/EPCTL01_W__ssdi_time_domain_results.mat")
print("="*60)
print("TIME DOMAIN FILE STRUCTURE")
print("="*60)
if time_file.exists():
    mat = loadmat(str(time_file), squeeze_me=True, struct_as_record=False)
    if 'results' in mat:
        results = mat['results']
        print(f"Results type: {type(results)}")
        if hasattr(results, 'dtype') and hasattr(results.dtype, 'names'):
            print(f"Field names: {results.dtype.names[:20] if results.dtype.names else 'None'}")
            # Check for dopto keys
            dopto_keys = [n for n in (results.dtype.names or []) if 'dopto' in n.lower()]
            print(f"Dopto keys found: {dopto_keys[:10]}")
            # Check for gopto keys
            gopto_keys = [n for n in (results.dtype.names or []) if 'gopto' in n.lower()]
            print(f"Gopto keys found: {gopto_keys[:10]}")
        else:
            print("Results is not a structured array")
            if isinstance(results, dict):
                print(f"Results keys: {list(results.keys())[:20]}")
else:
    print(f"File not found: {time_file}")

# Check spectral file
spectral_file = Path("code/results/ssdi_results/spectral/EPCTL01_N1_ssdi_spectral_domain_results.mat")
print("\n" + "="*60)
print("SPECTRAL FILE STRUCTURE")
print("="*60)
if spectral_file.exists():
    mat = loadmat(str(spectral_file), squeeze_me=True, struct_as_record=False)
    if 'results' in mat:
        results = mat['results']
        print(f"Results type: {type(results)}")
        if hasattr(results, 'dtype') and hasattr(results.dtype.names):
            print(f"Field names (first 30): {results.dtype.names[:30] if results.dtype.names else 'None'}")
            # Check for dopto keys
            dopto_keys = [n for n in (results.dtype.names or []) if 'dopto' in n.lower()]
            print(f"Dopto keys found: {dopto_keys[:10]}")
            # Check for Lopto keys
            lopto_keys = [n for n in (results.dtype.names or []) if 'lopto' in n.lower()]
            print(f"Lopto keys found: {lopto_keys[:10]}")
            
            # Check structure of one Lopto field
            if lopto_keys:
                first_lopto = lopto_keys[0]
                print(f"\nExamining {first_lopto}:")
                lopto_val = getattr(results, first_lopto) if hasattr(results, first_lopto) else results[first_lopto][0, 0]
                print(f"  Type: {type(lopto_val)}")
                if isinstance(lopto_val, np.ndarray):
                    print(f"  Shape: {lopto_val.shape}")
                    print(f"  Dtype: {lopto_val.dtype}")
                    if lopto_val.dtype == object:
                        print(f"  First element type: {type(lopto_val.flat[0])}")
                        if len(lopto_val) > 0:
                            first_elem = lopto_val.flat[0]
                            print(f"  First element: {first_elem}")
                            if isinstance(first_elem, np.ndarray):
                                print(f"    Shape: {first_elem.shape}")
                elif isinstance(lopto_val, (list, tuple)):
                    print(f"  Length: {len(lopto_val)}")
                    if len(lopto_val) > 0:
                        print(f"  First element type: {type(lopto_val[0])}")
        else:
            print("Results is not a structured array")
            if isinstance(results, dict):
                print(f"Results keys: {list(results.keys())[:20]}")
else:
    print(f"File not found: {spectral_file}")


