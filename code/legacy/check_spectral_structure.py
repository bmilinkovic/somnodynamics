#!/usr/bin/env python3
"""
Check the structure of spectral SSDI results to see if they contain frequency-specific dopto data
"""
import os
from pathlib import Path

# Try to import scipy, but provide instructions if not available
try:
    import scipy.io
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("scipy not available. Install with: pip install scipy")

def examine_spectral_file(filepath):
    """Examine a spectral .mat file structure"""
    if not SCIPY_AVAILABLE:
        print(f"Cannot examine {filepath} - scipy not available")
        return None
    
    print(f"\n{'='*60}")
    print(f"Examining: {os.path.basename(filepath)}")
    print(f"{'='*60}")
    
    try:
        data = scipy.io.loadmat(filepath, simplify_cells=True)
        print(f"Top-level keys: {list(data.keys())}")
        
        if 'results' in data:
            results = data['results']
            print(f"\nResults type: {type(results)}")
            
            if isinstance(results, dict):
                print(f"Results keys ({len(results)} total):")
                # Look for frequency-related keys
                all_keys = list(results.keys())
                freq_keys = [k for k in all_keys if any(b in str(k).lower() for b in ['freq', 'band', 'delta', 'theta', 'alpha', 'beta', 'gamma', 'sigma'])]
                dopto_keys = [k for k in all_keys if 'dopto' in str(k).lower()]
                
                print(f"  Frequency-related keys: {freq_keys[:10]}")
                print(f"  Dopto-related keys: {dopto_keys[:10]}")
                
                # Check structure of dopto keys
                if dopto_keys:
                    print(f"\nExamining first dopto key: {dopto_keys[0]}")
                    dopto_val = results[dopto_keys[0]]
                    print(f"  Type: {type(dopto_val)}")
                    if hasattr(dopto_val, 'shape'):
                        print(f"  Shape: {dopto_val.shape}")
                    if hasattr(dopto_val, 'dtype'):
                        print(f"  Dtype: {dopto_val.dtype}")
                
                # Check if there are nested frequency structures
                print(f"\nFirst 20 keys in results:")
                for key in all_keys[:20]:
                    val = results[key]
                    val_type = type(val).__name__
                    if hasattr(val, 'shape'):
                        print(f"  {key}: {val_type}, shape={val.shape}")
                    elif isinstance(val, dict):
                        print(f"  {key}: {val_type}, {len(val)} sub-keys")
                    else:
                        print(f"  {key}: {val_type}")
            else:
                print(f"Results is not a dict, it's a {type(results)}")
        else:
            print("No 'results' key found")
            print(f"Available keys: {list(data.keys())}")
            
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
    
    return None

if __name__ == "__main__":
    spectral_dir = Path("code/results/ssdi_results/spectral")
    
    if not spectral_dir.exists():
        print(f"Directory not found: {spectral_dir}")
    else:
        # Get first file
        files = list(spectral_dir.glob("*.mat"))
        if files:
            examine_spectral_file(files[0])
        else:
            print(f"No .mat files found in {spectral_dir}")


