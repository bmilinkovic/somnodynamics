# Canonical Anatomical Mapping (DMT-derived)

SomnoDynamics uses the canonical anatomical7 coarse-graining directly extracted from DMT methodology.

## Source of truth

Mapping source in DMT:
- `/Users/borjan/code/dmt_project/code/paper_style_mvgc_anatomical7/export_source_data_anatomical7.py`

Copied mapping artifact in SomnoDynamics:
- `code/coarse_graining/anatomical7_mapping_from_dmt.csv`

## Macro-region order

1. Frontal
2. Sensorimotor
3. Parietal
4. Temporal
5. Occipital
6. Cingulate
7. Limbic

## Validation script

Run:

```bash
python3 code/coarse_graining/validate_mapping.py
```

Outputs:
- `results/coarse_grained/anatomical7_mapping_counts.csv`
- `results/coarse_grained/anatomical7_mapping_summary.txt`

## Rule

The mapping must not be guessed or manually altered without explicit methodological justification and versioned documentation.
