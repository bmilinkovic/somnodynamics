# SomnoDynamics Validation

## Stage Validation Entry Point

Run:

```bash
python3 code/run_all/validate_pipeline_outputs.py
```

Outputs:
- `results/validation_report.csv`
- `results/validation_summary.txt`

## What it checks

- migration plan exists
- canonical anatomical mapping file exists
- mapping summary exists
- core analysis output tables (sensor power, GC, SSDI, cross-analysis)
- consolidated figure manifest

## Strict mode

```bash
python3 code/run_all/validate_pipeline_outputs.py --strict
```

Returns non-zero exit code when required artifacts are missing.

## Legacy comparison

Run:

```bash
python3 code/run_all/compare_with_legacy.py
```

Outputs:
- `results/legacy_comparison_report.csv`
- `results/legacy_comparison_summary.txt`

This compares new-vs-legacy artifact availability and basic subject/stage coverage where possible.
