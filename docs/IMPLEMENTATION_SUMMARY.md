# SomnoDynamics Implementation Summary

## Implemented

- Modular repository structure with legacy isolation under `code/legacy`
- Centralized configuration in `code/run_all/somnodynamics_config.toml`
- Sensor preprocessing pipeline integration
- Sensor absolute/relative power analysis with repeated-measures statistics
- Source reconstruction pipeline: anatomical46 reconstruction followed by anatomical7 coarse-graining
- Legacy Yeo7 source reconstruction wrapper kept for compatibility only
- Canonical anatomical mapping import/validation from DMT
- DMT-style stage-based GC module with MVGC2 vendored locally and repeated-measures multi-stage statistics
- DMT-style stage-based SSDI module with local `ssdi-1` integration and repeated-measures summary statistics
- GC-SSDI cross-analysis tables and figures
- Pipeline orchestrator for per-stage and full runs
- Consolidated figure build stage and manifest export
- Output validation stage with required/optional artifact checks

## Output conventions

- Machine-readable tables in CSV
- Figures in PNG and SVG
- Human-readable summaries in TXT/Markdown
- Standardized stage labels in all downstream exports:
  - `Wake`, `N1`, `N2`, `N3`, `REM`

## Known constraints

- MVGC2 is vendored in-repo under `code/external/MVGC2`; SSDI still depends on available toolbox install
- Full numeric validation requires complete dataset execution and parity benchmarking against legacy products

## Primary entrypoint

```bash
python3 code/run_all/run_pipeline.py all
```
