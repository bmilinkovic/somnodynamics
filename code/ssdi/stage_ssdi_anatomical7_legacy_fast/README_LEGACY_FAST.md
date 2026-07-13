# Legacy Fast SSDI Implementation (Preserved)

This directory preserves the earlier SomnoDynamics stage SSDI implementation that ran unusually fast.

It is retained for reproducibility/audit only and is **not** the authoritative SSDI pipeline.

Use the corrected DMT-faithful pipeline in:
- `code/ssdi/stage_ssdi_anatomical7`

Legacy path kept unchanged for comparison:
- single-stage optimization flow
- reduced default optimizer settings (`nrunsp=30`, `nitersp=3000`)
- reduced diagnostics/output schema versus DMT
