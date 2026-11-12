# DK68 Python dynamical dependence migration plan

Audit on 12 November 2026: the current DK68 branch fits 256 Hz, two-second
source epochs with pooled z-scoring, LRT model selection up to order 30, and
LWR regression. Its SSDI stage searches macro dimensions 1–10, with 100 proxy
starts, 10,000 proxy and 10,000 spectral iterations, a 256-bin frequency grid,
and inclusive delta 1–4, theta 5–8, alpha 9–12, sigma 13–16, and beta 17–30 Hz
bands. It expects 143 subject-stage files, writes stage folders, and the Python
summary and plotting scripts consume `*_ssdi_spectral_domain_results.mat` with
`<band>_dopto_m<m>` and `<band>_gopto_m<m>` fields. The DMT reference path named
in AGENTS.md is unavailable on this machine; the repository's existing DMT
parity audits and DK68 MATLAB implementation document the established method.

Implementation:

1. Vendor a pinned copy of the standalone emergence-engine multiscale package
   under `code/external`, keeping its licence and provenance.
2. Fit source recordings locally with the packaged full VAR LRT/LWR routine,
   save compact `A, V` model files, and transfer only those models to the HPC.
   The Python DK68 runner then invokes the packaged proxy-plus-canonical
   spectral optimisation with the exact DK68 grid.
3. Make runs resumable and support one file per Slurm array task. Keep the
   legacy MATLAB files intact as a reference.
4. Update the HPC environment and submitters so DK68 uses Python models
   through summary/figures, with a separate aggregation job after the
   array succeeds. Document the setup and resource expectations.
5. Verify the full local fit on real source data, a small DD run from a saved
   model, output compatibility, shell syntax, and clean package import.

## Implemented workflow

The vendored package is `code/external/emergence-multiscale/`. Its
`UPSTREAM.md` records the precise source revision. The local fitting script
is `fit_python_var_dk68.py` with serial launcher `fit_all_python_var_dk68.sh`
in the standalone DK68 `ssdi/` directory. It reads `data` or `source_ts`
from the original DK68 `.mat` files, normalises each region across samples
and epochs with the MATLAB branch's pooled sample standard deviation, and
uses the unchanged Emergence Engine `fit_var_auto` routine for LRT model
selection (maximum 30) and LWR fitting. It saves the selected `A, V`,
model order, candidate order diagnostics, labels and source provenance.
The HPC runner is `run_python_dd_dk68.py`; it loads these compact models.
Every dimension × band block runs the DDx proxy, clustering, then canonical
spectral DD. Output files include the existing summary scripts' `dopto` and
`gopto` fields, plus model and optimisation metadata.

The default HPC submitter `hpc/submit_ssdi_dk68_lrt30.sh` is now a Slurm
array over 143 subject-stage inputs. `hpc/submit_ssdi_dk68_python_summary.sh`
builds group tables and figures after all array tasks succeed. Python files
go to a new output root, so archived MATLAB outputs are untouched. Model and
dimension–band checkpoints are held under that output root and keyed by source
model metadata, analysis settings and engine revision, allowing interrupted
array tasks to resume without repeating completed blocks.

The transfer-function and CAK calculations in Emergence Engine use exact
block elimination of the VAR companion equations. A small-model parity test
checks them against the dense state-space calculations. This matters for DK68:
the dense companion state can reach 2,040 coordinates at model order 30.

## Verification and limits

The runner completed a synthetic subject-stage recording end to end. The
existing SomnoDynamics summary reader extracted its `dd_best` and run count,
and the figure reader extracted its structure-fragmentation index. The
vendored package's companion parity and time-changing tests pass. A synthetic
interruption test confirmed that a finished block and VAR fit are reused while
the missing block is recomputed.
The shell scripts pass `bash -n`. The external drive is mounted locally;
the full 143-file fitting loop has begun. A real DK68 source file (EPCTL06
N1) fitted at order 15 and produced a 549 KB compact model. A short DD smoke
run read that model successfully. The random-start streams
need not be identical to MATLAB's even with the same restart count, so
compare optimised DD distributions and solutions rather than expecting
bitwise agreement.
