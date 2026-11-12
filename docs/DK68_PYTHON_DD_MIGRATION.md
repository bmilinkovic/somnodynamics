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
2. Add a Python DK68 runner that reads the existing source recon `.mat` files,
   reproduces pooled z-scoring and VAR LRT/LWR selection, then invokes the
   packaged proxy-plus-canonical spectral optimisation with the exact DK68
   grid. Write compatible result fields, selected model metadata and status.
3. Make runs resumable and support one file per Slurm array task. Keep the
   legacy MATLAB files intact as a reference.
4. Update the HPC environment and submitters so DK68 uses Python from source
   files through summary/figures, with a separate aggregation job after the
   array succeeds. Document the setup and resource expectations.
5. Verify CLI help, a small synthetic end-to-end run, output compatibility,
   shell syntax, and clean package import. The full 143-file computation cannot
   run locally because the external data volume is not mounted.

## Implemented workflow

The vendored package is `code/external/emergence-multiscale/`. Its
`UPSTREAM.md` records the precise source revision. The runner is
`code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/run_python_dd_dk68.py`.
It reads `data` or `source_ts` from the original DK68 `.mat` files,
normalises each region across samples and epochs with the MATLAB branch's
pooled sample standard deviation, selects a VAR order with LRT (maximum 30)
using LWR fitting, and passes the resulting `A, V` to Emergence Engine.
Every dimension × band block runs the DDx proxy, clustering, then canonical
spectral DD. Output files include the existing summary scripts' `dopto` and
`gopto` fields, plus model and optimisation metadata.

The default HPC submitter `hpc/submit_ssdi_dk68_lrt30.sh` is now a Slurm
array over 143 subject-stage inputs. `hpc/submit_ssdi_dk68_python_summary.sh`
builds group tables and figures after all array tasks succeed. Python files
go to a new output root, so archived MATLAB outputs are untouched. Model and
dimension–band checkpoints are held under that output root and keyed by source
file metadata, analysis settings and engine revision, allowing interrupted
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
The shell scripts pass `bash -n`. The real 143-file dataset is unavailable
on this machine because the external data volume is not mounted; no real
DK68 values or MATLAB parity claim is made here. The random-start streams
need not be identical to MATLAB's even with the same restart count, so
compare optimised DD distributions and solutions rather than expecting
bitwise agreement.
