#!/usr/bin/env python3
"""Fit one or more DK68 source recordings with Emergence Engine's full LWR VAR.

This uses multiscale.fit_var_auto unchanged: LRT order selection up to 30 and
the original LWR regression at the selected order. Each compact output keeps
the selected A and V, model-order diagnostics, labels, and source provenance.
"""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import numpy as np

from run_python_dd_dk68 import (REPO, atomic_npz, discover, read_source,
                                STAGES)
from multiscale import fit_var_auto


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    data_root = Path(os.environ.get("SOMNODYNAMICS_DATA_ROOT",
                                   "/Volumes/ex_data/data_somnoemergence"))
    parser.add_argument("--source-root", type=Path, default=data_root /
                        "results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk68")
    parser.add_argument("--model-root", type=Path, default=data_root /
                        "results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz/model_inputs")
    parser.add_argument("--task-index", type=int, help="Zero-based subject-stage index")
    parser.add_argument("--expected-count", type=int, default=143)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--pmax", type=int, default=30)
    parser.add_argument("--criterion", choices=("AIC", "BIC", "HQC", "LRT"), default="LRT")
    parser.add_argument("--regmode", choices=("LWR", "OLS"), default="LWR")
    return parser.parse_args()


def fit_one(subject: str, stage: str, source: Path, model_root: Path,
            *, pmax: int, criterion: str, regmode: str,
            overwrite: bool = False) -> Path:
    destination = model_root / stage / f"{subject}_var_model.npz"
    if destination.exists() and not overwrite:
        print(f"[{subject} {stage}] already fitted: {destination}", flush=True)
        return destination
    started = time.monotonic()
    print(f"[{subject} {stage}] loading {source}", flush=True)
    signals, labels = read_source(source, 68)
    print(f"[{subject} {stage}] full {criterion}/{regmode} fit, pmax={pmax}, "
          f"shape={signals.shape}", flush=True)
    selected = fit_var_auto(signals, pmax=pmax, criterion=criterion,
                            regmode=regmode)
    if selected.order < 1:
        raise ValueError(f"{subject} {stage}: LRT selected order zero")
    model = selected.fit
    if not np.isfinite(model.coefficients).all() or not np.isfinite(model.covariance).all():
        raise ValueError(f"{subject} {stage}: non-finite fitted parameters")
    atomic_npz(
        destination,
        A=model.coefficients, V=model.covariance,
        order=np.asarray(selected.order),
        labels=np.asarray(labels, dtype=str),
        subject=np.asarray(subject), stage=np.asarray(stage),
        fs=np.asarray(256.0), pmax=np.asarray(pmax),
        criterion=np.asarray(criterion), regmode=np.asarray(regmode),
        aic_order=np.asarray(selected.candidates.aic),
        bic_order=np.asarray(selected.candidates.bic),
        hqc_order=np.asarray(selected.candidates.hqc),
        lrt_order=np.asarray(selected.candidates.lrt),
        aic_curve=selected.candidates.aic_curve,
        bic_curve=selected.candidates.bic_curve,
        hqc_curve=selected.candidates.hqc_curve,
        source_file=np.asarray(str(source)),
        source_size=np.asarray(source.stat().st_size),
        source_mtime_ns=np.asarray(source.stat().st_mtime_ns),
        samples=np.asarray(signals.shape[1]),
        trials=np.asarray(signals.shape[2]),
    )
    print(f"[{subject} {stage}] saved order={selected.order} "
          f"in {time.monotonic()-started:.1f}s: {destination}", flush=True)
    return destination


def main() -> int:
    args = parse_args()
    jobs = discover(args.source_root)
    if len(jobs) != args.expected_count:
        raise SystemExit(f"expected {args.expected_count} source files, found {len(jobs)}")
    if args.task_index is not None:
        if not 0 <= args.task_index < len(jobs):
            raise SystemExit(f"task index must be 0..{len(jobs)-1}")
        jobs = [jobs[args.task_index]]
    print(f"Fitting {len(jobs)} DK68 subject-stage file(s); stages={STAGES}",
          flush=True)
    for index, (subject, stage, source) in enumerate(jobs, 1):
        print(f"file {index}/{len(jobs)}", flush=True)
        fit_one(subject, stage, source, args.model_root, pmax=args.pmax,
                criterion=args.criterion, regmode=args.regmode,
                overwrite=args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
