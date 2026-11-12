#!/usr/bin/env python3
"""Fit and optimise DK68 spectral DD from the existing sleep source epochs.

Each invocation handles one subject-stage file (or a small selected set).
Results are written in the stage/subject MATLAB layout used by the existing
SomnoDynamics summaries and figures. The numerical engine is the vendored
Emergence Engine multiscale package; this file only handles study inputs,
parameters, progress, and output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat


REPO = Path(__file__).resolve().parents[3]
VENDORED_SRC = REPO / "code/external/emergence-multiscale/src"
if VENDORED_SRC.is_dir():
    sys.path.insert(0, str(VENDORED_SRC))

from multiscale import fit_var_auto, run_grid_parallel  # noqa: E402
from multiscale.analysis import similarity_matrix  # noqa: E402


STAGES = ("W_", "N1", "N2", "N3", "R_")
STAGE_LABELS = {"W_": "Wake", "N1": "N1", "N2": "N2", "N3": "N3", "R_": "REM"}
BANDS = {"delta": (1, 4), "theta": (5, 8), "alpha": (9, 12),
         "sigma": (13, 16), "beta": (17, 30)}
SOURCE_PATTERN = re.compile(r"^(EPCTL\d{2})_(.+)_source_time_series_68\.mat$")


def stage_from_token(token: str) -> str | None:
    token = token.upper()
    if token.startswith("W"):
        return "W_"
    if token.startswith("R"):
        return "R_"
    return token if token in ("N1", "N2", "N3") else None


def discover(source_root: Path) -> list[tuple[str, str, Path]]:
    """Find and uniquely identify subject-stage source files."""
    jobs = {}
    for path in source_root.glob("EPCTL*/mat/*_source_time_series_68.mat"):
        match = SOURCE_PATTERN.fullmatch(path.name)
        if match is None:
            continue
        subject, token = match.groups()
        stage = stage_from_token(token)
        if stage is None:
            continue
        key = (subject, stage)
        if key in jobs:
            raise ValueError(f"duplicate source file for {subject} {stage}")
        jobs[key] = path
    return [(subject, stage, jobs[(subject, stage)])
            for subject, stage in sorted(jobs, key=lambda k: (STAGES.index(k[1]), k[0]))]


def pooled_zscore(data: np.ndarray, epsilon: float = 1e-12) -> np.ndarray:
    """MATLAB model-input branch's pooled mean and sample standard deviation."""
    if data.ndim == 2:
        data = data[:, :, None]
    if data.ndim != 3 or not np.isfinite(data).all():
        raise ValueError("source data must be finite (regions, samples, epochs)")
    flat = data.reshape(data.shape[0], -1)
    mean = flat.mean(axis=1, keepdims=True)
    std = flat.std(axis=1, ddof=1, keepdims=True)
    std[~np.isfinite(std) | (std < epsilon)] = 1
    return (data - mean[:, None, :]) / std[:, None, :]


def read_source(path: Path, required_regions: int) -> tuple[np.ndarray, list[str]]:
    mat = loadmat(path, simplify_cells=True)
    key = "source_ts" if "source_ts" in mat else "data"
    if key not in mat:
        raise KeyError(f"{path}: expected source_ts or data")
    data = np.asarray(mat[key], dtype=float)
    if data.ndim == 2:
        data = data[:, :, None]
    if data.ndim != 3 or data.shape[0] != required_regions:
        raise ValueError(f"{path}: expected {required_regions} regions, got {data.shape}")
    labels = [str(label).strip() for label in np.ravel(mat.get("labels", []))]
    if len(labels) != required_regions:
        labels = [f"DK{i:02d}" for i in range(1, required_regions + 1)]
    return pooled_zscore(data), labels


def atomic_savemat(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".mat", delete=False) as tmp:
        temp_path = Path(tmp.name)
    try:
        savemat(temp_path, {"results": result}, do_compression=True)
        os.replace(temp_path, path)
    finally:
        temp_path.unlink(missing_ok=True)


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, suffix=".json",
                                     encoding="utf-8", delete=False) as tmp:
        json.dump(payload, tmp, indent=2)
        temp_path = Path(tmp.name)
    try:
        os.replace(temp_path, path)
    finally:
        temp_path.unlink(missing_ok=True)


def atomic_npz(path: Path, **arrays) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".npz", delete=False) as tmp:
        temp_path = Path(tmp.name)
    try:
        np.savez_compressed(temp_path, **arrays)
        os.replace(temp_path, path)
    finally:
        temp_path.unlink(missing_ok=True)


def checkpoint_directory(subject: str, stage: str, source: Path,
                         args: argparse.Namespace, bands: dict) -> Path:
    """Separate checkpoints when data, settings, or the engine change."""
    stat = source.stat()
    upstream = (REPO / "code/external/emergence-multiscale/UPSTREAM.md").read_text()
    settings = {
        "source": str(source.resolve()), "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns, "engine": upstream,
        "regions": args.required_regions, "fs": args.fs, "fres": args.fres,
        "pmax": args.pmax, "criterion": args.criterion,
        "regmode": args.regmode, "dimensions": args.dimensions,
        "bands": bands, "runs": args.runs, "pre_iters": args.pre_iters,
        "final_iters": args.final_iters, "seed": args.seed,
    }
    digest = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()[:16]
    return args.output_root / "checkpoints" / f"{subject}_{stage}" / digest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    data_root = Path(os.environ.get("SOMNODYNAMICS_DATA_ROOT", REPO)).expanduser()
    parser.add_argument("--source-root", type=Path, default=data_root /
                        "results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk68")
    parser.add_argument("--output-root", type=Path, default=data_root /
                        "results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz")
    parser.add_argument("--task-index", type=int, help="Zero-based index in the sorted subject-stage manifest")
    parser.add_argument("--expected-count", type=int, help="Require this many source files before running")
    parser.add_argument("--subject", help="Run only one EPCTL subject")
    parser.add_argument("--stage", choices=STAGES, help="Run only one sleep stage")
    parser.add_argument("--max-files", type=int, help="Limit selected files for a trial run")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--required-regions", type=int, default=68)
    parser.add_argument("--fs", type=float, default=256)
    parser.add_argument("--fres", type=int, default=256)
    parser.add_argument("--pmax", type=int, default=30)
    parser.add_argument("--criterion", choices=("AIC", "BIC", "HQC", "LRT"), default="LRT")
    parser.add_argument("--regmode", choices=("LWR", "OLS"), default="LWR")
    parser.add_argument("--dimensions", type=int, nargs="+", default=list(range(1, 11)))
    parser.add_argument("--band", action="append", metavar="NAME:LOW:HIGH",
                        help="Override default bands; repeat for each band")
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--pre-iters", type=int, default=10000)
    parser.add_argument("--final-iters", type=int, default=10000)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def selected_bands(specs: list[str] | None) -> dict[str, tuple[float, float]]:
    if not specs:
        return BANDS.copy()
    bands = {}
    for spec in specs:
        name, low, high = spec.split(":")
        bounds = (float(low), float(high))
        if (not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in bands or
                not 0 <= bounds[0] <= bounds[1]):
            raise ValueError(f"invalid or duplicate band: {spec}")
        bands[name] = bounds
    return bands


def run_one(subject: str, stage: str, source: Path, args: argparse.Namespace,
            bands: dict[str, tuple[float, float]]) -> Path:
    destination = args.output_root / stage / f"{subject}_ssdi_spectral_domain_results.mat"
    status_path = args.output_root / "logs/status" / f"{subject}_{stage}.json"
    if destination.exists() and not args.overwrite:
        print(f"[{subject} {stage}] already complete: {destination}", flush=True)
        return destination
    started = time.monotonic()
    try:
        checkpoint_dir = checkpoint_directory(subject, stage, source, args, bands)
        if args.overwrite:
            shutil.rmtree(checkpoint_dir, ignore_errors=True)
            destination.unlink(missing_ok=True)
        print(f"[{subject} {stage}] loading {source}", flush=True)
        signals, labels = read_source(source, args.required_regions)
        if max(args.dimensions) >= signals.shape[0] or min(args.dimensions) < 1:
            raise ValueError("dimensions must be from 1 to regions-1")
        model_cache = checkpoint_dir / "model.npz"
        if model_cache.exists():
            with np.load(model_cache, allow_pickle=False) as saved:
                a, v, order = saved["A"], saved["V"], int(saved["order"])
            print(f"[{subject} {stage}] reusing VAR order={order}", flush=True)
        else:
            print(f"[{subject} {stage}] selecting VAR order (pmax={args.pmax}, "
                  f"{args.criterion}/{args.regmode})", flush=True)
            selected = fit_var_auto(signals, pmax=args.pmax, criterion=args.criterion,
                                    regmode=args.regmode)
            order = selected.order
            if order < 1:
                raise ValueError("selected VAR order zero; DD requires dynamics")
            a, v = selected.fit.coefficients, selected.fit.covariance
            atomic_npz(model_cache, A=a, V=v, order=np.asarray(order))
        pairs = [(m, band) for m in args.dimensions for band in bands]
        paths = {(m, band): checkpoint_dir / f"{band}_m{m}.npz"
                 for m, band in pairs}
        completed = {pair for pair, path in paths.items() if path.exists()}
        pending = [pair for pair in pairs if pair not in completed]
        print(f"[{subject} {stage}] order={order}; "
              f"{len(completed)}/{len(pairs)} DD blocks already complete", flush=True)

        def progress(done, total, m, band, block):
            atomic_npz(paths[(m, band)],
                       values=block["final"].values,
                       bases=block["final_original_bases"],
                       pre_values=block["pre"].values,
                       representatives=block["representatives"],
                       cluster_sizes=block["cluster_sizes"])
            print(f"[{subject} {stage}] DD {len(completed)+done}/{len(pairs)} "
                  f"m={m} band={band} "
                  f"best={float(np.min(block['final'].values)):.8g}", flush=True)

        if pending:
            run_grid_parallel(
                a, v, args.fs, args.fres, args.dimensions, bands,
                runs=args.runs,
                seed=args.seed + sum(map(ord, subject + stage)), workers=args.workers,
                pre_iters=args.pre_iters, final_iters=args.final_iters,
                pairs=pending, progress=progress,
            )
        result = {
            "meta": {
                "subject_id": subject, "stage": stage,
                "stage_label": STAGE_LABELS[stage],
                "labels": np.asarray(labels, dtype=object),
                "n_regions": signals.shape[0], "morder": order,
                "fs": args.fs, "fres": args.fres, "input_file": str(source),
                "model_source": f"python_source_recon_var_{args.criterion.lower()}_{args.regmode.lower()}",
                "nrunsp": args.runs, "nitersp": args.pre_iters,
                "niterso": args.final_iters,
                "engine": "emergence-multiscale",
                "timestamp": datetime.now(timezone.utc).isoformat(),
            },
            "model": {"A": a, "V": v, "morder": order},
        }
        for m, band in pairs:
            with np.load(paths[(m, band)], allow_pickle=False) as saved:
                values = saved["values"]
                bases = saved["bases"]
                pre_values = saved["pre_values"]
                representatives = saved["representatives"]
                cluster_sizes = saved["cluster_sizes"]
            if not np.isfinite(values).all():
                raise ValueError(f"non-finite DD values for {band} m={m}")
            result[f"{band}_dopto_m{m}"] = values
            result[f"{band}_Lopto_m{m}"] = bases
            result[f"{band}_gopto_m{m}"] = similarity_matrix(bases)
            result[f"{band}_doptp_m{m}"] = pre_values
            result[f"{band}_preopt_m{m}"] = {
                "uidx": representatives + 1,
                "usiz": cluster_sizes,
                "nrunso": len(representatives),
            }
        atomic_savemat(destination, result)
        atomic_json(status_path, {
            "subject": subject, "stage": stage, "success": True,
            "source_file": str(source), "output_file": str(destination),
            "model_order": order, "blocks": len(pairs),
            "seconds": time.monotonic() - started,
        })
        print(f"[{subject} {stage}] finished in {time.monotonic() - started:.1f}s",
              flush=True)
        return destination
    except Exception as exc:
        atomic_json(status_path, {
            "subject": subject, "stage": stage, "success": False,
            "source_file": str(source), "error": repr(exc),
            "seconds": time.monotonic() - started,
        })
        raise


def main() -> int:
    args = parse_args()
    jobs = discover(args.source_root)
    if args.expected_count is not None and len(jobs) != args.expected_count:
        raise SystemExit(f"expected {args.expected_count} DK68 files; found {len(jobs)} "
                         f"in {args.source_root}")
    if not jobs:
        raise SystemExit(f"no DK68 source files in {args.source_root}")
    if args.task_index is not None:
        if not 0 <= args.task_index < len(jobs):
            raise SystemExit(f"task index must be 0..{len(jobs)-1}")
        jobs = [jobs[args.task_index]]
    if args.subject:
        jobs = [job for job in jobs if job[0] == args.subject]
    if args.stage:
        jobs = [job for job in jobs if job[1] == args.stage]
    if args.max_files is not None:
        jobs = jobs[:args.max_files]
    if not jobs:
        raise SystemExit("no jobs match the requested selection")
    bands = selected_bands(args.band)
    print(f"DK68 Python DD: {len(jobs)} selected files, {len(bands)} bands, "
          f"dimensions={args.dimensions}, workers={args.workers}", flush=True)
    for index, (subject, stage, source) in enumerate(jobs, 1):
        print(f"file {index}/{len(jobs)}: {subject} {stage}", flush=True)
        run_one(subject, stage, source, args, bands)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
