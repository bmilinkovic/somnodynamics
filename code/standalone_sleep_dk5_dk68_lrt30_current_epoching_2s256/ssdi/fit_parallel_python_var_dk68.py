#!/usr/bin/env python3
"""Fit DK68 VAR models with a small, RAM-aware pool of full multiscale fits.

Each child runs fit_python_var_dk68.py unchanged. Smaller recordings can fit
together; larger recordings run alone to leave sufficient memory. A completed
model is skipped, so the queue can be restarted after interruption.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

from run_python_dd_dk68 import discover


def alive(pid: int | None) -> bool:
    if pid is None:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    state = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)],
                           capture_output=True, text=True, check=False)
    return state.returncode == 0 and bool(state.stdout.strip()) and "Z" not in state.stdout


def weight(size: int) -> int:
    """Reserve more of the pool for source files that need more working RAM."""
    mb = size / 1_000_000
    return 1 if mb <= 100 else 2 if mb <= 200 else 3


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    data_root = Path(os.environ.get("SOMNODYNAMICS_DATA_ROOT",
                            "/Volumes/ex_data/data_somnoemergence"))
    parser.add_argument("--source-root", type=Path, default=data_root /
                        "results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk68")
    parser.add_argument("--model-root", type=Path, default=data_root /
                        "results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz/model_inputs")
    parser.add_argument("--log-root", type=Path, default=data_root /
                        "results/ssdi/stage_ssdi_dk68_lrt30_current_epoching_2s_256hz/logs/python_var_parallel")
    parser.add_argument("--max-workers", type=int, default=3)
    parser.add_argument("--external-index", type=int, default=0,
                        help="Task already running in a previous launcher")
    parser.add_argument("--external-pid", type=int,
                        help="PID of that existing fit; its memory occupies one slot")
    parser.add_argument("--expected-count", type=int, default=143)
    args = parser.parse_args()
    if not 1 <= args.max_workers <= 4:
        parser.error("--max-workers must be 1..4")

    jobs = discover(args.source_root)
    if len(jobs) != args.expected_count:
        raise SystemExit(f"expected {args.expected_count} source files, found {len(jobs)}")
    args.log_root.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name("fit_python_var_dk68.py")
    env = os.environ.copy()
    for name in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                 "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[name] = "1"

    pending = sorted(range(len(jobs)), key=lambda i: (jobs[i][2].stat().st_size, i))
    running: dict[int, tuple[subprocess.Popen[bytes], int]] = {}
    failures: list[int] = []
    while pending or running or alive(args.external_pid):
        for index, (process, _) in list(running.items()):
            code = process.poll()
            if code is None:
                continue
            del running[index]
            subject, stage, _ = jobs[index]
            print(f"{subject} {stage}: {'saved' if code == 0 else f'FAILED ({code})'}; "
                  f"{len(pending)} queued", flush=True)
            if code:
                failures.append(index)

        external = alive(args.external_pid)
        used = sum(reserved for _, reserved in running.values())
        if external:
            used += weight(jobs[args.external_index][2].stat().st_size)
        for index in list(pending):
            subject, stage, source = jobs[index]
            output = args.model_root / stage / f"{subject}_var_model.npz"
            if index == args.external_index and external:
                continue
            if output.exists():
                pending.remove(index)
                continue
            reserved = weight(source.stat().st_size)
            if used + reserved > args.max_workers:
                continue
            pending.remove(index)
            used += reserved
            log = args.log_root / f"task_{index:03d}_{subject}_{stage}.log"
            command = [sys.executable, "-u", str(script), "--task-index", str(index),
                       "--expected-count", str(args.expected_count),
                       "--source-root", str(args.source_root),
                       "--model-root", str(args.model_root)]
            with log.open("a") as stream:
                process = subprocess.Popen(command, stdout=stream,
                                           stderr=subprocess.STDOUT, env=env)
            running[index] = (process, reserved)
            print(f"Started {subject} {stage} ({source.stat().st_size/1e6:.0f} MB); "
                  f"PID {process.pid}; pool {used}/{args.max_workers}; "
                  f"{len(pending)} queued", flush=True)
        if pending or running or external:
            time.sleep(5)

    saved = sum((args.model_root / stage / f"{subject}_var_model.npz").is_file()
                for subject, stage, _ in jobs)
    print(f"Saved {saved}/{args.expected_count} DK68 models; "
          f"failed task indices: {failures}", flush=True)
    return 0 if saved == args.expected_count and not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
