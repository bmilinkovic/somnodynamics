#!/usr/bin/env python3
"""SomnoDynamics stage runner.

This script provides stage-specific and full-pipeline execution entrypoints.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
from config_utils import resolve_repo_path  # noqa: E402


def run_cmd(cmd: list[str]) -> int:
    print("[SomnoDynamics]", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=REPO_ROOT)
    return int(proc.returncode)


def run_python(script_rel: str, args: list[str]) -> int:
    script = REPO_ROOT / script_rel
    return run_cmd([sys.executable, str(script), *args])


def run_matlab(function_call: str, matlab_bin: str) -> int:
    return run_cmd([matlab_bin, "-batch", function_call])


def matlab_quote(text: str) -> str:
    return text.replace("'", "''")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SomnoDynamics pipeline stages")
    parser.add_argument(
        "stage",
        choices=[
            "preprocessing",
            "sensor_power",
            "source_recon",
            "coarse_graining",
            "gc",
            "ssdi",
            "figures",
            "compare_legacy",
            "validate",
            "cross_analysis",
            "all",
        ],
    )
    parser.add_argument(
        "--matlab-bin",
        default="/Applications/MATLAB_R2023b.app/bin/matlab",
        help="MATLAB executable used for GC/SSDI stages.",
    )
    parser.add_argument(
        "--skip-matlab",
        action="store_true",
        help="Skip MATLAB-based stages even when selected.",
    )
    parser.add_argument(
        "--n-jobs-subjects",
        type=int,
        default=1,
        help="Parallel workers across subjects for compatible Python stages.",
    )
    parser.add_argument(
        "--existing-preprocessed-dir",
        type=Path,
        default=None,
        help="Optional existing preprocessed epochs directory for preprocessing import mode.",
    )
    parser.add_argument(
        "--gc-source-root",
        type=Path,
        default=resolve_repo_path("results/source_recon/anatomical7"),
        help="Source reconstruction directory used by GC stage.",
    )
    parser.add_argument(
        "--gc-output-root",
        type=Path,
        default=resolve_repo_path("results/gc/stage_mvgc2_anatomical7_256hz"),
        help="GC output directory (defaults to 256 Hz-labeled results root).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    plan = []
    if args.stage in {"preprocessing", "all"}:
        preproc_args = ["--n-jobs", str(args.n_jobs_subjects)]
        if args.existing_preprocessed_dir is not None:
            preproc_args.extend(["--existing-preprocessed-dir", str(args.existing_preprocessed_dir)])
        plan.append(("python", "code/preprocessing/preprocess_sensor_space.py", preproc_args))
    if args.stage in {"sensor_power", "all"}:
        plan.append(("python", "code/sensor_power/run_sensor_power_analysis.py", ["--n-jobs", str(args.n_jobs_subjects)]))
    if args.stage in {"source_recon", "all"}:
        plan.append(("python", "code/src_recon/run_source_reconstruction_anatomical46.py", []))
    if args.stage in {"coarse_graining", "all"}:
        plan.append(("python", "code/coarse_graining/validate_mapping.py", []))
    if args.stage in {"gc", "all"}:
        gc_source = matlab_quote(str(args.gc_source_root))
        gc_output = matlab_quote(str(args.gc_output_root))
        gc_call = (
            "addpath('code/gc'); "
            f"run_gc_sleep_stage_pipeline('mode','stage_mvgc2','source_root','{gc_source}','output_root','{gc_output}');"
        )
        plan.append(("matlab", gc_call, []))
    if args.stage in {"ssdi", "all"}:
        plan.append((
            "matlab",
            "addpath('code/ssdi'); run_ssdi_sleep_stage_pipeline('mode','stage_ssdi');",
            [],
        ))
        plan.append((
            "python",
            "code/ssdi/summarize_ssdi_results.py",
            [
                "--spectral-dir",
                "results/ssdi/stage_ssdi_anatomical7",
                "--out-dir",
                "results/ssdi",
            ],
        ))
    if args.stage in {"figures", "all"}:
        plan.append(("python", "code/gc/plot_gc_stage_figures.py", []))
        plan.append(("python", "code/ssdi/plot_ssdi_stage_figures.py", []))
        plan.append(("python", "code/figures/build_publication_figures.py", []))
    if args.stage in {"cross_analysis", "all"}:
        plan.append(("python", "code/cross_analysis/run_gc_ssdi_cross_analysis.py", []))
    if args.stage in {"compare_legacy", "all"}:
        plan.append(("python", "code/run_all/compare_with_legacy.py", []))
    if args.stage in {"validate", "all"}:
        plan.append(("python", "code/run_all/validate_pipeline_outputs.py", []))

    for kind, target, extra in plan:
        if kind == "python":
            rc = run_python(target, extra)
        else:
            if args.skip_matlab:
                print(f"[SomnoDynamics] skipping MATLAB stage: {target}")
                continue
            rc = run_matlab(target, args.matlab_bin)

        if rc != 0:
            print(f"[SomnoDynamics] stage failed ({target}) with exit code {rc}")
            return rc

    print("[SomnoDynamics] completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
