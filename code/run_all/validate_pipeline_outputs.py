#!/usr/bin/env python3
"""Validate core SomnoDynamics outputs and report missing artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_utils import load_config, resolve_repo_path  # noqa: E402


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Validate SomnoDynamics outputs")
    parser.add_argument("--strict", action="store_true", help="Exit non-zero when required artifacts are missing")
    parser.add_argument("--report", type=Path, default=resolve_repo_path("results/validation_report.csv"))
    parser.add_argument("--summary", type=Path, default=resolve_repo_path("results/validation_summary.txt"))
    return parser.parse_args()


def checks() -> list[tuple[str, Path, bool]]:
    req = [
        ("migration_plan", resolve_repo_path("MIGRATION_PLAN_sleep_to_dmt_style.md"), True),
        ("mapping_csv", resolve_repo_path("code/coarse_graining/anatomical7_mapping_from_dmt.csv"), True),
        ("coarse_mapping_summary", resolve_repo_path("results/coarse_grained/anatomical7_mapping_summary.txt"), True),
        ("sensor_power_long", resolve_repo_path("results/sensor_power/sensor_bandpower_long.csv"), False),
        ("gc_omnibus_long", resolve_repo_path("results/gc/stage_mvgc2_anatomical7/stats/omnibus_long.csv"), False),
        ("gc_pairwise_long", resolve_repo_path("results/gc/stage_mvgc2_anatomical7/stats/pairwise_long.csv"), False),
        ("ssdi_long", resolve_repo_path("results/ssdi/ssdi_subject_long.csv"), False),
        ("cross_corr", resolve_repo_path("results/cross_analysis/gc_ssdi_correlations.csv"), False),
        ("fig_manifest", resolve_repo_path("results/figures/figure_manifest.csv"), False),
    ]
    return req


def main() -> int:
    args = parse_args()

    rows = []
    missing_required = 0
    for check_id, path, required in checks():
        exists = path.exists()
        rows.append(
            {
                "check_id": check_id,
                "path": str(path),
                "required": required,
                "exists": exists,
            }
        )
        if required and not exists:
            missing_required += 1

    df = pd.DataFrame(rows)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.report, index=False)

    with args.summary.open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics Validation Summary\n")
        handle.write("==============================\n")
        handle.write(f"Total checks: {len(df)}\n")
        handle.write(f"Required missing: {missing_required}\n")
        handle.write("\nMissing artifacts:\n")
        missing = df[~df["exists"]]
        if missing.empty:
            handle.write("- None\n")
        else:
            for _, row in missing.iterrows():
                req = "required" if row["required"] else "optional"
                handle.write(f"- [{req}] {row['check_id']}: {row['path']}\n")

    print(f"[SomnoDynamics] validation report -> {args.report}")
    print(f"[SomnoDynamics] validation summary -> {args.summary}")

    if args.strict and missing_required > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
