#!/usr/bin/env python3
"""Collect and normalize publication-facing figures for SomnoDynamics."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import load_config, resolve_repo_path  # noqa: E402


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Build consolidated figure set")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["figures_root"]),
    )
    return parser.parse_args()


def stage_sources() -> dict[str, list[Path]]:
    return {
        "sensor_power": [resolve_repo_path("results/sensor_power/figures")],
        "gc_stats": [
            resolve_repo_path("results/gc/stage_mvgc2_anatomical7/figures"),
            resolve_repo_path("results/stats/figures_gc"),
        ],
        "ssdi": [resolve_repo_path("results/ssdi/figures")],
        "cross_analysis": [resolve_repo_path("results/cross_analysis/figures")],
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, str]] = []

    for stage, dirs in stage_sources().items():
        stage_out = args.out_dir / stage
        stage_out.mkdir(parents=True, exist_ok=True)

        for src_dir in dirs:
            if not src_dir.exists():
                continue
            for fig in sorted(src_dir.glob("*")):
                if not fig.is_file():
                    continue
                if fig.suffix.lower() not in {".png", ".svg", ".pdf"}:
                    continue

                dst = stage_out / fig.name
                shutil.copy2(fig, dst)
                manifest_rows.append(
                    {
                        "stage": stage,
                        "source": str(fig),
                        "copied_to": str(dst),
                        "format": fig.suffix.lower().lstrip("."),
                    }
                )

    manifest = pd.DataFrame(manifest_rows)
    manifest_csv = args.out_dir / "figure_manifest.csv"
    manifest.to_csv(manifest_csv, index=False)

    summary = args.out_dir / "figure_summary.txt"
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics Figure Build Summary\n")
        handle.write("===============================\n")
        handle.write(f"Output root: {args.out_dir}\n")
        handle.write(f"Total figures copied: {len(manifest_rows)}\n")
        by_stage = manifest.groupby("stage").size().to_dict() if not manifest.empty else {}
        for stage, n in by_stage.items():
            handle.write(f"- {stage}: {n}\n")

    print(f"[SomnoDynamics] consolidated figures -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
