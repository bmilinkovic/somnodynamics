#!/usr/bin/env python3
"""Validate and summarize the DMT-derived canonical anatomical7 mapping."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import resolve_repo_path  # noqa: E402


MACRO_ORDER = [
    "Frontal",
    "Sensorimotor",
    "Parietal",
    "Temporal",
    "Occipital",
    "Cingulate",
    "Limbic",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate anatomical7 coarse-graining mapping")
    parser.add_argument(
        "--mapping-csv",
        type=Path,
        default=resolve_repo_path("code/coarse_graining/anatomical7_mapping_from_dmt.csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=resolve_repo_path("results/coarse_grained"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.mapping_csv)
    required = {"label", "macro_region", "rationale"}
    if not required.issubset(df.columns):
        raise ValueError(f"Mapping file missing required columns: {required}")

    unknown = sorted(set(df["macro_region"]) - set(MACRO_ORDER))
    if unknown:
        raise ValueError(f"Unexpected macro regions found: {unknown}")

    counts = (
        df.groupby("macro_region", as_index=False)
        .size()
        .rename(columns={"size": "n_labels"})
    )
    counts["macro_region"] = pd.Categorical(counts["macro_region"], categories=MACRO_ORDER, ordered=True)
    counts = counts.sort_values("macro_region")

    counts_csv = args.out_dir / "anatomical7_mapping_counts.csv"
    counts.to_csv(counts_csv, index=False)

    summary_path = args.out_dir / "anatomical7_mapping_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics Canonical Anatomical7 Mapping (from DMT)\n")
        handle.write("====================================================\n")
        handle.write(f"Source CSV: {args.mapping_csv}\n")
        handle.write(f"Total labels: {len(df)}\n")
        handle.write("Counts by macro region:\n")
        for _, row in counts.iterrows():
            handle.write(f"- {row['macro_region']}: {row['n_labels']}\n")

    print(f"[SomnoDynamics] mapping validated -> {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
