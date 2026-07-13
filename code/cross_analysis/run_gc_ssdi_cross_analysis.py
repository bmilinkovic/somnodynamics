#!/usr/bin/env python3
"""Cross-analysis between GC and SSDI stage-resolved outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import load_config, resolve_repo_path, stage_display, stage_order  # noqa: E402


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Run SomnoDynamics GC-SSDI cross analysis")
    parser.add_argument(
        "--gc-long",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["stats_root"]) / "gc_stage_node_values.csv",
    )
    parser.add_argument(
        "--ssdi-long",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["ssdi_root"]) / "ssdi_subject_long.csv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["cross_analysis_root"]),
    )
    return parser.parse_args()


def build_summary_tables(gc_long: pd.DataFrame, ssdi_long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    gc_summary = (
        gc_long.groupby(["subject_id", "stage_internal", "metric"], as_index=False)["value"]
        .mean()
        .pivot_table(index=["subject_id", "stage_internal"], columns="metric", values="value")
        .reset_index()
    )

    ssdi_summary = (
        ssdi_long.groupby(["subject_id", "stage_internal", "band", "m"], as_index=False)["dd_median"]
        .mean()
    )
    return gc_summary, ssdi_summary


def run_correlations(gc_summary: pd.DataFrame, ssdi_summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    merged = ssdi_summary.merge(gc_summary, on=["subject_id", "stage_internal"], how="inner")

    for (band, m), sub in merged.groupby(["band", "m"]):
        for metric in ["gc_in_time", "gc_out_time"]:
            if metric not in sub.columns:
                continue
            x = sub[metric].to_numpy(dtype=float)
            y = sub["dd_median"].to_numpy(dtype=float)
            mask = np.isfinite(x) & np.isfinite(y)
            if np.sum(mask) < 5:
                continue
            rho, p = spearmanr(x[mask], y[mask])
            rows.append(
                {
                    "band": band,
                    "m": int(m),
                    "gc_metric": metric,
                    "n": int(np.sum(mask)),
                    "spearman_rho": float(rho),
                    "spearman_p": float(p),
                }
            )

    return pd.DataFrame(rows)


def make_figures(gc_summary: pd.DataFrame, ssdi_summary: pd.DataFrame, corr_df: pd.DataFrame, out_dir: Path) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    merged = ssdi_summary.merge(gc_summary, on=["subject_id", "stage_internal"], how="inner")

    if corr_df.empty:
        return
    top = corr_df.sort_values("spearman_p").head(4)

    for _, row in top.iterrows():
        band = row["band"]
        m = int(row["m"])
        metric = row["gc_metric"]

        sub = merged[(merged["band"] == band) & (merged["m"] == m)]
        if metric not in sub.columns:
            continue

        fig, ax = plt.subplots(figsize=(6.5, 5.2))
        for stage in stage_order():
            pts = sub[sub["stage_internal"] == stage]
            if pts.empty:
                continue
            ax.scatter(pts[metric], pts["dd_median"], label=stage_display(stage), alpha=0.8)

        ax.set_xlabel(metric)
        ax.set_ylabel("SSDI dd_median")
        ax.set_title(f"GC vs SSDI | {band} m={m}\nrho={row['spearman_rho']:.3f}, p={row['spearman_p']:.3g}")
        ax.legend(frameon=False)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        stem = f"gc_ssdi_{metric}_{band}_m{m}"
        fig.savefig(fig_dir / f"{stem}.png", dpi=300)
        fig.savefig(fig_dir / f"{stem}.svg")
        plt.close(fig)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    gc_long = pd.read_csv(args.gc_long)
    ssdi_long = pd.read_csv(args.ssdi_long)

    gc_summary, ssdi_summary = build_summary_tables(gc_long, ssdi_long)
    gc_summary.to_csv(args.out_dir / "gc_summary_subject_stage.csv", index=False)
    ssdi_summary.to_csv(args.out_dir / "ssdi_summary_subject_stage_band_m.csv", index=False)

    corr_df = run_correlations(gc_summary, ssdi_summary)
    corr_df.to_csv(args.out_dir / "gc_ssdi_correlations.csv", index=False)

    make_figures(gc_summary, ssdi_summary, corr_df, args.out_dir)

    with (args.out_dir / "cross_analysis_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics GC-SSDI Cross Analysis\n")
        handle.write("===================================\n")
        handle.write(f"GC long: {args.gc_long}\n")
        handle.write(f"SSDI long: {args.ssdi_long}\n")
        handle.write(f"Correlation rows: {len(corr_df)}\n")
        handle.write("Method: Spearman correlation at subject-stage level.\n")

    print(f"[SomnoDynamics] cross-analysis complete -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
