#!/usr/bin/env python3
"""Repeated-measures statistics for GC outputs across sleep stages."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import friedmanchisquare, wilcoxon

try:
    from statsmodels.stats.multitest import multipletests
except Exception:  # pragma: no cover
    multipletests = None

import matplotlib.pyplot as plt

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import stage_display, stage_order  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="GC repeated-measures statistics")
    parser.add_argument("--gc-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def mat_struct_to_dict(obj: Any) -> Any:
    from scipy.io.matlab.mio5_params import mat_struct  # type: ignore

    if isinstance(obj, mat_struct):
        out = {}
        for name in obj._fieldnames:
            out[name] = mat_struct_to_dict(getattr(obj, name))
        return out
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        return [mat_struct_to_dict(v) for v in obj.tolist()]
    return obj


def load_gc_struct(path: Path, var_name: str) -> dict[str, Any]:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if var_name not in mat:
        raise KeyError(f"{var_name} not found in {path}")
    data = mat_struct_to_dict(mat[var_name])
    if not isinstance(data, dict):
        raise ValueError(f"Expected dict-like MATLAB struct for {var_name}")
    return data


def to_node_array(v: Any) -> np.ndarray:
    arr = np.asarray(v, dtype=float).squeeze()
    if arr.ndim == 0:
        arr = np.array([float(arr)])
    return arr.astype(float)


def build_time_long(gc_in: dict[str, Any], gc_out: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    for subject, in_struct in gc_in.items():
        if subject not in gc_out:
            continue
        out_struct = gc_out[subject]
        if not isinstance(in_struct, dict) or not isinstance(out_struct, dict):
            continue

        for stage in stage_order():
            if stage not in in_struct or stage not in out_struct:
                continue

            in_arr = to_node_array(in_struct[stage])
            out_arr = to_node_array(out_struct[stage])
            n_nodes = min(in_arr.size, out_arr.size)

            for node_idx in range(n_nodes):
                rows.append(
                    {
                        "subject_id": subject,
                        "stage_internal": stage,
                        "stage_label": stage_display(stage),
                        "metric": "gc_in_time",
                        "node": node_idx + 1,
                        "value": float(in_arr[node_idx]),
                    }
                )
                rows.append(
                    {
                        "subject_id": subject,
                        "stage_internal": stage,
                        "stage_label": stage_display(stage),
                        "metric": "gc_out_time",
                        "node": node_idx + 1,
                        "value": float(out_arr[node_idx]),
                    }
                )

    return pd.DataFrame(rows)


def repeated_stats(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    f_rows: list[dict[str, Any]] = []
    p_rows: list[dict[str, Any]] = []

    for metric in sorted(df["metric"].unique()):
        mdf = df[df["metric"] == metric]
        metric_pair_idx: list[int] = []

        for node in sorted(mdf["node"].unique()):
            sub = mdf[mdf["node"] == node]
            pivot = (
                sub.pivot_table(index="subject_id", columns="stage_internal", values="value", aggfunc="mean")
                .reindex(columns=stage_order())
                .dropna()
            )
            if pivot.shape[0] < 3:
                continue

            values = [pivot[s].to_numpy() for s in stage_order()]
            stat, p = friedmanchisquare(*values)
            f_rows.append(
                {
                    "metric": metric,
                    "node": int(node),
                    "n_subjects": int(pivot.shape[0]),
                    "friedman_stat": float(stat),
                    "friedman_p": float(p),
                }
            )

            for a, b in itertools.combinations(stage_order(), 2):
                try:
                    w_stat, w_p = wilcoxon(pivot[a], pivot[b], zero_method="wilcox")
                except ValueError:
                    w_stat, w_p = np.nan, np.nan
                p_rows.append(
                    {
                        "metric": metric,
                        "node": int(node),
                        "stage_a": a,
                        "stage_b": b,
                        "n_subjects": int(pivot.shape[0]),
                        "wilcoxon_stat": float(w_stat) if np.isfinite(w_stat) else np.nan,
                        "wilcoxon_p": float(w_p) if np.isfinite(w_p) else np.nan,
                        "mean_diff_b_minus_a": float(np.nanmean(pivot[b] - pivot[a])),
                    }
                )
                metric_pair_idx.append(len(p_rows) - 1)

        if multipletests is not None and metric_pair_idx:
            pvals = np.array([p_rows[i]["wilcoxon_p"] for i in metric_pair_idx], dtype=float)
            valid = np.isfinite(pvals)
            qvals = np.full_like(pvals, np.nan, dtype=float)
            reject = np.zeros_like(pvals, dtype=bool)
            if np.any(valid):
                reject_v, qvals_v, _, _ = multipletests(pvals[valid], alpha=0.05, method="fdr_bh")
                qvals[valid] = qvals_v
                reject[valid] = reject_v
            for local_i, row_i in enumerate(metric_pair_idx):
                p_rows[row_i]["wilcoxon_q_bh"] = float(qvals[local_i]) if np.isfinite(qvals[local_i]) else np.nan
                p_rows[row_i]["sig_q_bh_0p05"] = bool(reject[local_i])
        else:
            for row_i in metric_pair_idx:
                p_rows[row_i]["wilcoxon_q_bh"] = np.nan
                p_rows[row_i]["sig_q_bh_0p05"] = False

    return pd.DataFrame(f_rows), pd.DataFrame(p_rows)


def make_figures(df: pd.DataFrame, out_dir: Path) -> None:
    fig_dir = out_dir / "figures_gc"
    fig_dir.mkdir(parents=True, exist_ok=True)

    for metric in sorted(df["metric"].unique()):
        mdf = df[df["metric"] == metric]
        summary = (
            mdf.groupby(["node", "stage_internal"], as_index=False)["value"]
            .mean()
            .pivot(index="node", columns="stage_internal", values="value")
            .reindex(columns=stage_order())
        )

        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.imshow(summary.to_numpy(), aspect="auto", cmap="viridis")
        ax.set_xticks(np.arange(len(stage_order())), [stage_display(s) for s in stage_order()])
        ax.set_yticks(np.arange(len(summary.index)), [str(int(n)) for n in summary.index])
        ax.set_xlabel("Sleep Stage")
        ax.set_ylabel("Node")
        ax.set_title(f"{metric}: mean across subjects")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("GC value")
        fig.tight_layout()
        fig.savefig(fig_dir / f"{metric}_stage_node_heatmap.png", dpi=300)
        fig.savefig(fig_dir / f"{metric}_stage_node_heatmap.svg")
        plt.close(fig)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    gc_in = load_gc_struct(args.gc_dir / "gcin_time_subject_data_struct.mat", "gc_in_time_data")
    gc_out = load_gc_struct(args.gc_dir / "gcout_time_subject_data_struct.mat", "gc_out_time_data")

    long_df = build_time_long(gc_in, gc_out)
    long_csv = args.out_dir / "gc_stage_node_values.csv"
    long_df.to_csv(long_csv, index=False)

    f_df, p_df = repeated_stats(long_df)
    f_df.to_csv(args.out_dir / "gc_friedman_stats.csv", index=False)
    p_df.to_csv(args.out_dir / "gc_pairwise_wilcoxon_stats.csv", index=False)

    make_figures(long_df, args.out_dir)

    with (args.out_dir / "gc_stats_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics GC Repeated-Measures Statistics\n")
        handle.write("==========================================\n")
        handle.write(f"Input: {args.gc_dir}\n")
        handle.write(f"Long table rows: {len(long_df)}\n")
        handle.write("Tests: Friedman + pairwise Wilcoxon with BH-FDR by metric.\n")

    print(f"[SomnoDynamics] GC stats complete -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
