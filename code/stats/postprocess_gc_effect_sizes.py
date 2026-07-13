#!/usr/bin/env python3
"""Post-process existing GC statistics with effect sizes and diagnostics."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Post-process GC stats: effect sizes + filtering + plots")
    parser.add_argument(
        "--stats-dir",
        type=Path,
        default=Path("results/gc/stage_mvgc2_anatomical7/stats"),
        help="Directory containing pairwise_long.csv and omnibus_long.csv",
    )
    parser.add_argument(
        "--filtered-dir",
        type=Path,
        default=Path("results/statistics_filtered"),
        help="Directory for filtered result tables and count summary",
    )
    parser.add_argument(
        "--effects-dir",
        type=Path,
        default=Path("results/statistics_effect_sizes"),
        help="Directory for effect size outputs and figures",
    )
    parser.add_argument(
        "--q-threshold",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--effect-threshold",
        type=float,
        default=0.3,
    )
    return parser.parse_args()


def export_plot(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def infer_stage_count(pairwise_df: pd.DataFrame) -> int:
    cols = [c for c in ("stage_a", "stage_b") if c in pairwise_df.columns]
    labels = set()
    for col in cols:
        labels.update(pairwise_df[col].dropna().astype(str).unique().tolist())
    return max(len(labels), 2)


def pairwise_effect_sizes(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["effect_size_method"] = "unknown"
    out["effect_size"] = np.nan

    if "rank_corr" in out.columns:
        out["effect_size"] = pd.to_numeric(out["rank_corr"], errors="coerce")
        out["effect_size_method"] = "rank_biserial_signed_rank"
        return out

    if {"W", "n_subjects_effective"}.issubset(out.columns):
        n = pd.to_numeric(out["n_subjects_effective"], errors="coerce").to_numpy(dtype=float)
        w = pd.to_numeric(out["W"], errors="coerce").to_numpy(dtype=float)
        rank_total = n * (n + 1.0) / 2.0
        out["effect_size"] = np.where(rank_total > 0, (2.0 * w / rank_total) - 1.0, np.nan)
        out["effect_size"] = np.clip(out["effect_size"], -1.0, 1.0)
        out["effect_size_method"] = "rank_biserial_signed_rank"
        return out

    if {"t", "n_subjects_effective"}.issubset(out.columns):
        t = pd.to_numeric(out["t"], errors="coerce").to_numpy(dtype=float)
        n = pd.to_numeric(out["n_subjects_effective"], errors="coerce").to_numpy(dtype=float)
        out["effect_size"] = np.where(n > 0, t / np.sqrt(n), np.nan)
        out["effect_size_method"] = "cohen_dz_from_t"
        return out

    if {"mean_diff", "sd_diff"}.issubset(out.columns):
        mean_diff = pd.to_numeric(out["mean_diff"], errors="coerce").to_numpy(dtype=float)
        sd_diff = pd.to_numeric(out["sd_diff"], errors="coerce").to_numpy(dtype=float)
        out["effect_size"] = np.where(sd_diff > 0, mean_diff / sd_diff, np.nan)
        out["effect_size_method"] = "cohen_dz"
        return out

    return out


def omnibus_effect_sizes(df: pd.DataFrame, n_stages: int) -> pd.DataFrame:
    out = df.copy()
    out["effect_size_method"] = "unknown"
    out["effect_size"] = np.nan

    if {"friedman_chi2", "n_subjects_effective"}.issubset(out.columns):
        chi2 = pd.to_numeric(out["friedman_chi2"], errors="coerce").to_numpy(dtype=float)
        n = pd.to_numeric(out["n_subjects_effective"], errors="coerce").to_numpy(dtype=float)
        denom = n * max(n_stages - 1, 1)
        kendall_w = np.where(denom > 0, chi2 / denom, np.nan)
        out["effect_size"] = np.clip(kendall_w, 0.0, 1.0)
        out["effect_size_method"] = "kendall_w_friedman"
        return out

    if {"F", "df_num", "df_den"}.issubset(out.columns):
        F = pd.to_numeric(out["F"], errors="coerce").to_numpy(dtype=float)
        df_num = pd.to_numeric(out["df_num"], errors="coerce").to_numpy(dtype=float)
        df_den = pd.to_numeric(out["df_den"], errors="coerce").to_numpy(dtype=float)
        out["effect_size"] = np.where((F >= 0) & (df_num > 0) & (df_den > 0), (F * df_num) / (F * df_num + df_den), np.nan)
        out["effect_size_method"] = "partial_eta_squared"
        return out

    return out


def filter_df(df: pd.DataFrame, q_threshold: float, effect_threshold: float) -> pd.DataFrame:
    q = pd.to_numeric(df["q_primary"], errors="coerce")
    e = pd.to_numeric(df["effect_size"], errors="coerce")
    return df[(q <= q_threshold) & (np.abs(e) > effect_threshold)].copy()


def summarize_counts(df: pd.DataFrame, q_threshold: float, effect_threshold: float) -> tuple[int, int, int]:
    q = pd.to_numeric(df["q_primary"], errors="coerce")
    e = pd.to_numeric(df["effect_size"], errors="coerce")
    total = int(len(df))
    sig = int((q <= q_threshold).sum())
    filt = int(((q <= q_threshold) & (np.abs(e) > effect_threshold)).sum())
    return total, sig, filt


def make_effect_plots(df_combined: pd.DataFrame, fig_dir: Path, q_threshold: float, effect_threshold: float) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    q = pd.to_numeric(df_combined["q_primary"], errors="coerce")
    e = pd.to_numeric(df_combined["effect_size"], errors="coerce")
    valid = np.isfinite(q) & np.isfinite(e)
    d = df_combined.loc[valid].copy()
    q = pd.to_numeric(d["q_primary"], errors="coerce")
    e = pd.to_numeric(d["effect_size"], errors="coerce")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(e, bins=60, color="#356f9f", edgecolor="white", alpha=0.9)
    ax.axvline(0.0, color="black", lw=1.0)
    ax.axvline(effect_threshold, color="#c0392b", lw=1.0, ls="--")
    ax.axvline(-effect_threshold, color="#c0392b", lw=1.0, ls="--")
    ax.set_xlabel("Effect size")
    ax.set_ylabel("Count")
    ax.set_title("GC effect-size distribution (all tests)")
    export_plot(fig, fig_dir / "effect_size_hist_all")

    sig_mask = q <= q_threshold
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(e[sig_mask], bins=50, color="#2a9d8f", edgecolor="white", alpha=0.9)
    ax.axvline(0.0, color="black", lw=1.0)
    ax.axvline(effect_threshold, color="#c0392b", lw=1.0, ls="--")
    ax.axvline(-effect_threshold, color="#c0392b", lw=1.0, ls="--")
    ax.set_xlabel("Effect size")
    ax.set_ylabel("Count")
    ax.set_title("GC effect-size distribution (q_primary <= 0.05)")
    export_plot(fig, fig_dir / "effect_size_hist_significant")

    y = -np.log10(np.clip(q.to_numpy(dtype=float), 1e-300, None))
    filtered_mask = (q <= q_threshold) & (np.abs(e) > effect_threshold)

    fig, ax = plt.subplots(figsize=(9.5, 6))
    pair_mask = d["test_scope"] == "pairwise"
    omni_mask = d["test_scope"] == "omnibus"
    ax.scatter(e[pair_mask], y[pair_mask], s=13, alpha=0.35, c="#6c757d", label="Pairwise")
    ax.scatter(e[omni_mask], y[omni_mask], s=16, alpha=0.45, c="#457b9d", label="Omnibus")
    ax.scatter(e[filtered_mask], y[filtered_mask], s=28, alpha=0.9, c="#d62828", label=f"Filtered (q<={q_threshold}, |ES|>{effect_threshold})")
    ax.axvline(effect_threshold, color="#c0392b", lw=1.0, ls="--")
    ax.axvline(-effect_threshold, color="#c0392b", lw=1.0, ls="--")
    ax.axhline(-np.log10(q_threshold), color="#34495e", lw=1.0, ls=":")
    ax.set_xlabel("Effect size")
    ax.set_ylabel("-log10(q_primary)")
    ax.set_title("GC effect size vs q-value")
    ax.legend(loc="upper left", frameon=False)
    export_plot(fig, fig_dir / "effect_size_vs_q_primary")


def build_summary_text(
    pair_total: int,
    pair_sig: int,
    pair_filt: int,
    omni_total: int,
    omni_sig: int,
    omni_filt: int,
    all_df: pd.DataFrame,
    q_threshold: float,
) -> str:
    total = pair_total + omni_total
    sig = pair_sig + omni_sig
    filt = pair_filt + omni_filt

    q = pd.to_numeric(all_df["q_primary"], errors="coerce")
    e = pd.to_numeric(all_df["effect_size"], errors="coerce")
    sig_es = np.abs(e[(q <= q_threshold) & np.isfinite(e)])
    n_sig = int(sig_es.shape[0])
    small_count = int((sig_es < 0.2).sum()) if n_sig > 0 else 0
    small_frac = (small_count / n_sig) if n_sig > 0 else np.nan

    lines = [
        "GC effect-size post-processing summary",
        "=====================================",
        "",
        "Effect-size definitions:",
        "- Pairwise (Wilcoxon signed-rank): rank correlation (rank-biserial) in [-1, 1]",
        "- Omnibus (Friedman): Kendall's W in [0, 1]",
        "",
        "Counts:",
        f"pairwise_total_tests: {pair_total}",
        f"pairwise_significant_tests_q<=0.05: {pair_sig}",
        f"pairwise_filtered_tests_q<=0.05_and_absES>0.3: {pair_filt}",
        f"omnibus_total_tests: {omni_total}",
        f"omnibus_significant_tests_q<=0.05: {omni_sig}",
        f"omnibus_filtered_tests_q<=0.05_and_absES>0.3: {omni_filt}",
        f"combined_total_tests: {total}",
        f"combined_significant_tests_q<=0.05: {sig}",
        f"combined_filtered_tests_q<=0.05_and_absES>0.3: {filt}",
        "",
    ]

    if n_sig > 0:
        lines.append(
            f"Among significant tests, {small_count}/{n_sig} ({small_frac:.1%}) have |effect size| < 0.2."
        )
    else:
        lines.append("No significant tests found at q_primary <= 0.05.")

    if n_sig > 0 and small_frac >= 0.5:
        lines.append("Interpretation: many statistically significant results have small effect sizes.")
    elif n_sig > 0:
        lines.append("Interpretation: a substantial subset of significant results are above small-effect range.")

    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    pairwise_path = args.stats_dir / "pairwise_long.csv"
    omnibus_path = args.stats_dir / "omnibus_long.csv"
    if not pairwise_path.exists() or not omnibus_path.exists():
        raise FileNotFoundError(f"Missing required stats files in {args.stats_dir}")

    pairwise = pd.read_csv(pairwise_path)
    omnibus = pd.read_csv(omnibus_path)
    n_stages = infer_stage_count(pairwise)

    pair_es = pairwise_effect_sizes(pairwise)
    omni_es = omnibus_effect_sizes(omnibus, n_stages=n_stages)
    pair_es["test_scope"] = "pairwise"
    omni_es["test_scope"] = "omnibus"

    args.filtered_dir.mkdir(parents=True, exist_ok=True)
    args.effects_dir.mkdir(parents=True, exist_ok=True)

    pair_es.to_csv(args.effects_dir / "gc_pairwise_effect_sizes.csv", index=False)
    omni_es.to_csv(args.effects_dir / "gc_omnibus_effect_sizes.csv", index=False)

    pair_filt = filter_df(pair_es, args.q_threshold, args.effect_threshold)
    omni_filt = filter_df(omni_es, args.q_threshold, args.effect_threshold)

    pair_filt.to_csv(args.filtered_dir / "gc_pairwise_filtered_q0p05_absES0p3.csv", index=False)
    omni_filt.to_csv(args.filtered_dir / "gc_omnibus_filtered_q0p05_absES0p3.csv", index=False)

    pair_total, pair_sig, pair_filt_n = summarize_counts(pair_es, args.q_threshold, args.effect_threshold)
    omni_total, omni_sig, omni_filt_n = summarize_counts(omni_es, args.q_threshold, args.effect_threshold)

    combined = pd.concat([pair_es, omni_es], axis=0, ignore_index=True)
    make_effect_plots(combined, args.effects_dir / "figures", args.q_threshold, args.effect_threshold)

    summary = build_summary_text(
        pair_total=pair_total,
        pair_sig=pair_sig,
        pair_filt=pair_filt_n,
        omni_total=omni_total,
        omni_sig=omni_sig,
        omni_filt=omni_filt_n,
        all_df=combined,
        q_threshold=args.q_threshold,
    )
    (args.filtered_dir / "gc_effect_size_filter_summary.txt").write_text(summary, encoding="utf-8")
    (args.effects_dir / "gc_effect_size_summary.txt").write_text(summary, encoding="utf-8")

    print(summary.strip())
    print(f"\nSaved effect-size tables -> {args.effects_dir}")
    print(f"Saved filtered tables -> {args.filtered_dir}")
    print(f"Saved figures -> {args.effects_dir / 'figures'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
