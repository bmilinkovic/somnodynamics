#!/usr/bin/env python3
"""Standalone DMT-style publication figures for DK5 stage-based GC outputs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "run_all"))
from config_utils import load_config, resolve_repo_path, stage_display  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    }
)

STAGE_BASE = "W_"
STAGE_ORDER = ["N1", "N2", "N3", "R_"]
BAND_ORDER = ["delta", "theta", "alpha", "sigma", "beta"]
BAND_SYMBOL = {
    "delta": "δ",
    "theta": "θ",
    "alpha": "α",
    "sigma": "σ",
    "beta": "β",
}
METRIC_ORDER = ["directed_src_in", "directed_src_out", "undirected_src"]
METRIC_TITLES = {
    "directed_src_in": "GC in",
    "directed_src_out": "GC out",
    "undirected_src": "Conditional MI",
}

FONT_TICK = 17
FONT_PANEL_TITLE = 20
FONT_SUPTITLE = 24
FONT_CBAR_LABEL = 20
FONT_CBAR_TICK = 17
FONT_CBAR_ENDPOINT = 17


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Build DMT-style GC figures for SomnoDynamics")
    parser.add_argument(
        "--stats-dir",
        type=Path,
        default=resolve_repo_path("results/gc/stage_mvgc2_dk5_lrt30_current_epoching_2s_256hz/stats"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=resolve_repo_path("results/gc/stage_mvgc2_dk5_lrt30_current_epoching_2s_256hz/figures"),
    )
    parser.add_argument(
        "--value-col",
        type=str,
        default="rank_corr",
        choices=["effect_mean", "effect_median", "rank_corr"],
        help="Pairwise statistic column to visualize.",
    )
    return parser.parse_args()


def cool_diverging(n: int = 256) -> np.ndarray:
    x = np.linspace(0, 1, n)
    anchors_x = np.array([0.00, 0.22, 0.50, 0.78, 1.00])
    anchors = np.array(
        [
            [0.06, 0.22, 0.41],
            [0.20, 0.56, 0.73],
            [0.95, 0.97, 0.98],
            [0.96, 0.67, 0.51],
            [0.71, 0.24, 0.24],
        ]
    )
    r = np.interp(x, anchors_x, anchors[:, 0])
    g = np.interp(x, anchors_x, anchors[:, 1])
    b = np.interp(x, anchors_x, anchors[:, 2])
    return np.column_stack([r, g, b])


def export_triplet(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")


def annotate_cbar_extremes(cbar: plt.colorbar, left_text: str, right_text: str) -> None:
    cbar.ax.text(0.0, -1.32, left_text, transform=cbar.ax.transAxes, ha="left", va="top", fontsize=FONT_CBAR_ENDPOINT)
    cbar.ax.text(1.0, -1.32, right_text, transform=cbar.ax.transAxes, ha="right", va="top", fontsize=FONT_CBAR_ENDPOINT)


def add_compact_cbar(
    fig: plt.Figure,
    mappable: plt.cm.ScalarMappable,
    label: str,
    left_text: str,
    right_text: str,
    *,
    y: float = 0.085,
    width: float = 0.56,
    height: float = 0.024,
) -> plt.colorbar:
    x = (1.0 - width) / 2.0
    cax = fig.add_axes([x, y, width, height])
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal")
    cb.set_label(label, fontsize=FONT_CBAR_LABEL, labelpad=9)
    cb.ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    cb.ax.tick_params(labelsize=FONT_CBAR_TICK)
    annotate_cbar_extremes(cb, left_text, right_text)
    return cb


def plot_roi_matrix(
    ax: plt.Axes,
    M: np.ndarray,
    Sig: np.ndarray,
    labels: list[str],
    title: str,
    vmax: float,
    show_ylabels: bool = True,
    show_xlabels: bool = True,
) -> None:
    Mplot = M.copy().astype(float)
    np.fill_diagonal(Mplot, 0.0)
    mask = np.eye(Mplot.shape[0], dtype=bool)
    m = np.ma.array(Mplot, mask=mask)

    cmap = plt.matplotlib.colors.ListedColormap(cool_diverging(256))
    cmap.set_bad("white")
    ax.imshow(m, cmap=cmap, vmin=-vmax, vmax=vmax, interpolation="nearest", origin="upper")
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels if show_xlabels else [], rotation=45, ha="right", fontsize=FONT_TICK)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels if show_ylabels else [], fontsize=FONT_TICK)
    ax.set_title(title, fontsize=FONT_PANEL_TITLE, weight="bold", pad=10)
    ax.tick_params(length=0)
    if show_xlabels:
        ax.set_xlabel("From region", fontsize=FONT_TICK, labelpad=9)
    if show_ylabels:
        ax.set_ylabel("To region", fontsize=FONT_TICK, labelpad=8)

    for k in np.arange(-0.5, len(labels), 1):
        ax.axhline(k, color=(0.86, 0.89, 0.93), lw=0.4)
        ax.axvline(k, color=(0.86, 0.89, 0.93), lw=0.4)

    rr, cc = np.where(Sig & ~np.eye(len(labels), dtype=bool))
    if len(rr) > 0:
        ax.scatter(cc, rr, s=18, facecolors="none", edgecolors="yellow", linewidths=0.8)

    for spine in ax.spines.values():
        spine.set_visible(False)


def matrix_from_long(
    df: pd.DataFrame,
    labels: list[str],
    value_col: str,
    sig_col: str,
) -> tuple[np.ndarray, np.ndarray]:
    mat = np.full((len(labels), len(labels)), np.nan, dtype=float)
    sig = np.zeros((len(labels), len(labels)), dtype=bool)
    li = {lab: i for i, lab in enumerate(labels)}
    for _, r in df.iterrows():
        t = str(r["target"])
        s = str(r["source"])
        if t not in li or s not in li:
            continue
        mat[li[t], li[s]] = float(r[value_col])
        sig[li[t], li[s]] = bool(r[sig_col])
    return mat, sig


def figure_time_roi_wake_contrasts(
    df: pd.DataFrame,
    labels: list[str],
    out_dir: Path,
    measure_family: str,
    stem_name: str,
    title_prefix: str,
    value_col: str,
    value_label: str,
) -> None:
    sub = df[(df["domain"] == "time") & (df["band"] == "all") & (df["measure_family"] == measure_family) & (df["stage_a"] == STAGE_BASE)]

    mats = []
    sigs = []
    for st in STAGE_ORDER:
        sdf = sub[sub["stage_b"] == st]
        M, S = matrix_from_long(sdf, labels, value_col, "sig_primary")
        mats.append(M)
        sigs.append(S)

    vmax = np.nanmax(np.abs(np.stack(mats)))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 12.1))
    fig.subplots_adjust(left=0.045, right=0.995, top=0.90, bottom=0.23, wspace=0.001, hspace=0.16)
    axes = axes.ravel()
    for i, st in enumerate(STAGE_ORDER):
        plot_roi_matrix(
            axes[i],
            mats[i],
            sigs[i],
            labels,
            f"{stage_display(st)} - Wake",
            vmax,
            show_ylabels=(i % 2 == 0),
            show_xlabels=(i >= 2),
        )

    sm = plt.cm.ScalarMappable(cmap=plt.matplotlib.colors.ListedColormap(cool_diverging(256)), norm=plt.Normalize(vmin=-vmax, vmax=vmax))
    sm.set_array([])
    _cb = add_compact_cbar(
        fig,
        sm,
        f"{value_label} (Stage - Wake)",
        "higher in Wake",
        "higher in stage",
        y=0.062,
        width=0.58,
        height=0.026,
    )
    fig.suptitle(title_prefix, fontsize=FONT_SUPTITLE, weight="bold", y=0.97)

    export_triplet(fig, out_dir / stem_name)
    plt.close(fig)


def figure_time_src_wake_contrasts(df: pd.DataFrame, labels: list[str], out_dir: Path, value_col: str, value_label: str) -> None:
    sub = df[(df["domain"] == "time") & (df["band"] == "all") & (df["stage_a"] == STAGE_BASE) & (df["measure_family"].isin(METRIC_ORDER))]

    fig, axes = plt.subplots(1, 3, figsize=(17.2, 8.2))
    fig.subplots_adjust(left=0.05, right=0.995, top=0.90, bottom=0.22, wspace=0.01)
    global_max = 0.0
    payload = []

    for metric in METRIC_ORDER:
        mdf = sub[sub["measure_family"] == metric]
        M = np.full((len(labels), len(STAGE_ORDER)), np.nan, dtype=float)
        Sig = np.zeros((len(labels), len(STAGE_ORDER)), dtype=bool)
        li = {lab: i for i, lab in enumerate(labels)}
        si = {s: i for i, s in enumerate(STAGE_ORDER)}

        for _, r in mdf.iterrows():
            t = str(r["target"])
            st = str(r["stage_b"])
            if t in li and st in si:
                M[li[t], si[st]] = float(r[value_col])
                Sig[li[t], si[st]] = bool(r["sig_primary"])

        payload.append((metric, M, Sig))
        if np.isfinite(M).any():
            global_max = max(global_max, float(np.nanmax(np.abs(M))))

    if global_max == 0:
        global_max = 1.0

    for ax, (metric, M, Sig) in zip(axes, payload):
        im = ax.imshow(M, aspect="auto", cmap=plt.matplotlib.colors.ListedColormap(cool_diverging(256)), vmin=-global_max, vmax=global_max)
        ax.set_xticks(range(len(STAGE_ORDER)))
        ax.set_xticklabels([stage_display(s) for s in STAGE_ORDER], fontsize=FONT_TICK)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels if metric == METRIC_ORDER[0] else [], fontsize=FONT_TICK)
        ax.set_title(METRIC_TITLES[metric], fontsize=FONT_PANEL_TITLE, weight="bold")
        rr, cc = np.where(Sig)
        if len(rr) > 0:
            ax.scatter(cc, rr, s=22, facecolors="none", edgecolors="yellow", linewidths=0.9)

    _cb = add_compact_cbar(
        fig,
        im,
        f"{value_label} (Stage - Wake)",
        "higher in Wake",
        "higher in stage",
        y=0.062,
        width=0.58,
        height=0.026,
    )
    fig.suptitle("Time-Domain Source-Level Wake Contrasts", fontsize=FONT_SUPTITLE, weight="bold")
    export_triplet(fig, out_dir / "time_src_wake_contrasts")
    plt.close(fig)


def figure_spectral_src_wake_contrasts(df: pd.DataFrame, labels: list[str], out_dir: Path, value_col: str, value_label: str) -> None:
    sub = df[(df["domain"] == "spectral") & (df["stage_a"] == STAGE_BASE) & (df["measure_family"].isin(METRIC_ORDER))]

    for st in STAGE_ORDER:
        sdf = sub[sub["stage_b"] == st]
        fig, axes = plt.subplots(1, 3, figsize=(17.2, 8.2))
        fig.subplots_adjust(left=0.05, right=0.995, top=0.90, bottom=0.22, wspace=0.01)

        global_max = 0.0
        payload = []
        for metric in METRIC_ORDER:
            mdf = sdf[sdf["measure_family"] == metric]
            M = np.full((len(labels), len(BAND_ORDER)), np.nan, dtype=float)
            Sig = np.zeros((len(labels), len(BAND_ORDER)), dtype=bool)
            li = {lab: i for i, lab in enumerate(labels)}
            bi = {b: i for i, b in enumerate(BAND_ORDER)}

            for _, r in mdf.iterrows():
                t = str(r["target"])
                b = str(r["band"])
                if t in li and b in bi:
                    M[li[t], bi[b]] = float(r[value_col])
                    Sig[li[t], bi[b]] = bool(r["sig_primary"])

            payload.append((metric, M, Sig))
            if np.isfinite(M).any():
                global_max = max(global_max, float(np.nanmax(np.abs(M))))

        if global_max == 0:
            global_max = 1.0

        for ax, (metric, M, Sig) in zip(axes, payload):
            im = ax.imshow(M, aspect="auto", cmap=plt.matplotlib.colors.ListedColormap(cool_diverging(256)), vmin=-global_max, vmax=global_max)
            ax.set_xticks(range(len(BAND_ORDER)))
            ax.set_xticklabels([BAND_SYMBOL[b] for b in BAND_ORDER], rotation=0, ha="center", fontsize=FONT_TICK + 1)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels if metric == METRIC_ORDER[0] else [], fontsize=FONT_TICK)
            ax.set_title(METRIC_TITLES[metric], fontsize=FONT_PANEL_TITLE, weight="bold", pad=10)
            rr, cc = np.where(Sig)
            if len(rr) > 0:
                ax.scatter(cc, rr, s=18, facecolors="none", edgecolors="yellow", linewidths=0.8)

        _cb = add_compact_cbar(
            fig,
            im,
            f"{value_label} (Stage - Wake)",
            "higher in Wake",
            f"higher in {stage_display(st)}",
            y=0.062,
            width=0.58,
            height=0.026,
        )
        fig.suptitle(
            f"Conditional Granger Causality / Conditional Mutual Information contrast: {stage_display(st)} - Wake",
            fontsize=FONT_SUPTITLE,
            weight="bold",
            y=0.98,
        )
        export_triplet(fig, out_dir / f"spectral_src_wake_contrast_{st}")
        plt.close(fig)


def figure_spectral_roi_wake_contrasts(
    df: pd.DataFrame,
    labels: list[str],
    out_dir: Path,
    measure_family: str,
    stem_prefix: str,
    value_col: str,
    value_label: str,
) -> None:
    sub = df[(df["domain"] == "spectral") & (df["stage_a"] == STAGE_BASE) & (df["measure_family"] == measure_family)]

    for st in STAGE_ORDER:
        sdf = sub[sub["stage_b"] == st]

        mats = []
        sigs = []
        for band in BAND_ORDER:
            bdf = sdf[sdf["band"] == band]
            M, S = matrix_from_long(bdf, labels, value_col, "sig_primary")
            mats.append(M)
            sigs.append(S)

        vmax = np.nanmax(np.abs(np.stack(mats)))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0

        fig, axes = plt.subplots(1, len(BAND_ORDER), figsize=(18.5, 6.4))
        fig.subplots_adjust(left=0.035, right=0.995, top=0.82, bottom=0.29, wspace=0.018)
        axes = axes.ravel()

        for i, band in enumerate(BAND_ORDER):
            plot_roi_matrix(
                axes[i],
                mats[i],
                sigs[i],
                labels,
                BAND_SYMBOL[band],
                vmax,
                show_ylabels=(i == 0),
                show_xlabels=True,
            )

        sm = plt.cm.ScalarMappable(cmap=plt.matplotlib.colors.ListedColormap(cool_diverging(256)), norm=plt.Normalize(vmin=-vmax, vmax=vmax))
        sm.set_array([])
        _cb = add_compact_cbar(
            fig,
            sm,
            f"{value_label} (Stage - Wake)",
            "higher in Wake",
            f"higher in {stage_display(st)}",
            y=0.062,
            width=0.58,
            height=0.026,
        )

        if measure_family == "undirected_roi":
            title_txt = f"Conditional Mutual Information contrast: {stage_display(st)} - Wake"
        else:
            title_txt = f"Pairwise Conditional Granger Causality contrast: {stage_display(st)} - Wake"
        fig.suptitle(title_txt, fontsize=FONT_SUPTITLE, weight="bold", y=0.97)
        export_triplet(fig, out_dir / f"{stem_prefix}_{st}")
        plt.close(fig)


def save_source_tables(pairwise_df: pd.DataFrame, out_dir: Path) -> None:
    src = out_dir / "source_data"
    src.mkdir(parents=True, exist_ok=True)

    for st in STAGE_ORDER:
        sdf = pairwise_df[(pairwise_df["stage_a"] == STAGE_BASE) & (pairwise_df["stage_b"] == st)]
        sdf.to_csv(src / f"wake_contrast_{st}_pairwise_long.csv", index=False)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    pairwise = pd.read_csv(args.stats_dir / "pairwise_long.csv")
    _omni = pd.read_csv(args.stats_dir / "omnibus_long.csv")
    if args.value_col not in pairwise.columns:
        raise ValueError(f"Requested --value-col '{args.value_col}' not found in pairwise_long.csv")
    value_label = {
        "effect_mean": "Mean difference",
        "effect_median": "Median difference",
        "rank_corr": "Rank correlation",
    }[args.value_col]

    dk5_order = ["Frontal", "Sensorimotor", "Parietal", "Temporal", "Occipital"]
    observed_labels = set(pairwise["target"].dropna().astype(str)) | set(pairwise["source"].dropna().astype(str))
    labels = [lab for lab in dk5_order if lab in observed_labels]
    labels.extend(sorted(observed_labels.difference(labels)))

    figure_time_roi_wake_contrasts(
        pairwise,
        labels,
        args.out_dir,
        measure_family="directed_roi",
        stem_name="time_directed_roi_wake_contrasts",
        title_prefix="Directed Time-Domain ROI Wake Contrasts",
        value_col=args.value_col,
        value_label=value_label,
    )
    figure_time_roi_wake_contrasts(
        pairwise,
        labels,
        args.out_dir,
        measure_family="undirected_roi",
        stem_name="time_undirected_roi_wake_contrasts",
        title_prefix="Conditional Mutual Information Time-Domain ROI Wake Contrasts",
        value_col=args.value_col,
        value_label=value_label,
    )
    figure_time_src_wake_contrasts(pairwise, labels, args.out_dir, args.value_col, value_label)
    figure_spectral_src_wake_contrasts(pairwise, labels, args.out_dir, args.value_col, value_label)
    figure_spectral_roi_wake_contrasts(
        pairwise,
        labels,
        args.out_dir,
        "directed_roi",
        "spectral_directed_roi_wake_contrast",
        args.value_col,
        value_label,
    )
    figure_spectral_roi_wake_contrasts(
        pairwise,
        labels,
        args.out_dir,
        "undirected_roi",
        "spectral_undirected_roi_wake_contrast",
        args.value_col,
        value_label,
    )
    save_source_tables(pairwise, args.out_dir)

    with (args.out_dir / "gc_figure_summary.txt").open("w", encoding="utf-8") as f:
        f.write("SomnoDynamics GC Figure Summary\n")
        f.write("=============================\n")
        f.write(f"Stats source: {args.stats_dir}\n")
        f.write(f"Output dir: {args.out_dir}\n")
        f.write("Style: DMT-like contrast heatmaps with significance overlays (yellow circles).\n")
        f.write("Contrasts shown: Stage - Wake for N1/N2/N3/REM.\n")
        f.write(f"Color values: {args.value_col}\n")
        if args.value_col == "rank_corr":
            f.write("Direction: rank_corr > 0 means Stage > Wake; rank_corr < 0 means Wake > Stage.\n")

    print(f"[SomnoDynamics] GC figures written -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
