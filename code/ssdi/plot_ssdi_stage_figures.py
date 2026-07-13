#!/usr/bin/env python3
"""DMT-parity SSDI figures for SomnoDynamics sleep-stage analysis."""

from __future__ import annotations

import argparse
import os
import re
from itertools import combinations
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
import seaborn as sns
from matplotlib import patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.patches import Polygon
from scipy.stats import norm
from scipy.stats import wilcoxon

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import load_config, resolve_repo_path, stage_display, stage_order  # noqa: E402

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

DEFAULT_BAND_ORDER = ["delta", "theta", "alpha", "sigma", "beta", "gamma"]
BAND_GREEK = {
    "delta": "δ",
    "theta": "θ",
    "alpha": "α",
    "sigma": "σ",
    "beta": "β",
    "gamma": "γ",
}
STAGE_ORDER = stage_order()
STAGE_PALETTE = {
    "W_": "#2EC4B6",
    "N1": "#F4A261",
    "N2": "#E76F51",
    "N3": "#7B2CBF",
    "R_": "#3A86FF",
}
WAKE_STAGE = "W_"


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Build DMT-parity SSDI stage figures")
    parser.add_argument(
        "--long-csv",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["ssdi_root"]) / "ssdi_subject_long.csv",
    )
    parser.add_argument(
        "--pairwise-csv",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["ssdi_root"]) / "ssdi_dd_median_pairwise_wilcoxon.csv",
        help="Legacy arg retained for compatibility; not required for DMT-parity figures.",
    )
    parser.add_argument(
        "--spectral-dir",
        type=Path,
        required=True,
        help="Stage SSDI results root containing stage subfolders with *_ssdi_spectral_domain_results.mat files.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["ssdi_root"]) / "figures",
    )
    parser.add_argument(
        "--pairwise-root",
        type=Path,
        default=None,
        help="Root directory for stage-pairwise downstream outputs.",
    )
    parser.add_argument(
        "--sf-threshold",
        type=float,
        default=float(os.environ.get("SSDI_SF_THRESHOLD", "1e-4")),
    )
    return parser.parse_args()


def export_triplet(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")


def bh_fdr(pvals: Iterable[float], alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(list(pvals), dtype=float)
    q = np.full_like(p, np.nan)
    sig = np.zeros_like(p, dtype=bool)
    ok = np.isfinite(p)
    if not np.any(ok):
        return q, sig
    pv = p[ok]
    order = np.argsort(pv)
    ps = pv[order]
    m = len(ps)
    r = np.arange(1, m + 1)
    qs = np.minimum.accumulate((ps * m / r)[::-1])[::-1]
    qu = np.empty_like(qs)
    qu[order] = np.clip(qs, 0, 1)
    q[ok] = qu
    sig[ok] = qu <= alpha
    return q, sig


def wilson_ci95(success: int, total: int) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    z = 1.959963984540054
    phat = success / total
    den = 1.0 + (z * z) / total
    ctr = (phat + (z * z) / (2.0 * total)) / den
    half = z * np.sqrt((phat * (1.0 - phat) + (z * z) / (4.0 * total)) / total) / den
    return float(ctr - half), float(ctr + half)


def star_label(q: float) -> str:
    if not np.isfinite(q) or q > 0.05:
        return ""
    if q < 1e-3:
        return "***"
    if q < 1e-2:
        return "**"
    return "*"


def canonical_band_order(df: pd.DataFrame) -> list[str]:
    observed = [b for b in DEFAULT_BAND_ORDER if b in set(df["band"].astype(str))]
    leftovers = sorted(set(df["band"].astype(str)) - set(observed))
    return observed + leftovers


def prep_long(df: pd.DataFrame, bands: list[str]) -> pd.DataFrame:
    out = df.copy()
    out = out[out["band"].isin(bands)].copy()
    out["band"] = pd.Categorical(out["band"], bands, ordered=True)
    out["band_greek"] = out["band"].astype(str).map(lambda b: BAND_GREEK.get(b, b))
    out["stage_internal"] = pd.Categorical(out["stage_internal"], STAGE_ORDER, ordered=True)
    out["stage_label"] = out["stage_internal"].astype(str).map(stage_display)
    out["dd_min"] = out["dd_best"].astype(float)
    return out


def parse_subject_id(name: str) -> str | None:
    m = re.match(r"^(EPCTL\d{2})_ssdi_spectral_domain_results\.mat$", name)
    if not m:
        return None
    return m.group(1)


def load_results(path: Path) -> dict:
    try:
        mat = sio.loadmat(path, simplify_cells=True)
    except Exception:
        return {}
    out = mat.get("results", {})
    return out if isinstance(out, dict) else {}


def upper_tri_values(M: np.ndarray) -> np.ndarray:
    idx = np.triu_indices_from(M, k=1)
    return np.asarray(M[idx], dtype=float)


def extract_sf_subject_long(
    spectral_dir: Path, bands: list[str], mvals: list[int], sf_threshold: float
) -> pd.DataFrame:
    rows: list[dict] = []
    for stage in STAGE_ORDER:
        stage_dir = spectral_dir / stage
        if not stage_dir.exists():
            continue
        for f in sorted(stage_dir.glob("*_ssdi_spectral_domain_results.mat")):
            sid = parse_subject_id(f.name)
            if sid is None:
                continue
            res = load_results(f)
            if not res:
                continue
            for band in bands:
                for m in mvals:
                    g = None
                    for key in (f"{band}_gopto_m{m}", f"{band}_goptp_m{m}"):
                        if key in res:
                            gt = np.asarray(res[key], dtype=float)
                            if gt.ndim == 2 and gt.shape[0] == gt.shape[1]:
                                g = gt
                                break
                    if g is None or g.shape[0] < 2:
                        continue
                    tri = upper_tri_values(g)
                    if tri.size == 0:
                        continue
                    rows.append(
                        {
                            "subject_id": sid,
                            "stage_internal": stage,
                            "stage_label": stage_display(stage),
                            "band": band,
                            "band_greek": BAND_GREEK.get(band, band),
                            "m": int(m),
                            "n_runs": int(g.shape[0]),
                            "n_pairs": int(tri.size),
                            "n_pairs_near_zero": int(np.sum(tri < sf_threshold)),
                            "sf_index": float(np.mean(tri < sf_threshold)),
                        }
                    )
    return pd.DataFrame(rows)


def paired_wake_stats(df: pd.DataFrame, value_col: str, q_col: str, sig_col: str) -> pd.DataFrame:
    rows: list[dict] = []
    bands = [str(b) for b in sorted(df["band"].unique())]
    mvals = sorted(df["m"].astype(int).unique())
    for stage in STAGE_ORDER:
        if stage == WAKE_STAGE:
            continue
        for band in bands:
            pvals: list[float] = []
            idxs: list[int] = []
            for m in mvals:
                a = df[
                    (df["stage_internal"] == WAKE_STAGE)
                    & (df["band"] == band)
                    & (df["m"] == m)
                ][["subject_id", value_col]].rename(columns={value_col: "x_a"})
                b = df[
                    (df["stage_internal"] == stage)
                    & (df["band"] == band)
                    & (df["m"] == m)
                ][["subject_id", value_col]].rename(columns={value_col: "x_b"})
                x = b.merge(a, on="subject_id")
                if len(x) < 3:
                    continue
                delta = x["x_b"].to_numpy(dtype=float) - x["x_a"].to_numpy(dtype=float)
                try:
                    stat, p = wilcoxon(
                        x["x_b"].to_numpy(dtype=float),
                        x["x_a"].to_numpy(dtype=float),
                        alternative="two-sided",
                        zero_method="wilcox",
                        mode="approx",
                    )
                except Exception:
                    stat, p = np.nan, np.nan
                rows.append(
                    {
                        "stage_a": WAKE_STAGE,
                        "stage_b": stage,
                        "band": band,
                        "m": int(m),
                        "n_subjects": int(len(x)),
                        "delta_mean": float(np.nanmean(delta)),
                        "delta_median": float(np.nanmedian(delta)),
                        "stat": float(stat),
                        "p": float(p),
                    }
                )
                pvals.append(p)
                idxs.append(len(rows) - 1)
            if pvals:
                q, sig = bh_fdr(pvals, alpha=0.05)
                for i, ridx in enumerate(idxs):
                    rows[ridx][q_col] = float(q[i])
                    rows[ridx][sig_col] = bool(sig[i])
    return pd.DataFrame(rows)


def compute_sf_pooled(sf_df: pd.DataFrame) -> pd.DataFrame:
    grp = (
        sf_df.groupby(["stage_internal", "band", "m"], as_index=False, observed=True)[["n_pairs", "n_pairs_near_zero"]]
        .sum()
    )
    grp["sf_index_pooled"] = grp["n_pairs_near_zero"] / grp["n_pairs"]
    grp["stage_label"] = grp["stage_internal"].map(stage_display)
    grp["band_greek"] = grp["band"].map(lambda b: BAND_GREEK.get(str(b), str(b)))
    return grp


def compute_pooled_wake_tests(sf_pooled: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for stage in STAGE_ORDER:
        if stage == WAKE_STAGE:
            continue
        for band in sorted(sf_pooled["band"].astype(str).unique()):
            pvals: list[float] = []
            idxs: list[int] = []
            for m in sorted(sf_pooled["m"].astype(int).unique()):
                wa = sf_pooled[
                    (sf_pooled["stage_internal"] == WAKE_STAGE)
                    & (sf_pooled["band"] == band)
                    & (sf_pooled["m"] == m)
                ]
                st = sf_pooled[
                    (sf_pooled["stage_internal"] == stage)
                    & (sf_pooled["band"] == band)
                    & (sf_pooled["m"] == m)
                ]
                if wa.empty or st.empty:
                    continue
                wa0 = wa.iloc[0]
                st0 = st.iloc[0]
                n1 = int(wa0["n_pairs"])
                n2 = int(st0["n_pairs"])
                x1 = int(wa0["n_pairs_near_zero"])
                x2 = int(st0["n_pairs_near_zero"])
                if n1 <= 0 or n2 <= 0:
                    continue
                p1 = x1 / n1
                p2 = x2 / n2
                p_pool = (x1 + x2) / (n1 + n2)
                se = np.sqrt(p_pool * (1.0 - p_pool) * (1.0 / n1 + 1.0 / n2))
                if not np.isfinite(se) or se <= 0:
                    z = np.nan
                    p = np.nan
                else:
                    z = (p2 - p1) / se
                    p = 2.0 * (1.0 - norm.cdf(abs(z)))
                rows.append(
                    {
                        "stage_a": WAKE_STAGE,
                        "stage_b": stage,
                        "band": band,
                        "m": int(m),
                        "delta": float(p2 - p1),
                        "z": float(z),
                        "p": float(p),
                    }
                )
                pvals.append(p)
                idxs.append(len(rows) - 1)
            if pvals:
                q, sig = bh_fdr(pvals, alpha=0.05)
                for i, ridx in enumerate(idxs):
                    rows[ridx]["q_bh"] = float(q[i])
                    rows[ridx]["sig_bh"] = bool(sig[i])
    return pd.DataFrame(rows)


def figure_violin_by_macro(long_df: pd.DataFrame, out_dir: Path, bands: list[str]) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    mvals = sorted(long_df["m"].unique())
    n_panels = len(mvals)
    ncols = 3
    nrows = int(np.ceil(n_panels / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(7.2 * ncols, 5.4 * nrows), sharey=True)
    axes = np.atleast_1d(axes).ravel()

    for i, m in enumerate(mvals):
        ax = axes[i]
        sub = long_df[long_df["m"] == m].copy()
        order = [BAND_GREEK.get(b, b) for b in bands]

        sns.violinplot(
            data=sub,
            x="band_greek",
            y="dd_min",
            hue="stage_internal",
            order=order,
            hue_order=STAGE_ORDER,
            palette=STAGE_PALETTE,
            density_norm="width",
            inner="quartile",
            cut=0,
            linewidth=0.8,
            width=0.82,
            dodge=True,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x="band_greek",
            y="dd_min",
            hue="stage_internal",
            order=order,
            hue_order=STAGE_ORDER,
            palette=STAGE_PALETTE,
            dodge=True,
            size=2.4,
            alpha=0.35,
            edgecolor="white",
            linewidth=0.2,
            ax=ax,
        )
        if ax.get_legend() is not None:
            ax.get_legend().remove()
        ax.set_title(f"{int(m)}-Macro", fontsize=18, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Dynamical Dependence" if i % ncols == 0 else "")
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=12)

    for j in range(n_panels, len(axes)):
        axes[j].axis("off")

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=STAGE_PALETTE[s], markersize=8, label=stage_display(s))
        for s in STAGE_ORDER
    ]
    fig.legend(handles=handles, labels=[stage_display(s) for s in STAGE_ORDER], loc="upper center", ncol=5, frameon=False)
    fig.suptitle("Canonical DD Min by Macro Size and Band (Sleep Stages)", fontsize=22, weight="bold", y=1.02)
    fig.tight_layout()

    export_triplet(fig, out_dir / "canonical_dd_min_violin_by_macro")
    plt.close(fig)


def figure_grand_average(long_df: pd.DataFrame, out_dir: Path, bands: list[str]) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    sub = long_df.copy()
    order = [BAND_GREEK.get(b, b) for b in bands]
    fig, ax = plt.subplots(figsize=(22, 8))
    sns.violinplot(
        data=sub,
        x="band_greek",
        y="dd_min",
        hue="stage_internal",
        order=order,
        hue_order=STAGE_ORDER,
        palette=STAGE_PALETTE,
        density_norm="width",
        inner="quartile",
        cut=0,
        linewidth=0.8,
        width=0.72,
        dodge=True,
        ax=ax,
    )
    sns.stripplot(
        data=sub,
        x="band_greek",
        y="dd_min",
        hue="stage_internal",
        order=order,
        hue_order=STAGE_ORDER,
        palette=STAGE_PALETTE,
        dodge=True,
        size=5.0,
        alpha=0.45,
        edgecolor="none",
        linewidth=0.0,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    by_label = {}
    for h, l in zip(handles, labels):
        if l in STAGE_ORDER and l not in by_label:
            by_label[l] = h
    leg_handles = [by_label[s] for s in STAGE_ORDER if s in by_label]
    leg_labels = [stage_display(s) for s in STAGE_ORDER if s in by_label]
    if leg_handles:
        ax.legend(
            leg_handles,
            leg_labels,
            title="Stage",
            loc="upper right",
            bbox_to_anchor=(1.0, 0.92),
            frameon=False,
            fontsize=16,
            title_fontsize=18,
        )

    ann_ax = fig.add_axes([0.84, 0.17, 0.05, 0.70])
    ann_ax.axis("off")
    dd_cmap = LinearSegmentedColormap.from_list("dd_arrow", ["#FFFFFF", "#F8E79A", "#F4A261"])
    grad = np.linspace(0.0, 1.0, 256, dtype=float).reshape(-1, 1)
    ann_ax.imshow(grad, aspect="auto", origin="lower", cmap=dd_cmap, extent=[0.40, 0.60, 0.00, 0.84])
    ann_ax.add_patch(
        Polygon(
            [[0.08, 0.84], [0.92, 0.84], [0.50, 1.00]],
            closed=True,
            facecolor="#F4A261",
            edgecolor="none",
        )
    )
    ann_ax.set_xlim(0.0, 1.0)
    ann_ax.set_ylim(0.0, 1.0)
    ann_ax.text(1.00, 0.98, "Scale\nintegrated", ha="left", va="top", fontsize=24, weight="bold", linespacing=0.95)
    ann_ax.text(1.00, 0.02, "Scale\nseparated", ha="left", va="bottom", fontsize=24, weight="bold", linespacing=0.95)

    ax.set_xlabel("Frequency band", fontsize=24)
    ax.set_ylabel("Dynamical Dependence", fontsize=24)
    ax.tick_params(axis="x", labelsize=24)
    ax.tick_params(axis="y", labelsize=20)
    fig.tight_layout(rect=[0, 0, 0.80, 1])

    export_triplet(fig, out_dir / "canonical_dd_min_grand_average")
    plt.close(fig)


def figure_sf_pooled_main(
    sf_pooled: pd.DataFrame,
    pooled_wake_stats: pd.DataFrame,
    out_dir: Path,
    bands: list[str],
    sf_threshold: float,
) -> None:
    mvals = sorted(sf_pooled["m"].astype(int).unique())
    fig, axes = plt.subplots(1, len(bands), figsize=(6.2 * len(bands), 9.0), sharey=True)
    if len(bands) == 1:
        axes = np.asarray([axes])

    non_wake_stages = [s for s in STAGE_ORDER if s != WAKE_STAGE]
    x_offsets = {st: off for st, off in zip(non_wake_stages, np.linspace(-0.22, 0.22, len(non_wake_stages)))}
    y_rows = {st: y for st, y in zip(non_wake_stages, np.linspace(0.97, 0.82, len(non_wake_stages)))}
    star_pe = [pe.withStroke(linewidth=2.4, foreground="white")]
    for i, band in enumerate(bands):
        ax = axes[i]
        for stage in STAGE_ORDER:
            sub = sf_pooled[(sf_pooled["band"] == band) & (sf_pooled["stage_internal"] == stage)].sort_values("m")
            if sub.empty:
                continue
            x = sub["m"].to_numpy(dtype=float)
            y = sub["sf_index_pooled"].to_numpy(dtype=float)
            lo = []
            hi = []
            for _, r in sub.iterrows():
                l, h = wilson_ci95(int(r["n_pairs_near_zero"]), int(r["n_pairs"]))
                lo.append(l)
                hi.append(h)
            lo = np.asarray(lo, dtype=float)
            hi = np.asarray(hi, dtype=float)
            ax.plot(x, y, "-o", color=STAGE_PALETTE[stage], linewidth=3.0, markersize=7, label=stage_display(stage) if i == 0 else None)
            ax.fill_between(x, lo, hi, color=STAGE_PALETTE[stage], alpha=0.14, linewidth=0)

        st = pooled_wake_stats[pooled_wake_stats["band"] == band]
        for stage in non_wake_stages:
            ss = st[st["stage_b"] == stage].set_index("m")
            for m in mvals:
                if m not in ss.index:
                    continue
                star = star_label(float(ss.loc[m, "q_bh"]))
                if not star:
                    continue
                ax.text(
                    m + x_offsets[stage],
                    y_rows[stage],
                    star,
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="top",
                    fontsize=18,
                    weight="bold",
                    color=STAGE_PALETTE[stage],
                    path_effects=star_pe,
                    bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.68},
                )
        ax.set_title(BAND_GREEK.get(band, band), fontsize=34, weight="bold")
        ax.set_xticks(mvals)
        ax.set_xlabel("Macro", fontsize=28)
        ax.tick_params(axis="x", labelsize=24)
        ax.tick_params(axis="y", labelsize=22)
        ax.grid(alpha=0.25, linewidth=0.6)
    axes[0].set_ylabel("SF index", fontsize=30)

    handles = [
        plt.Line2D([], [], color=STAGE_PALETTE[s], marker="o", linestyle="-", linewidth=3, label=stage_display(s))
        for s in STAGE_ORDER
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.47, 1.02), ncol=5, fontsize=18, frameon=False)

    ann_ax = fig.add_axes([0.86, 0.14, 0.05, 0.74])
    ann_ax.axis("off")
    sf_cmap = LinearSegmentedColormap.from_list("sf_arrow", ["#FFFFFF", "#F8E79A", "#F4A261"])
    grad = np.linspace(0.0, 1.0, 256, dtype=float).reshape(-1, 1)
    ann_ax.imshow(grad, aspect="auto", origin="lower", cmap=sf_cmap, extent=[0.40, 0.60, 0.00, 0.84])
    ann_ax.add_patch(
        Polygon(
            [[0.08, 0.84], [0.92, 0.84], [0.50, 1.00]],
            closed=True,
            facecolor="#F4A261",
            edgecolor="none",
        )
    )
    ann_ax.set_xlim(0.0, 1.0)
    ann_ax.set_ylim(0.0, 1.0)
    ann_ax.text(1.10, 0.98, "Structured", ha="left", va="top", fontsize=26, weight="bold")
    ann_ax.text(1.10, 0.02, "Fragmented", ha="left", va="bottom", fontsize=26, weight="bold")

    fig.suptitle(f"Structure-Fragmentation Pooled Index (threshold < {sf_threshold:g})", fontsize=22, y=1.05, weight="bold")
    fig.subplots_adjust(left=0.05, right=0.84, top=0.86, bottom=0.22, wspace=0.14)
    export_triplet(fig, out_dir / "sf_index_pooled_main")
    plt.close(fig)


def figure_sf_stage_composite(
    sf_df: pd.DataFrame, sf_wake_stats: pd.DataFrame, out_dir: Path, bands: list[str]
) -> None:
    mvals = sorted(sf_df["m"].astype(int).unique())
    fig, axes = plt.subplots(1, len(bands), figsize=(6.0 * len(bands), 8.5), sharey=True)
    if len(bands) == 1:
        axes = np.asarray([axes])
    non_wake_stages = [s for s in STAGE_ORDER if s != WAKE_STAGE]
    x_offsets = {st: off for st, off in zip(non_wake_stages, np.linspace(-0.22, 0.22, len(non_wake_stages)))}
    y_rows = {st: y for st, y in zip(non_wake_stages, np.linspace(0.97, 0.82, len(non_wake_stages)))}
    star_pe = [pe.withStroke(linewidth=2.2, foreground="white")]
    for j, band in enumerate(bands):
        ax = axes[j]
        sub = sf_df[sf_df["band"] == band]
        if sub.empty:
            ax.axis("off")
            continue
        for stage in STAGE_ORDER:
            stage_sub = sub[sub["stage_internal"] == stage]
            means = []
            sems = []
            for m in mvals:
                vals = stage_sub[stage_sub["m"] == m]["sf_index"].to_numpy(dtype=float)
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    means.append(np.nan)
                    sems.append(np.nan)
                    continue
                means.append(float(np.mean(vals)))
                sems.append(float(np.std(vals, ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0)
            means = np.asarray(means, dtype=float)
            sems = np.asarray(sems, dtype=float)
            ax.plot(mvals, means, "-o", color=STAGE_PALETTE[stage], linewidth=3.2, markersize=6.5, label=stage_display(stage) if j == 0 else None)
            ax.fill_between(mvals, means - sems, means + sems, color=STAGE_PALETTE[stage], alpha=0.18, linewidth=0)

        st = sf_wake_stats[sf_wake_stats["band"] == band]
        for stage in non_wake_stages:
            ss = st[st["stage_b"] == stage].set_index("m")
            for m in mvals:
                if m not in ss.index:
                    continue
                q = float(ss.loc[m, "q_bh"])
                star = star_label(q)
                if not star:
                    continue
                ax.text(
                    m + x_offsets[stage],
                    y_rows[stage],
                    star,
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="top",
                    fontsize=20,
                    color=STAGE_PALETTE[stage],
                    weight="bold",
                    path_effects=star_pe,
                    bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.68},
                )
        ax.set_title(BAND_GREEK.get(band, band), fontsize=40, weight="bold")
        ax.set_xlabel("macro", fontsize=30)
        ax.tick_params(axis="x", labelsize=24)
        ax.tick_params(axis="y", labelsize=24)
        ax.grid(alpha=0.25, axis="y")

    axes[0].set_ylabel("SF index", fontsize=32)
    handles = [plt.Line2D([], [], color=STAGE_PALETTE[s], marker="o", linestyle="-", label=stage_display(s)) for s in STAGE_ORDER]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=20, bbox_to_anchor=(0.5, 1.02))

    ann_ax = fig.add_axes([0.86, 0.14, 0.05, 0.74])
    ann_ax.axis("off")
    sf_cmap = LinearSegmentedColormap.from_list("sf_intuition_arrow", ["#FFFFFF", "#F8E79A", "#F4A261"])
    grad = np.linspace(0.0, 1.0, 256, dtype=float).reshape(-1, 1)
    ann_ax.imshow(grad, aspect="auto", origin="lower", cmap=sf_cmap, extent=[0.40, 0.60, 0.00, 0.84])
    ann_ax.add_patch(
        Polygon(
            [[0.08, 0.84], [0.92, 0.84], [0.50, 1.00]],
            closed=True,
            facecolor="#F4A261",
            edgecolor="none",
        )
    )
    ann_ax.set_xlim(0.0, 1.0)
    ann_ax.set_ylim(0.0, 1.0)
    ann_ax.text(1.10, 0.98, "Structured", ha="left", va="top", fontsize=26, weight="bold")
    ann_ax.text(1.10, 0.02, "Fragmented", ha="left", va="bottom", fontsize=26, weight="bold")

    fig.tight_layout(rect=[0, 0, 0.86, 0.90])
    export_triplet(fig, out_dir / "SF_index_stage_band_composite")
    export_triplet(fig, out_dir / "SF_index_4x5_timewindow_band_composite")
    plt.close(fig)


def figure_dd_stage_composite(
    long_df: pd.DataFrame, dd_wake_grand: pd.DataFrame, out_dir: Path, bands: list[str]
) -> None:
    dd_avg = (
        long_df.groupby(["subject_id", "stage_internal", "band", "band_greek"], as_index=False, observed=True)["dd_min"]
        .mean()
    )
    x_base = np.arange(len(bands), dtype=float)
    width = 0.15
    offsets = {st: (i - 2) * width for i, st in enumerate(STAGE_ORDER)}
    fig, ax = plt.subplots(1, 1, figsize=(24, 9))
    rng = np.random.default_rng(0)

    y_global = [np.inf, -np.inf]
    for stage in STAGE_ORDER:
        for bi, band in enumerate(bands):
            vals = dd_avg[(dd_avg["stage_internal"] == stage) & (dd_avg["band"] == band)]["dd_min"].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            x = x_base[bi] + offsets[stage]
            vp = ax.violinplot(dataset=[vals], positions=[x], widths=0.95 * width, showmeans=False, showmedians=False, showextrema=False)
            for body in vp["bodies"]:
                body.set_facecolor(STAGE_PALETTE[stage])
                body.set_edgecolor("none")
                body.set_alpha(0.55 if stage == WAKE_STAGE else 0.42)
            jitter = (rng.random(vals.size) - 0.5) * (0.28 * width)
            ax.scatter(np.full(vals.size, x) + jitter, vals, s=8, color="black", alpha=0.25, zorder=3)
            y_global[0] = min(y_global[0], float(np.min(vals)))
            y_global[1] = max(y_global[1], float(np.max(vals)))

    st = dd_wake_grand[dd_wake_grand["stage_a"] == WAKE_STAGE]
    yr = (y_global[1] - y_global[0] + 1e-12) if np.isfinite(y_global[0]) else 1.0
    y_star = y_global[1] + 0.18 * yr
    for bi, band in enumerate(bands):
        sb = st[st["band"] == band]
        if sb.empty:
            continue
        for stage in [s for s in STAGE_ORDER if s != WAKE_STAGE]:
            rr = sb[sb["stage_b"] == stage]
            if rr.empty:
                continue
            star = star_label(float(rr.iloc[0].get("q_bh_pair", np.nan)))
            if not star:
                continue
            x = x_base[bi] + offsets[stage]
            ax.text(x, y_star, star, ha="center", va="bottom", fontsize=26, weight="bold", color=STAGE_PALETTE[stage])

    if np.isfinite(y_global[0]) and np.isfinite(y_global[1]):
        pad_low = 0.08 * (y_global[1] - y_global[0] + 1e-12)
        pad_high = 0.35 * (y_global[1] - y_global[0] + 1e-12)
        ax.set_ylim(y_global[0] - pad_low, y_global[1] + pad_high)

    ax.set_xticks(x_base)
    ax.set_xticklabels([BAND_GREEK.get(b, b) for b in bands], fontsize=23)
    ax.tick_params(axis="y", labelsize=20)
    ax.set_xlabel("Frequency band", fontsize=24)
    ax.set_ylabel("DD", fontsize=24, weight="bold")
    ax.grid(alpha=0.25, axis="y")
    legend_handles = [Patch(facecolor=STAGE_PALETTE[s], edgecolor="none", alpha=0.55 if s == WAKE_STAGE else 0.42, label=stage_display(s)) for s in STAGE_ORDER]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=13, ncol=1)

    ann_ax = fig.add_axes([0.86, 0.14, 0.05, 0.74])
    ann_ax.axis("off")
    dd_cmap = LinearSegmentedColormap.from_list("dd_intuition_arrow", ["#FFFFFF", "#F8E79A", "#F4A261"])
    grad = np.linspace(0.0, 1.0, 256, dtype=float).reshape(-1, 1)
    ann_ax.imshow(grad, aspect="auto", origin="lower", cmap=dd_cmap, extent=[0.40, 0.60, 0.00, 0.84])
    ann_ax.add_patch(
        Polygon(
            [[0.08, 0.84], [0.92, 0.84], [0.50, 1.00]],
            closed=True,
            facecolor="#F4A261",
            edgecolor="none",
        )
    )
    ann_ax.set_xlim(0.0, 1.0)
    ann_ax.set_ylim(0.0, 1.0)
    ann_ax.text(1.10, 0.98, "Scale\nintegrated", ha="left", va="top", fontsize=24, weight="bold", linespacing=0.95)
    ann_ax.text(1.10, 0.02, "Scale\nseparated", ha="left", va="bottom", fontsize=24, weight="bold", linespacing=0.95)

    fig.suptitle("DD by frequency band, Wake + sleep stages", fontsize=24, weight="bold")
    fig.tight_layout(rect=[0, 0, 0.84, 0.95])
    export_triplet(fig, out_dir / "grand_average_DD_composite_stages")
    export_triplet(fig, out_dir / "grand_average_DD_composite_timewindows")
    plt.close(fig)


def build_state_long(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    out = df[["subject_id", "stage_internal", "band", "band_greek", "m", value_col]].copy()
    out["state"] = out["stage_internal"].map(stage_display)
    out["m"] = out["m"].astype(int)
    out = (
        out.groupby(["subject_id", "state", "band", "band_greek", "m"], as_index=False, observed=True)[value_col]
        .mean()
    )
    return out


def compute_state_pairwise(
    state_df: pd.DataFrame,
    value_col: str,
    detailed_q_col: str,
    grand_q_col: str,
    detailed_family: str,
    grand_family: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    states = [stage_display(s) for s in STAGE_ORDER]
    pairs = list(combinations(states, 2))
    mvals = sorted(state_df["m"].astype(int).unique())
    bands = sorted(state_df["band"].astype(str).unique())
    rows = []
    grand_rows = []

    for state_a, state_b in pairs:
        for m in mvals:
            pvals = []
            idxs = []
            for band in bands:
                a = state_df[(state_df["state"] == state_a) & (state_df["band"] == band) & (state_df["m"] == m)][["subject_id", value_col]].rename(columns={value_col: "x_a"})
                b = state_df[(state_df["state"] == state_b) & (state_df["band"] == band) & (state_df["m"] == m)][["subject_id", value_col]].rename(columns={value_col: "x_b"})
                x = b.merge(a, on="subject_id")
                if len(x) < 3:
                    continue
                delta = x["x_b"].to_numpy(dtype=float) - x["x_a"].to_numpy(dtype=float)
                try:
                    stat, p = wilcoxon(
                        x["x_b"].to_numpy(dtype=float),
                        x["x_a"].to_numpy(dtype=float),
                        alternative="two-sided",
                        zero_method="wilcox",
                        mode="approx",
                    )
                except Exception:
                    stat, p = np.nan, np.nan
                rows.append(
                    {
                        "state_a": state_a,
                        "state_b": state_b,
                        "comparison": f"{state_a}_vs_{state_b}",
                        "band": band,
                        "m": int(m),
                        "n_subjects": int(len(x)),
                        "delta_mean": float(np.nanmean(delta)),
                        "delta_median": float(np.nanmedian(delta)),
                        "stat": float(stat),
                        "p": float(p),
                        "correction_family": detailed_family,
                    }
                )
                pvals.append(p)
                idxs.append(len(rows) - 1)
            if pvals:
                q, sig = bh_fdr(pvals, alpha=0.05)
                for i, ridx in enumerate(idxs):
                    rows[ridx][detailed_q_col] = float(q[i])
                    rows[ridx]["sig_detailed"] = bool(sig[i])

        avg = (
            state_df.groupby(["subject_id", "state", "band"], as_index=False, observed=True)[value_col]
            .mean()
            .rename(columns={value_col: "x"})
        )
        pvals = []
        idxs = []
        for band in bands:
            a = avg[(avg["state"] == state_a) & (avg["band"] == band)][["subject_id", "x"]].rename(columns={"x": "x_a"})
            b = avg[(avg["state"] == state_b) & (avg["band"] == band)][["subject_id", "x"]].rename(columns={"x": "x_b"})
            x = b.merge(a, on="subject_id")
            if len(x) < 3:
                continue
            delta = x["x_b"].to_numpy(dtype=float) - x["x_a"].to_numpy(dtype=float)
            try:
                stat, p = wilcoxon(
                    x["x_b"].to_numpy(dtype=float),
                    x["x_a"].to_numpy(dtype=float),
                    alternative="two-sided",
                    zero_method="wilcox",
                    mode="approx",
                )
            except Exception:
                stat, p = np.nan, np.nan
            grand_rows.append(
                {
                    "state_a": state_a,
                    "state_b": state_b,
                    "comparison": f"{state_a}_vs_{state_b}",
                    "band": band,
                    "n_subjects": int(len(x)),
                    "delta_mean": float(np.nanmean(delta)),
                    "delta_median": float(np.nanmedian(delta)),
                    "stat": float(stat),
                    "p": float(p),
                    "correction_family": grand_family,
                }
            )
            pvals.append(p)
            idxs.append(len(grand_rows) - 1)
        if pvals:
            q, sig = bh_fdr(pvals, alpha=0.05)
            for i, ridx in enumerate(idxs):
                grand_rows[ridx][grand_q_col] = float(q[i])
                grand_rows[ridx]["sig_grand"] = bool(sig[i])

    return pd.DataFrame(rows), pd.DataFrame(grand_rows)


def plot_state_trajectory(
    state_long: pd.DataFrame,
    value_col: str,
    out_stem: str,
    pairwise_fig_root: Path,
    bands: list[str],
    main_color: str,
    subject_color: str,
) -> None:
    agg = (
        state_long.groupby(["subject_id", "state", "band", "band_greek"], as_index=False, observed=True)[value_col]
        .mean()
    )
    states = [stage_display(s) for s in STAGE_ORDER]
    fig, axes = plt.subplots(1, len(bands), figsize=(6.2 * len(bands), 8.2), sharey=True)
    if len(bands) == 1:
        axes = np.asarray([axes])
    x_vals = np.arange(len(states))
    y_global = [np.inf, -np.inf]

    for bi, band in enumerate(bands):
        ax = axes[bi]
        sub = agg[agg["band"] == band]
        if sub.empty:
            ax.axis("off")
            continue
        piv = sub.pivot(index="subject_id", columns="state", values=value_col).reindex(columns=states)
        for _, row in piv.iterrows():
            y = row.to_numpy(dtype=float)
            mask = np.isfinite(y)
            if np.sum(mask) >= 2:
                ax.plot(x_vals[mask], y[mask], color=subject_color, alpha=0.25, linewidth=0.8)
        means = []
        sems = []
        for s in states:
            vals = sub[sub["state"] == s][value_col].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                means.append(np.nan)
                sems.append(np.nan)
            else:
                means.append(float(np.mean(vals)))
                sems.append(float(np.std(vals, ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0)
                y_global[0] = min(y_global[0], float(np.min(vals)))
                y_global[1] = max(y_global[1], float(np.max(vals)))
        means = np.asarray(means, dtype=float)
        sems = np.asarray(sems, dtype=float)
        ax.plot(x_vals, means, "-o", color=main_color, linewidth=4.0, markersize=7.0)
        ax.fill_between(x_vals, means - sems, means + sems, color=main_color, alpha=0.18, linewidth=0)
        ax.set_title(BAND_GREEK.get(band, band), fontsize=30, weight="bold")
        ax.set_xticks(x_vals)
        ax.set_xticklabels(states, rotation=25, ha="right", fontsize=22)
        ax.tick_params(axis="y", labelsize=20)
        ax.grid(alpha=0.25, axis="y")
        if bi != 0:
            ax.tick_params(axis="y", labelleft=False)

    if np.isfinite(y_global[0]) and np.isfinite(y_global[1]):
        pad = 0.08 * (y_global[1] - y_global[0] + 1e-12)
        for ax in axes:
            ax.set_ylim(y_global[0] - pad, y_global[1] + pad)
    fig.tight_layout()
    export_triplet(fig, pairwise_fig_root / out_stem)
    plt.close(fig)


def pairwise_matrices_for_band(grand_df: pd.DataFrame, band: str, q_col: str) -> tuple[np.ndarray, np.ndarray]:
    states = [stage_display(s) for s in STAGE_ORDER]
    n = len(states)
    delta = np.full((n, n), np.nan, dtype=float)
    qmat = np.full((n, n), np.nan, dtype=float)
    sub = grand_df[grand_df["band"] == band]
    idx = {s: i for i, s in enumerate(states)}
    for _, r in sub.iterrows():
        a = r["state_a"]
        b = r["state_b"]
        if a not in idx or b not in idx:
            continue
        ia = idx[a]
        ib = idx[b]
        d = float(r["delta_mean"])
        q = float(r[q_col]) if q_col in r and np.isfinite(r[q_col]) else np.nan
        delta[ia, ib] = d
        delta[ib, ia] = -d
        qmat[ia, ib] = q
        qmat[ib, ia] = q
    return delta, qmat


def plot_pairwise_heatmaps(
    grand_df: pd.DataFrame,
    out_stem: str,
    pairwise_fig_root: Path,
    bands: list[str],
    q_col: str = "q_bh_pair",
    cbar_label: str = "Mean stage difference (column - row)",
) -> None:
    states = [stage_display(s) for s in STAGE_ORDER]
    fig, axes = plt.subplots(1, len(bands), figsize=(7.0 * len(bands), 7.8), squeeze=False)
    axes = axes[0]

    vmax = 0.0
    mats = {}
    for band in bands:
        d, q = pairwise_matrices_for_band(grand_df, band, q_col=q_col)
        mats[band] = (d, q)
        if np.any(np.isfinite(d)):
            vmax = max(vmax, float(np.nanmax(np.abs(d))))
    if vmax <= 0:
        vmax = 1.0

    im_ref = None
    for i, band in enumerate(bands):
        ax = axes[i]
        d, q = mats[band]
        # Show upper triangle only (including no diagonal), hide lower triangle.
        tri_mask = np.tril(np.ones_like(d, dtype=bool), k=0)
        d_upper = np.ma.array(d, mask=tri_mask | ~np.isfinite(d))
        cmap = plt.cm.get_cmap("RdBu_r").copy()
        cmap.set_bad(color="white")
        im = ax.imshow(
            d_upper,
            cmap=cmap,
            vmin=-vmax,
            vmax=vmax,
            interpolation="nearest",
        )
        im_ref = im
        ax.set_title(BAND_GREEK.get(band, band), fontsize=30, weight="bold")
        ax.set_xticks(np.arange(len(states)))
        ax.set_yticks(np.arange(len(states)))
        ax.set_xticklabels(states, rotation=45, ha="right", fontsize=20)
        ax.set_yticklabels(states, fontsize=20)
        ax.grid(False)
        for r in range(len(states)):
            for c in range(len(states)):
                if r == c or not np.isfinite(d[r, c]) or c <= r:
                    continue
                star = star_label(float(q[r, c])) if np.isfinite(q[r, c]) else ""
                ax.text(c, r, f"{d[r, c]:+.3f}\n{star}", ha="center", va="center", fontsize=14, color="black")
        ax.set_xlim(-0.5, len(states) - 0.5)
        ax.set_ylim(len(states) - 0.5, -0.5)
        if i != 0:
            ax.tick_params(axis="y", labelleft=False)

    if im_ref is not None:
        cax = fig.add_axes([0.885, 0.24, 0.020, 0.66])
        cb = fig.colorbar(im_ref, cax=cax)
        cb.set_label(cbar_label, fontsize=24)
        cb.ax.tick_params(labelsize=20)
    fig.tight_layout(rect=[0, 0, 0.86, 1])
    export_triplet(fig, pairwise_fig_root / out_stem)
    plt.close(fig)


def compose_pairwise_overview_and_cleanup(pairwise_fig_root: Path) -> None:
    left_panels = [
        ("dd_stage_trajectory_by_band", "DD across stages"),
        ("sf_stage_trajectory_by_band", "SF across stages"),
    ]
    right_panels = [
        ("dd_pairwise_stage_heatmap_by_band", "change in DD"),
        ("sf_pairwise_stage_heatmap_by_band", "change in SF"),
    ]

    def _load_images(panels: list[tuple[str, str]]) -> list[np.ndarray]:
        imgs = []
        for stem, _ in panels:
            p = pairwise_fig_root / f"{stem}.png"
            if not p.exists():
                return []
            imgs.append(plt.imread(p))
        return imgs

    def _crop_whitespace(img: np.ndarray) -> np.ndarray:
        if img.ndim == 3 and img.shape[2] == 4:
            mask = img[..., 3] > 0.01
        else:
            rgb = img[..., :3] if img.ndim == 3 else img
            mask = np.any(rgb < 0.995, axis=-1) if img.ndim == 3 else (rgb < 0.995)
        if not np.any(mask):
            return img
        ys, xs = np.where(mask)
        y0, y1 = int(np.min(ys)), int(np.max(ys)) + 1
        x0, x1 = int(np.min(xs)), int(np.max(xs)) + 1
        pad = max(8, int(0.01 * max(img.shape[0], img.shape[1])))
        y0 = max(0, y0 - pad)
        y1 = min(img.shape[0], y1 + pad)
        x0 = max(0, x0 - pad)
        x1 = min(img.shape[1], x1 + pad)
        return img[y0:y1, x0:x1]

    def _render_vertical(images: list[np.ndarray], panels: list[tuple[str, str]], out_stem: str, figsize: tuple[float, float]) -> None:
        fig, axes = plt.subplots(2, 1, figsize=figsize)
        for ax, img, (_, title) in zip(np.atleast_1d(axes), images, panels):
            ax.imshow(img)
            ax.set_title(title, fontsize=28, weight="bold", pad=8)
            ax.axis("off")
        fig.tight_layout(rect=[0, 0, 1, 1], h_pad=0.0)
        export_triplet(fig, pairwise_fig_root / out_stem)
        plt.close(fig)

    left_images = [_crop_whitespace(im) for im in _load_images(left_panels)]
    right_images = [_crop_whitespace(im) for im in _load_images(right_panels)]
    if left_images and right_images:
        _render_vertical(left_images, left_panels, "stage_pairwise_overview_left_column", (20, 12.5))
        _render_vertical(right_images, right_panels, "stage_pairwise_overview_right_column", (24, 13.5))
        _render_vertical(left_images, left_panels, "timewindow_pairwise_overview_left_column", (20, 12.5))
        _render_vertical(right_images, right_panels, "timewindow_pairwise_overview_right_column", (24, 13.5))

    for stem, _ in left_panels + right_panels:
        for ext in ("png", "svg", "pdf"):
            p = pairwise_fig_root / f"{stem}.{ext}"
            if p.exists():
                p.unlink()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pairwise_root = args.pairwise_root or (args.out_dir.parent / "downstream_stage_pairwise")
    pairwise_stats_root = pairwise_root / "stats"
    pairwise_fig_root = pairwise_root / "figures"
    pairwise_rep_root = pairwise_root / "reports"
    pairwise_stats_root.mkdir(parents=True, exist_ok=True)
    pairwise_fig_root.mkdir(parents=True, exist_ok=True)
    pairwise_rep_root.mkdir(parents=True, exist_ok=True)

    long_df_raw = pd.read_csv(args.long_csv)
    bands = canonical_band_order(long_df_raw)
    long_df = prep_long(long_df_raw, bands=bands)
    mvals = sorted(long_df["m"].astype(int).unique())

    figure_violin_by_macro(long_df, args.out_dir, bands=bands)
    figure_grand_average(long_df, args.out_dir, bands=bands)

    sf_df = extract_sf_subject_long(args.spectral_dir, bands=bands, mvals=mvals, sf_threshold=args.sf_threshold)
    if sf_df.empty:
        raise RuntimeError(f"No SF rows from spectral dir: {args.spectral_dir}")
    sf_df["band"] = pd.Categorical(sf_df["band"], bands, ordered=True)
    sf_wake = paired_wake_stats(sf_df, value_col="sf_index", q_col="q_bh", sig_col="sig_bh")
    sf_pooled = compute_sf_pooled(sf_df)
    sf_pooled_wake = compute_pooled_wake_tests(sf_pooled)
    figure_sf_pooled_main(sf_pooled, sf_pooled_wake, args.out_dir, bands=bands, sf_threshold=args.sf_threshold)
    figure_sf_stage_composite(sf_df, sf_wake, args.out_dir, bands=bands)

    dd_state = build_state_long(long_df, value_col="dd_min")
    sf_state = build_state_long(sf_df, value_col="sf_index")
    dd_pair, dd_pair_grand = compute_state_pairwise(
        dd_state,
        value_col="dd_min",
        detailed_q_col="q_bh_pair_m",
        grand_q_col="q_bh_pair",
        detailed_family="BH within comparison x m across bands",
        grand_family="BH within comparison across bands (grand over m)",
    )
    sf_pair, sf_pair_grand = compute_state_pairwise(
        sf_state,
        value_col="sf_index",
        detailed_q_col="q_bh_pair_band",
        grand_q_col="q_bh_pair",
        detailed_family="BH within comparison x band across m",
        grand_family="BH within comparison across bands (grand over m)",
    )

    dd_wake_grand = dd_pair_grand[dd_pair_grand["state_a"] == stage_display(WAKE_STAGE)].copy()
    dd_wake_grand["stage_a"] = WAKE_STAGE
    dd_wake_grand["stage_b"] = dd_wake_grand["state_b"].map({stage_display(s): s for s in STAGE_ORDER})
    figure_dd_stage_composite(long_df, dd_wake_grand, args.out_dir, bands=bands)

    dd_state.to_csv(pairwise_stats_root / "dd_stage_long.csv", index=False)
    sf_state.to_csv(pairwise_stats_root / "sf_stage_long.csv", index=False)
    dd_pair.to_csv(pairwise_stats_root / "dd_stage_pairwise_stats.csv", index=False)
    dd_pair_grand.to_csv(pairwise_stats_root / "dd_stage_pairwise_grand_stats.csv", index=False)
    sf_pair.to_csv(pairwise_stats_root / "sf_stage_pairwise_stats.csv", index=False)
    sf_pair_grand.to_csv(pairwise_stats_root / "sf_stage_pairwise_grand_stats.csv", index=False)
    sf_df.to_csv(args.out_dir.parent / "sf_index_subject_long.csv", index=False)
    sf_wake.to_csv(args.out_dir.parent / "sf_index_wake_pairwise_stats.csv", index=False)
    sf_pooled.to_csv(args.out_dir.parent / "sf_index_pooled_counts.csv", index=False)
    sf_pooled_wake.to_csv(args.out_dir.parent / "sf_index_pooled_wake_stats.csv", index=False)

    plot_state_trajectory(
        state_long=dd_state,
        value_col="dd_min",
        out_stem="dd_stage_trajectory_by_band",
        pairwise_fig_root=pairwise_fig_root,
        bands=bands,
        main_color="#D97706",
        subject_color="#E7C9A8",
    )
    plot_state_trajectory(
        state_long=sf_state,
        value_col="sf_index",
        out_stem="sf_stage_trajectory_by_band",
        pairwise_fig_root=pairwise_fig_root,
        bands=bands,
        main_color="#F59E0B",
        subject_color="#F6D7A4",
    )
    plot_pairwise_heatmaps(
        grand_df=dd_pair_grand,
        out_stem="dd_pairwise_stage_heatmap_by_band",
        pairwise_fig_root=pairwise_fig_root,
        bands=bands,
        q_col="q_bh_pair",
        cbar_label="DD difference",
    )
    plot_pairwise_heatmaps(
        grand_df=sf_pair_grand,
        out_stem="sf_pairwise_stage_heatmap_by_band",
        pairwise_fig_root=pairwise_fig_root,
        bands=bands,
        q_col="q_bh_pair",
        cbar_label="SF difference",
    )
    compose_pairwise_overview_and_cleanup(pairwise_fig_root)

    with (pairwise_rep_root / "stage_pairwise_summary.txt").open("w", encoding="utf-8") as f:
        f.write("SomnoDynamics stage-pairwise summary\n")
        f.write("===================================\n")
        f.write("Test: paired Wilcoxon signed-rank (two-sided)\n")
        f.write("DD grand correction: BH within stage-pair across bands\n")
        f.write("SF grand correction: BH within stage-pair across bands\n")
        f.write(f"DD pairwise rows: {len(dd_pair_grand)}\n")
        f.write(f"SF pairwise rows: {len(sf_pair_grand)}\n")

    with (args.out_dir / "ssdi_figure_summary.txt").open("w", encoding="utf-8") as f:
        f.write("SomnoDynamics SSDI DMT-parity Figure Summary\n")
        f.write("===========================================\n")
        f.write(f"Long table: {args.long_csv}\n")
        f.write(f"Spectral source: {args.spectral_dir}\n")
        f.write(f"Bands: {', '.join(bands)}\n")
        f.write(f"Stages: {', '.join(stage_display(s) for s in STAGE_ORDER)}\n")
        f.write("Outputs include canonical DD summaries, pooled SF, stage composites, and stage pairwise overview panels.\n")

    print(f"[SomnoDynamics] SSDI DMT-parity figures -> {args.out_dir}")
    print(f"[SomnoDynamics] Stage pairwise outputs -> {pairwise_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
