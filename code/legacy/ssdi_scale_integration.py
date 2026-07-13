#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ssdi_scale_integration.py

Broadband SSDI scale-integration analysis.

Local index per scale:
    SI_local(n)  = 1 - d(M_n, M_{n+1})

Global index per scale:
    SI_global(n) = 1 - mean_m d(M_n, M_m), m != n

We:
  - compute these per subject, stage, and scale n
  - plot bar+scatter per scale (bars = mean per stage, points = subjects)
  - run within-subject Wilcoxon tests across stages for each scale
  - FDR-correct p-values (BH)
  - draw significance brackets ONLY on the local plot
"""

from pathlib import Path
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple
import itertools

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.io import loadmat
from scipy.linalg import subspace_angles
from scipy.stats import wilcoxon

try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("[WARN] statsmodels not found; FDR correction will be skipped.")

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------

TIME_DIR = Path(
    "/Users/borjan/code/healthy-sleep-project/code/results/ssdi_results/time"
)

TIME_PATTERN = "EPCTL{subject}_{stage}_ssdi_time_domain_results.mat"

STAGES = ["W_", "N1", "N2", "N3", "R_"]
MACRO_SIZES = list(range(2, 10))  # 2..9

def lopto_cell_index(n: int) -> int:
    return n - 1

FIG_DIR = TIME_DIR / "figs_scale_integration"

plt.style.use("seaborn-v0_8-whitegrid")


# ---------------------------------------------------------------------
# MATLAB helpers
# ---------------------------------------------------------------------
def parse_subject_stage(filename: str) -> Tuple[str, str]:
    m = re.match(r"EPCTL(\d+)_(\w+)_ssdi_.*\.mat", filename)
    if not m:
        raise ValueError(f"Cannot parse subject/stage from {filename}")
    subj = f"{int(m.group(1)):02d}"
    stage = m.group(2)
    return subj, stage


def load_results_mat(path: Path) -> Dict[str, np.ndarray]:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "results" not in mat:
        raise KeyError(f"'results' not found in {path}")
    res = mat["results"]

    if not hasattr(res, "dtype") or res.dtype.names is None:
        return _matlab_struct_to_dict(res)

    res_dict = {}
    for name in res.dtype.names:
        field_val = res[name]
        if isinstance(field_val, np.ndarray):
            if field_val.ndim == 0:
                res_dict[name] = field_val.item()
            elif field_val.size == 1:
                res_dict[name] = field_val.flat[0]
            else:
                res_dict[name] = field_val
        else:
            res_dict[name] = field_val
    return res_dict


def _matlab_struct_to_dict(matobj) -> Dict:
    from scipy.io.matlab.mio5_params import mat_struct

    if isinstance(matobj, mat_struct):
        out = {}
        for f in matobj._fieldnames:
            out[f] = _matlab_struct_to_dict(getattr(matobj, f))
        return out
    elif isinstance(matobj, np.ndarray) and matobj.dtype == object:
        return [_matlab_struct_to_dict(e) for e in matobj]
    else:
        return matobj


# ---------------------------------------------------------------------
# Subspace helpers
# ---------------------------------------------------------------------
def extract_subspace(L_data, n: int) -> np.ndarray:
    if isinstance(L_data, dict):
        if n not in L_data:
            raise IndexError(f"Macro size n={n} not in L_dict.")
        L = np.asarray(L_data[n])
    elif isinstance(L_data, (list, tuple, np.ndarray)):
        idx = lopto_cell_index(n)
        if not (0 <= idx < len(L_data)):
            raise IndexError(f"Macro size n={n} not available in L_list.")
        L = np.asarray(L_data[idx])
    else:
        L = np.asarray(L_data)

    if L.ndim == 0:
        raise ValueError("L cell contains a scalar, not a matrix.")
    elif L.ndim == 1:
        L = L.reshape(-1, 1)
    elif L.ndim > 2:
        L = L.reshape(L.shape[0], -1)

    if L.shape[1] < n:
        return L
    return L[:, :n]


def subspace_distance(U: np.ndarray, V: np.ndarray) -> float:
    U = np.asarray(U)
    V = np.asarray(V)
    angles = subspace_angles(U, V)
    return float(np.sqrt(np.mean(np.sin(angles) ** 2)))


# ---------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------
@dataclass
class ScaleIntegrationRecord:
    stage: str
    subject: str
    si_local: float
    si_global: float
    min_scale: int
    max_scale: int


@dataclass
class LocalScaleRecord:
    stage: str
    subject: str
    scale: int
    si_local_n: float


@dataclass
class GlobalScaleRecord:
    stage: str
    subject: str
    scale: int
    si_global_n: float


# ---------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------
def compute_scale_integration() -> Tuple[
    List[ScaleIntegrationRecord],
    List[LocalScaleRecord],
    List[GlobalScaleRecord],
]:
    agg_records: List[ScaleIntegrationRecord] = []
    local_scale_records: List[LocalScaleRecord] = []
    global_scale_records: List[GlobalScaleRecord] = []

    files = sorted(TIME_DIR.glob("EPCTL*_ssdi_time_domain_results.mat"))
    if not files:
        raise FileNotFoundError(f"No time-domain SSDI files found in {TIME_DIR}")

    for path in files:
        subj, stage = parse_subject_stage(path.name)
        if stage not in STAGES:
            continue

        print(f"Processing subject {subj}, stage {stage}")
        res = load_results_mat(path)

        L_data = None
        for name in ["Lopto_m", "lopto_m", "L_m", "Lopto"]:
            if name in res:
                L_data = res[name]
                if isinstance(L_data, np.ndarray):
                    if L_data.ndim > 1 and L_data.shape[0] == 1 and L_data.shape[1] == 1:
                        L_data = L_data[0, 0]
                    elif L_data.ndim == 0:
                        L_data = L_data.item()
                break

        if L_data is None:
            print(f"  [WARN] L-matrix field not found in {path}")
            continue

        subspaces: Dict[int, np.ndarray] = {}

        if isinstance(L_data, (list, np.ndarray)):
            available_sizes = list(range(2, len(L_data) + 2))
        elif isinstance(L_data, dict):
            available_sizes = [n for n in MACRO_SIZES if n in L_data]
        else:
            available_sizes = MACRO_SIZES

        for n in MACRO_SIZES:
            if n not in available_sizes:
                continue
            try:
                U_n = extract_subspace(L_data, n)
                subspaces[n] = U_n
            except Exception as e:
                print(f"    [WARN] Could not extract subspace for n={n}: {e}")
                continue

        if len(subspaces) < 2:
            print(f"  [WARN] Need at least 2 subspaces to compute distances, got {len(subspaces)}.")
            continue

        scales = sorted(subspaces.keys())
        local_dists = []
        pair_dists = []

        # local per-scale
        for i in range(len(scales) - 1):
            n1 = scales[i]
            n2 = scales[i + 1]
            d = subspace_distance(subspaces[n1], subspaces[n2])
            local_dists.append(d)
            local_scale_records.append(
                LocalScaleRecord(
                    stage=stage,
                    subject=subj,
                    scale=n1,
                    si_local_n=1.0 - float(d),
                )
            )

        # global per-scale
        for i, n1 in enumerate(scales):
            d_list = []
            for j, n2 in enumerate(scales):
                if n1 == n2:
                    continue
                d = subspace_distance(subspaces[n1], subspaces[n2])
                d_list.append(d)
                if j > i:
                    pair_dists.append(d)
            if d_list:
                global_scale_records.append(
                    GlobalScaleRecord(
                        stage=stage,
                        subject=subj,
                        scale=n1,
                        si_global_n=1.0 - float(np.mean(d_list)),
                    )
                )

        if not local_dists or not pair_dists:
            continue

        si_local = 1.0 - float(np.mean(local_dists))
        si_global = 1.0 - float(np.mean(pair_dists))

        agg_records.append(
            ScaleIntegrationRecord(
                stage=stage,
                subject=subj,
                si_local=si_local,
                si_global=si_global,
                min_scale=min(scales),
                max_scale=max(scales),
            )
        )

    return agg_records, local_scale_records, global_scale_records


# ---------------------------------------------------------------------
# Stats (Wilcoxon + FDR)
# ---------------------------------------------------------------------
def run_stats(
    local_scale_records: List[LocalScaleRecord],
    global_scale_records: List[GlobalScaleRecord],
    outdir: Path,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    outdir.mkdir(parents=True, exist_ok=True)

    def records_to_df(records, value_field: str) -> pd.DataFrame:
        rows = []
        for r in records:
            val = getattr(r, value_field)
            if np.isfinite(val):
                rows.append(
                    {
                        "Subject": r.subject,
                        "Stage": r.stage,
                        "Scale": r.scale,
                        "Value": float(val),
                    }
                )
        if not rows:
            return pd.DataFrame(columns=["Subject", "Stage", "Scale", "Value"])
        df_ = pd.DataFrame(rows)
        df_["Stage"] = pd.Categorical(df_["Stage"], categories=STAGES, ordered=True)
        return df_

    df_local = records_to_df(local_scale_records, "si_local_n")
    df_global = records_to_df(global_scale_records, "si_global_n")

    def do_wilcoxon_table(df: pd.DataFrame, metric_name: str) -> pd.DataFrame:
        if df.empty:
            print(f"[WARN] No data for {metric_name}, skipping stats.")
            return pd.DataFrame()

        results = []
        unique_scales = sorted(df["Scale"].unique())

        for scale in unique_scales:
            df_s = df[df["Scale"] == scale]
            if df_s.empty:
                continue

            for st1, st2 in itertools.combinations(STAGES, 2):
                d1 = df_s[df_s["Stage"] == st1][["Subject", "Value"]]
                d2 = df_s[df_s["Stage"] == st2][["Subject", "Value"]]

                if d1.empty or d2.empty:
                    continue

                merged = d1.merge(d2, on="Subject", suffixes=("_" + st1, "_" + st2))
                if len(merged) < 5:
                    continue

                vals1 = merged["Value_" + st1].values
                vals2 = merged["Value_" + st2].values

                try:
                    stat, p = wilcoxon(vals1, vals2, alternative="two-sided")
                except ValueError:
                    stat, p = np.nan, np.nan

                results.append(
                    {
                        "metric": metric_name,
                        "scale": scale,
                        "stage1": st1,
                        "stage2": st2,
                        "n": len(merged),
                        "median_" + st1: float(np.median(vals1)),
                        "median_" + st2: float(np.median(vals2)),
                        "stat": float(stat) if np.isfinite(stat) else np.nan,
                        "p_raw": float(p) if np.isfinite(p) else np.nan,
                    }
                )

        if not results:
            return pd.DataFrame()

        res_df = pd.DataFrame(results)

        if HAS_STATSMODELS:
            mask = res_df["p_raw"].notna()
            if mask.any():
                rej, p_corr, _, _ = multipletests(
                    res_df.loc[mask, "p_raw"].values, method="fdr_bh"
                )
                res_df.loc[mask, "p_fdr"] = p_corr
                res_df.loc[mask, "reject_fdr"] = rej
            else:
                res_df["p_fdr"] = np.nan
                res_df["reject_fdr"] = False
        else:
            res_df["p_fdr"] = np.nan
            res_df["reject_fdr"] = False

        return res_df

    stats_local = do_wilcoxon_table(df_local, "local")
    stats_global = do_wilcoxon_table(df_global, "global")

    if not stats_local.empty:
        stats_local.to_csv(outdir / "stats_scale_integration_local.csv", index=False)
        print(f"Saved stats_scale_integration_local.csv in {outdir}")

    if not stats_global.empty:
        stats_global.to_csv(outdir / "stats_scale_integration_global.csv", index=False)
        print(f"Saved stats_scale_integration_global.csv in {outdir}")

    return stats_local, stats_global


# ---------------------------------------------------------------------
# Plotting (with sig brackets for local only)
# ---------------------------------------------------------------------
def _p_to_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "n.s."


def plot_scale_integration(
    agg_records: List[ScaleIntegrationRecord],
    local_scale_records: List[LocalScaleRecord],
    global_scale_records: List[GlobalScaleRecord],
    stats_local: pd.DataFrame,
    stats_global: pd.DataFrame,
    outdir: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    stages = STAGES
    unique_scales = sorted({r.scale for r in local_scale_records + global_scale_records})
    if not unique_scales:
        print("[WARN] No scales found in per-scale records.")
        return

    scale_palette = sns.color_palette("viridis", n_colors=len(unique_scales))
    scale_to_color = {s: scale_palette[i] for i, s in enumerate(unique_scales)}

    def _bar_scatter_scale(
        scale_records: List,
        title: str,
        ylabel: str,
        filename: str,
        value_attr: str,
        stats_df: pd.DataFrame = None,
        add_brackets: bool = False,
    ) -> None:
        if not scale_records:
            print(f"[WARN] No data for {filename}")
            return

        plot_data = []
        for r in scale_records:
            val = getattr(r, value_attr)
            if np.isfinite(val):
                plot_data.append(
                    {
                        "Stage": r.stage,
                        "Subject": r.subject,
                        "Scale": r.scale,
                        "Value": float(val),
                    }
                )

        if not plot_data:
            print(f"[WARN] All values NaN for {filename}")
            return

        df = pd.DataFrame(plot_data)
        df["Stage"] = pd.Categorical(df["Stage"], categories=stages, ordered=True)

        fig, ax = plt.subplots(figsize=(10, 5))

        x_stage = np.arange(len(stages))
        n_scales = len(unique_scales)
        total_width = 0.8
        bar_width = total_width / n_scales
        offset_start = -total_width / 2 + bar_width / 2

        # Bars
        for si, st in enumerate(stages):
            for k, scale in enumerate(unique_scales):
                sub = df[(df["Stage"] == st) & (df["Scale"] == scale)]
                if sub.empty:
                    continue
                mean_val = sub["Value"].mean()
                x_pos = x_stage[si] + offset_start + k * bar_width
                ax.bar(
                    x_pos,
                    mean_val,
                    width=bar_width * 0.9,
                    color=scale_to_color[scale],
                    edgecolor="black",
                    alpha=0.7,
                    zorder=1,
                )

        # Scatter
        max_y = -np.inf
        for si, st in enumerate(stages):
            for k, scale in enumerate(unique_scales):
                sub = df[(df["Stage"] == st) & (df["Scale"] == scale)]
                if sub.empty:
                    continue
                base_x = x_stage[si] + offset_start + k * bar_width
                xs = base_x + (np.random.rand(len(sub)) - 0.5) * (bar_width * 0.6)
                ys = sub["Value"].values
                ax.scatter(
                    xs,
                    ys,
                    c=[scale_to_color[scale]] * len(sub),
                    edgecolor="black",
                    linewidth=0.4,
                    s=25,
                    alpha=0.9,
                    zorder=5,
                )
                max_y = max(max_y, float(ys.max()))

        # y-limits slightly above max for brackets
        if np.isfinite(max_y):
            ax.set_ylim(-0.05, min(1.1, max_y + 0.15))
        else:
            ax.set_ylim(-0.05, 1.05)

        ax.set_xticks(x_stage)
        ax.set_xticklabels(stages)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        # Legend for scales
        legend_elements = [
            plt.Line2D(
                [0],
                [0],
                marker="s",
                color="w",
                markerfacecolor=scale_to_color[s],
                markeredgecolor="black",
                markeredgewidth=0.5,
                markersize=8,
                label=f"{s} → {s+1}",
            )
            for s in unique_scales
        ]
        ax.legend(
            handles=legend_elements,
            title="Scale integration",
            loc="upper right",
            frameon=True,
            fontsize=8,
            title_fontsize=9,
        )

        ax.text(
            0.01,
            0.97,
            "1 = strong scale-integration\n0 = no overlap between scales\n"
            "Bars = mean per (stage, scale)\nPoints = subjects, coloured by scale n",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        )

        ax.axhline(0.5, color="k", linestyle="--", linewidth=0.7, alpha=0.5)
        ax.grid(axis="y", linestyle="--", alpha=0.5)

        # -------- brackets for significant pairs (only when requested) --------
        if add_brackets and stats_df is not None and not stats_df.empty:
            # only FDR-significant rows
            if "reject_fdr" in stats_df.columns:
                sig_rows = stats_df[stats_df["reject_fdr"] == True]
            else:
                sig_rows = stats_df[stats_df["p_raw"] < 0.05]

            # offset per scale to stack brackets
            h_step = 0.03
            for scale in unique_scales:
                rows_s = sig_rows[sig_rows["scale"] == scale]
                if rows_s.empty:
                    continue

                # determine base height just above all data for this scale
                # use max of points for that scale across stages
                sub_scale = df[df["Scale"] == scale]
                if sub_scale.empty:
                    continue
                y_base = float(sub_scale["Value"].max()) + 0.03
                for idx, row in rows_s.iterrows():
                    st1 = row["stage1"]
                    st2 = row["stage2"]
                    p_val = row.get("p_fdr", row["p_raw"])

                    stars = _p_to_stars(p_val)
                    if stars == "n.s.":
                        continue

                    k = unique_scales.index(scale)
                    x1 = x_stage[stages.index(st1)] + offset_start + k * bar_width
                    x2 = x_stage[stages.index(st2)] + offset_start + k * bar_width

                    y = y_base + h_step * idx  # crude stacking; works given few sigs
                    ax.plot([x1, x1, x2, x2], [y, y + h_step, y + h_step, y], color="k", lw=0.8)
                    ax.text(
                        (x1 + x2) / 2,
                        y + h_step,
                        stars,
                        ha="center",
                        va="bottom",
                        fontsize=8,
                    )

            # extend ylim a bit if needed
            ylim = ax.get_ylim()
            ax.set_ylim(ylim[0], max(ylim[1], max_y + 0.2))

        fig.tight_layout()
        fig.savefig(outdir / filename, dpi=300)

    # Local SI(n) with brackets
    _bar_scatter_scale(
        local_scale_records,
        title="Local scale-integration per scale\n(adjacent macro scales)",
        ylabel=r"$SI_{\mathrm{local}}(n) = 1 - d(M_n, M_{n+1})$",
        filename="scale_integration_local_per_scale.png",
        value_attr="si_local_n",
        stats_df=stats_local,
        add_brackets=True,
    )

    # Global SI(n) WITHOUT brackets
    _bar_scatter_scale(
        global_scale_records,
        title="Global scale-integration per scale\n(all macro scale pairs)",
        ylabel=r"$SI_{\mathrm{global}}(n) = 1 - \langle d(M_n, M_m) \rangle_{m\neq n}$",
        filename="scale_integration_global_per_scale.png",
        value_attr="si_global_n",
        stats_df=None,
        add_brackets=False,
    )


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    agg_records, local_scale_records, global_scale_records = compute_scale_integration()
    print(f"Computed aggregate scale-integration for {len(agg_records)} subject–stage combinations.")
    print(f"Local per-scale records:  {len(local_scale_records)}")
    print(f"Global per-scale records: {len(global_scale_records)}")

    stats_local, stats_global = run_stats(local_scale_records, global_scale_records, FIG_DIR)
    plot_scale_integration(
        agg_records,
        local_scale_records,
        global_scale_records,
        stats_local,
        stats_global,
        FIG_DIR,
    )
    print(f"Saved scale-integration figures and stats to: {FIG_DIR}")


if __name__ == "__main__":
    main()
