#!/usr/bin/env python3
"""Analyze 3-macro beta contributions across subjects and sleep stages."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.io import loadmat
from scipy.stats import friedmanchisquare
from scipy.stats import wilcoxon

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import load_config, resolve_repo_path, stage_display, stage_order  # noqa: E402

STAGES = stage_order()
BANDS = ["delta", "theta", "alpha", "sigma", "beta", "gamma"]
REGION_ORDER = ["Frontal", "Sensorimotor", "Parietal", "Temporal", "Occipital", "Cingulate", "Limbic"]
BAND_DISPLAY = {
    "delta": "δ",
    "theta": "θ",
    "alpha": "α",
    "sigma": "σ",
    "beta": "β",
    "gamma": "γ",
}
SUBJECT_RE = re.compile(r"^(EPCTL\d{2})_ssdi_spectral_domain_results\.mat$")


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    p = argparse.ArgumentParser(description="Analyze region beta contributions for m=3 SSDI models")
    p.add_argument(
        "--spectral-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["ssdi_root"]) / "stage_ssdi_anatomical7",
        help="SSDI stage root with stage subfolders containing *_ssdi_spectral_domain_results.mat",
    )
    p.add_argument(
        "--out-root",
        type=Path,
        default=None,
        help="Output root for beta analysis (default: <spectral-dir>/beta_macro3_analysis)",
    )
    p.add_argument("--m", type=int, default=3, help="Macro dimension m (default: 3)")
    p.add_argument("--alpha", type=float, default=0.05, help="Alpha for significance calls")
    return p.parse_args()


def mat_struct_to_dict(obj: Any) -> Any:
    from scipy.io.matlab import mat_struct  # type: ignore

    if isinstance(obj, mat_struct):
        out: dict[str, Any] = {}
        for name in obj._fieldnames:
            out[name] = mat_struct_to_dict(getattr(obj, name))
        return out
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        return [mat_struct_to_dict(v) for v in obj.tolist()]
    return obj


def load_results(path: Path) -> dict[str, Any]:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "results" not in mat:
        raise KeyError(f"'results' not found in {path}")
    res = mat_struct_to_dict(mat["results"])
    if not isinstance(res, dict):
        raise ValueError(f"'results' is not dict-like in {path}")
    return res


def fdr_adjust(pvals: list[float], method: str) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    sig = np.zeros_like(p, dtype=bool)
    ok = np.isfinite(p)
    if not np.any(ok):
        return q, sig
    pv = p[ok]
    order = np.argsort(pv)
    ps = pv[order]
    m = len(ps)
    r = np.arange(1, m + 1, dtype=float)
    factor = 1.0
    if method.lower() == "by":
        factor = float(np.sum(1.0 / r))
    qs = np.minimum.accumulate((ps * m * factor / r)[::-1])[::-1]
    qu = np.empty_like(qs)
    qu[order] = np.clip(qs, 0.0, 1.0)
    q[ok] = qu
    return q, sig


def wilcoxon_vs_null(values: np.ndarray, null_value: float) -> dict[str, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    out = {
        "n_subjects": float(x.size),
        "median_beta": np.nan,
        "mean_beta": np.nan,
        "median_diff": np.nan,
        "mean_diff": np.nan,
        "p_two_sided": np.nan,
        "p_greater": np.nan,
        "stat_two_sided": np.nan,
        "stat_greater": np.nan,
        "effect_r_abs": np.nan,
        "prop_above_null": np.nan,
    }
    if x.size == 0:
        return out
    d = x - null_value
    out["median_beta"] = float(np.nanmedian(x))
    out["mean_beta"] = float(np.nanmean(x))
    out["median_diff"] = float(np.nanmedian(d))
    out["mean_diff"] = float(np.nanmean(d))
    out["prop_above_null"] = float(np.mean(d > 0.0))
    nz = d[np.abs(d) > 0.0]
    if nz.size < 3:
        return out
    try:
        two = wilcoxon(d, alternative="two-sided", zero_method="wilcox", method="approx")
        out["p_two_sided"] = float(two.pvalue)
        out["stat_two_sided"] = float(two.statistic)
        z = getattr(two, "zstatistic", None)
        if z is not None and np.isfinite(z):
            out["effect_r_abs"] = float(abs(z) / math.sqrt(float(nz.size)))
    except Exception:
        pass
    try:
        gt = wilcoxon(d, alternative="greater", zero_method="wilcox", method="approx")
        out["p_greater"] = float(gt.pvalue)
        out["stat_greater"] = float(gt.statistic)
    except Exception:
        pass
    return out


def stage_name(stage_internal: str) -> str:
    return stage_display(stage_internal)


def extract_beta_rows(spectral_dir: Path, m: int) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for st in STAGES:
        for path in sorted((spectral_dir / st).glob("*_ssdi_spectral_domain_results.mat")):
            msub = SUBJECT_RE.match(path.name)
            if msub is None:
                continue
            subject = msub.group(1)
            res = load_results(path)
            meta = res.get("meta", {})
            labels = np.asarray(meta.get("labels", []), dtype=object).reshape(-1)
            if labels.size == 0:
                labels = np.asarray(
                    ["Frontal", "Sensorimotor", "Parietal", "Temporal", "Occipital", "Cingulate", "Limbic"],
                    dtype=object,
                )
            labels = np.asarray([str(x) for x in labels], dtype=object)
            n_regions = int(labels.size)
            null_beta = float(m) / float(n_regions) if n_regions > 0 else np.nan
            for band in BANDS:
                key = f"{band}_beta_m{m}"
                if key not in res:
                    continue
                arr = np.asarray(res[key], dtype=float)
                if arr.ndim == 1:
                    arr = arr.reshape(-1, 1)
                if arr.ndim != 2:
                    continue
                if arr.shape[0] != n_regions:
                    continue
                n_runs = int(arr.shape[1])
                best = arr[:, 0]
                for idx, region in enumerate(labels):
                    rows.append(
                        {
                            "subject_id": subject,
                            "stage_internal": st,
                            "stage_display": stage_name(st),
                            "band": band,
                            "m": int(m),
                            "region": str(region),
                            "n_regions": int(n_regions),
                            "n_runs": int(n_runs),
                            "beta_null": float(null_beta),
                            "beta_best": float(best[idx]),
                        }
                    )
    return pd.DataFrame(rows)


def run_stage_collapsed_tests(df: pd.DataFrame, alpha: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Stage-level (band-averaged per subject) is the primary "condition" analysis.
    agg = (
        df.groupby(["subject_id", "stage_internal", "stage_display", "region", "m", "n_regions", "beta_null"], as_index=False, observed=True)[
            "beta_best"
        ]
        .mean()
        .rename(columns={"beta_best": "beta_stage_avg"})
    )
    out_rows: list[dict[str, Any]] = []
    for st in STAGES:
        sdf = agg[agg["stage_internal"] == st]
        regions = sorted(sdf["region"].astype(str).unique().tolist())
        pvals: list[float] = []
        idxs: list[int] = []
        for rg in regions:
            r = sdf[sdf["region"] == rg]
            vals = r["beta_stage_avg"].to_numpy(dtype=float)
            null_value = float(r["beta_null"].iloc[0]) if len(r) else np.nan
            stats = wilcoxon_vs_null(vals, null_value)
            out_rows.append(
                {
                    "stage_internal": st,
                    "stage_display": stage_name(st),
                    "region": rg,
                    "m": int(r["m"].iloc[0]) if len(r) else np.nan,
                    "n_regions": int(r["n_regions"].iloc[0]) if len(r) else np.nan,
                    "beta_null": null_value,
                    **stats,
                }
            )
            pvals.append(stats["p_two_sided"])
            idxs.append(len(out_rows) - 1)
        q_by, _ = fdr_adjust(pvals, method="by")
        for j, ridx in enumerate(idxs):
            out_rows[ridx]["q_by_within_stage"] = float(q_by[j]) if np.isfinite(q_by[j]) else np.nan
            out_rows[ridx]["sig_by_within_stage"] = bool(np.isfinite(q_by[j]) and q_by[j] <= alpha)
            out_rows[ridx]["contributes_above_null"] = bool(
                np.isfinite(out_rows[ridx]["q_by_within_stage"])
                and out_rows[ridx]["q_by_within_stage"] <= alpha
                and np.isfinite(out_rows[ridx]["median_diff"])
                and out_rows[ridx]["median_diff"] > 0.0
            )
    stage_region = pd.DataFrame(out_rows)
    q_global, _ = fdr_adjust(stage_region["p_two_sided"].tolist(), method="by")
    stage_region["q_by_global"] = q_global
    stage_region["sig_by_global"] = np.isfinite(q_global) & (q_global <= alpha)
    stage_region["contributes_above_null_global"] = (
        stage_region["sig_by_global"] & np.isfinite(stage_region["median_diff"]) & (stage_region["median_diff"] > 0.0)
    )
    return agg, stage_region


def run_band_specific_tests(df: pd.DataFrame, alpha: float) -> pd.DataFrame:
    out_rows: list[dict[str, Any]] = []
    for st in STAGES:
        for band in BANDS:
            sdf = df[(df["stage_internal"] == st) & (df["band"] == band)]
            regions = sorted(sdf["region"].astype(str).unique().tolist())
            pvals: list[float] = []
            idxs: list[int] = []
            for rg in regions:
                r = sdf[sdf["region"] == rg]
                vals = r["beta_best"].to_numpy(dtype=float)
                null_value = float(r["beta_null"].iloc[0]) if len(r) else np.nan
                stats = wilcoxon_vs_null(vals, null_value)
                out_rows.append(
                    {
                        "stage_internal": st,
                        "stage_display": stage_name(st),
                        "band": band,
                        "region": rg,
                        "m": int(r["m"].iloc[0]) if len(r) else np.nan,
                        "n_regions": int(r["n_regions"].iloc[0]) if len(r) else np.nan,
                        "beta_null": null_value,
                        **stats,
                    }
                )
                pvals.append(stats["p_two_sided"])
                idxs.append(len(out_rows) - 1)
            q_by, _ = fdr_adjust(pvals, method="by")
            for j, ridx in enumerate(idxs):
                out_rows[ridx]["q_by_within_stage_band"] = float(q_by[j]) if np.isfinite(q_by[j]) else np.nan
                out_rows[ridx]["sig_by_within_stage_band"] = bool(np.isfinite(q_by[j]) and q_by[j] <= alpha)
                out_rows[ridx]["contributes_above_null"] = bool(
                    np.isfinite(out_rows[ridx]["q_by_within_stage_band"])
                    and out_rows[ridx]["q_by_within_stage_band"] <= alpha
                    and np.isfinite(out_rows[ridx]["median_diff"])
                    and out_rows[ridx]["median_diff"] > 0.0
                )
    band_region = pd.DataFrame(out_rows)
    q_global, _ = fdr_adjust(band_region["p_two_sided"].tolist(), method="by")
    band_region["q_by_global"] = q_global
    band_region["sig_by_global"] = np.isfinite(q_global) & (q_global <= alpha)
    band_region["contributes_above_null_global"] = (
        band_region["sig_by_global"] & np.isfinite(band_region["median_diff"]) & (band_region["median_diff"] > 0.0)
    )
    return band_region


def run_friedman_omnibus(stage_avg: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for st in STAGES:
        sdf = stage_avg[stage_avg["stage_internal"] == st]
        if sdf.empty:
            continue
        wide = sdf.pivot_table(index="subject_id", columns="region", values="beta_stage_avg", aggfunc="mean")
        wide = wide.dropna(axis=0, how="any")
        if wide.shape[0] < 3 or wide.shape[1] < 3:
            rows.append(
                {
                    "stage_internal": st,
                    "stage_display": stage_name(st),
                    "n_subjects": int(wide.shape[0]),
                    "n_regions": int(wide.shape[1]),
                    "friedman_stat": np.nan,
                    "friedman_p": np.nan,
                }
            )
            continue
        vals = [wide[c].to_numpy(dtype=float) for c in wide.columns]
        stat, p = friedmanchisquare(*vals)
        rows.append(
            {
                "stage_internal": st,
                "stage_display": stage_name(st),
                "n_subjects": int(wide.shape[0]),
                "n_regions": int(wide.shape[1]),
                "friedman_stat": float(stat),
                "friedman_p": float(p),
            }
        )
    omni = pd.DataFrame(rows)
    if len(omni) > 0:
        q, _ = fdr_adjust(omni["friedman_p"].tolist(), method="by")
        omni["q_by_across_stages"] = q
        omni["sig_by_across_stages"] = np.isfinite(q) & (q <= 0.05)
    return omni


def plot_stage_region_heatmap(stage_region: pd.DataFrame, out_path_stem: Path, alpha: float) -> None:
    if stage_region.empty:
        return
    order_stage = [stage_name(s) for s in STAGES]
    order_region = REGION_ORDER
    mat = (
        stage_region.pivot_table(index="stage_display", columns="region", values="median_diff", aggfunc="first")
        .reindex(index=order_stage, columns=order_region)
    )
    sig = (
        stage_region.pivot_table(index="stage_display", columns="region", values="contributes_above_null", aggfunc="first")
        .reindex(index=order_stage, columns=order_region)
        .fillna(False)
    )
    fig, ax = plt.subplots(figsize=(11, 4.8))
    vmax = float(np.nanmax(np.abs(mat.to_numpy(dtype=float)))) if np.isfinite(np.nanmax(np.abs(mat.to_numpy(dtype=float)))) else 1.0
    vmax = max(vmax, 1e-6)
    sns.heatmap(
        mat,
        ax=ax,
        cmap="coolwarm",
        center=0.0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Median(beta - 3/7)"},
    )
    for i, st in enumerate(mat.index):
        for j, rg in enumerate(mat.columns):
            if bool(sig.loc[st, rg]):
                ax.text(j + 0.5, i + 0.5, "*", ha="center", va="center", color="black", fontsize=13, fontweight="bold")
    ax.set_title(f"m=3 beta contributions by stage (star = BY q<{alpha:g} and beta>3/7)")
    ax.set_xlabel("Region")
    ax.set_ylabel("Stage")
    fig.tight_layout()
    fig.savefig(out_path_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_frequency_specific_stage_panels(band_region: pd.DataFrame, out_path_stem: Path, alpha: float) -> None:
    if band_region.empty:
        return
    order_stage = [stage_name(s) for s in STAGES]
    vmax = np.nanmax(np.abs(band_region["median_diff"].to_numpy(dtype=float)))
    if not np.isfinite(vmax):
        vmax = 1.0
    vmax = float(max(vmax, 1e-6))

    fig, axes = plt.subplots(1, len(order_stage), figsize=(22, 5.8), sharey=True)
    cbar_ax = fig.add_axes([0.92, 0.18, 0.012, 0.66])
    for i, st in enumerate(STAGES):
        ax = axes[i]
        sdf = band_region[band_region["stage_internal"] == st]
        mat = (
            sdf.pivot_table(index="band", columns="region", values="median_diff", aggfunc="first")
            .reindex(index=BANDS, columns=REGION_ORDER)
        )
        sig = (
            sdf.pivot_table(index="band", columns="region", values="contributes_above_null", aggfunc="first")
            .reindex(index=BANDS, columns=REGION_ORDER)
            .fillna(False)
        )
        sns.heatmap(
            mat,
            ax=ax,
            cmap="coolwarm",
            center=0.0,
            vmin=-vmax,
            vmax=vmax,
            linewidths=0.4,
            linecolor="white",
            cbar=(i == len(order_stage) - 1),
            cbar_ax=cbar_ax if i == len(order_stage) - 1 else None,
            cbar_kws={"label": "Median(beta - 3/7)"},
        )
        for r in range(mat.shape[0]):
            for c in range(mat.shape[1]):
                if bool(sig.iloc[r, c]):
                    ax.text(c + 0.5, r + 0.5, "*", ha="center", va="center", color="black", fontsize=11, fontweight="bold")
        ax.set_title(order_stage[i], fontsize=13, fontweight="bold")
        ax.set_xlabel("Region", fontsize=10)
        if i == 0:
            ax.set_ylabel("Band", fontsize=10)
        else:
            ax.set_ylabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)
        ax.set_yticklabels([BAND_DISPLAY.get(b, b) for b in BANDS], rotation=0, fontsize=11)
    fig.suptitle(f"m=3 frequency-specific beta contributions (star = BY q<{alpha:g} and beta>3/7)", fontsize=14, y=0.98)
    fig.tight_layout(rect=[0.01, 0.01, 0.91, 0.95])
    fig.savefig(out_path_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_path_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_report(
    stage_region: pd.DataFrame,
    band_region: pd.DataFrame,
    friedman: pd.DataFrame,
    out_txt: Path,
    alpha: float,
    m: int,
) -> None:
    with out_txt.open("w", encoding="utf-8") as f:
        f.write("SomnoDynamics beta contribution analysis (m=3)\n")
        f.write("============================================\n\n")
        f.write(f"Primary null: beta == m/n_regions = {m}/7 = {m/7:.6f}\n")
        f.write("Primary test per region/condition: one-sample Wilcoxon signed-rank (two-sided), evaluated against null.\n")
        f.write("Correction: Benjamini-Yekutieli (BY), within each condition family; global BY also reported.\n\n")
        f.write("Stage-level omnibus (band-averaged beta across subjects): Friedman test across regions\n")
        if friedman.empty:
            f.write("- No omnibus rows available.\n")
        else:
            for _, r in friedman.iterrows():
                f.write(
                    f"- {r['stage_display']}: n={int(r['n_subjects'])}, chi2={r['friedman_stat']:.4g}, "
                    f"p={r['friedman_p']:.4g}, q_BY={r.get('q_by_across_stages', np.nan):.4g}\n"
                )
        f.write("\nSignificant contributing regions by stage (band-averaged)\n")
        for st in STAGES:
            sdf = stage_region[(stage_region["stage_internal"] == st) & (stage_region["contributes_above_null"])]
            sdf = sdf.sort_values("median_diff", ascending=False)
            if len(sdf) == 0:
                f.write(f"- {stage_name(st)}: none (BY q<{alpha:g})\n")
                continue
            parts = [f"{r['region']} (median_diff={r['median_diff']:.4f}, q={r['q_by_within_stage']:.4g})" for _, r in sdf.iterrows()]
            f.write(f"- {stage_name(st)}: " + "; ".join(parts) + "\n")
        f.write("\nBand-specific rows are available in stats CSV (stage x band x region).\n")
        f.write(f"Total band-specific tests: {len(band_region)}\n")
        sig_band = band_region[band_region["contributes_above_null"]]
        f.write(f"Band-specific significant rows (BY q<{alpha:g} within stage-band): {len(sig_band)}\n")
        if len(sig_band) > 0:
            f.write("Top band-specific positives:\n")
            top = sig_band.sort_values("median_diff", ascending=False).head(15)
            for _, r in top.iterrows():
                f.write(
                    f"- {r['stage_display']} | {BAND_DISPLAY.get(str(r['band']), str(r['band']))} | {r['region']}: "
                    f"median_diff={r['median_diff']:.4f}, q={r['q_by_within_stage_band']:.4g}\n"
                )


def main() -> int:
    args = parse_args()
    out_root = args.out_root if args.out_root is not None else (args.spectral_dir / "beta_macro3_analysis")
    stats_dir = out_root / "stats"
    fig_dir = out_root / "figures"
    rep_dir = out_root / "reports"
    stats_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    rep_dir.mkdir(parents=True, exist_ok=True)

    beta_long = extract_beta_rows(args.spectral_dir, args.m)
    if beta_long.empty:
        raise RuntimeError(f"No beta rows extracted from {args.spectral_dir}")

    stage_avg, stage_region = run_stage_collapsed_tests(beta_long, alpha=args.alpha)
    band_region = run_band_specific_tests(beta_long, alpha=args.alpha)
    friedman = run_friedman_omnibus(stage_avg)

    beta_long.to_csv(stats_dir / "beta_m3_subject_stage_band_region_long.csv", index=False)
    stage_avg.to_csv(stats_dir / "beta_m3_subject_stage_region_bandavg_long.csv", index=False)
    stage_region.to_csv(stats_dir / "beta_m3_stage_region_wilcoxon.csv", index=False)
    band_region.to_csv(stats_dir / "beta_m3_stage_band_region_wilcoxon.csv", index=False)
    friedman.to_csv(stats_dir / "beta_m3_stage_region_friedman.csv", index=False)

    plot_stage_region_heatmap(stage_region, fig_dir / "beta_m3_stage_region_median_diff_heatmap", alpha=args.alpha)
    plot_frequency_specific_stage_panels(
        band_region,
        fig_dir / "beta_m3_stage_band_region_median_diff_composite",
        alpha=args.alpha,
    )
    write_report(stage_region, band_region, friedman, rep_dir / "beta_m3_analysis_summary.txt", args.alpha, args.m)

    print(f"[SomnoDynamics] beta_m3 analysis complete -> {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
