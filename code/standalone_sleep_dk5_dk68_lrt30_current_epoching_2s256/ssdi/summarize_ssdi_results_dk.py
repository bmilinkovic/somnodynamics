#!/usr/bin/env python3
"""Summarize stage-based SSDI outputs and run repeated-measures stats."""

from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import friedmanchisquare, wilcoxon

try:
    from statsmodels.stats.multitest import multipletests
except Exception:  # pragma: no cover
    multipletests = None

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "run_all"))
from config_utils import load_config, resolve_repo_path, stage_display, stage_order  # noqa: E402

BANDS = ["delta", "theta", "alpha", "sigma", "beta"]


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Summarize SSDI outputs")
    parser.add_argument(
        "--spectral-dir",
        type=Path,
        default=resolve_repo_path("results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=resolve_repo_path("results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz"),
    )
    return parser.parse_args()


def mat_struct_to_dict(obj: Any) -> Any:
    from scipy.io.matlab import mat_struct  # type: ignore

    if isinstance(obj, mat_struct):
        out = {}
        for name in obj._fieldnames:
            out[name] = mat_struct_to_dict(getattr(obj, name))
        return out
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        return [mat_struct_to_dict(v) for v in obj.tolist()]
    return obj


def load_results(path: Path) -> dict[str, Any]:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "results" not in mat:
        raise KeyError(f"'results' missing in {path}")
    res = mat_struct_to_dict(mat["results"])
    if not isinstance(res, dict):
        raise ValueError(f"Results struct in {path} is not dict-like")
    return res


def parse_file(path: Path) -> tuple[str, str] | None:
    # New DMT-faithful layout: <root>/<stage>/<subject>_ssdi_spectral_domain_results.mat
    if path.parent.name in stage_order():
        m = re.match(r"^(EPCTL\d{2})_ssdi_spectral_domain_results\.mat$", path.name)
        if m:
            return m.group(1), path.parent.name

    # Legacy fast layout: <root>/spectral/<subject>_<stage>_ssdi_spectral_domain_results.mat
    m = re.match(r"^(EPCTL\d{2})_(W_|N1|N2|N3|R_)_ssdi_spectral_domain_results\.mat$", path.name)
    if m:
        return m.group(1), m.group(2)
    return None


def extract_rows(spectral_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    files = sorted(spectral_dir.glob("*/*_ssdi_spectral_domain_results.mat"))
    if not files:
        files = sorted(spectral_dir.glob("EPCTL*_ssdi_spectral_domain_results.mat"))

    for file in files:
        parsed = parse_file(file)
        if parsed is None:
            continue
        subject_id, stage = parsed
        res = load_results(file)

        for band in BANDS:
            for m in range(1, 11):
                key = f"{band}_dopto_m{m}"
                if key not in res:
                    continue
                vals = np.asarray(res[key], dtype=float).ravel()
                if vals.size == 0:
                    continue
                rows.append(
                    {
                        "subject_id": subject_id,
                        "stage_internal": stage,
                        "stage_label": stage_display(stage),
                        "band": band,
                        "m": m,
                        "dd_best": float(np.nanmin(vals)),
                        "dd_median": float(np.nanmedian(vals)),
                        "dd_mean": float(np.nanmean(vals)),
                        "n_runs": int(np.sum(np.isfinite(vals))),
                        "nrunsp_cfg": float(np.asarray(res.get("meta", {}).get("nrunsp", np.nan)).ravel()[0]) if isinstance(res.get("meta"), dict) and "nrunsp" in res["meta"] else np.nan,
                        "nitersp_cfg": float(np.asarray(res.get("meta", {}).get("nitersp", np.nan)).ravel()[0]) if isinstance(res.get("meta"), dict) and "nitersp" in res["meta"] else np.nan,
                        "niterso_cfg": float(np.asarray(res.get("meta", {}).get("niterso", np.nan)).ravel()[0]) if isinstance(res.get("meta"), dict) and "niterso" in res["meta"] else np.nan,
                    }
                )

    return pd.DataFrame(rows)


def repeated_stats(df: pd.DataFrame, metric: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    f_rows: list[dict[str, Any]] = []
    p_rows: list[dict[str, Any]] = []

    for band in sorted(df["band"].unique()):
        bdf = df[df["band"] == band]
        band_pair_idx: list[int] = []

        for m in sorted(bdf["m"].unique()):
            sub = bdf[bdf["m"] == m]
            pivot = (
                sub.pivot_table(index="subject_id", columns="stage_internal", values=metric, aggfunc="mean")
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
                    "band": band,
                    "m": int(m),
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
                        "band": band,
                        "m": int(m),
                        "stage_a": a,
                        "stage_b": b,
                        "n_subjects": int(pivot.shape[0]),
                        "wilcoxon_stat": float(w_stat) if np.isfinite(w_stat) else np.nan,
                        "wilcoxon_p": float(w_p) if np.isfinite(w_p) else np.nan,
                        "mean_diff_b_minus_a": float(np.nanmean(pivot[b] - pivot[a])),
                    }
                )
                band_pair_idx.append(len(p_rows) - 1)

        if multipletests is not None and band_pair_idx:
            pvals = np.array([p_rows[i]["wilcoxon_p"] for i in band_pair_idx], dtype=float)
            valid = np.isfinite(pvals)
            qvals = np.full_like(pvals, np.nan, dtype=float)
            reject = np.zeros_like(pvals, dtype=bool)
            if np.any(valid):
                reject_v, qvals_v, _, _ = multipletests(pvals[valid], alpha=0.05, method="fdr_bh")
                qvals[valid] = qvals_v
                reject[valid] = reject_v
            for local_i, row_i in enumerate(band_pair_idx):
                p_rows[row_i]["wilcoxon_q_bh"] = float(qvals[local_i]) if np.isfinite(qvals[local_i]) else np.nan
                p_rows[row_i]["sig_q_bh_0p05"] = bool(reject[local_i])
        else:
            for row_i in band_pair_idx:
                p_rows[row_i]["wilcoxon_q_bh"] = np.nan
                p_rows[row_i]["sig_q_bh_0p05"] = False

    return pd.DataFrame(f_rows), pd.DataFrame(p_rows)


def make_figures(df: pd.DataFrame, out_dir: Path) -> None:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    for band in sorted(df["band"].unique()):
        bdf = df[df["band"] == band]
        stage_mean = (
            bdf.groupby(["m", "stage_internal"], as_index=False)["dd_median"]
            .mean()
            .pivot(index="m", columns="stage_internal", values="dd_median")
            .reindex(columns=stage_order())
        )

        fig, ax = plt.subplots(figsize=(9, 5))
        for stage in stage_order():
            if stage not in stage_mean.columns:
                continue
            ax.plot(stage_mean.index, stage_mean[stage], marker="o", label=stage_display(stage))

        ax.set_title(f"SSDI median DD across macro sizes | {band}")
        ax.set_xlabel("Macroscopic size (m)")
        ax.set_ylabel("DD median")
        ax.legend(frameon=False)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(fig_dir / f"ssdi_dd_median_{band}.png", dpi=300)
        fig.savefig(fig_dir / f"ssdi_dd_median_{band}.svg")
        plt.close(fig)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    long_df = extract_rows(args.spectral_dir)
    if long_df.empty:
        raise SystemExit(f"No SSDI spectral outputs found in {args.spectral_dir}")

    long_df.to_csv(args.out_dir / "ssdi_subject_long.csv", index=False)

    for metric in ["dd_median", "dd_best", "dd_mean"]:
        f_df, p_df = repeated_stats(long_df, metric)
        f_df.to_csv(args.out_dir / f"ssdi_{metric}_friedman.csv", index=False)
        p_df.to_csv(args.out_dir / f"ssdi_{metric}_pairwise_wilcoxon.csv", index=False)

    make_figures(long_df, args.out_dir)

    # Generate publication-style SSDI figures (DMT-like composites).
    plot_script = Path(__file__).resolve().parent / "plot_ssdi_stage_figures.py"
    if plot_script.exists():
        env = dict(os.environ)
        env.setdefault("MPLCONFIGDIR", "/tmp/mpl")
        subprocess.run(
            [
                sys.executable,
                str(plot_script),
                "--long-csv",
                str(args.out_dir / "ssdi_subject_long.csv"),
                "--spectral-dir",
                str(args.spectral_dir),
                "--out-dir",
                str(args.out_dir / "figures"),
            ],
            check=True,
            env=env,
        )

    with (args.out_dir / "ssdi_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics SSDI Summary\n")
        handle.write("==========================\n")
        handle.write(f"Input spectral dir: {args.spectral_dir}\n")
        handle.write(f"Rows: {len(long_df)}\n")
        handle.write("Stats: Friedman + pairwise Wilcoxon with BH-FDR by band.\n")

    print(f"[SomnoDynamics] SSDI summary complete -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
