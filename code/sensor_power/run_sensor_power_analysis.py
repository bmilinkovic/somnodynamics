#!/usr/bin/env python3
"""Sensor-space spectral analysis with DMT-style metric computation.

This stage computes:
- bandpower metrics with Welch settings aligned to DMT scripts
- GFP time series and GFP-PSD metrics
- repeated-measures statistics across sleep stages
- publication-ready topomap figures
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon

try:
    from statsmodels.stats.multitest import multipletests
except Exception:  # pragma: no cover
    multipletests = None

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import load_config, resolve_repo_path, stage_display, stage_order  # noqa: E402

try:
    import mne
    MNE_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover
    mne = None  # type: ignore
    MNE_IMPORT_ERROR = exc


BANDS = {
    "delta": (1.0, 4.0),
    "theta": (5.0, 8.0),
    "alpha": (9.0, 12.0),
    "sigma": (13.0, 16.0),
    "beta": (17.0, 30.0),
    "gamma": (31.0, 44.0),
}
BAND_NAMES = list(BANDS.keys())
FMIN_ALL = 1.0
FMAX_ALL = 44.0
WELCH_WIN_SECONDS = 2.0
WELCH_WINDOW = "hann"
EPS = 1e-20


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Run sensor-level power analysis")
    parser.add_argument(
        "--sensor-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["preprocessed_sensor_space"]),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["sensor_power_root"]),
    )
    parser.add_argument("--n-jobs", type=int, default=1, help="Parallel workers across epoch files.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def parse_epochs_file(path: Path) -> tuple[str, str] | None:
    m = re.match(r"^(EPCTL\d{2})_(W_|N1|N2|N3|R_)_epochs\.fif$", path.name)
    if not m:
        return None
    return m.group(1), m.group(2)


def _collect_epoch_files(sensor_dir: Path) -> list[Path]:
    files = sorted(sensor_dir.glob("EPCTL*/EPCTL*_epochs.fif"))
    if not files:
        files = sorted(sensor_dir.glob("EPCTL*_epochs.fif"))
    return files


def _bandpower_sum(psd: np.ndarray, freqs: np.ndarray, lo: float, hi: float) -> np.ndarray:
    idx = np.logical_and(freqs >= lo, freqs <= hi)
    if not np.any(idx):
        return np.zeros(psd.shape[:-1], dtype=float)
    return psd[..., idx].sum(axis=-1)


def _welch_nfft_novlp(sfreq: float, n_times: int) -> tuple[int, int]:
    target = int(round(sfreq * WELCH_WIN_SECONDS))
    n_fft = max(8, min(target, n_times))
    n_overlap = max(0, n_fft // 2)
    return n_fft, n_overlap


def _process_epochs_file(file: Path, out_dir: Path, overwrite: bool) -> list[dict]:
    parsed = parse_epochs_file(file)
    if parsed is None:
        return []
    subject_id, stage = parsed
    stage_label = stage_display(stage)

    subj_stage_dir = out_dir / "per_subject_stage" / subject_id
    subj_stage_dir.mkdir(parents=True, exist_ok=True)
    epoch_csv = subj_stage_dir / f"{subject_id}_{stage}_bandpower.csv"
    chan_csv = subj_stage_dir / f"{subject_id}_{stage}_bandpower_channels.csv"
    abs_npy = subj_stage_dir / f"{subject_id}_{stage}_bandpower_abs.npy"
    rel_npy = subj_stage_dir / f"{subject_id}_{stage}_bandpower_rel.npy"

    gfp_dir = out_dir / "gfp" / subject_id
    gfp_dir.mkdir(parents=True, exist_ok=True)
    gfp_ts_npy = gfp_dir / f"{subject_id}_{stage}_gfp_ts.npy"
    gfp_npy = gfp_dir / f"{subject_id}_{stage}_gfp.npy"
    gfp_csv = gfp_dir / f"{subject_id}_{stage}_gfp.csv"
    gfp_psd_npy = gfp_dir / f"{subject_id}_{stage}_gfp_psd.npy"
    gfp_psd_csv = gfp_dir / f"{subject_id}_{stage}_gfp_psd.csv"

    if chan_csv.exists() and not overwrite:
        return pd.read_csv(chan_csv).to_dict(orient="records")

    epochs = mne.read_epochs(file, preload=True, verbose="ERROR")
    ep_data = epochs.get_data()  # (ep, ch, t)
    sfreq = float(epochs.info["sfreq"])
    n_ep, n_ch, n_times = ep_data.shape
    n_fft, n_overlap = _welch_nfft_novlp(sfreq, n_times)

    psd_ep = []
    for ep in ep_data:
        psd, freqs = mne.time_frequency.psd_array_welch(
            ep,
            sfreq=sfreq,
            fmin=FMIN_ALL,
            fmax=FMAX_ALL,
            n_fft=n_fft,
            n_overlap=n_overlap,
            window=WELCH_WINDOW,
            average="mean",
            verbose="ERROR",
        )
        psd_ep.append(psd)
    psd_ep = np.stack(psd_ep, axis=0)  # (ep, ch, freq)

    abs_bp_ep_ch: dict[str, np.ndarray] = {}
    abs_bp_ep = np.empty((n_ep, len(BAND_NAMES)), dtype=float)
    for i, band in enumerate(BAND_NAMES):
        lo, hi = BANDS[band]
        bp = _bandpower_sum(psd_ep, freqs, lo, hi)  # (ep, ch)
        abs_bp_ep_ch[band] = bp
        abs_bp_ep[:, i] = bp.mean(axis=1)  # epoch x mean-channel

    total_pow_ep = abs_bp_ep.sum(axis=1, keepdims=True) + EPS
    rel_bp_ep = abs_bp_ep / total_pow_ep
    log_bp_ep = 10.0 * np.log(abs_bp_ep + EPS)

    abs_ch = {f"abs_{b}": abs_bp_ep_ch[b].mean(axis=0) for b in BAND_NAMES}
    tot_ch = sum(abs_ch[f"abs_{b}"] for b in BAND_NAMES) + EPS
    rel_ch = {f"rel_{b}": abs_ch[f"abs_{b}"] / tot_ch for b in BAND_NAMES}
    log_ch = {f"log_{b}": 10.0 * np.log(abs_ch[f"abs_{b}"] + EPS) for b in BAND_NAMES}

    ch_names = list(epochs.ch_names)
    channel_table = {
        "subject_id": np.repeat(subject_id, n_ch),
        "stage_internal": np.repeat(stage, n_ch),
        "stage_label": np.repeat(stage_label, n_ch),
        "channel": ch_names,
    }
    channel_table.update(abs_ch)
    channel_table.update(rel_ch)
    channel_table.update(log_ch)
    chan_df = pd.DataFrame(channel_table)
    chan_df.to_csv(chan_csv, index=False)

    np.save(abs_npy, abs_bp_ep.astype(np.float32))
    np.save(rel_npy, rel_bp_ep.astype(np.float32))
    epoch_df = pd.DataFrame(
        np.hstack([np.arange(1, n_ep + 1)[:, None], abs_bp_ep, rel_bp_ep, log_bp_ep]),
        columns=["epoch"] + [f"abs_{b}" for b in BAND_NAMES] + [f"rel_{b}" for b in BAND_NAMES] + [f"log_{b}" for b in BAND_NAMES],
    )
    epoch_df.insert(0, "subject_id", subject_id)
    epoch_df.insert(1, "stage_internal", stage)
    epoch_df.insert(2, "stage_label", stage_label)
    epoch_df.to_csv(epoch_csv, index=False)

    gfp_ts = ep_data.std(axis=1, ddof=0)  # (ep, time)
    gfp_mean = gfp_ts.mean(axis=1)  # (ep,)
    np.save(gfp_ts_npy, gfp_ts.astype(np.float32))
    np.save(gfp_npy, gfp_mean.astype(np.float32))
    pd.DataFrame({"epoch": np.arange(1, n_ep + 1), "gfp": gfp_mean}).to_csv(gfp_csv, index=False)

    gfp_psd_ep = []
    for ep in gfp_ts:
        psd, freqs_g = mne.time_frequency.psd_array_welch(
            ep[None, :],
            sfreq=sfreq,
            fmin=FMIN_ALL,
            fmax=FMAX_ALL,
            n_fft=n_fft,
            n_overlap=n_overlap,
            window=WELCH_WINDOW,
            average="mean",
            verbose="ERROR",
        )
        gfp_psd_ep.append(psd[0])
    gfp_psd_ep = np.stack(gfp_psd_ep, axis=0)
    gfp_psd_mean = gfp_psd_ep.mean(axis=0)
    np.save(gfp_psd_npy, np.vstack([freqs_g, gfp_psd_mean]).astype(np.float32))
    pd.DataFrame({"freq_hz": freqs_g, "psd": gfp_psd_mean}).to_csv(gfp_psd_csv, index=False)

    out_rows: list[dict] = []
    for _, row in chan_df.iterrows():
        for band in BAND_NAMES:
            out_rows.append(
                {
                    "subject_id": row["subject_id"],
                    "stage_internal": row["stage_internal"],
                    "stage_label": row["stage_label"],
                    "channel": row["channel"],
                    "band": band,
                    "abs_power": float(row[f"abs_{band}"]),
                    "rel_power": float(row[f"rel_{band}"]),
                    "log_power": float(row[f"log_{band}"]),
                }
            )
    return out_rows


def load_power_records(sensor_dir: Path, out_dir: Path, n_jobs: int, overwrite: bool) -> tuple[pd.DataFrame, mne.Info]:
    records: list[dict] = []
    epoch_files = _collect_epoch_files(sensor_dir)
    if not epoch_files:
        raise FileNotFoundError(f"No epochs files found in {sensor_dir}")

    info_ref = mne.read_epochs(epoch_files[0], preload=False, verbose="ERROR").info

    if n_jobs <= 1:
        for file in epoch_files:
            records.extend(_process_epochs_file(file, out_dir, overwrite))
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {
                pool.submit(_process_epochs_file, file, out_dir, overwrite): file
                for file in epoch_files
            }
            for fut in as_completed(futures):
                file = futures[fut]
                try:
                    rows = fut.result()
                except Exception as exc:
                    raise RuntimeError(f"Failed to process {file}: {exc}") from exc
                records.extend(rows)

    if not records:
        raise RuntimeError("No valid epochs files were parsed.")
    return pd.DataFrame(records), info_ref


def repeated_measures_stats(df: pd.DataFrame, value_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    friedman_rows: list[dict] = []
    pair_rows: list[dict] = []

    for band in [b for b in BAND_NAMES if b in set(df["band"].unique())]:
        band_df = df[df["band"] == band]
        band_pair_idx = []

        for channel in sorted(band_df["channel"].unique()):
            sub = band_df[band_df["channel"] == channel]
            pivot = (
                sub.pivot_table(index="subject_id", columns="stage_internal", values=value_col, aggfunc="mean")
                .reindex(columns=stage_order())
                .dropna()
            )
            if pivot.shape[0] < 3:
                continue

            values = [pivot[s].to_numpy() for s in stage_order()]
            stat, p = friedmanchisquare(*values)
            friedman_rows.append(
                {
                    "measure": value_col,
                    "band": band,
                    "channel": channel,
                    "n_subjects": int(pivot.shape[0]),
                    "friedman_stat": float(stat),
                    "friedman_p": float(p),
                }
            )

            for s1, s2 in itertools.combinations(stage_order(), 2):
                try:
                    w_stat, w_p = wilcoxon(pivot[s1], pivot[s2], zero_method="wilcox")
                except ValueError:
                    w_stat, w_p = np.nan, np.nan
                row = {
                    "measure": value_col,
                    "band": band,
                    "channel": channel,
                    "stage_a": s1,
                    "stage_b": s2,
                    "n_subjects": int(pivot.shape[0]),
                    "wilcoxon_stat": float(w_stat) if np.isfinite(w_stat) else np.nan,
                    "wilcoxon_p": float(w_p) if np.isfinite(w_p) else np.nan,
                    "mean_diff_b_minus_a": float(np.nanmean(pivot[s2] - pivot[s1])),
                }
                pair_rows.append(row)
                band_pair_idx.append(len(pair_rows) - 1)

        if multipletests is not None and band_pair_idx:
            pvals = np.array([pair_rows[i]["wilcoxon_p"] for i in band_pair_idx], dtype=float)
            valid = np.isfinite(pvals)
            qvals = np.full_like(pvals, np.nan, dtype=float)
            reject = np.zeros_like(pvals, dtype=bool)
            if np.any(valid):
                reject_v, qvals_v, _, _ = multipletests(pvals[valid], alpha=0.05, method="fdr_bh")
                qvals[valid] = qvals_v
                reject[valid] = reject_v
            for loc, i in enumerate(band_pair_idx):
                pair_rows[i]["wilcoxon_q_bh"] = float(qvals[loc]) if np.isfinite(qvals[loc]) else np.nan
                pair_rows[i]["sig_q_bh_0p05"] = bool(reject[loc])
        else:
            for i in band_pair_idx:
                pair_rows[i]["wilcoxon_q_bh"] = np.nan
                pair_rows[i]["sig_q_bh_0p05"] = False

    return pd.DataFrame(friedman_rows), pd.DataFrame(pair_rows)


def make_topomap_figures(df: pd.DataFrame, info: mne.Info, out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    ch_order = info["ch_names"]
    for measure in ["abs_power", "rel_power"]:
        for band in [b for b in BAND_NAMES if b in set(df["band"].unique())]:
            fig, axes = plt.subplots(1, len(stage_order()), figsize=(4.0 * len(stage_order()), 3.8), constrained_layout=True)

            band_df = df[df["band"] == band]
            all_vals = []
            stage_vals: dict[str, np.ndarray] = {}
            for stage in stage_order():
                stage_df = (
                    band_df[band_df["stage_internal"] == stage]
                    .groupby("channel", as_index=False)[measure]
                    .mean()
                )
                stage_df["channel"] = pd.Categorical(stage_df["channel"], categories=ch_order, ordered=True)
                stage_df = stage_df.sort_values("channel")
                vals = stage_df[measure].to_numpy(dtype=float)
                stage_vals[stage] = vals
                all_vals.append(vals)

            all_stack = np.concatenate(all_vals)
            vmin = float(np.nanpercentile(all_stack, 5))
            vmax = float(np.nanpercentile(all_stack, 95))
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
                vmin, vmax = float(np.nanmin(all_stack)), float(np.nanmax(all_stack))

            for ax, stage in zip(np.atleast_1d(axes), stage_order()):
                vals = stage_vals[stage]
                im, _ = mne.viz.plot_topomap(
                    vals,
                    info,
                    axes=ax,
                    show=False,
                    contours=4,
                    cmap="viridis",
                    vlim=(vmin, vmax),
                )
                ax.set_title(stage_display(stage))

            cbar = fig.colorbar(im, ax=np.atleast_1d(axes).ravel().tolist(), shrink=0.8)
            cbar.set_label(measure)
            fig.suptitle(f"{measure} | {band}")

            png = fig_dir / f"{measure}_{band}_topomap.png"
            svg = fig_dir / f"{measure}_{band}_topomap.svg"
            fig.savefig(png, dpi=300)
            fig.savefig(svg)
            plt.close(fig)


def main() -> int:
    args = parse_args()

    if MNE_IMPORT_ERROR is not None:
        raise SystemExit(
            "MNE is required for sensor power analysis. Install dependencies first "
            f"(import error: {MNE_IMPORT_ERROR})."
        )
    if args.n_jobs < 1:
        raise SystemExit("--n-jobs must be >= 1")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    power_df, info = load_power_records(args.sensor_dir, args.out_dir, args.n_jobs, args.overwrite)

    power_csv = args.out_dir / "sensor_bandpower_long.csv"
    power_df.to_csv(power_csv, index=False)

    f_abs, p_abs = repeated_measures_stats(power_df, "abs_power")
    f_rel, p_rel = repeated_measures_stats(power_df, "rel_power")
    f_log, p_log = repeated_measures_stats(power_df, "log_power")

    f_abs.to_csv(args.out_dir / "stats_abs_friedman.csv", index=False)
    p_abs.to_csv(args.out_dir / "stats_abs_pairwise_wilcoxon.csv", index=False)
    f_rel.to_csv(args.out_dir / "stats_rel_friedman.csv", index=False)
    p_rel.to_csv(args.out_dir / "stats_rel_pairwise_wilcoxon.csv", index=False)
    f_log.to_csv(args.out_dir / "stats_log_friedman.csv", index=False)
    p_log.to_csv(args.out_dir / "stats_log_pairwise_wilcoxon.csv", index=False)

    make_topomap_figures(power_df, info, args.out_dir)

    summary = args.out_dir / "summary.txt"
    with summary.open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics Sensor Power Analysis\n")
        handle.write("================================\n")
        handle.write("Computation parity with DMT scripts:\n")
        handle.write("- 02_01_compute_global_power.py\n")
        handle.write("- 02_02_compute_gfp_to_psd.py\n")
        handle.write("- 03_01_compute_bandpower.py\n")
        handle.write("- 03_02_bandpower_to_topoplot.py\n\n")
        handle.write(f"Input: {args.sensor_dir}\n")
        handle.write(f"Rows: {len(power_df)}\n")
        handle.write(f"Output table: {power_csv}\n")
        handle.write("Statistics: Friedman + pairwise Wilcoxon (BH-FDR within each band).\n")

    print(f"[SomnoDynamics] sensor power complete -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
