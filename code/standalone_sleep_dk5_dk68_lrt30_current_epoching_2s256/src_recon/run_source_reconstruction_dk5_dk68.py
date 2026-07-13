#!/usr/bin/env python3
"""Standalone source reconstruction for DK5 and DK68 outputs.

Spec-aligned behavior:
- LCMV-NAI on fsaverage with 3-layer BEM and template trans.
- Data covariance from continuous EDF (`mne.compute_raw_covariance`).
- Noise covariance from ad-hoc white diagonal (`std=1e-5`).
- `weight_norm='nai'`, `pick_ori='max-power'`, `reg=0.05`.
- Source windows: 2-second non-overlapping epochs (current established pipeline), then polyphase resample to target sfreq (default 256 Hz).
- Save DK5 and DK68 `.mat` files with keys `data` and `labels`.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import savemat
from scipy.signal import resample_poly

import mne
from mne.beamformer import apply_lcmv_epochs, make_lcmv

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "run_all"))
from config_utils import canonical_stage, load_config, resolve_repo_path, stage_order  # noqa: E402


EXCLUDE_CHANNELS = {
    "ChEMG1",
    "ChEMG2",
    "RLEG-",
    "RLEG+",
    "LLEG-",
    "LLEG+",
    "EOG1",
    "EOG2",
    "ECG1",
    "ECG2",
    "SO1",
    "SO2",
    "ZY1",
    "ZY2",
}

DK5_PARCELS = {
    "Frontal": [
        "superiorfrontal", "rostralmiddlefrontal", "caudalmiddlefrontal", "parsopercularis",
        "parstriangularis", "parsorbitalis", "lateralorbitofrontal", "medialorbitofrontal", "frontalpole",
    ],
    "Sensorimotor": ["precentral", "postcentral", "paracentral"],
    "Parietal": ["superiorparietal", "inferiorparietal", "supramarginal", "precuneus"],
    "Temporal": ["superiortemporal", "middletemporal", "inferiortemporal", "transversetemporal", "bankssts", "temporalpole", "fusiform"],
    "Occipital": ["lateraloccipital", "lingual", "cuneus", "pericalcarine"],
}


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    p = argparse.ArgumentParser(description="Standalone DK5/DK68 source reconstruction")
    p.add_argument(
        "--sensor-dir",
        type=Path,
        default=resolve_repo_path("results/preprocessing/preprocessed_sensor_space_artifact_removed_256hz"),
    )
    p.add_argument("--raw-root", type=Path, default=Path(cfg["raw_data_root"]))
    p.add_argument(
        "--cov-source",
        choices=["epochs", "raw"],
        default="epochs",
        help="Use preprocessed epochs for LCMV data covariance by default; raw mode follows the slower raw EDF covariance path.",
    )
    p.add_argument("--montage-file", type=Path, default=Path(cfg["montage_file"]))
    p.add_argument("--subjects-dir", type=Path, default=Path("/Users/borjan/mne_data/MNE-fsaverage-data"))
    p.add_argument("--out-root", type=Path, default=resolve_repo_path("results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz"))
    p.add_argument("--subject", action="append", default=[])
    p.add_argument("--target-sfreq", type=float, default=256.0)
    p.add_argument("--epoch-seconds", type=float, default=2.0)
    p.add_argument("--min-epochs", type=int, default=15)
    p.add_argument("--l-freq", type=float, default=1.0)
    p.add_argument("--h-freq", type=float, default=45.0)
    p.add_argument("--notch-freqs", type=float, nargs="+", default=[50.0, 100.0])
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def load_montage(montage_file: Path) -> mne.channels.DigMontage:
    try:
        df = pd.read_csv(montage_file, sep="\t", header=None)
    except Exception:
        df = pd.read_csv(montage_file, sep=r"\s+", header=None, engine="python")
    names = df[1].astype(str).str.upper().tolist()
    xyz = df[[2, 3, 4]].to_numpy(dtype=float)
    xyz[:, [0, 1]] = xyz[:, [1, 0]]
    return mne.channels.make_dig_montage(ch_pos=dict(zip(names, xyz)), coord_frame="head")


def parse_epochs_file(path: Path) -> tuple[str, str] | None:
    m = re.match(r"^(EPCTL\d{2})_(W|W_|N1|N2|N3|R|R_)_epochs\.fif$", path.name)
    if not m:
        return None
    return m.group(1), canonical_stage(m.group(2))


def discover_subjects(sensor_dir: Path, include: list[str]) -> list[str]:
    include_set = {s.strip() for s in include if s.strip()}
    subjects = sorted(p.name for p in sensor_dir.glob("EPCTL*") if p.is_dir() and re.match(r"EPCTL\d{2}$", p.name))
    if not subjects:
        subjects = sorted({m.group(1) for p in sensor_dir.glob("EPCTL*_epochs.fif") if (m := re.match(r"^(EPCTL\d{2})_", p.name))})
    if include_set:
        subjects = [s for s in subjects if s in include_set]
    return subjects


def collect_stage_files(sensor_dir: Path, subject_id: str) -> dict[str, Path]:
    out: dict[str, Path] = {}
    nested = sorted((sensor_dir / subject_id).glob(f"{subject_id}_*_epochs.fif"))
    flat = sorted(sensor_dir.glob(f"{subject_id}_*_epochs.fif")) if not nested else []
    for f in nested + flat:
        parsed = parse_epochs_file(f)
        if parsed is None:
            continue
        _, stage = parsed
        out[stage] = f
    return out


def preprocess_raw_for_cov(raw: mne.io.BaseRaw, montage: mne.channels.DigMontage, ref_ch_names: list[str], args: argparse.Namespace) -> mne.io.BaseRaw:
    raw.notch_filter(freqs=args.notch_freqs, method="spectrum_fit", n_jobs=1, verbose="ERROR")
    raw.filter(args.l_freq, args.h_freq, fir_design="firwin", n_jobs=1, verbose="ERROR")

    raw.rename_channels({ch: ch.replace("-Ref", "").upper() for ch in raw.ch_names if "-Ref" in ch})
    raw.set_montage(montage, on_missing="ignore")

    bad_locs = [ch["ch_name"] for ch in raw.info["chs"] if not np.isfinite(ch["loc"][:3]).all()]
    drop = [ch for ch in raw.ch_names if ch in EXCLUDE_CHANNELS or ch in bad_locs]
    if drop:
        raw.drop_channels(drop)

    keep = [ch for ch in ref_ch_names if ch in raw.ch_names]
    if len(keep) < 10:
        raise RuntimeError(f"Too few overlap channels between raw and epochs info ({len(keep)}).")
    raw.pick(keep)

    raw.set_eeg_reference("average", projection=True)
    raw.apply_proj()
    return raw


def decimate_polyphase(data: np.ndarray, sfreq_in: float, sfreq_out: float) -> np.ndarray:
    if np.isclose(sfreq_in, sfreq_out):
        return data
    if sfreq_out > sfreq_in:
        raise ValueError(f"Upsampling not supported in this pipeline: {sfreq_in} -> {sfreq_out}")
    ratio = sfreq_in / sfreq_out
    q = int(round(ratio))
    if not np.isclose(ratio, q):
        raise ValueError(f"Non-integer downsample ratio: {sfreq_in} -> {sfreq_out}")
    return resample_poly(data, up=1, down=q, axis=1).astype(np.float32, copy=False)


def split_nonoverlap(label_ts: np.ndarray, sfreq: float, epoch_seconds: float) -> np.ndarray:
    samples = int(round(sfreq * epoch_seconds))
    nseg = label_ts.shape[1] // samples
    if nseg < 1:
        return np.empty((label_ts.shape[0], samples, 0), dtype=np.float32)
    trimmed = label_ts[:, : nseg * samples]
    return trimmed.reshape(label_ts.shape[0], nseg, samples).transpose(0, 2, 1).astype(np.float32, copy=False)


def build_dk_labels(subjects_dir: Path):
    labels = mne.read_labels_from_annot("fsaverage", parc="aparc", subjects_dir=str(subjects_dir))
    labels = [lb for lb in labels if not lb.name.startswith("unknown") and not lb.name.startswith("???")]
    labels_sorted = sorted(labels, key=lambda x: x.name)

    by_name = {lb.name: lb for lb in labels_sorted}
    dk5_labels = []
    for macro, parcels in DK5_PARCELS.items():
        macro_label = None
        for hemi in ("lh", "rh"):
            for parcel in parcels:
                key = f"{parcel}-{hemi}"
                if key not in by_name:
                    continue
                macro_label = by_name[key] if macro_label is None else (macro_label + by_name[key])
        if macro_label is None:
            raise RuntimeError(f"Could not build DK5 macro label for {macro}")
        macro_label.name = macro
        dk5_labels.append(macro_label)

    return labels_sorted, dk5_labels


def compute_epoch_data_cov(stage_files: dict[str, Path], batch_size: int) -> tuple[mne.Covariance, mne.Info]:
    info_ref: mne.Info | None = None
    sum_x: np.ndarray | None = None
    sum_xx: np.ndarray | None = None
    n_samples = 0

    for stage in stage_order():
        stage_file = stage_files.get(stage)
        if stage_file is None:
            continue
        epochs = mne.read_epochs(stage_file, preload=False, verbose="ERROR")
        if info_ref is None:
            info_ref = epochs.info.copy()
        for start in range(0, len(epochs), max(1, batch_size)):
            stop = min(start + max(1, batch_size), len(epochs))
            data = epochs[start:stop].get_data()
            x = np.transpose(data, (1, 0, 2)).reshape(data.shape[1], -1).T.astype(np.float64, copy=False)
            if sum_x is None:
                sum_x = np.zeros(x.shape[1], dtype=np.float64)
                sum_xx = np.zeros((x.shape[1], x.shape[1]), dtype=np.float64)
            sum_x += x.sum(axis=0)
            sum_xx += x.T @ x
            n_samples += x.shape[0]

    if info_ref is None or sum_x is None or sum_xx is None or n_samples < 2:
        raise RuntimeError("No sensor samples available for epoch covariance")

    cov_mat = (sum_xx - np.outer(sum_x, sum_x) / float(n_samples)) / float(max(n_samples - 1, 1))
    data_cov = mne.Covariance(
        data=cov_mat,
        names=info_ref["ch_names"],
        bads=info_ref["bads"],
        projs=info_ref["projs"],
        nfree=max(n_samples - 1, 1),
    )
    return data_cov, info_ref


def build_common_filter(subject_id: str, stage_files: dict[str, Path], src, bem, args: argparse.Namespace):
    first_stage = next(iter(stage_files.values()))
    e0 = mne.read_epochs(first_stage, preload=False, verbose="ERROR")
    info_ref = e0.info.copy()

    if args.cov_source == "epochs":
        data_cov, info_ref = compute_epoch_data_cov(stage_files, args.batch_size)
    else:
        raw_file = args.raw_root / subject_id / f"{subject_id}.edf"
        if not raw_file.exists():
            raise FileNotFoundError(f"Missing raw EDF for covariance: {raw_file}")

        montage = load_montage(args.montage_file)
        raw = mne.io.read_raw_edf(raw_file, preload=True, verbose="ERROR")
        raw = preprocess_raw_for_cov(raw, montage, info_ref["ch_names"], args)
        data_cov = mne.compute_raw_covariance(raw, method="empirical", rank="info", verbose="ERROR")
        raw.close()

    noise_cov = mne.make_ad_hoc_cov(info_ref, std=1e-5)
    fwd = mne.make_forward_solution(info_ref, trans="fsaverage", src=src, bem=bem, eeg=True, mindist=5.0)

    filters = make_lcmv(
        info_ref,
        fwd,
        data_cov,
        reg=0.05,
        noise_cov=noise_cov,
        pick_ori="max-power",
        weight_norm="nai",
        rank="info",
    )
    return filters


def reconstruct_stage(stage_file: Path, filters, labels68, labels5, src, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, float]:
    epochs = mne.read_epochs(stage_file, preload=False, verbose="ERROR")
    sfreq_in = float(epochs.info["sfreq"])

    p68: list[np.ndarray] = []
    p5: list[np.ndarray] = []

    for start in range(0, len(epochs), max(1, args.batch_size)):
        stop = min(start + max(1, args.batch_size), len(epochs))
        batch = epochs[start:stop]
        stcs = apply_lcmv_epochs(batch, filters, return_generator=True, verbose="ERROR")
        for stc in stcs:
            ts68 = np.asarray(mne.extract_label_time_course(stc, labels68, src, mode="mean_flip", return_generator=False, verbose="ERROR"), dtype=np.float32)
            ts5 = np.asarray(mne.extract_label_time_course(stc, labels5, src, mode="mean_flip", return_generator=False, verbose="ERROR"), dtype=np.float32)

            seg68 = split_nonoverlap(ts68, sfreq_in, args.epoch_seconds)
            seg5 = split_nonoverlap(ts5, sfreq_in, args.epoch_seconds)
            if seg68.shape[2] == 0:
                continue

            if not np.isclose(sfreq_in, args.target_sfreq):
                seg68 = np.stack([decimate_polyphase(seg68[:, :, i], sfreq_in, args.target_sfreq) for i in range(seg68.shape[2])], axis=2)
                seg5 = np.stack([decimate_polyphase(seg5[:, :, i], sfreq_in, args.target_sfreq) for i in range(seg5.shape[2])], axis=2)

            p68.append(seg68)
            p5.append(seg5)

    if not p68:
        return np.empty((68, 0, 0), dtype=np.float32), np.empty((5, 0, 0), dtype=np.float32), sfreq_in

    dk68 = np.concatenate(p68, axis=2)
    dk5 = np.concatenate(p5, axis=2)
    return dk68, dk5, sfreq_in


def save_stage_outputs(subject: str, stage: str, dk68: np.ndarray, dk5: np.ndarray, labels68: list[str], labels5: list[str], out_root: Path) -> tuple[Path, Path]:
    out68 = out_root / "dk68" / subject / "mat"
    out5 = out_root / "dk5" / subject / "mat"
    out68.mkdir(parents=True, exist_ok=True)
    out5.mkdir(parents=True, exist_ok=True)

    f68 = out68 / f"{subject}_{stage}_source_time_series_68.mat"
    f5 = out5 / f"{subject}_{stage}_source_time_series.mat"

    savemat(f68, {"data": dk68, "labels": np.asarray(labels68, dtype=object)})
    savemat(f5, {"data": dk5, "labels": np.asarray(labels5, dtype=object)})
    return f68, f5


def main() -> int:
    args = parse_args()
    if args.target_sfreq < 256.0:
        raise SystemExit("--target-sfreq must be >= 256 Hz")

    subjects = discover_subjects(args.sensor_dir, args.subject)
    if not subjects:
        raise SystemExit(f"No subjects found in {args.sensor_dir}")

    src = mne.setup_source_space("fsaverage", spacing="ico4", subjects_dir=str(args.subjects_dir), add_dist=False)
    bem = mne.make_bem_solution(mne.make_bem_model(subject="fsaverage", ico=4, subjects_dir=str(args.subjects_dir)))
    labels68_obj, labels5_obj = build_dk_labels(args.subjects_dir)
    labels68 = [lb.name for lb in labels68_obj]
    labels5 = [lb.name for lb in labels5_obj]

    rows = []
    for subject in subjects:
        stage_files = collect_stage_files(args.sensor_dir, subject)
        if not stage_files:
            print(f"[WARN] no stage files for {subject}")
            continue

        filters = build_common_filter(subject, stage_files, src, bem, args)

        for stage in stage_order():
            stage_file = stage_files.get(stage)
            if stage_file is None:
                continue
            f5 = args.out_root / "dk5" / subject / "mat" / f"{subject}_{stage}_source_time_series.mat"
            f68 = args.out_root / "dk68" / subject / "mat" / f"{subject}_{stage}_source_time_series_68.mat"
            if not args.overwrite and f5.exists() and f68.exists():
                print(f"[skip] {subject} {stage}")
                continue

            dk68, dk5, sfreq_in = reconstruct_stage(stage_file, filters, labels68_obj, labels5_obj, src, args)
            if dk5.shape[2] < args.min_epochs:
                print(f"[WARN] {subject} {stage}: {dk5.shape[2]} epochs < min {args.min_epochs}; skipping save")
                continue

            out68, out5 = save_stage_outputs(subject, stage, dk68, dk5, labels68, labels5, args.out_root)
            rows.append(
                {
                    "subject_id": subject,
                    "stage": stage,
                    "sfreq_in": sfreq_in,
                    "sfreq_out": args.target_sfreq,
                    "epoch_seconds": args.epoch_seconds,
                    "dk5_shape": str(tuple(dk5.shape)),
                    "dk68_shape": str(tuple(dk68.shape)),
                    "dk5_file": str(out5),
                    "dk68_file": str(out68),
                }
            )
            print(f"[ok] {subject} {stage}: dk5={dk5.shape}, dk68={dk68.shape}")

    args.out_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out_root / "source_recon_summary.csv", index=False)
    print(f"[done] summary -> {args.out_root / 'source_recon_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
