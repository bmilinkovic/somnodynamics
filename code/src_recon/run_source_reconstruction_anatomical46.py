#!/usr/bin/env python3
"""Source reconstruction: anatomical46 first, then anatomical7 coarse-graining.

Implemented behavior (memory-safe + methodologically aligned):
- Use all stage data to build one common LCMV filter per subject.
- Reconstruct source time series from sensor epochs first (10s source chunks).
- Segment reconstructed source signals into non-overlapping 2-second trials.
- Enforce strict expected trial counts from raw sleep-stage durations.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.io import savemat
    SCIPY_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover
    savemat = None
    SCIPY_IMPORT_ERROR = exc

try:
    import mne
    from mne.beamformer import apply_lcmv_epochs, make_lcmv
    MNE_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover
    mne = None
    apply_lcmv_epochs = None
    make_lcmv = None
    MNE_IMPORT_ERROR = exc

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import canonical_stage, load_config, resolve_repo_path, stage_display, stage_order  # noqa: E402


MACRO_ORDER = [
    "Frontal",
    "Sensorimotor",
    "Parietal",
    "Temporal",
    "Occipital",
    "Cingulate",
    "Limbic",
]


def is_unknown_label(name: str) -> bool:
    s = str(name).strip()
    return s in {"???-lh", "???-rh"} or s.startswith("???")


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Run anatomical46->anatomical7 source reconstruction")
    parser.add_argument(
        "--sensor-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["preprocessed_sensor_space"]),
        help="Directory containing preprocessed sensor epochs.",
    )
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=Path(cfg["raw_data_root"]),
        help="Raw root containing subject TXT stage-annotation files.",
    )
    parser.add_argument(
        "--out-fine-root",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["source_recon_anatomical46"]),
    )
    parser.add_argument(
        "--out-coarse-root",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["source_recon_anatomical7"]),
    )
    parser.add_argument(
        "--mapping-csv",
        type=Path,
        default=resolve_repo_path("code/coarse_graining/anatomical7_mapping_from_dmt.csv"),
    )
    parser.add_argument("--subjects-dir", type=Path, default=Path("/Users/borjan/mne_data/MNE-fsaverage-data"))
    parser.add_argument("--subject", action="append", default=[])
    parser.add_argument("--batch-size", type=int, default=8, help="Sensor-epoch batch size for source extraction.")
    parser.add_argument("--save-csv", action="store_true", help="Also save large CSV time-series exports.")
    parser.add_argument(
        "--target-sfreq",
        type=float,
        default=256.0,
        help="Target source sampling frequency (Hz). Must be >= 256 for legacy-compatible downstream runs.",
    )
    parser.add_argument("--target-epoch-seconds", type=float, default=2.0, help="Target source epoch length (seconds).")
    parser.add_argument(
        "--no-strict-count-check",
        action="store_true",
        help="Allow mismatches between expected and reconstructed 2-second trial counts.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def parse_epochs_file(path: Path) -> tuple[str, str] | None:
    m = re.match(r"^(EPCTL\d{2})_(W|W_|N1|N2|N3|R|R_)_epochs\.fif$", path.name)
    if not m:
        return None
    return m.group(1), canonical_stage(m.group(2))


def discover_subjects(sensor_dir: Path, include: list[str]) -> list[str]:
    include_set = {s.strip() for s in include if s.strip()}
    subjects = sorted(
        p.name for p in sensor_dir.glob("EPCTL*") if p.is_dir() and re.match(r"EPCTL\d{2}$", p.name)
    )
    if not subjects:
        subjects = sorted(
            {
                m.group(1)
                for p in sensor_dir.glob("EPCTL*_epochs.fif")
                if (m := re.match(r"^(EPCTL\d{2})_", p.name))
            }
        )
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


def read_epochs_resilient(path: Path, preload: bool = False) -> mne.Epochs:
    try:
        return mne.read_epochs(path, preload=preload, verbose="ERROR")
    except OSError:
        cache_dir = resolve_repo_path("data/interim/epoch_read_cache")
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = cache_dir / path.name
        shutil.copy2(path.resolve(), cache_file)
        return mne.read_epochs(cache_file, preload=preload, verbose="ERROR")


def build_mapping(mapping_csv: Path) -> dict[str, str]:
    df = pd.read_csv(mapping_csv)
    return dict(zip(df["label"].astype(str), df["macro_region"].astype(str)))


def aggregate_anatomical7(fine_ts: np.ndarray, labels: list[str], label_to_macro: dict[str, str]) -> tuple[np.ndarray, dict[str, int]]:
    coverage = {m: 0 for m in MACRO_ORDER}
    macro = np.full((len(MACRO_ORDER), fine_ts.shape[1], fine_ts.shape[2]), np.nan, dtype=float)

    for i, macro_name in enumerate(MACRO_ORDER):
        idx = [j for j, lab in enumerate(labels) if label_to_macro.get(lab) == macro_name]
        coverage[macro_name] = len(idx)
        if not idx:
            raise ValueError(f"No labels mapped to macro-region: {macro_name}")
        macro[i, :, :] = np.nanmean(fine_ts[idx, :, :], axis=0)

    return macro, coverage


def save_timeseries_csv(path: Path, ts: np.ndarray, labels: list[str]) -> None:
    flat = ts.transpose(2, 1, 0).reshape(-1, ts.shape[0])
    pd.DataFrame(flat, columns=labels).to_csv(path, index=False)


def load_expected_windows(raw_root: Path, subject_id: str, epoch_seconds: float) -> dict[str, dict[str, float]]:
    txt = raw_root / subject_id / f"{subject_id}.txt"
    if not txt.exists():
        raise FileNotFoundError(f"Missing annotation TXT for {subject_id}: {txt}")

    df = pd.read_csv(txt, sep="\t", header=None, names=["label", "start_time", "duration"])
    out = {
        stage: {"duration_s": 0.0, "row_count": 0, "expected_windows": 0}
        for stage in stage_order()
    }
    ignored = 0
    for _, row in df.iterrows():
        raw_label = str(row["label"]).strip()
        try:
            stage = canonical_stage(raw_label)
        except ValueError:
            ignored += 1
            continue
        out[stage]["duration_s"] += float(row["duration"])
        out[stage]["row_count"] += 1

    for stage in stage_order():
        out[stage]["expected_windows"] = int(out[stage]["duration_s"] // epoch_seconds)
    out["_ignored_rows"] = {"count": ignored}
    return out


def build_common_filter_from_epochs(
    stage_files: dict[str, Path],
    src,
    bem,
    batch_size: int,
):
    info_ref: mne.Info | None = None
    sum_x: np.ndarray | None = None
    sum_xx: np.ndarray | None = None
    n_samples = 0

    for stage in stage_order():
        file = stage_files.get(stage)
        if file is None:
            continue
        epochs = read_epochs_resilient(file, preload=False)
        if info_ref is None:
            info_ref = epochs.info.copy()

        for start in range(0, len(epochs), batch_size):
            stop = min(start + batch_size, len(epochs))
            batch = epochs[start:stop]
            data = batch.get_data()  # (ep, ch, t)
            x = np.transpose(data, (1, 0, 2)).reshape(data.shape[1], -1).T.astype(np.float64, copy=False)
            if sum_x is None:
                sum_x = np.zeros(x.shape[1], dtype=np.float64)
                sum_xx = np.zeros((x.shape[1], x.shape[1]), dtype=np.float64)
            sum_x += x.sum(axis=0)
            sum_xx += x.T @ x
            n_samples += x.shape[0]

    if info_ref is None or sum_x is None or sum_xx is None or n_samples < 2:
        raise RuntimeError("No sensor samples available to build common filter")

    cov_mat = (sum_xx - np.outer(sum_x, sum_x) / float(n_samples)) / float(max(n_samples - 1, 1))
    data_cov = mne.Covariance(
        data=cov_mat,
        names=info_ref["ch_names"],
        bads=info_ref["bads"],
        projs=info_ref["projs"],
        nfree=max(n_samples - 1, 1),
    )

    fwd = mne.make_forward_solution(info_ref, trans="fsaverage", src=src, bem=bem, eeg=True, mindist=5.0)
    noise_cov = mne.make_ad_hoc_cov(info_ref, std=1e-5)
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


def decimate_signal(data: np.ndarray, sfreq_in: float, sfreq_out: float) -> np.ndarray:
    ratio = sfreq_in / sfreq_out
    if not np.isclose(ratio, round(ratio)):
        raise ValueError(f"Non-integer resample ratio not supported: {sfreq_in} -> {sfreq_out}")
    step = int(round(ratio))
    if step <= 0:
        raise ValueError(f"Invalid decimation step from {sfreq_in} -> {sfreq_out}")
    if step == 1:
        return data
    return data[:, ::step]


def segment_single_epoch(label_ts: np.ndarray, sfreq: float, epoch_seconds: float) -> np.ndarray:
    samples = int(round(sfreq * epoch_seconds))
    if samples < 1:
        raise ValueError(f"Invalid target epoch samples: {samples}")
    n_seg = label_ts.shape[1] // samples
    if n_seg < 1:
        return np.empty((label_ts.shape[0], samples, 0), dtype=np.float32)
    trimmed = label_ts[:, : n_seg * samples]
    return trimmed.reshape(label_ts.shape[0], n_seg, samples).transpose(0, 2, 1).astype(np.float32, copy=False)


def extract_stage_fine_ts(
    stage_file: Path,
    filters,
    labels_all,
    src,
    target_sfreq: float,
    target_epoch_seconds: float,
    batch_size: int,
) -> tuple[np.ndarray, int, float]:
    epochs = read_epochs_resilient(stage_file, preload=False)
    sfreq_in = float(epochs.info["sfreq"])
    pieces: list[np.ndarray] = []

    for start in range(0, len(epochs), batch_size):
        stop = min(start + batch_size, len(epochs))
        batch = epochs[start:stop]
        stc_gen = apply_lcmv_epochs(batch, filters, return_generator=True, verbose="ERROR")
        for stc in stc_gen:
            label_ts = mne.extract_label_time_course(
                stc,
                labels_all,
                src,
                mode="mean_flip",
                return_generator=False,
                verbose="ERROR",
            )
            label_ts = np.asarray(label_ts, dtype=np.float32)
            label_ts = decimate_signal(label_ts, sfreq_in, target_sfreq)
            seg = segment_single_epoch(label_ts, target_sfreq, target_epoch_seconds)
            if seg.shape[2] > 0:
                pieces.append(seg)

    if not pieces:
        raise RuntimeError(f"No source segments generated from {stage_file}")

    fine_ts = np.concatenate(pieces, axis=2)
    return fine_ts, len(epochs), sfreq_in


def run_subject(
    subject_id: str,
    args: argparse.Namespace,
    src,
    bem,
    labels_all,
    label_names: list[str],
    label_to_macro: dict[str, str],
) -> list[dict[str, object]]:
    stage_files = collect_stage_files(args.sensor_dir, subject_id)
    if not stage_files:
        print(f"[WARN] no epochs found for {subject_id}")
        return []

    expected = load_expected_windows(args.raw_root, subject_id, float(args.target_epoch_seconds))
    strict = not args.no_strict_count_check

    if strict:
        for stage in stage_order():
            expected_windows = int(expected[stage]["expected_windows"])
            if expected_windows > 0 and stage not in stage_files:
                raise RuntimeError(
                    f"{subject_id} {stage}: expected {expected_windows} windows from raw durations, "
                    "but no preprocessed sensor epochs file was found."
                )

    filters = build_common_filter_from_epochs(stage_files, src, bem, batch_size=max(1, args.batch_size))

    rows: list[dict[str, object]] = []
    for stage in stage_order():
        stage_file = stage_files.get(stage)
        if stage_file is None:
            continue

        fine_mat = args.out_fine_root / subject_id / "mat" / f"{subject_id}_{stage}_source_time_series_46.mat"
        fine_csv = args.out_fine_root / subject_id / "csv" / f"{subject_id}_{stage}_source_time_series_46.csv"
        coarse_mat = args.out_coarse_root / subject_id / "mat" / f"{subject_id}_{stage}_source_time_series.mat"
        coarse_csv = args.out_coarse_root / subject_id / "csv" / f"{subject_id}_{stage}_source_time_series.csv"

        if fine_mat.exists() and coarse_mat.exists() and not args.overwrite:
            print(f"[skip] {subject_id} {stage}")
            continue

        fine_ts, input_epochs, input_sfreq = extract_stage_fine_ts(
            stage_file,
            filters,
            labels_all,
            src,
            target_sfreq=float(args.target_sfreq),
            target_epoch_seconds=float(args.target_epoch_seconds),
            batch_size=max(1, args.batch_size),
        )

        expected_windows = int(expected[stage]["expected_windows"])
        actual_windows = int(fine_ts.shape[2])
        if strict and actual_windows != expected_windows:
            raise RuntimeError(
                f"{subject_id} {stage}: expected {expected_windows} non-overlap windows "
                f"from raw durations, got {actual_windows}."
            )

        coarse_ts, coverage = aggregate_anatomical7(fine_ts, label_names, label_to_macro)
        coarse_ts = coarse_ts.astype(np.float32, copy=False)

        fine_mat.parent.mkdir(parents=True, exist_ok=True)
        fine_csv.parent.mkdir(parents=True, exist_ok=True)
        coarse_mat.parent.mkdir(parents=True, exist_ok=True)
        coarse_csv.parent.mkdir(parents=True, exist_ok=True)

        savemat(fine_mat, {"source_ts": fine_ts, "labels": np.asarray(label_names, dtype=object), "atlas": "hcpmmp1_combined_46"})
        if args.save_csv:
            save_timeseries_csv(fine_csv, fine_ts, label_names)
        savemat(coarse_mat, {"source_ts": coarse_ts, "labels": np.asarray(MACRO_ORDER, dtype=object), "atlas": "anatomical7"})
        if args.save_csv:
            save_timeseries_csv(coarse_csv, coarse_ts, MACRO_ORDER)

        rows.append(
            {
                "subject_id": subject_id,
                "stage_internal": stage,
                "stage_label": stage_display(stage),
                "n_fine_labels": int(fine_ts.shape[0]),
                "n_macro_labels": int(coarse_ts.shape[0]),
                "n_timepoints": int(fine_ts.shape[1]),
                "n_epochs": int(fine_ts.shape[2]),
                "expected_epochs_from_raw": expected_windows,
                "raw_stage_duration_s": float(expected[stage]["duration_s"]),
                "raw_stage_row_count": int(expected[stage]["row_count"]),
                "source_sfreq_hz": float(args.target_sfreq),
                "source_epoch_seconds": float(args.target_epoch_seconds),
                "input_sensor_epochs": int(input_epochs),
                "input_sensor_sfreq_hz": float(input_sfreq),
                "source_input_epochs_file": str(stage_file),
                "filter_type": "common_across_stages",
                "count_check_passed": bool(actual_windows == expected_windows),
                "fine_mat": str(fine_mat),
                "coarse_mat": str(coarse_mat),
                **{f"mapped_{k}": int(v) for k, v in coverage.items()},
            }
        )
        print(
            f"[ok] {subject_id} {stage}: fine={fine_ts.shape} coarse={coarse_ts.shape} "
            f"expected={expected_windows}"
        )

    return rows


def main() -> int:
    args = parse_args()

    if MNE_IMPORT_ERROR is not None:
        raise SystemExit(f"MNE is required for source reconstruction (error: {MNE_IMPORT_ERROR})")
    if SCIPY_IMPORT_ERROR is not None:
        raise SystemExit(f"scipy is required for source reconstruction (error: {SCIPY_IMPORT_ERROR})")
    if not args.sensor_dir.exists():
        raise SystemExit(f"Sensor epochs directory not found: {args.sensor_dir}")
    if not args.raw_root.exists():
        raise SystemExit(f"Raw root not found: {args.raw_root}")
    if float(args.target_sfreq) < 256.0:
        raise SystemExit("--target-sfreq must be >= 256 Hz (no lower than preprocessing standard).")

    fs_dir = args.subjects_dir / "fsaverage"
    if not fs_dir.exists():
        raise SystemExit(
            f"Missing fsaverage at {fs_dir}. "
            "Install/download fsaverage into --subjects-dir before running source reconstruction."
        )

    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    mpl_cfg = resolve_repo_path("data/interim/mplconfig")
    mpl_cfg.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cfg))
    mne.set_log_level("WARNING")

    subjects = discover_subjects(args.sensor_dir, args.subject)
    if not subjects:
        raise SystemExit(f"No subjects found in {args.sensor_dir}")

    args.out_fine_root.mkdir(parents=True, exist_ok=True)
    args.out_coarse_root.mkdir(parents=True, exist_ok=True)

    label_to_macro = build_mapping(args.mapping_csv)
    print(f"[SomnoDynamics] source input epochs root: {args.sensor_dir}")
    print(f"[SomnoDynamics] raw annotation root: {args.raw_root}")

    src = mne.setup_source_space("fsaverage", spacing="ico4", subjects_dir=str(args.subjects_dir), add_dist=False)
    bem_model = mne.make_bem_model(subject="fsaverage", ico=4, subjects_dir=str(args.subjects_dir))
    bem = mne.make_bem_solution(bem_model)

    labels_all_raw = mne.read_labels_from_annot(
        "fsaverage",
        parc="HCPMMP1_combined",
        subjects_dir=str(args.subjects_dir),
    )
    labels_all = [lbl for lbl in labels_all_raw if not is_unknown_label(lbl.name)]
    label_names = [str(l.name).strip() for l in labels_all]

    missing = sorted(set(label_names) - set(label_to_macro.keys()))
    if missing:
        raise SystemExit(
            "Mapping file does not cover all surviving HCP labels. "
            f"Missing: {missing}"
        )
    print(
        f"[SomnoDynamics] HCP labels: total={len(labels_all_raw)} "
        f"filtered={len(labels_all)}"
    )

    all_rows: list[dict[str, object]] = []
    for subject in subjects:
        all_rows.extend(run_subject(subject, args, src, bem, labels_all, label_names, label_to_macro))

    summary_csv = args.out_coarse_root / "source_reconstruction_summary.csv"
    pd.DataFrame(all_rows).to_csv(summary_csv, index=False)
    print(f"[SomnoDynamics] source reconstruction summary -> {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
