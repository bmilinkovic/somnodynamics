#!/usr/bin/env python3
"""Sensor-space preprocessing for SomnoDynamics.

Reads raw EDF/TXT files from the configured dataset path and writes stage-specific
MNE Epochs files under `results/preprocessing/preprocessed_sensor_space`.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import canonical_stage, load_config, resolve_repo_path, stage_order  # noqa: E402

try:
    import mne
    MNE_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover
    mne = None  # type: ignore
    MNE_IMPORT_ERROR = exc


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

LEGACY_STAGE_FILENAME = {
    "W_": "W",
    "N1": "N1",
    "N2": "N2",
    "N3": "N3",
    "R_": "R",
}


@dataclass
class SubjectSummary:
    subject_id: str
    stage: str
    n_epochs: int
    out_file: str


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Preprocess sensor-space sleep EEG")
    parser.add_argument("--raw-root", type=Path, default=Path(cfg["raw_data_root"]))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["preprocessed_sensor_space"]),
    )
    parser.add_argument("--montage-file", type=Path, default=Path(cfg["montage_file"]))
    parser.add_argument("--subject", action="append", default=[])
    parser.add_argument("--epoch-seconds", type=float, default=10.0)
    parser.add_argument("--sfreq", type=float, default=256.0)
    parser.add_argument("--l-freq", type=float, default=1.0)
    parser.add_argument("--h-freq", type=float, default=45.0)
    parser.add_argument(
        "--notch-freqs",
        type=float,
        nargs="+",
        default=[50.0, 100.0],
        help="Notch frequencies in Hz.",
    )
    parser.add_argument(
        "--filter-n-jobs",
        type=int,
        default=1,
        help="Worker threads used inside MNE filtering/resampling for one subject.",
    )
    parser.add_argument(
        "--memmap-dir",
        type=Path,
        default=resolve_repo_path("data/interim/mne_memmap"),
        help="Directory used for disk-backed MNE memmap buffers.",
    )
    parser.add_argument(
        "--preload-mode",
        choices=["memory", "memmap"],
        default="memory",
        help="Raw EDF preload mode. 'memory' is faster; 'memmap' minimizes RAM spikes.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=resolve_repo_path("data/interim/mne_cache"),
        help="MNE cache directory.",
    )
    parser.add_argument(
        "--memmap-min-size",
        type=str,
        default="1M",
        help="MNE memmap threshold passed to mne.set_memmap_min_size (e.g., 1M, 10M).",
    )
    parser.add_argument("--n-jobs", type=int, default=1, help="Parallel workers across subjects.")
    parser.add_argument(
        "--existing-preprocessed-dir",
        type=Path,
        default=None,
        help="Import existing preprocessed epochs (EPCTLxx_stage_epochs.fif) instead of reprocessing raw EDF.",
    )
    parser.add_argument(
        "--existing-mode",
        choices=["symlink", "copy"],
        default="symlink",
        help="Import mode used with --existing-preprocessed-dir.",
    )
    parser.add_argument(
        "--artifact-matrix-dir",
        type=Path,
        default=None,
        help=(
            "Optional directory containing per-subject artifact matrices "
            "(e.g., EPCTLxx_artndxn.mat). When provided, 30 s rows flagged "
            "as bad are excluded before stage epoching."
        ),
    )
    parser.add_argument(
        "--artifact-missing-policy",
        choices=["error", "warn_skip"],
        default="error",
        help=(
            "Behavior when --artifact-matrix-dir is set but a subject artifact "
            "matrix is missing."
        ),
    )
    parser.add_argument(
        "--legacy-compatible",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Match legacy preprocessing output behavior: flat output directory, "
            "W/R stage filenames, and cleaned.csv export."
        ),
    )
    parser.add_argument(
        "--write-cleaned-csv",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write <subject>_cleaned.csv as in the legacy preprocessing script.",
    )
    parser.add_argument(
        "--cleaned-csv-chunk-samples",
        type=int,
        default=20000,
        help="Chunk size for streaming cleaned.csv export without loading all samples at once.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def discover_subjects(raw_root: Path, include_subjects: Iterable[str]) -> List[str]:
    include = {s.strip() for s in include_subjects if s.strip()}
    subjects = sorted(
        p.name for p in raw_root.glob("EPCTL*") if p.is_dir() and re.match(r"EPCTL\d{2}$", p.name)
    )
    if include:
        subjects = [s for s in subjects if s in include]
    return subjects


def load_montage(montage_file: Path) -> mne.channels.DigMontage:
    if not montage_file.exists():
        raise FileNotFoundError(f"Configured montage file not found: {montage_file}")

    try:
        df = pd.read_csv(montage_file, sep="\t", header=None)
    except Exception:
        df = pd.read_csv(montage_file, sep=r"\s+", header=None, engine="python")

    names = df[1].astype(str).str.upper().tolist()
    xyz = df[[2, 3, 4]].to_numpy(dtype=float)
    xyz[:, [0, 1]] = xyz[:, [1, 0]]
    return mne.channels.make_dig_montage(ch_pos=dict(zip(names, xyz)), coord_frame="head")


def time_to_sample(raw: mne.io.BaseRaw, tsec: float) -> int:
    return int(tsec * raw.info["sfreq"])


def normalize_subject_id(text: str) -> str | None:
    m = re.search(r"(EPCTL\d{2})", str(text).upper())
    if m is None:
        return None
    return m.group(1)


def build_artifact_index(artifact_dir: Path) -> dict[str, Path]:
    index: dict[str, Path] = {}
    for path in sorted(artifact_dir.glob("*.mat")):
        if path.name.startswith("._"):
            continue
        subject_id = normalize_subject_id(path.name)
        if subject_id is None:
            continue
        prev = index.get(subject_id)
        if prev is None:
            index[subject_id] = path
            continue
        if "artndxn" in path.name.lower() and "artndxn" not in prev.name.lower():
            index[subject_id] = path
    return index


def load_artifact_matrix(mat_file: Path) -> np.ndarray:
    h5_error: Exception | None = None
    try:
        import h5py  # type: ignore

        with h5py.File(mat_file, "r") as handle:
            if "artndxn" in handle:
                arr = np.asarray(handle["artndxn"])
            else:
                keys = list(handle.keys())
                if not keys:
                    raise ValueError(f"{mat_file} has no datasets")
                arr = np.asarray(handle[keys[0]])
            arr = np.asarray(arr).squeeze()
            if arr.ndim == 1:
                arr = arr[:, np.newaxis]
            if arr.ndim != 2:
                raise ValueError(f"artifact matrix in {mat_file} is not 2D: shape={arr.shape}")
            return arr
    except Exception as exc:
        h5_error = exc

    try:
        import scipy.io as sio  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "Artifact-matrix mode requires h5py or scipy to read .mat files."
        ) from exc

    try:
        mat = sio.loadmat(mat_file)
    except Exception as exc:
        raise RuntimeError(
            f"Could not load artifact matrix {mat_file} "
            f"(h5py error: {h5_error}; scipy error: {exc})"
        ) from exc

    keys = [k for k in mat.keys() if not k.startswith("__")]
    if not keys:
        raise ValueError(f"{mat_file} contains no user variables")

    arr = np.asarray(mat[keys[0]]).squeeze()
    if arr.ndim == 1:
        arr = arr[:, np.newaxis]
    if arr.ndim != 2:
        raise ValueError(f"artifact matrix in {mat_file} is not 2D: shape={arr.shape}")
    return arr


def artifact_bad_row_mask(subject_id: str, annot_rows: int, artifact_matrix: np.ndarray, artifact_file: Path) -> np.ndarray:
    d0, d1 = artifact_matrix.shape
    if d0 == annot_rows and d1 != annot_rows:
        epoch_axis = 0
    elif d1 == annot_rows and d0 != annot_rows:
        epoch_axis = 1
    else:
        # In this dataset, epochs are the larger axis (e.g., ~900x83).
        epoch_axis = 0 if d0 >= d1 else 1

    bad_rows = np.any(artifact_matrix == 0, axis=1 if epoch_axis == 0 else 0)
    bad_rows = np.asarray(bad_rows, dtype=bool)
    if bad_rows.size < annot_rows:
        pad = np.zeros(annot_rows - bad_rows.size, dtype=bool)
        bad_rows = np.concatenate([bad_rows, pad], axis=0)
        print(
            f"[WARN] {subject_id} artifact rows ({bad_rows.size - pad.size}) < "
            f"annotation rows ({annot_rows}) in {artifact_file.name}; padding with clean rows."
        )
    elif bad_rows.size > annot_rows:
        print(
            f"[WARN] {subject_id} artifact rows ({bad_rows.size}) > "
            f"annotation rows ({annot_rows}) in {artifact_file.name}; truncating extras."
        )
        bad_rows = bad_rows[:annot_rows]
    return bad_rows


def get_subject_artifact_mask(subject_id: str, annot_df: pd.DataFrame, args: argparse.Namespace) -> np.ndarray | None:
    if args.artifact_matrix_dir is None:
        return None

    artifact_map = getattr(args, "artifact_map", {})
    artifact_file_str = None
    if isinstance(artifact_map, dict):
        artifact_file_str = artifact_map.get(subject_id)
    artifact_file = Path(artifact_file_str) if artifact_file_str else None

    if artifact_file is None or not artifact_file.exists():
        msg = (
            f"Artifact matrix not found for {subject_id} under {args.artifact_matrix_dir}. "
            "Expected names like EPCTLxx_artndxn.mat (also supports known legacy typos)."
        )
        if args.artifact_missing_policy == "warn_skip":
            print(f"[WARN] {msg} Proceeding without artifact-row rejection.")
            return None
        raise FileNotFoundError(msg)

    artifact = load_artifact_matrix(artifact_file)
    bad_row_mask = artifact_bad_row_mask(
        subject_id=subject_id,
        annot_rows=len(annot_df),
        artifact_matrix=artifact,
        artifact_file=artifact_file,
    )
    print(
        f"[{subject_id}] artifact rejection rows: {int(bad_row_mask.sum())} / {len(bad_row_mask)} "
        f"(file: {artifact_file.name})"
    )
    return bad_row_mask


def annotation_to_events(
    raw: mne.io.BaseRaw,
    annot_df: pd.DataFrame,
    epoch_seconds: float,
    drop_row_mask: np.ndarray | None = None,
) -> Dict[str, np.ndarray]:
    out: Dict[str, list[list[int]]] = {s: [] for s in stage_order()}
    epoch_samples = int(epoch_seconds * raw.info["sfreq"])

    for idx, row in annot_df.iterrows():
        row_idx = int(idx)
        if drop_row_mask is not None and row_idx < drop_row_mask.size and drop_row_mask[row_idx]:
            continue
        raw_label = str(row["label"]).strip()
        try:
            stage = canonical_stage(raw_label)
        except ValueError:
            continue

        onset = time_to_sample(raw, float(row["start_time"]))
        duration_samples = int(float(row["duration"]) * raw.info["sfreq"])
        n_epochs = duration_samples // epoch_samples

        stage_id = stage_order().index(stage) + 1
        for i in range(n_epochs):
            out[stage].append([onset + i * epoch_samples, 0, stage_id])

    return {k: np.asarray(v, dtype=int) for k, v in out.items() if v}


def stage_to_output_label(stage: str, args: argparse.Namespace) -> str:
    if not args.legacy_compatible:
        return stage
    return LEGACY_STAGE_FILENAME.get(stage, stage.replace("_", ""))


def stage_output_file(subject_id: str, stage: str, args: argparse.Namespace) -> Path:
    stage_label = stage_to_output_label(stage, args)
    if args.legacy_compatible:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        return args.output_dir / f"{subject_id}_{stage_label}_epochs.fif"
    out_subject_dir = args.output_dir / subject_id
    out_subject_dir.mkdir(parents=True, exist_ok=True)
    return out_subject_dir / f"{subject_id}_{stage_label}_epochs.fif"


def save_cleaned_csv_streaming(raw: mne.io.BaseRaw, csv_path: Path, chunk_samples: int) -> None:
    if chunk_samples < 1:
        raise ValueError("--cleaned-csv-chunk-samples must be >= 1")

    mode = "w"
    first_chunk = True
    for start in range(0, raw.n_times, chunk_samples):
        stop = min(start + chunk_samples, raw.n_times)
        data, times = raw[:, start:stop]
        # Match legacy format (time as index, channels as columns), but streamed.
        frame = pd.DataFrame(data=data.T, index=times, columns=raw.ch_names)
        frame.to_csv(csv_path, mode=mode, header=first_chunk)
        mode = "a"
        first_chunk = False


def preprocess_subject(subject_id: str, args: argparse.Namespace, montage: mne.channels.DigMontage | None) -> List[SubjectSummary]:
    subj_dir = args.raw_root / subject_id
    edf = subj_dir / f"{subject_id}.edf"
    txt = subj_dir / f"{subject_id}.txt"
    if not edf.exists() or not txt.exists():
        print(f"[WARN] skipping {subject_id}: missing EDF/TXT")
        return []

    memmap_file: Path | None = None
    if args.preload_mode == "memmap":
        args.memmap_dir.mkdir(parents=True, exist_ok=True)
        memmap_file = args.memmap_dir / f"{subject_id}_raw_preload.dat"
        if memmap_file.exists():
            memmap_file.unlink()
        print(f"[{subject_id}] loading EDF -> memmap: {memmap_file}")
        raw = mne.io.read_raw_edf(str(edf), preload=str(memmap_file), verbose="ERROR")
    else:
        print(f"[{subject_id}] loading EDF -> memory")
        raw = mne.io.read_raw_edf(str(edf), preload=True, verbose="ERROR")
    print(f"[{subject_id}] notch_filter {args.notch_freqs} Hz")
    raw.notch_filter(freqs=args.notch_freqs, method="spectrum_fit", n_jobs=args.filter_n_jobs, verbose="ERROR")
    print(f"[{subject_id}] bandpass {args.l_freq}-{args.h_freq} Hz")
    raw.filter(args.l_freq, args.h_freq, fir_design="firwin", n_jobs=args.filter_n_jobs, verbose="ERROR")
    print(f"[{subject_id}] resample -> {args.sfreq} Hz")
    raw.resample(args.sfreq, npad="auto", n_jobs=args.filter_n_jobs, verbose="ERROR")

    raw.rename_channels(
        {
            ch: ch.replace("-Ref", "").upper()
            for ch in raw.ch_names
            if "-Ref" in ch
        }
    )
    if montage is not None:
        raw.set_montage(montage, on_missing="ignore")

    bad_locs = [
        ch["ch_name"]
        for ch in raw.info["chs"]
        if not np.isfinite(ch["loc"][:3]).all()
    ]
    drop = [ch for ch in raw.ch_names if ch in EXCLUDE_CHANNELS or ch in bad_locs]
    if drop:
        raw.drop_channels(drop)

    raw.set_eeg_reference("average", projection=True)
    raw.apply_proj()

    if args.write_cleaned_csv:
        if args.legacy_compatible:
            cleaned_csv = args.output_dir / f"{subject_id}_cleaned.csv"
        else:
            out_subject_dir = args.output_dir / subject_id
            out_subject_dir.mkdir(parents=True, exist_ok=True)
            cleaned_csv = out_subject_dir / f"{subject_id}_cleaned.csv"
        if cleaned_csv.exists() and not args.overwrite:
            print(f"[skip] {cleaned_csv}")
        else:
            print(f"[{subject_id}] writing cleaned csv: {cleaned_csv}")
            save_cleaned_csv_streaming(raw, cleaned_csv, args.cleaned_csv_chunk_samples)
            print(f"[ok] cleaned csv saved: {cleaned_csv}")

    print(f"[{subject_id}] reading stage annotations")
    annot = pd.read_csv(txt, sep="\t", header=None, names=["label", "start_time", "duration"])
    artifact_mask = get_subject_artifact_mask(subject_id, annot, args)
    events_by_stage = annotation_to_events(raw, annot, args.epoch_seconds, drop_row_mask=artifact_mask)
    print(f"[{subject_id}] stages with events: {sorted(events_by_stage.keys())}")

    summaries: list[SubjectSummary] = []
    for stage, events in events_by_stage.items():
        out_file = stage_output_file(subject_id, stage, args)
        if out_file.exists() and not args.overwrite:
            print(f"[skip] {out_file}")
            continue

        event_id = {stage: int(stage_order().index(stage) + 1)}
        epochs = mne.Epochs(
            raw,
            events,
            event_id=event_id,
            tmin=0.0,
            tmax=args.epoch_seconds,
            baseline=None,
            preload=False,
            verbose="ERROR",
        )
        epochs.save(str(out_file), overwrite=True, verbose="ERROR")

        summaries.append(
            SubjectSummary(
                subject_id=subject_id,
                stage=stage_to_output_label(stage, args),
                n_epochs=len(epochs),
                out_file=str(out_file),
            )
        )
        print(f"[ok] {subject_id} {stage}: {len(epochs)} epochs")

    raw.close()
    try:
        if memmap_file is not None and memmap_file.exists():
            memmap_file.unlink()
    except Exception:
        pass

    return summaries


def preprocess_subject_worker(subject_id: str, worker_args: dict[str, object]) -> list[dict[str, object]]:
    args = argparse.Namespace(**worker_args)
    montage = load_montage(args.montage_file)
    rows = preprocess_subject(subject_id, args, montage)
    return [r.__dict__ for r in rows]


def parse_existing_epoch(path: Path) -> tuple[str, str] | None:
    m = re.match(r"^(EPCTL\d{2})_(W|W_|N1|N2|N3|R|R_)_epochs\.fif$", path.name)
    if not m:
        return None
    subject_id = m.group(1)
    stage = canonical_stage(m.group(2))
    return subject_id, stage


def import_existing_preprocessed(args: argparse.Namespace) -> list[SubjectSummary]:
    if args.existing_preprocessed_dir is None or not args.existing_preprocessed_dir.exists():
        raise FileNotFoundError(f"existing preprocessed dir not found: {args.existing_preprocessed_dir}")

    include = {s.strip() for s in args.subject if str(s).strip()}
    epoch_files = sorted(args.existing_preprocessed_dir.glob("EPCTL*_epochs.fif"))
    if not epoch_files:
        raise SystemExit(f"No epochs files found in {args.existing_preprocessed_dir}")

    out_rows: list[SubjectSummary] = []
    for src_file in epoch_files:
        parsed = parse_existing_epoch(src_file)
        if parsed is None:
            continue
        subject_id, stage = parsed
        if include and subject_id not in include:
            continue

        dst_file = stage_output_file(subject_id, stage, args)
        if dst_file.exists():
            if not args.overwrite:
                print(f"[skip] {dst_file}")
                continue
            if dst_file.is_symlink() or dst_file.is_file():
                dst_file.unlink()

        if args.existing_mode == "copy":
            shutil.copy2(src_file, dst_file)
        else:
            dst_file.symlink_to(src_file.resolve())

        n_epochs = -1
        try:
            epochs = mne.read_epochs(str(dst_file), preload=False, verbose="ERROR")
            n_epochs = len(epochs)
        except Exception:
            pass

        out_rows.append(
            SubjectSummary(
                subject_id=subject_id,
                stage=stage_to_output_label(stage, args),
                n_epochs=n_epochs,
                out_file=str(dst_file),
            )
        )
        print(f"[ok] imported {subject_id} {stage} -> {dst_file}")

    if not out_rows:
        raise SystemExit("No matching epochs were imported.")
    return out_rows


def main() -> int:
    args = parse_args()

    if MNE_IMPORT_ERROR is not None:
        raise SystemExit(
            "MNE is required for preprocessing. Install dependencies first "
            f"(import error: {MNE_IMPORT_ERROR})."
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.memmap_dir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    if args.n_jobs < 1:
        raise SystemExit("--n-jobs must be >= 1")
    if args.filter_n_jobs < 1:
        raise SystemExit("--filter-n-jobs must be >= 1")
    if args.sfreq < 256.0:
        raise SystemExit("--sfreq must be >= 256 Hz (no lower than legacy standard).")
    if args.artifact_matrix_dir is not None and args.existing_preprocessed_dir is not None:
        raise SystemExit("--artifact-matrix-dir cannot be combined with --existing-preprocessed-dir")
    if args.artifact_matrix_dir is not None:
        if not args.artifact_matrix_dir.exists():
            raise SystemExit(f"--artifact-matrix-dir not found: {args.artifact_matrix_dir}")
        artifact_index = build_artifact_index(args.artifact_matrix_dir)
        if not artifact_index:
            raise SystemExit(f"No artifact .mat files found under {args.artifact_matrix_dir}")
        args.artifact_map = {subject_id: str(path) for subject_id, path in artifact_index.items()}
    else:
        args.artifact_map = {}

    mne.set_cache_dir(str(args.cache_dir))
    mne.set_memmap_min_size(args.memmap_min_size)

    all_rows: list[SubjectSummary] = []
    if args.existing_preprocessed_dir is not None:
        print(f"[SomnoDynamics] importing existing preprocessed epochs from: {args.existing_preprocessed_dir}")
        all_rows.extend(import_existing_preprocessed(args))
    else:
        subjects = discover_subjects(args.raw_root, args.subject)
        if not subjects:
            raise SystemExit(f"No subjects found in {args.raw_root}")

        print(f"[SomnoDynamics] using montage file: {args.montage_file}")
        if args.n_jobs == 1:
            montage = load_montage(args.montage_file)
            for subject in subjects:
                all_rows.extend(preprocess_subject(subject, args, montage))
        else:
            print(f"[SomnoDynamics] preprocessing in parallel mode across subjects (n_jobs={args.n_jobs})")
            worker_args = {
                "raw_root": args.raw_root,
                "output_dir": args.output_dir,
                "montage_file": args.montage_file,
                "subject": [],
                "epoch_seconds": args.epoch_seconds,
                "sfreq": args.sfreq,
                "l_freq": args.l_freq,
                "h_freq": args.h_freq,
                "notch_freqs": args.notch_freqs,
                "filter_n_jobs": args.filter_n_jobs,
                "memmap_dir": args.memmap_dir,
                "cache_dir": args.cache_dir,
                "memmap_min_size": args.memmap_min_size,
                "preload_mode": args.preload_mode,
                "n_jobs": 1,
                "existing_preprocessed_dir": None,
                "existing_mode": args.existing_mode,
                "artifact_matrix_dir": args.artifact_matrix_dir,
                "artifact_missing_policy": args.artifact_missing_policy,
                "artifact_map": args.artifact_map,
                "legacy_compatible": args.legacy_compatible,
                "write_cleaned_csv": args.write_cleaned_csv,
                "cleaned_csv_chunk_samples": args.cleaned_csv_chunk_samples,
                "overwrite": args.overwrite,
            }
            failures: list[str] = []
            with ProcessPoolExecutor(max_workers=args.n_jobs) as pool:
                futures = {pool.submit(preprocess_subject_worker, subject, worker_args): subject for subject in subjects}
                for fut in as_completed(futures):
                    subject = futures[fut]
                    try:
                        rows = fut.result()
                    except Exception as exc:
                        failures.append(subject)
                        print(f"[ERROR] subject {subject} failed: {exc}")
                        continue
                    all_rows.extend(SubjectSummary(**row) for row in rows)
            if failures:
                failed = ", ".join(sorted(failures))
                raise SystemExit(f"Preprocessing failed for subjects: {failed}")

    summary_csv = args.output_dir / "preprocessing_summary.csv"
    pd.DataFrame([r.__dict__ for r in all_rows]).to_csv(summary_csv, index=False)
    print(f"[SomnoDynamics] preprocessing summary: {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
