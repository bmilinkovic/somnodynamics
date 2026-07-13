#!/usr/bin/env python3
"""
Yeo-7 source reconstruction for the healthy sleep project.

This script mirrors the existing LCMV pipeline and reuses the 7-region Yeo
logic used in dmt_project:
  - fsaverage forward model
  - LCMV beamformer (NAI, max-power)
  - deterministic 7-region ordering
  - per-subject/per-stage ROI time series export
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

try:
    from scipy.io import savemat
    SCIPY_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover - environment dependent
    savemat = None
    SCIPY_IMPORT_ERROR = exc

try:
    import pandas as pd
    PANDAS_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover - environment dependent
    pd = None
    PANDAS_IMPORT_ERROR = exc

try:
    import mne
    from mne.beamformer import apply_lcmv_epochs, make_lcmv
    MNE_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover - environment dependent
    mne = None
    apply_lcmv_epochs = None
    make_lcmv = None
    MNE_IMPORT_ERROR = exc


LOGGER = logging.getLogger("source_reconstruction_yeo7")

SLEEP_STAGE_MAPPING: Dict[str, int] = {"W": 1, "N1": 2, "N2": 3, "N3": 4, "R": 5}

EXCLUDE_CHANNELS: List[str] = [
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
]

# Deterministic label order, matching dmt_project 7-region mapping.
REGION_DEFINITIONS: Sequence[Tuple[str, Sequence[str]]] = (
    ("Visual", ("7Networks_1-lh", "7Networks_1-rh")),
    ("SomMot", ("7Networks_2-lh", "7Networks_2-rh")),
    ("DorsAttn", ("7Networks_3-lh", "7Networks_3-rh")),
    ("SalVentAttn", ("7Networks_4-lh", "7Networks_4-rh")),
    ("Limbic", ("7Networks_5-lh", "7Networks_5-rh")),
    ("FrontoParietal", ("7Networks_6-lh", "7Networks_6-rh")),
    ("Default", ("7Networks_7-lh", "7Networks_7-rh")),
)

DEFAULT_APHNY_DATA_DIR = Path(
    os.environ.get(
        "APHNY_DATA_DIR",
        "/Users/borjan/Library/CloudStorage/OneDrive-ICM/APHNY data",
    )
)
DEFAULT_RAW_DIR = Path(os.environ.get("APHNY_RAW_DIR", str(DEFAULT_APHNY_DATA_DIR / "raw")))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Yeo-7 source reconstruction")
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help="Directory with subject folders containing EDF/TXT files.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_APHNY_DATA_DIR,
        help="Data root containing montage file.",
    )
    parser.add_argument(
        "--montage-file",
        type=Path,
        default=None,
        help="Optional explicit montage file path (.pos). Overrides --data-dir lookup.",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=Path(
            "/Users/borjan/code/healthy-sleep-project/code/results/source_recon/yeo7"
        ),
        help="Output root for Yeo-7 source time series.",
    )
    parser.add_argument(
        "--subjects-dir",
        type=Path,
        default=Path("/Users/borjan/mne_data/MNE-fsaverage-data"),
        help="Directory containing fsaverage files.",
    )
    parser.add_argument("--sfreq", type=float, default=256.0)
    parser.add_argument("--fmin", type=float, default=1.0)
    parser.add_argument("--fmax", type=float, default=45.0)
    parser.add_argument("--notch", nargs="+", type=float, default=[50.0, 100.0])
    parser.add_argument("--epoch-seconds", type=float, default=2.0)
    parser.add_argument("--lcmv-reg", type=float, default=0.05)
    parser.add_argument("--mindist-mm", type=float, default=5.0)
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def build_yeo7_labels(subjects_dir: Path) -> List[mne.Label]:
    labels_all = mne.read_labels_from_annot(
        "fsaverage",
        parc="Yeo2011_7Networks_N1000",
        subjects_dir=str(subjects_dir),
    )

    labels_7: List[mne.Label] = []
    for region_name, network_labels in REGION_DEFINITIONS:
        region_labels = [
            label
            for label in labels_all
            if any(label.name == net_label for net_label in network_labels)
        ]
        if len(region_labels) != len(network_labels):
            raise RuntimeError(
                f"Could not resolve all Yeo labels for {region_name}. "
                f"Expected {list(network_labels)}, got {[l.name for l in region_labels]}"
            )
        combined_label = region_labels[0]
        for label in region_labels[1:]:
            combined_label += label
        combined_label.name = region_name
        labels_7.append(combined_label)

    return labels_7


def load_sleep_stages(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None)
    df.columns = ["label", "start_time", "duration"]
    return df


def load_montage(montage_file: Path | None) -> mne.channels.DigMontage:
    candidates: List[Path] = []
    if montage_file is not None:
        candidates.append(montage_file)
    env_montage = os.environ.get("APHNY_MONTAGE_FILE", "").strip()
    if env_montage:
        candidates.insert(0, Path(env_montage))
    project_fallback = Path(__file__).resolve().parents[1] / "Co-registered_average_positions.pos"
    candidates.append(project_fallback)

    last_error: Exception | None = None
    selected_local_copy: Path | None = None
    for candidate in candidates:
        if not candidate.exists():
            continue
        try:
            selected_local_copy = copy_with_retries(candidate)
            break
        except Exception as exc:  # pragma: no cover - filesystem dependent
            last_error = exc
            continue

    if selected_local_copy is None:
        raise FileNotFoundError(
            "Unable to access montage file. Tried: "
            f"{[str(path) for path in candidates]}. "
            "If using OneDrive, set the file to 'Always keep on this device' "
            "or export APHNY_MONTAGE_FILE to a local copy."
        ) from last_error

    montage_data = read_montage_table(selected_local_copy)
    electrode_names = montage_data[1].astype(str).str.upper().tolist()
    coordinates = montage_data[[2, 3, 4]].to_numpy(dtype=float)
    coordinates[:, [0, 1]] = coordinates[:, [1, 0]]
    return mne.channels.make_dig_montage(
        ch_pos=dict(zip(electrode_names, coordinates)),
        coord_frame="head",
    )


def copy_with_retries(source_path: Path, retries: int = 5, delay_seconds: float = 2.0) -> Path:
    temp_path = Path(tempfile.gettempdir()) / (
        f"aphny_montage_cache_{int(time.time() * 1_000_000)}_{source_path.name}"
    )
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            shutil.copyfile(source_path, temp_path)
            if file_is_all_zero_bytes(temp_path):
                raise RuntimeError(
                    f"{source_path} appears to be a cloud placeholder (all-zero bytes)."
                )
            return temp_path
        except Exception as exc:  # pragma: no cover - filesystem dependent
            last_error = exc
            LOGGER.warning(
                "Failed to copy montage file (attempt %d/%d): %s",
                attempt,
                retries,
                exc,
            )
            if attempt < retries:
                time.sleep(delay_seconds)
    raise RuntimeError(f"Could not copy montage file from {source_path}") from last_error


def file_is_all_zero_bytes(path: Path, sample_size: int = 512) -> bool:
    with path.open("rb") as handle:
        sample = handle.read(sample_size)
    return bool(sample) and all(byte == 0 for byte in sample)


def read_montage_table(path: Path):
    try:
        return pd.read_csv(path, sep="\t", header=None)
    except Exception:
        # Some files are whitespace-delimited rather than tab-delimited.
        return pd.read_csv(path, sep=r"\s+", header=None, engine="python")


def time_to_sample(raw: mne.io.BaseRaw, time_in_seconds: float) -> int:
    return int(time_in_seconds * raw.info["sfreq"])


def create_epochs_by_stage(
    raw_eeg: mne.io.BaseRaw,
    sleep_stages_df: pd.DataFrame,
    epoch_seconds: float,
) -> Dict[str, mne.Epochs]:
    epochs_by_stage: Dict[str, mne.Epochs] = {}
    for stage_label, events_df in sleep_stages_df.groupby("label"):
        if stage_label not in SLEEP_STAGE_MAPPING:
            continue

        events: List[List[int]] = []
        for _, row in events_df.iterrows():
            onset_sample = time_to_sample(raw_eeg, float(row["start_time"]))
            event_id = SLEEP_STAGE_MAPPING[stage_label]
            events.append([onset_sample, 0, event_id])

        if not events:
            continue

        events_np = np.asarray(events, dtype=int)
        epochs = mne.Epochs(
            raw_eeg,
            events_np,
            event_id=SLEEP_STAGE_MAPPING[stage_label],
            tmin=0,
            tmax=epoch_seconds,
            baseline=None,
            preload=True,
        )
        epochs_by_stage[stage_label] = epochs
        LOGGER.info("Created %d epochs for stage %s", len(epochs), stage_label)
    return epochs_by_stage


def preprocess_subject(
    edf_file: Path,
    sleep_stage_file: Path,
    montage: mne.channels.DigMontage,
    sfreq: float,
    fmin: float,
    fmax: float,
    notch: Sequence[float],
    epoch_seconds: float,
) -> Tuple[mne.io.BaseRaw, Dict[str, mne.Epochs]]:
    validate_edf_access(edf_file)
    raw = mne.io.read_raw_edf(str(edf_file), preload=True)
    raw.notch_filter(freqs=list(notch), method="spectrum_fit")
    raw.filter(fmin, fmax, fir_design="firwin")
    raw.resample(sfreq, npad="auto")

    raw.rename_channels(
        {
            ch_name: ch_name.replace("-Ref", "").upper()
            for ch_name in raw.ch_names
            if "-Ref" in ch_name
        }
    )
    raw.set_montage(montage, on_missing="ignore")

    bad_channels = [
        ch["ch_name"] for ch in raw.info["chs"] if not np.isfinite(ch["loc"][:3]).all()
    ]
    raw_eeg_cleaned = raw.copy().drop_channels(bad_channels + EXCLUDE_CHANNELS)
    raw_eeg_cleaned.set_eeg_reference("average", projection=True)
    raw_eeg_cleaned.apply_proj()

    sleep_stages_df = load_sleep_stages(sleep_stage_file)
    epochs_by_stage = create_epochs_by_stage(
        raw_eeg_cleaned, sleep_stages_df, epoch_seconds
    )
    return raw_eeg_cleaned, epochs_by_stage


def validate_edf_access(edf_file: Path) -> None:
    try:
        with edf_file.open("rb") as handle:
            sample = handle.read(256)
    except Exception as exc:
        raise RuntimeError(
            f"Cannot access EDF file: {edf_file}. "
            "If this is in OneDrive CloudStorage, move data to a local folder "
            "or mark files as 'Always keep on this device'."
        ) from exc
    if not sample or all(byte == 0 for byte in sample):
        raise RuntimeError(
            f"EDF file appears unavailable/offline: {edf_file}. "
            "The file header is empty/all-zero, which is typical of a cloud placeholder."
        )


def save_stage_outputs(
    out_mat_file: Path,
    out_csv_file: Path,
    source_ts: np.ndarray,
    labels: Sequence[str],
    sfreq: float,
    epoch_seconds: float,
    inverse_params: Dict[str, object],
) -> None:
    out_mat_file.parent.mkdir(parents=True, exist_ok=True)
    out_csv_file.parent.mkdir(parents=True, exist_ok=True)

    reshaped = source_ts.transpose(2, 1, 0).reshape(-1, source_ts.shape[0])
    pd.DataFrame(reshaped, columns=list(labels)).to_csv(out_csv_file, index=False)

    metadata = {
        "atlas": "yeo7",
        "label_order": np.asarray(labels, dtype=object),
        "sampling_rate_hz": float(sfreq),
        "epoch_length_seconds": float(epoch_seconds),
        "inverse_params": inverse_params,
    }
    savemat(
        out_mat_file,
        {
            "source_ts": source_ts,
            "labels": np.asarray(labels, dtype=object),
            "atlas": "yeo7",
            "metadata": metadata,
        },
    )


def process_subject(
    subject_dir: Path,
    args: argparse.Namespace,
    montage: mne.channels.DigMontage,
    src: mne.SourceSpaces,
    bem_solution: mne.bem.ConductorModel,
    labels_7: Sequence[mne.Label],
) -> None:
    subject = subject_dir.name
    edf_file = subject_dir / f"{subject}.edf"
    sleep_stage_file = subject_dir / f"{subject}.txt"
    if not edf_file.exists() or not sleep_stage_file.exists():
        LOGGER.warning("Missing EDF or TXT for %s, skipping", subject)
        return

    LOGGER.info("Processing subject %s", subject)
    raw_eeg, epochs_by_stage = preprocess_subject(
        edf_file=edf_file,
        sleep_stage_file=sleep_stage_file,
        montage=montage,
        sfreq=args.sfreq,
        fmin=args.fmin,
        fmax=args.fmax,
        notch=args.notch,
        epoch_seconds=args.epoch_seconds,
    )
    if not epochs_by_stage:
        LOGGER.warning("No valid stage epochs found for %s", subject)
        return

    noise_cov = mne.compute_raw_covariance(raw_eeg, method="empirical")
    fwd = mne.make_forward_solution(
        raw_eeg.info,
        trans="fsaverage",
        src=src,
        bem=bem_solution,
        eeg=True,
        mindist=args.mindist_mm,
    )
    inverse_operator = make_lcmv(
        raw_eeg.info,
        fwd,
        noise_cov,
        reg=args.lcmv_reg,
        pick_ori="max-power",
        weight_norm="nai",
    )

    ordered_labels = [label.name for label in labels_7]
    inverse_params = {
        "method": "LCMV",
        "weight_norm": "nai",
        "pick_ori": "max-power",
        "reg": float(args.lcmv_reg),
        "mindist_mm": float(args.mindist_mm),
        "src_spacing": "ico4",
    }

    for stage_label, epochs in epochs_by_stage.items():
        stcs = apply_lcmv_epochs(epochs, inverse_operator)
        label_ts = mne.extract_label_time_course(
            stcs,
            labels_7,
            fwd["src"],
            mode="mean_flip",
            return_generator=False,
            allow_empty=True,
        )
        source_ts = np.transpose(np.asarray(label_ts), (1, 2, 0))
        if source_ts.shape[0] != len(ordered_labels):
            raise RuntimeError(
                f"Unexpected region count for {subject} {stage_label}: "
                f"{source_ts.shape[0]} (expected {len(ordered_labels)})"
            )

        subject_out = args.out_root / subject
        mat_file = (
            subject_out
            / "mat"
            / f"{subject}_{stage_label}_source_time_series.mat"
        )
        csv_file = (
            subject_out
            / "csv"
            / f"{subject}_{stage_label}_source_time_series.csv"
        )
        save_stage_outputs(
            out_mat_file=mat_file,
            out_csv_file=csv_file,
            source_ts=source_ts,
            labels=ordered_labels,
            sfreq=args.sfreq,
            epoch_seconds=args.epoch_seconds,
            inverse_params=inverse_params,
        )
        LOGGER.info(
            "Saved %s %s source time series: %s",
            subject,
            stage_label,
            mat_file,
        )


def main() -> None:
    args = parse_args()
    if MNE_IMPORT_ERROR is not None:
        raise ModuleNotFoundError(
            "MNE-Python is required for source reconstruction. "
            "Install dependencies first."
        ) from MNE_IMPORT_ERROR
    if PANDAS_IMPORT_ERROR is not None:
        raise ModuleNotFoundError(
            "pandas is required for source reconstruction. "
            "Install dependencies first."
        ) from PANDAS_IMPORT_ERROR
    if SCIPY_IMPORT_ERROR is not None:
        raise ModuleNotFoundError(
            "scipy is required for source reconstruction. "
            "Install dependencies first."
        ) from SCIPY_IMPORT_ERROR

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    args.out_root.mkdir(parents=True, exist_ok=True)
    mne.datasets.fetch_fsaverage(subjects_dir=str(args.subjects_dir))
    src = mne.setup_source_space(
        "fsaverage",
        spacing="ico4",
        subjects_dir=str(args.subjects_dir),
        add_dist=False,
    )
    bem_model = mne.make_bem_model(
        subject="fsaverage",
        ico=4,
        subjects_dir=str(args.subjects_dir),
    )
    bem_solution = mne.make_bem_solution(bem_model)

    labels_7 = build_yeo7_labels(args.subjects_dir)
    label_names = [label.name for label in labels_7]
    LOGGER.info("Using Yeo-7 label order: %s", label_names)
    savemat(
        args.out_root / "yeo7_label_order.mat",
        {"labels": np.asarray(label_names, dtype=object), "atlas": "yeo7"},
    )

    montage_file = args.montage_file or (args.data_dir / "Co-registered_average_positions.pos")
    montage = load_montage(montage_file)

    if not args.raw_dir.exists():
        raise FileNotFoundError(f"Raw directory does not exist: {args.raw_dir}")
    try:
        subject_dirs = sorted(
            [path for path in args.raw_dir.iterdir() if path.is_dir()],
            key=lambda path: path.name,
        )
    except Exception as exc:
        raise RuntimeError(
            f"Cannot read raw directory: {args.raw_dir}. "
            "Use --raw-dir with a local, accessible folder."
        ) from exc
    for subject_dir in subject_dirs:
        process_subject(
            subject_dir=subject_dir,
            args=args,
            montage=montage,
            src=src,
            bem_solution=bem_solution,
            labels_7=labels_7,
        )

    LOGGER.info("Finished Yeo-7 source reconstruction. Output root: %s", args.out_root)


if __name__ == "__main__":
    main()
