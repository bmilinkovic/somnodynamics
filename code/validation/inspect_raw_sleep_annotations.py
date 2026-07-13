#!/usr/bin/env python3
"""Inspect raw sleep annotations and expected 2-second window counts.

This validation script reads raw EDF/TXT files from the configured dataset root
and writes stage-duration and EDF-annotation summaries for audit purposes.
"""

from __future__ import annotations

import argparse
import os
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import canonical_stage, load_config, resolve_repo_path, stage_order  # noqa: E402


os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
try:
    import mne
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"MNE is required for this validation script (error: {exc})")


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Inspect raw sleep annotations and stage durations.")
    parser.add_argument("--raw-root", type=Path, default=Path(cfg["raw_data_root"]))
    parser.add_argument("--subject", action="append", default=[])
    parser.add_argument(
        "--stage-summary-csv",
        type=Path,
        default=resolve_repo_path("results/stats/validation/raw_stage_duration_summary.csv"),
    )
    parser.add_argument(
        "--edf-annotation-csv",
        type=Path,
        default=resolve_repo_path("results/stats/validation/raw_edf_annotation_summary.csv"),
    )
    return parser.parse_args()


def discover_subjects(raw_root: Path, include_subjects: Iterable[str]) -> list[str]:
    include = {s.strip() for s in include_subjects if s.strip()}
    subjects = sorted(
        p.name for p in raw_root.glob("EPCTL*") if p.is_dir() and re.match(r"EPCTL\d{2}$", p.name)
    )
    if include:
        subjects = [s for s in subjects if s in include]
    return subjects


def summarize_subject(raw_root: Path, subject_id: str) -> tuple[list[dict], list[dict]]:
    subj_dir = raw_root / subject_id
    txt = subj_dir / f"{subject_id}.txt"
    edf = subj_dir / f"{subject_id}.edf"
    if not txt.exists() or not edf.exists():
        return [], []

    annot = pd.read_csv(txt, sep="\t", header=None, names=["label", "start_time", "duration"])
    rows: list[dict] = []
    ignored_labels = Counter()
    for _, row in annot.iterrows():
        raw_stage = str(row["label"]).strip()
        try:
            stage = canonical_stage(raw_stage)
        except ValueError:
            ignored_labels[raw_stage] += 1
            continue
        rows.append(
            {
                "subject_id": subject_id,
                "stage_internal": stage,
                "duration_s": float(row["duration"]),
            }
        )

    if rows:
        stage_df = pd.DataFrame(rows)
        grouped = (
            stage_df.groupby(["subject_id", "stage_internal"], as_index=False)
            .agg(row_count=("duration_s", "size"), total_duration_s=("duration_s", "sum"))
        )
        grouped["expected_nonoverlap_2s_windows"] = (grouped["total_duration_s"] // 2.0).astype(int)
    else:
        grouped = pd.DataFrame(columns=["subject_id", "stage_internal", "row_count", "total_duration_s", "expected_nonoverlap_2s_windows"])

    existing_stages = set(grouped["stage_internal"].tolist()) if not grouped.empty else set()
    missing_rows = []
    for stage in stage_order():
        if stage in existing_stages:
            continue
        missing_rows.append(
            {
                "subject_id": subject_id,
                "stage_internal": stage,
                "row_count": 0,
                "total_duration_s": 0.0,
                "expected_nonoverlap_2s_windows": 0,
            }
        )
    if missing_rows:
        grouped = pd.concat([grouped, pd.DataFrame(missing_rows)], ignore_index=True)

    grouped = grouped.sort_values(["subject_id", "stage_internal"]).reset_index(drop=True)
    grouped["ignored_nonstage_rows"] = int(sum(ignored_labels.values()))
    grouped["ignored_nonstage_labels"] = ",".join(sorted(ignored_labels.keys()))

    raw = mne.io.read_raw_edf(str(edf), preload=False, verbose="ERROR")
    edf_rows = [
        {
            "subject_id": subject_id,
            "sfreq_hz": float(raw.info["sfreq"]),
            "n_channels": int(len(raw.ch_names)),
            "edf_annotation_count": int(len(raw.annotations)),
            "annotation_description": str(desc),
            "annotation_description_count": int(count),
        }
        for desc, count in sorted(Counter(raw.annotations.description).items())
    ]
    return grouped.to_dict(orient="records"), edf_rows


def main() -> int:
    args = parse_args()
    if not args.raw_root.exists():
        raise SystemExit(f"Raw root not found: {args.raw_root}")

    subjects = discover_subjects(args.raw_root, args.subject)
    if not subjects:
        raise SystemExit(f"No subjects found under {args.raw_root}")

    stage_rows: list[dict] = []
    edf_rows: list[dict] = []
    for subject_id in subjects:
        s_rows, e_rows = summarize_subject(args.raw_root, subject_id)
        stage_rows.extend(s_rows)
        edf_rows.extend(e_rows)

    stage_df = pd.DataFrame(stage_rows).sort_values(["subject_id", "stage_internal"]).reset_index(drop=True)
    edf_df = pd.DataFrame(edf_rows).sort_values(["subject_id", "annotation_description"]).reset_index(drop=True)

    args.stage_summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.edf_annotation_csv.parent.mkdir(parents=True, exist_ok=True)
    stage_df.to_csv(args.stage_summary_csv, index=False)
    edf_df.to_csv(args.edf_annotation_csv, index=False)

    print(f"[ok] stage summary -> {args.stage_summary_csv}")
    print(f"[ok] EDF annotation summary -> {args.edf_annotation_csv}")
    print(f"[ok] subjects processed: {len(subjects)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
