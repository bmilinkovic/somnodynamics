#!/usr/bin/env python3
"""Compare legacy vs current source epoch counts against raw-stage durations.

Outputs a machine-readable CSV with:
- raw stage durations from TXT annotations
- expected non-overlapping 2-second window counts
- legacy source epoch counts (old trusted outputs)
- current source epoch counts (SomnoDynamics outputs)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

try:
    from scipy.io import loadmat
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"scipy is required for this validation script (error: {exc})")

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import canonical_stage, load_config, resolve_repo_path, stage_order  # noqa: E402


OLD_STAGE_ALIAS = {"W_": "W", "N1": "N1", "N2": "N2", "N3": "N3", "R_": "R"}


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Compare old vs new source epoch counts.")
    parser.add_argument("--raw-root", type=Path, default=Path(cfg["raw_data_root"]))
    parser.add_argument(
        "--old-source-root",
        type=Path,
        default=Path("/Volumes/dataSets/healthy-sleep-project-datasets/data/preprocessed_source_reconstructed"),
    )
    parser.add_argument(
        "--new-source-root",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["source_recon_anatomical7"]),
    )
    parser.add_argument("--subject", action="append", default=[])
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=resolve_repo_path("results/stats/validation/source_epoch_count_comparison.csv"),
    )
    return parser.parse_args()


def discover_subjects(raw_root: Path, include_subjects: list[str]) -> list[str]:
    include = {s.strip() for s in include_subjects if s.strip()}
    subjects = sorted(
        p.name for p in raw_root.glob("EPCTL*") if p.is_dir() and re.match(r"EPCTL\d{2}$", p.name)
    )
    if include:
        subjects = [s for s in subjects if s in include]
    return subjects


def read_source_shape(path: Path) -> tuple[int, int, int] | None:
    if not path.exists():
        return None
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    key = "source_ts" if "source_ts" in mat else ("data" if "data" in mat else None)
    if key is None:
        return None
    shape = tuple(int(x) for x in mat[key].shape)
    if len(shape) != 3:
        return None
    return shape


def load_stage_totals(txt_file: Path) -> tuple[dict[str, float], dict[str, int], int]:
    df = pd.read_csv(txt_file, sep="\t", header=None, names=["label", "start_time", "duration"])
    stage_rows: dict[str, int] = {s: 0 for s in stage_order()}
    stage_duration: dict[str, float] = {s: 0.0 for s in stage_order()}
    ignored = 0
    for _, row in df.iterrows():
        raw_label = str(row["label"]).strip()
        try:
            stage = canonical_stage(raw_label)
        except ValueError:
            ignored += 1
            continue
        stage_rows[stage] += 1
        stage_duration[stage] += float(row["duration"])
    return stage_duration, stage_rows, ignored


def compare_subject(subject_id: str, args: argparse.Namespace) -> list[dict]:
    txt = args.raw_root / subject_id / f"{subject_id}.txt"
    if not txt.exists():
        return []
    duration_s, row_count, ignored_rows = load_stage_totals(txt)

    rows: list[dict] = []
    for stage in stage_order():
        expected_2s = int(duration_s[stage] // 2.0)

        old_stage = OLD_STAGE_ALIAS[stage]
        old_mat = args.old_source_root / subject_id / "mat" / f"{subject_id}_{old_stage}_source_time_series.mat"
        new_mat = args.new_source_root / subject_id / "mat" / f"{subject_id}_{stage}_source_time_series.mat"

        old_shape = read_source_shape(old_mat)
        new_shape = read_source_shape(new_mat)
        old_epochs = old_shape[2] if old_shape is not None else None
        new_epochs = new_shape[2] if new_shape is not None else None

        rows.append(
            {
                "subject_id": subject_id,
                "stage_internal": stage,
                "raw_stage_duration_s": float(duration_s[stage]),
                "raw_stage_row_count": int(row_count[stage]),
                "expected_nonoverlap_2s_windows": expected_2s,
                "legacy_old_epochs": old_epochs,
                "legacy_old_n_regions": old_shape[0] if old_shape is not None else None,
                "legacy_old_n_timepoints": old_shape[1] if old_shape is not None else None,
                "current_new_epochs": new_epochs,
                "current_new_n_regions": new_shape[0] if new_shape is not None else None,
                "current_new_n_timepoints": new_shape[1] if new_shape is not None else None,
                "legacy_matches_row_count": bool(old_epochs == row_count[stage]) if old_epochs is not None else None,
                "current_matches_expected_2s_windows": bool(new_epochs == expected_2s) if new_epochs is not None else None,
                "ratio_new_over_old_epochs": (float(new_epochs) / float(old_epochs)) if (new_epochs and old_epochs) else None,
                "legacy_mat_file": str(old_mat),
                "current_mat_file": str(new_mat),
                "ignored_nonstage_rows": int(ignored_rows),
            }
        )
    return rows


def main() -> int:
    args = parse_args()
    if not args.raw_root.exists():
        raise SystemExit(f"Raw root not found: {args.raw_root}")
    if not args.old_source_root.exists():
        raise SystemExit(f"Old source root not found: {args.old_source_root}")
    if not args.new_source_root.exists():
        raise SystemExit(f"New source root not found: {args.new_source_root}")

    subjects = discover_subjects(args.raw_root, args.subject)
    if not subjects:
        raise SystemExit(f"No subjects found under {args.raw_root}")

    rows: list[dict] = []
    for subject_id in subjects:
        rows.extend(compare_subject(subject_id, args))

    df = pd.DataFrame(rows).sort_values(["subject_id", "stage_internal"]).reset_index(drop=True)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    print(f"[ok] wrote {len(df)} rows -> {args.out_csv}")
    if "ratio_new_over_old_epochs" in df.columns:
        ratios = df["ratio_new_over_old_epochs"].dropna()
        if not ratios.empty:
            print(
                f"[ok] ratio_new_over_old_epochs min={ratios.min():.3f} "
                f"max={ratios.max():.3f} median={ratios.median():.3f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
