#!/usr/bin/env python3
"""Compare SomnoDynamics outputs with legacy outputs where available."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config_utils import external_results_root, resolve_repo_path  # noqa: E402

try:
    from scipy.io import loadmat
except Exception:  # pragma: no cover
    loadmat = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare new outputs against legacy outputs")
    parser.add_argument("--repo-root", type=Path, default=resolve_repo_path("."))
    parser.add_argument("--report-csv", type=Path, default=resolve_repo_path("results/legacy_comparison_report.csv"))
    parser.add_argument("--summary-txt", type=Path, default=resolve_repo_path("results/legacy_comparison_summary.txt"))
    return parser.parse_args()


def _norm(p: Path, root: Path) -> Path:
    return p if p.is_absolute() else root / p


def _list_stage_files(root: Path, pattern: str) -> set[tuple[str, str]]:
    out: set[tuple[str, str]] = set()
    for f in root.glob(pattern):
        # New layout: <root>/<stage>/<subject>_ssdi_spectral_domain_results.mat
        if f.parent.name in {"W_", "N1", "N2", "N3", "R_"}:
            m = re.match(r"^(EPCTL\d{2})_ssdi_spectral_domain_results\\.mat$", f.name)
            if m:
                out.add((m.group(1), f.parent.name))
                continue

        # Legacy/new-fast layout: <subject>_<stage>_...
        m = re.match(r"^(EPCTL\d{2})_(W_|N1|N2|N3|R_)", f.name)
        if m:
            out.add((m.group(1), m.group(2)))
    return out


def _safe_load_mat(path: Path, var: str) -> Any:
    if loadmat is None or not path.exists():
        return None
    try:
        m = loadmat(path, squeeze_me=True, struct_as_record=False)
    except Exception:
        return None
    return m.get(var)


def compare_gc(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    new_dir = external_results_root() / "gc"
    legacy_candidates = [
        root / "code/legacy/results/new_results_GC",
    ]
    legacy_dir = next((d for d in legacy_candidates if d.exists()), None)

    keys = [
        "gcin_time_subject_data_struct.mat",
        "gcout_time_subject_data_struct.mat",
        "gcin_freq_subject_data_struct.mat",
        "gcout_freq_subject_data_struct.mat",
    ]
    for k in keys:
        rows.append(
            {
                "domain": "gc",
                "artifact": k,
                "new_exists": (new_dir / k).exists(),
                "legacy_exists": (legacy_dir / k).exists() if legacy_dir else False,
                "match": (new_dir / k).exists() and ((legacy_dir / k).exists() if legacy_dir else False),
                "details": str(legacy_dir) if legacy_dir else "legacy dir not found",
            }
        )
    return rows


def compare_ssdi(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    new_spec = external_results_root() / "ssdi/stage_ssdi_anatomical7"
    legacy_candidates = [
        root / "code/legacy/results/ssdi_results/spectral",
    ]
    legacy_spec = next((d for d in legacy_candidates if d.exists()), None)

    new_set = _list_stage_files(new_spec, "*/*_ssdi_spectral_domain_results.mat") if new_spec.exists() else set()
    if not new_set and new_spec.exists():
        new_set = _list_stage_files(new_spec, "EPCTL*_ssdi_spectral_domain_results.mat")
    old_set = _list_stage_files(legacy_spec, "EPCTL*_ssdi_spectral*.mat") if legacy_spec else set()

    rows.append(
        {
            "domain": "ssdi",
            "artifact": "subject_stage_coverage_spectral",
            "new_exists": new_spec.exists(),
            "legacy_exists": bool(legacy_spec and legacy_spec.exists()),
            "match": len(new_set.intersection(old_set)) > 0,
            "details": f"new={len(new_set)} legacy={len(old_set)} overlap={len(new_set.intersection(old_set))}",
        }
    )
    return rows


def compare_sensor_power(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    new_csv = external_results_root() / "sensor_power/sensor_bandpower_long.csv"
    legacy_candidates = [
        root / "code/legacy/results/sleep-results/bandpower",
    ]
    legacy_dir = next((d for d in legacy_candidates if d.exists()), None)

    n_rows = 0
    if new_csv.exists():
        try:
            n_rows = len(pd.read_csv(new_csv))
        except Exception:
            n_rows = -1

    rows.append(
        {
            "domain": "sensor_power",
            "artifact": "sensor_bandpower_long.csv",
            "new_exists": new_csv.exists(),
            "legacy_exists": bool(legacy_dir and legacy_dir.exists()),
            "match": new_csv.exists() and bool(legacy_dir and legacy_dir.exists()),
            "details": f"new_rows={n_rows} legacy_dir={legacy_dir}",
        }
    )
    return rows


def compare_cross_analysis(root: Path) -> list[dict[str, Any]]:
    out = external_results_root() / "cross_analysis/gc_ssdi_correlations.csv"
    return [
        {
            "domain": "cross_analysis",
            "artifact": "gc_ssdi_correlations.csv",
            "new_exists": out.exists(),
            "legacy_exists": False,
            "match": out.exists(),
            "details": "No direct legacy equivalent; existence check only.",
        }
    ]


def main() -> int:
    args = parse_args()
    root = args.repo_root

    rows: list[dict[str, Any]] = []
    rows.extend(compare_sensor_power(root))
    rows.extend(compare_gc(root))
    rows.extend(compare_ssdi(root))
    rows.extend(compare_cross_analysis(root))

    df = pd.DataFrame(rows)

    report_csv = _norm(args.report_csv, root)
    summary_txt = _norm(args.summary_txt, root)
    report_csv.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(report_csv, index=False)

    with summary_txt.open("w", encoding="utf-8") as handle:
        handle.write("SomnoDynamics vs Legacy Comparison Summary\n")
        handle.write("========================================\n")
        handle.write(f"Checks: {len(df)}\n")
        handle.write(f"Matched checks: {int(df['match'].sum())}\n")
        handle.write("\nDetails:\n")
        for _, row in df.iterrows():
            handle.write(
                f"- [{row['domain']}] {row['artifact']}: new={row['new_exists']} "
                f"legacy={row['legacy_exists']} match={row['match']} | {row['details']}\n"
            )

    print(f"[SomnoDynamics] legacy comparison report -> {report_csv}")
    print(f"[SomnoDynamics] legacy comparison summary -> {summary_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
