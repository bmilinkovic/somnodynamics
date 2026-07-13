#!/usr/bin/env python3
"""
Validation utilities for Yeo-7 source reconstruction outputs.

Checks:
1) Shape / finite values.
2) Stable label ordering across files.
3) Basic amplitude sanity.
4) Optional parity comparison with dmt_project output.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np

try:
    from scipy.io import loadmat
    SCIPY_IMPORT_ERROR = None
except ImportError as exc:  # pragma: no cover - environment dependent
    loadmat = None
    SCIPY_IMPORT_ERROR = exc


EXPECTED_LABEL_COUNT = 7
UPPER_TRIANGLE_CORR_THRESHOLD = 0.80
FROBENIUS_THRESHOLD = 1.50


@dataclass
class FileCheckResult:
    path: Path
    ok: bool
    message: str
    labels: List[str]
    shape: Tuple[int, int, int]
    median_abs: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate Yeo-7 source reconstruction")
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path(
            "/Users/borjan/code/healthy-sleep-project/code/results/source_recon/yeo7"
        ),
        help="Root directory containing subject/mat/*.mat files.",
    )
    parser.add_argument(
        "--sleep-file",
        type=Path,
        default=None,
        help="Optional sleep-project .mat file for parity check.",
    )
    parser.add_argument(
        "--dmt-file",
        type=Path,
        default=None,
        help="Optional dmt_project .mat file for parity check.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return non-zero exit code on warnings/failures.",
    )
    return parser.parse_args()


def normalize_labels(raw_labels: object) -> List[str]:
    labels_array = np.asarray(raw_labels).squeeze()
    if labels_array.ndim == 0:
        return [str(labels_array.item()).strip()]
    labels: List[str] = []
    for value in labels_array.tolist():
        if isinstance(value, bytes):
            labels.append(value.decode("utf-8").strip())
        else:
            labels.append(str(value).strip())
    return labels


def load_source_mat(path: Path) -> Tuple[np.ndarray, List[str]]:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "source_ts" in mat:
        source_ts = np.asarray(mat["source_ts"])
    elif "data" in mat:
        source_ts = np.asarray(mat["data"])
    else:
        raise KeyError(f"{path} does not contain 'source_ts' or 'data'.")

    if "labels" not in mat:
        raise KeyError(f"{path} does not contain 'labels'.")
    labels = normalize_labels(mat["labels"])
    return source_ts, labels


def corr_matrix_from_source(source_ts: np.ndarray) -> np.ndarray:
    if source_ts.ndim != 3:
        raise ValueError(f"Expected 3D source_ts, got shape {source_ts.shape}")
    flattened = source_ts.reshape(source_ts.shape[0], -1)
    flattened = flattened - flattened.mean(axis=1, keepdims=True)
    std = flattened.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    flattened = flattened / std
    return np.corrcoef(flattened)


def scan_source_files(source_root: Path) -> List[Path]:
    return sorted(source_root.glob("*/mat/*_source_time_series.mat"))


def run_sanity_checks(source_root: Path) -> Tuple[List[FileCheckResult], int]:
    mat_files = scan_source_files(source_root)
    if not mat_files:
        raise FileNotFoundError(f"No source .mat files found under {source_root}")

    results: List[FileCheckResult] = []
    reference_labels: List[str] | None = None
    warning_count = 0

    for path in mat_files:
        try:
            source_ts, labels = load_source_mat(path)
            if source_ts.ndim != 3:
                results.append(
                    FileCheckResult(
                        path=path,
                        ok=False,
                        message=f"invalid shape {source_ts.shape}",
                        labels=labels,
                        shape=tuple(source_ts.shape),
                        median_abs=float("nan"),
                    )
                )
                warning_count += 1
                continue

            if source_ts.shape[0] != EXPECTED_LABEL_COUNT:
                warning_count += 1
                ok = False
                message = (
                    f"expected {EXPECTED_LABEL_COUNT} regions, got {source_ts.shape[0]}"
                )
            else:
                ok = True
                message = "ok"

            if len(labels) != EXPECTED_LABEL_COUNT:
                ok = False
                warning_count += 1
                message = f"expected {EXPECTED_LABEL_COUNT} labels, got {len(labels)}"

            if reference_labels is None:
                reference_labels = labels
            elif labels != reference_labels:
                ok = False
                warning_count += 1
                message = "label ordering mismatch vs first file"

            if not np.isfinite(source_ts).all():
                ok = False
                warning_count += 1
                message = "contains NaN/Inf values"

            median_abs = float(np.median(np.abs(source_ts)))
            if median_abs <= 0:
                ok = False
                warning_count += 1
                message = "zero-amplitude median"

            results.append(
                FileCheckResult(
                    path=path,
                    ok=ok,
                    message=message,
                    labels=labels,
                    shape=tuple(source_ts.shape),
                    median_abs=median_abs,
                )
            )
        except Exception as exc:  # pylint: disable=broad-except
            warning_count += 1
            results.append(
                FileCheckResult(
                    path=path,
                    ok=False,
                    message=f"error: {exc}",
                    labels=[],
                    shape=(0, 0, 0),
                    median_abs=float("nan"),
                )
            )

    medians = np.array([r.median_abs for r in results if np.isfinite(r.median_abs)])
    if medians.size >= 3:
        log_medians = np.log10(np.maximum(medians, np.finfo(float).eps))
        std = np.std(log_medians)
        if std > 0:
            zscores = (log_medians - np.mean(log_medians)) / std
            outliers = int(np.sum(np.abs(zscores) > 3.0))
            warning_count += outliers

    return results, warning_count


def parity_check(sleep_file: Path, dmt_file: Path) -> Tuple[int, List[str]]:
    messages: List[str] = []
    warnings = 0

    sleep_ts, sleep_labels = load_source_mat(sleep_file)
    dmt_ts, dmt_labels = load_source_mat(dmt_file)

    if len(sleep_labels) != EXPECTED_LABEL_COUNT or len(dmt_labels) != EXPECTED_LABEL_COUNT:
        warnings += 1
        messages.append("WARN: label count is not 7 in one or both files.")

    if set(sleep_labels) != set(dmt_labels):
        warnings += 1
        messages.append("WARN: sleep/dmt label sets do not match.")
        return warnings, messages

    common_order = sleep_labels
    dmt_reindex = [dmt_labels.index(label) for label in common_order]
    dmt_ts = dmt_ts[dmt_reindex, :, :]

    sleep_corr = corr_matrix_from_source(sleep_ts)
    dmt_corr = corr_matrix_from_source(dmt_ts)

    frob = float(np.linalg.norm(sleep_corr - dmt_corr, ord="fro"))
    iu = np.triu_indices(EXPECTED_LABEL_COUNT, k=1)
    upper_corr = float(np.corrcoef(sleep_corr[iu], dmt_corr[iu])[0, 1])

    messages.append(f"Parity labels: PASS ({EXPECTED_LABEL_COUNT} labels)")
    messages.append(f"Parity Frobenius diff: {frob:.4f}")
    messages.append(f"Parity upper-triangle corr: {upper_corr:.4f}")

    if frob > FROBENIUS_THRESHOLD:
        warnings += 1
        messages.append(
            f"WARN: Frobenius difference {frob:.4f} exceeds threshold {FROBENIUS_THRESHOLD:.2f}."
        )
    else:
        messages.append("Parity Frobenius status: PASS")

    if not np.isfinite(upper_corr) or upper_corr < UPPER_TRIANGLE_CORR_THRESHOLD:
        warnings += 1
        messages.append(
            "WARN: upper-triangle correlation "
            f"{upper_corr:.4f} below threshold {UPPER_TRIANGLE_CORR_THRESHOLD:.2f}."
        )
    else:
        messages.append("Parity upper-triangle correlation status: PASS")

    return warnings, messages


def main() -> int:
    args = parse_args()
    if SCIPY_IMPORT_ERROR is not None:
        raise ModuleNotFoundError(
            "scipy is required for validation. Install dependencies first."
        ) from SCIPY_IMPORT_ERROR

    results, warning_count = run_sanity_checks(args.source_root)

    print("=" * 72)
    print("YEO7 SANITY CHECKS")
    print("=" * 72)
    for result in results:
        status = "PASS" if result.ok else "WARN"
        print(
            f"{status} | {result.path.name} | shape={result.shape} | "
            f"median|x|={result.median_abs:.6g} | {result.message}"
        )

    print("-" * 72)
    print(f"Files checked: {len(results)}")
    print(f"Warnings: {warning_count}")

    if args.sleep_file and args.dmt_file:
        parity_warnings, messages = parity_check(args.sleep_file, args.dmt_file)
        warning_count += parity_warnings
        print("=" * 72)
        print("CROSS-PROJECT PARITY CHECK")
        print("=" * 72)
        for message in messages:
            print(message)
        print(f"Parity warnings: {parity_warnings}")

    print("=" * 72)
    print("LEGACY PIPELINE STATUS")
    print("=" * 72)
    print("Legacy scripts were not modified by this validation tool.")

    if args.strict and warning_count > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
