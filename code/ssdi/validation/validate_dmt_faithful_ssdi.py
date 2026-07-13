#!/usr/bin/env python3
"""Validate DMT-faithful SSDI stage outputs and emit a markdown summary."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.io import loadmat

STAGES = ["W_", "N1", "N2", "N3", "R_"]
BANDS = ["delta", "theta", "alpha", "sigma", "beta", "gamma"]
MVALS = [2, 3, 4, 5, 6]


@dataclass
class KeyCheck:
    file: Path
    missing: list[str]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Validate SSDI stage outputs")
    p.add_argument(
        "--subject-measure-root",
        type=Path,
        default=Path("results/gc/stage_mvgc2_anatomical7/data_intermediate/subject_measures"),
    )
    p.add_argument(
        "--ssdi-root",
        type=Path,
        default=Path("results/ssdi/stage_ssdi_anatomical7"),
    )
    p.add_argument(
        "--status-csv",
        type=Path,
        default=Path("results/ssdi/stage_ssdi_anatomical7/logs/stage1_ssdi_status.csv"),
    )
    p.add_argument(
        "--dmt-reference-file",
        type=Path,
        default=Path("/Users/borjan/code/dmt_project/results/ssdi_analysis/rest/P2M63_ssdi_spectral_domain_results.mat"),
    )
    p.add_argument(
        "--out-markdown",
        type=Path,
        default=Path("docs/ssdi_validation_summary.md"),
    )
    return p.parse_args()


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
        raise ValueError(f"results is not dict-like in {path}")
    return res


def count_inputs(root: Path) -> dict[str, int]:
    return {
        st: len(list((root / st).glob("*_subject_measures.mat")))
        for st in STAGES
    }


def count_outputs(root: Path) -> dict[str, int]:
    return {
        st: len(list((root / st).glob("*_ssdi_spectral_domain_results.mat")))
        for st in STAGES
    }


def check_required_keys(result_file: Path) -> KeyCheck:
    res = load_results(result_file)
    missing: list[str] = []
    for band in BANDS:
        for m in MVALS:
            req = [
                f"{band}_preopt_m{m}",
                f"{band}_Loptp_m{m}",
                f"{band}_dopto_m{m}",
                f"{band}_diag_m{m}",
                f"{band}_node_weights_m{m}",
                f"{band}_beta_m{m}",
            ]
            for key in req:
                if key not in res:
                    missing.append(key)
    return KeyCheck(file=result_file, missing=missing)


def find_stage_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for st in STAGES:
        files.extend(sorted((root / st).glob("*_ssdi_spectral_domain_results.mat")))
    return sorted(files)


def shape_summary(res: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    keys = [
        "delta_dopto_m2",
        "delta_doptp_m2",
        "delta_Lopto_m2",
        "delta_Loptp_m2",
        "delta_node_weights_m2",
        "delta_beta_m2",
        "gamma_dopto_m6",
        "gamma_Lopto_m6",
    ]
    for k in keys:
        if k not in res:
            out[k] = None
            continue
        arr = np.asarray(res[k])
        out[k] = {
            "shape": list(arr.shape),
            "dtype": str(arr.dtype),
        }
    return out


def load_runtime_table(status_csv: Path) -> pd.DataFrame:
    if not status_csv.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(status_csv)
    except Exception:
        return pd.DataFrame()


def pick_representative(files: list[Path]) -> Path | None:
    for st in ["N2", "N3", "N1", "W_", "R_"]:
        for f in files:
            if f.parent.name == st:
                return f
    return files[0] if files else None


def parse_subject_stage(path: Path) -> tuple[str, str]:
    stage = path.parent.name
    m = re.match(r"^(EPCTL\d{2})_ssdi_spectral_domain_results\.mat$", path.name)
    subject = m.group(1) if m else path.stem
    return subject, stage


def main() -> int:
    args = parse_args()

    inputs = count_inputs(args.subject_measure_root)
    outputs = count_outputs(args.ssdi_root)
    files = find_stage_files(args.ssdi_root)

    missing_checks: list[KeyCheck] = []
    for f in files:
        kc = check_required_keys(f)
        if kc.missing:
            missing_checks.append(kc)

    runtime_df = load_runtime_table(args.status_csv)

    rep_file = pick_representative(files)
    rep_subject = ""
    rep_stage = ""
    rep_shapes: dict[str, Any] = {}
    rep_meta: dict[str, Any] = {}
    if rep_file is not None:
        rep_subject, rep_stage = parse_subject_stage(rep_file)
        rep = load_results(rep_file)
        rep_shapes = shape_summary(rep)
        meta = rep.get("meta", {}) if isinstance(rep, dict) else {}
        if isinstance(meta, dict):
            def scalar_or_none(name: str):
                if name not in meta:
                    return None
                val = np.asarray(meta[name]).ravel()
                if val.size == 0:
                    return None
                return val[0]

            rep_meta = {
                "n_regions": int(np.asarray(meta.get("n_regions", np.nan)).ravel()[0]) if "n_regions" in meta else None,
                "nrunsp": int(np.asarray(meta.get("nrunsp", np.nan)).ravel()[0]) if "nrunsp" in meta else None,
                "nitersp": int(np.asarray(meta.get("nitersp", np.nan)).ravel()[0]) if "nitersp" in meta else None,
                "niterso": int(np.asarray(meta.get("niterso", np.nan)).ravel()[0]) if "niterso" in meta else None,
                "stage": str(meta.get("stage", rep_stage)),
                "pre_optimiser": str(scalar_or_none("pre_optimiser")) if scalar_or_none("pre_optimiser") is not None else None,
                "pre_optimiser_objective": str(scalar_or_none("pre_optimiser_objective")) if scalar_or_none("pre_optimiser_objective") is not None else None,
                "final_optimiser": str(scalar_or_none("final_optimiser")) if scalar_or_none("final_optimiser") is not None else None,
                "final_optimiser_objective": str(scalar_or_none("final_optimiser_objective")) if scalar_or_none("final_optimiser_objective") is not None else None,
            }

    dmt_summary: dict[str, Any] = {
        "reference_exists": args.dmt_reference_file.exists(),
    }
    if args.dmt_reference_file.exists() and rep_file is not None:
        dmt_res = load_results(args.dmt_reference_file)
        rep_res = load_results(rep_file)
        rep_keys = set(rep_res.keys())
        dmt_keys = set(dmt_res.keys())
        dmt_summary.update(
            {
                "reference_file": str(args.dmt_reference_file),
                "somno_key_count": len(rep_keys),
                "dmt_key_count": len(dmt_keys),
                "shared_keys": len(rep_keys.intersection(dmt_keys)),
                "somno_only_keys": len(rep_keys - dmt_keys),
                "dmt_only_keys": len(dmt_keys - rep_keys),
            }
        )

    total_expected_cases = sum(inputs.values())
    total_output_cases = sum(outputs.values())

    rt_info = {
        "n_status_rows": int(len(runtime_df)),
        "runtime_sec_median": float(runtime_df["runtime_sec"].median()) if "runtime_sec" in runtime_df else np.nan,
        "runtime_sec_min": float(runtime_df["runtime_sec"].min()) if "runtime_sec" in runtime_df else np.nan,
        "runtime_sec_max": float(runtime_df["runtime_sec"].max()) if "runtime_sec" in runtime_df else np.nan,
    }

    md_lines = []
    md_lines.append("# SSDI Validation Summary (DMT-Faithful Stage Pipeline)")
    md_lines.append("")
    md_lines.append(f"Generated: {dt.datetime.now().isoformat(timespec='seconds')}")
    md_lines.append("")
    md_lines.append("## Coverage")
    md_lines.append(f"- Input subject-stage cases found: **{total_expected_cases}**")
    md_lines.append(f"- Output subject-stage cases found: **{total_output_cases}**")
    md_lines.append("- Stage input counts:")
    for st in STAGES:
        md_lines.append(f"  - {st}: {inputs[st]}")
    md_lines.append("- Stage output counts:")
    for st in STAGES:
        md_lines.append(f"  - {st}: {outputs[st]}")

    md_lines.append("")
    md_lines.append("## Optimizer Validation")
    md_lines.append(f"- Macro sizes per case: {MVALS} (count={len(MVALS)})")
    md_lines.append(f"- Bands per case: {BANDS} (count={len(BANDS)})")
    md_lines.append("- Expected optimization blocks per subject-stage case:")
    md_lines.append(f"  - preopt calls: {len(MVALS) * len(BANDS)}")
    md_lines.append(f"  - second-stage calls: {len(MVALS) * len(BANDS)}")
    if rep_meta:
        md_lines.append(
            f"- Representative config from output meta: nrunsp={rep_meta.get('nrunsp')}, "
            f"nitersp={rep_meta.get('nitersp')}, niterso={rep_meta.get('niterso')}"
        )
        if rep_meta.get("final_optimiser_objective") is not None:
            md_lines.append(
                f"- Final objective: `{rep_meta.get('final_optimiser_objective')}` "
                f"(optimiser `{rep_meta.get('final_optimiser')}`)"
            )
        if rep_meta.get("pre_optimiser_objective") is not None:
            md_lines.append(
                f"- Pre objective: `{rep_meta.get('pre_optimiser_objective')}` "
                f"(optimiser `{rep_meta.get('pre_optimiser')}`)"
            )

    md_lines.append("")
    md_lines.append("## Runtime")
    md_lines.append(f"- Status rows: {rt_info['n_status_rows']}")
    if np.isfinite(rt_info["runtime_sec_median"]):
        md_lines.append(
            f"- Runtime seconds (min/median/max): "
            f"{rt_info['runtime_sec_min']:.2f} / {rt_info['runtime_sec_median']:.2f} / {rt_info['runtime_sec_max']:.2f}"
        )

    md_lines.append("")
    md_lines.append("## Structural Checks")
    md_lines.append(f"- Files with missing required DMT-style key families: {len(missing_checks)}")
    if missing_checks:
        first = missing_checks[0]
        md_lines.append(f"- Example missing file: `{first.file}`")
        md_lines.append(f"- Missing keys (first 20): `{first.missing[:20]}`")
    else:
        md_lines.append("- All inspected outputs contained required key families for all band/m pairs.")

    if rep_file is not None:
        md_lines.append("")
        md_lines.append("## Representative Output")
        md_lines.append(f"- File: `{rep_file}`")
        md_lines.append(f"- Subject: `{rep_subject}`")
        md_lines.append(f"- Stage: `{rep_stage}`")
        md_lines.append("- Shape summary:")
        for k, v in rep_shapes.items():
            md_lines.append(f"  - {k}: {json.dumps(v)}")

    md_lines.append("")
    md_lines.append("## DMT Structure Comparison")
    if dmt_summary.get("reference_exists"):
        md_lines.append(f"- Reference file: `{dmt_summary.get('reference_file')}`")
        md_lines.append(f"- Somno key count: {dmt_summary.get('somno_key_count')}")
        md_lines.append(f"- DMT key count: {dmt_summary.get('dmt_key_count')}")
        md_lines.append(f"- Shared keys: {dmt_summary.get('shared_keys')}")
        md_lines.append(f"- Somno-only keys: {dmt_summary.get('somno_only_keys')}")
        md_lines.append(f"- DMT-only keys: {dmt_summary.get('dmt_only_keys')}")
        md_lines.append("- Note: available DMT reference output appears from an older result set; key-count differences are expected.")
    else:
        md_lines.append("- DMT reference output file not found; structural comparison skipped.")

    args.out_markdown.parent.mkdir(parents=True, exist_ok=True)
    args.out_markdown.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"Wrote validation summary: {args.out_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
