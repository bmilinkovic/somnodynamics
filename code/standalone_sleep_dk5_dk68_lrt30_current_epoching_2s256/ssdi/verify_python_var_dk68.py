#!/usr/bin/env python3
"""Check every compact DK68 VAR model before transferring it to the cluster."""

from __future__ import annotations

import argparse
import os
import zipfile
from pathlib import Path

import numpy as np

from run_python_dd_dk68 import discover, discover_models


def main() -> int:
    data_root = Path(os.environ.get("SOMNODYNAMICS_DATA_ROOT",
                            "/Volumes/ex_data/data_somnoemergence"))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=data_root /
                        "results/source_recon/dk5_dk68_lcmv_nai_current_epoching_2s_256hz/dk68")
    parser.add_argument("--model-root", type=Path, default=data_root /
                        "results/ssdi/stage_ssdi_dk68_python_lrt30_current_epoching_2s_256hz/model_inputs")
    parser.add_argument("--expected-count", type=int, default=143)
    args = parser.parse_args()

    sources = {(subject, stage): path for subject, stage, path in discover(args.source_root)}
    models = {(subject, stage): path for subject, stage, path in discover_models(args.model_root)}
    if len(sources) != args.expected_count:
        raise SystemExit(f"expected {args.expected_count} source files, found {len(sources)}")
    missing = sorted(set(sources) - set(models))
    extra = sorted(set(models) - set(sources))
    errors = []
    for subject, stage in sorted(set(sources) & set(models)):
        source, path = sources[(subject, stage)], models[(subject, stage)]
        try:
            with np.load(path, allow_pickle=False) as saved:
                required = ("A", "V", "order", "labels", "subject", "stage", "fs",
                            "pmax", "criterion", "regmode", "source_file", "source_size",
                            "source_mtime_ns", "samples", "trials")
                absent = set(required) - set(saved.files)
                if absent:
                    raise ValueError(f"missing fields: {sorted(absent)}")
                order = int(saved["order"])
                a, v = saved["A"], saved["V"]
                if not 1 <= order <= 30 or a.shape != (68, 68, order) or v.shape != (68, 68):
                    raise ValueError(f"invalid VAR dimensions: A={a.shape}, V={v.shape}, p={order}")
                if not np.isfinite(a).all() or not np.isfinite(v).all():
                    raise ValueError("non-finite A or V")
                if not np.allclose(v, v.T, atol=1e-8, rtol=1e-8):
                    raise ValueError("V is not symmetric")
                np.linalg.cholesky(v)
                if saved["labels"].shape != (68,):
                    raise ValueError("expected 68 region labels")
                metadata = (str(saved["subject"]), str(saved["stage"]), float(saved["fs"]),
                            int(saved["pmax"]), str(saved["criterion"]), str(saved["regmode"]))
                if metadata != (subject, stage, 256.0, 30, "LRT", "LWR"):
                    raise ValueError(f"unexpected fitting metadata: {metadata}")
                stat = source.stat()
                if (str(saved["source_file"]) != str(source)
                        or int(saved["source_size"]) != stat.st_size
                        or int(saved["source_mtime_ns"]) != stat.st_mtime_ns):
                    raise ValueError("source provenance differs from the current source file")
                if int(saved["samples"]) != 512 or int(saved["trials"]) < 1:
                    raise ValueError("invalid sample or trial count")
        except (OSError, EOFError, ValueError, KeyError, zipfile.BadZipFile,
                np.linalg.LinAlgError) as exc:
            errors.append(f"{subject} {stage}: {exc}")
    print(f"Verified {len(models) - len(extra) - len(errors)}/{args.expected_count} DK68 models")
    if missing:
        print(f"Missing {len(missing)}; first ten: {missing[:10]}")
    if extra:
        print(f"Unexpected {len(extra)}: {extra}")
    for error in errors:
        print(f"Invalid: {error}")
    return 1 if missing or extra or errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
