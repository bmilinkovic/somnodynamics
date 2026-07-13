#!/usr/bin/env python3
"""
Atlas-aware wrapper for ssdi_cross_frequency_tests.py.

This reuses the existing cross-frequency SSDI script without changing its
legacy defaults, by patching path globals at runtime.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run cross-frequency SSDI for an atlas")
    parser.add_argument(
        "--atlas",
        choices=["legacy", "yeo7"],
        default="yeo7",
    )
    parser.add_argument(
        "--spectral-dir",
        type=Path,
        default=None,
        help="Override spectral SSDI directory.",
    )
    parser.add_argument(
        "--time-dir",
        type=Path,
        default=None,
        help="Override time-domain SSDI directory.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=None,
        help="Override figure output directory.",
    )
    return parser.parse_args()


def default_dirs(atlas: str) -> tuple[Path, Path]:
    if atlas == "legacy":
        base = Path("/Users/borjan/code/healthy-sleep-project/code/results/ssdi_results")
    else:
        base = Path("/Users/borjan/code/healthy-sleep-project/code/results/yeo7/ssdi_results")
    return base / "spectral", base / "time"


def load_module(module_path: Path):
    spec = importlib.util.spec_from_file_location("ssdi_cross_frequency_tests", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> int:
    args = parse_args()
    module_path = Path("/Users/borjan/code/healthy-sleep-project/code/ssdi_cross_frequency_tests.py")
    module = load_module(module_path)

    spectral_default, time_default = default_dirs(args.atlas)
    spectral_dir = args.spectral_dir or spectral_default
    time_dir = args.time_dir or time_default
    fig_dir = args.fig_dir or (spectral_dir / "figs_crossfreq")

    module.SPECTRAL_DIR = spectral_dir
    module.TIME_DIR = time_dir
    module.FIG_DIR = fig_dir
    module.main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
