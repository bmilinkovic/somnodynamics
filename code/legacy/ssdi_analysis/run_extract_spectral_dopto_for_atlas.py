#!/usr/bin/env python3
"""
Atlas-aware wrapper for extract_spectral_dopto.py.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract spectral dopto vectors for an atlas")
    parser.add_argument(
        "--atlas",
        choices=["legacy", "yeo7"],
        default="yeo7",
    )
    parser.add_argument(
        "--spectral-dir",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
    )
    return parser.parse_args()


def default_dirs(atlas: str) -> tuple[Path, Path]:
    if atlas == "legacy":
        spectral = Path("/Users/borjan/code/healthy-sleep-project/code/results/ssdi_results/spectral")
        output = Path("/Users/borjan/code/healthy-sleep-project/code/ssdi_analysis/data")
    else:
        spectral = Path("/Users/borjan/code/healthy-sleep-project/code/results/yeo7/ssdi_results/spectral")
        output = Path("/Users/borjan/code/healthy-sleep-project/code/ssdi_analysis/data/yeo7")
    return spectral, output


def load_module(module_path: Path):
    spec = importlib.util.spec_from_file_location("extract_spectral_dopto", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> int:
    args = parse_args()
    module_path = Path("/Users/borjan/code/healthy-sleep-project/code/extract_spectral_dopto.py")
    module = load_module(module_path)

    spectral_default, output_default = default_dirs(args.atlas)
    module.SPECTRAL_DIR = args.spectral_dir or spectral_default
    module.OUTPUT_DIR = args.output_dir or output_default
    module.main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
