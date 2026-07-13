#!/usr/bin/env python3
"""Legacy Yeo-7 source reconstruction wrapper.

SomnoDynamics primary source branch is anatomical46 -> anatomical7.
This wrapper is preserved only for compatibility.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run_all"))
from config_utils import load_config, resolve_repo_path  # noqa: E402


def parse_args() -> argparse.Namespace:
    cfg = load_config()
    parser = argparse.ArgumentParser(description="Run legacy Yeo-7 source reconstruction")
    parser.add_argument("--raw-root", type=Path, default=Path(cfg["raw_data_root"]))
    parser.add_argument("--data-dir", type=Path, default=Path(cfg["raw_data_root"]))
    parser.add_argument("--montage-file", type=Path, default=Path(cfg["montage_file"]))
    parser.add_argument(
        "--out-root",
        type=Path,
        default=resolve_repo_path(cfg["paths"]["source_recon_anatomical7"]) / "legacy_yeo7",
    )
    parser.add_argument("--subjects-dir", type=Path, default=Path("/Users/borjan/mne_data/MNE-fsaverage-data"))
    parser.add_argument("--sfreq", type=float, default=256.0)
    parser.add_argument("--epoch-seconds", type=float, default=2.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    legacy_script = resolve_repo_path("code/legacy/source/source_reconstruction_yeo7.py")
    if not legacy_script.exists():
        raise FileNotFoundError(f"Legacy script not found: {legacy_script}")

    cmd = [
        sys.executable,
        str(legacy_script),
        "--raw-dir",
        str(args.raw_root),
        "--data-dir",
        str(args.data_dir),
        "--out-root",
        str(args.out_root),
        "--subjects-dir",
        str(args.subjects_dir),
        "--sfreq",
        str(args.sfreq),
        "--epoch-seconds",
        str(args.epoch_seconds),
    ]
    if args.montage_file is not None:
        cmd.extend(["--montage-file", str(args.montage_file)])

    print("[SomnoDynamics] invoking legacy source reconstruction:")
    print(" ", " ".join(cmd))
    proc = subprocess.run(cmd)
    return int(proc.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
