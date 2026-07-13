#!/usr/bin/env python3
"""Shared configuration utilities for SomnoDynamics."""

from __future__ import annotations

from functools import lru_cache
import os
from pathlib import Path
from typing import Any, Dict

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def config_path() -> Path:
    return Path(__file__).resolve().with_name("somnodynamics_config.toml")


@lru_cache(maxsize=1)
def load_config() -> Dict[str, Any]:
    with config_path().open("rb") as handle:
        cfg = tomllib.load(handle)
    cfg["_repo_root"] = str(repo_root())
    cfg["_external_data_root"] = str(external_data_root(cfg))
    cfg["_external_results_root"] = str(external_results_root(cfg))
    return cfg


def external_data_root(cfg: Dict[str, Any] | None = None) -> Path:
    if cfg is None:
        with config_path().open("rb") as handle:
            cfg = tomllib.load(handle)
    return Path(os.environ.get("SOMNODYNAMICS_DATA_ROOT", cfg["external_data_root"])).expanduser()


def external_results_root(cfg: Dict[str, Any] | None = None) -> Path:
    if cfg is None:
        with config_path().open("rb") as handle:
            cfg = tomllib.load(handle)
    results_root = Path(str(cfg.get("paths", {}).get("results_root", "results")))
    if results_root.is_absolute():
        return results_root
    return external_data_root(cfg) / results_root


def resolve_repo_path(path_value: str) -> Path:
    p = Path(path_value)
    if p.is_absolute():
        return p
    if p.parts and p.parts[0] == "results":
        return external_results_root() / Path(*p.parts[1:])
    return repo_root() / p


def stage_alias_map() -> Dict[str, str]:
    return {str(k): str(v) for k, v in load_config()["stage_label_aliases"].items()}


def stage_display_map() -> Dict[str, str]:
    return {str(k): str(v) for k, v in load_config()["stage_labels"].items()}


def canonical_stage(stage: str) -> str:
    aliases = stage_alias_map()
    s = str(stage).strip()
    if s in aliases:
        return aliases[s]
    s_no = s.rstrip("_")
    if s_no in aliases:
        return aliases[s_no]
    raise ValueError(f"Unknown sleep stage label: {stage}")


def stage_display(stage: str) -> str:
    return stage_display_map()[canonical_stage(stage)]


def stage_order() -> list[str]:
    return ["W_", "N1", "N2", "N3", "R_"]


if __name__ == "__main__":
    cfg = load_config()
    print(f"Loaded config for: {cfg['project_name']}")
    print(f"Raw data root: {cfg['raw_data_root']}")
