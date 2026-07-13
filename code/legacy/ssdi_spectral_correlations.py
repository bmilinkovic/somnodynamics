#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ssdi_spectral_correlations.py

Correlation analysis between:
  - band power (P)
  - band-limited dynamical dependence (DD)
  - band-limited structure–function index (SF)

All DD and SF come from spectral SSDI .mat files:
  DD:  first element of <band>_dopto_mn
  SF:  proportion of near-zero entries in <band>_gopto_mn

Band power is loaded from a separate band_power.npz, expected structure:
  np.load(...) -> keys are stages ("W_", "N1", ...), each value is an
  object storing a nested dict:
     power[stage][subject][band] = float

You only need to adapt the CONFIG paths and make sure band_power.npz
is prepared as described.

Author: Borjan Milinkovic
"""

from pathlib import Path
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.stats import spearmanr

# Threshold for "near-zero" in the structure–function index
EPSILON_SF = 0.1

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------

SPECTRAL_DIR = Path(
    "/Users/borjan/code/healthy-sleep-project/code/results/ssdi_results/_spectral"
)

BAND_POWER_FILE = SPECTRAL_DIR.parent / "band_power.npz"  # adjust path if needed

BANDS = ["delta", "theta", "alpha", "sigma", "beta", "gamma"]
STAGES = ["W_", "N1", "N2", "N3", "R_"]
MACRO_SIZES = list(range(2, 10))

FIG_DIR = SPECTRAL_DIR / "figs_correlations"
plt.style.use("seaborn-v0_8-whitegrid")


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def _matlab_struct_to_dict(matobj) -> Dict:
    """
    Recursively convert a MATLAB struct (as returned by scipy.io.loadmat)
    into a nested Python dict. Non-struct entries are returned as-is.

    Assumes loadmat(..., squeeze_me=True, struct_as_record=False).
    """
    from scipy.io.matlab.mio5_params import mat_struct

    if isinstance(matobj, mat_struct):
        out = {}
        for f in matobj._fieldnames:
            elem = getattr(matobj, f)
            out[f] = _matlab_struct_to_dict(elem)
        return out
    elif isinstance(matobj, np.ndarray) and matobj.dtype == object:
        # cell array -> list, with recursive conversion
        return [_matlab_struct_to_dict(e) for e in matobj]
    else:
        return matobj


def structure_function_index(D: np.ndarray, epsilon: float = EPSILON_SF) -> float:
    """
    Compute structure–function (dynamical organisation) index from
    optima distance matrix D (shape [N, N]) as the proportion of
    near-zero entries in the upper triangle (excluding diagonal).
    """
    D = np.asarray(D)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D must be a square matrix.")
    iu = np.triu_indices_from(D, k=1)
    vals = D[iu]
    if vals.size == 0:
        return np.nan
    return float(np.mean(vals < epsilon))


def parse_subject_stage(filename: str) -> Tuple[str, str]:
    m = re.match(r"EPCTL(\d+)_(\w+)_ssdi_.*\.mat", filename)
    if not m:
        raise ValueError(f"Cannot parse subject/stage from {filename}")
    subj = f"{int(m.group(1)):02d}"
    stage = m.group(2)
    return subj, stage


def load_results_mat(path: Path) -> Dict[str, np.ndarray]:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    res_struct = mat["results"]
    return _matlab_struct_to_dict(res_struct)


def load_band_power(path: Path) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Load band power nested dict:

    power[stage][subject][band] = float
    """
    npz = np.load(path, allow_pickle=True)
    power = {}
    for stage in npz.files:
        power[stage] = npz[stage].item()
    return power


@dataclass
class CorrRecord:
    stage: str
    subject: str
    n: int
    band: str
    power: float
    dd: float
    sf: float


# ---------------------------------------------------------------------
# Build long table
# ---------------------------------------------------------------------
def build_records() -> List[CorrRecord]:
    power = load_band_power(BAND_POWER_FILE)
    spectral_files = sorted(SPECTRAL_DIR.glob("EPCTL*_ssdi_spectral_results.mat"))

    records: List[CorrRecord] = []

    for path in spectral_files:
        subj, stage = parse_subject_stage(path.name)
        if stage not in STAGES:
            continue
        if stage not in power or subj not in power[stage]:
            # no power info for this subject/stage
            continue

        res = load_results_mat(path)

        for n in MACRO_SIZES:
            for band in BANDS:
                dopto_key = f"{band}_dopto_m{n}"
                gopto_key = f"{band}_gopto_m{n}"

                if dopto_key not in res or gopto_key not in res:
                    continue

                dopto = np.asarray(res[dopto_key]).ravel()
                dd = float(dopto[0])

                gopto = np.asarray(res[gopto_key])
                sf = structure_function_index(gopto)

                try:
                    P = float(power[stage][subj][band])
                except KeyError:
                    P = np.nan

                records.append(
                    CorrRecord(
                        stage=stage,
                        subject=subj,
                        n=n,
                        band=band,
                        power=P,
                        dd=dd,
                        sf=sf,
                    )
                )

    return records


# ---------------------------------------------------------------------
# Correlation + plotting
# ---------------------------------------------------------------------
def corr_by_band_stage(
    records: List[CorrRecord],
    x_attr: str,
    y_attr: str,
) -> Dict[Tuple[str, str], Tuple[float, float, int]]:
    """
    Compute Spearman correlation between attributes x_attr and y_attr
    for each (band, stage).
    """
    out: Dict[Tuple[str, str], Tuple[float, float, int]] = {}
    for band in BANDS:
        for stage in STAGES:
            xs = []
            ys = []
            for r in records:
                if r.band == band and r.stage == stage:
                    x = getattr(r, x_attr)
                    y = getattr(r, y_attr)
                    if np.isfinite(x) and np.isfinite(y):
                        xs.append(x)
                        ys.append(y)
            xs = np.array(xs)
            ys = np.array(ys)
            if len(xs) < 4:
                out[(band, stage)] = (np.nan, np.nan, len(xs))
            else:
                rho, p = spearmanr(xs, ys)
                out[(band, stage)] = (float(rho), float(p), len(xs))
    return out


def plot_corr_bars(
    corr: Dict[Tuple[str, str], Tuple[float, float, int]],
    title: str,
    outpath: Path,
) -> None:
    bands = BANDS
    stages = STAGES
    nb = len(bands)
    ns = len(stages)

    R = np.full((nb, ns), np.nan)
    for bi, band in enumerate(bands):
        for si, stage in enumerate(stages):
            R[bi, si] = corr.get((band, stage), (np.nan, np.nan, 0))[0]

    x = np.arange(nb)
    width = 0.8 / ns
    colors = plt.cm.tab10(np.linspace(0, 1, ns))

    fig, ax = plt.subplots(figsize=(10, 5))
    for si, stage in enumerate(stages):
        offsets = x - 0.4 + width / 2 + si * width
        ax.bar(
            offsets,
            R[:, si],
            width=width,
            label=stage,
            color=colors[si],
            alpha=0.85,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(bands)
    ax.axhline(0, color="k", linestyle="--", lw=1)
    ax.set_ylabel("Spearman ρ")
    ax.set_title(title)
    ax.legend(title="Stage", loc="lower right")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("Building records (P, DD, SF) ...")
    records = build_records()
    print(f"Total records: {len(records)}")

    # Power vs DD
    corr_P_DD = corr_by_band_stage(records, "power", "dd")
    plot_corr_bars(
        corr_P_DD,
        title="Band power vs DD (Spearman ρ)",
        outpath=FIG_DIR / "corr_P_vs_DD.png",
    )

    # Power vs SF
    corr_P_SF = corr_by_band_stage(records, "power", "sf")
    plot_corr_bars(
        corr_P_SF,
        title="Band power vs SF index (Spearman ρ)",
        outpath=FIG_DIR / "corr_P_vs_SF.png",
    )

    # DD vs SF
    corr_DD_SF = corr_by_band_stage(records, "dd", "sf")
    plot_corr_bars(
        corr_DD_SF,
        title="DD vs SF index (Spearman ρ)",
        outpath=FIG_DIR / "corr_DD_vs_SF.png",
    )

    print(f"Saved correlation figures to: {FIG_DIR}")


if __name__ == "__main__":
    main()
