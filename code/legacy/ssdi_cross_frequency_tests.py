#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ssdi_spectral_crossfreq.py

Cross-frequency SSDI analysis on MATLAB .mat files.

1) Approximate additivity:
   For each subject, stage, and macro size n:
     - broadband DD from time-domain file (first element of <broad_prefix>_dopto_mn)
     - band-wise DD from spectral file (first element of <band>_dopto_mn)
   Compare DD_broadband vs c * sum_bands DD_band, where c is a global scale factor.

2) Macro alignment across bands:
   For each subject, stage, and macro size n:
     - extract macroscopic subspace from each band-specific L-matrix
       (field: <band>_Lopto_m; take first n columns of the cell for that scale)
     - compute subspace distances between all band pairs using principal angles.

3) Band-level contributions and single-band synergy (for a target macro size):
   - contribution_b = (c * DD_b) / DD_broad
   - synergy_single = DD_broad - max_b (c * DD_b)

Dynamical organisation (structure–function index) is always computed from the
optima distance matrices <freq>_gopto_mn, as proportion of near-zero entries.

Author: Borjan Milinkovic (template by ChatGPT)
"""

from pathlib import Path
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.linalg import subspace_angles

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------

# Directory containing the spectral SSDI results (.mat)
SPECTRAL_DIR = Path(
    "/Users/borjan/code/healthy-sleep-project/code/results/ssdi_results/spectral"
)

# Directory + filename pattern for the *time-domain* broadband SSDI results.
TIME_DIR = Path(
    "/Users/borjan/code/healthy-sleep-project/code/results/ssdi_results"
) / "time"

TIME_PATTERN = "EPCTL{subject}_{stage}_ssdi_time_domain_results.mat"

# Time-domain files use 'dopto_m{n}' without prefix
BROADBAND_PREFIX = ""  # empty string means no prefix

# Frequency bands in the spectral SSDI .mat files
BANDS = ["delta", "theta", "alpha", "sigma", "beta", "gamma"]

# Sleep stages (as they appear in filenames)
STAGES = ["W_", "N1", "N2", "N3", "R_"]

# Macroscopic sizes to analyse
MACRO_SIZES = list(range(2, 10))  # 2..9

# For the band-contribution / synergy analysis, choose one representative n
TARGET_N = 6  # <-- change if you want another macro size

# L-matrices are in a 1x8 cell; we assume cell index = n-1 for n in 2..9
def lopto_cell_index(n: int) -> int:
    return n - 1

# Threshold for "near-zero" in the structure–function index
EPSILON_SF = 0.1

# Figure output directory (created under SPECTRAL_DIR)
FIG_DIR = SPECTRAL_DIR / "figs_crossfreq"

# Matplotlib style
plt.style.use("seaborn-v0_8-whitegrid")


# ---------------------------------------------------------------------
# Helpers for loading MATLAB structs
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


def load_results_mat(path: Path) -> Dict[str, np.ndarray]:
    """
    Load 'results' struct from a .mat file and convert to a flat dict
    mapping field name -> numpy array / list.

    Handles both structured arrays (MATLAB structs) and regular dicts.
    """
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    if "results" not in mat:
        raise KeyError(f"'results' not found in {path}")
    res_struct = mat["results"]

    # If it's already a dict-like from the recursive conversion, use that
    if not hasattr(res_struct, "dtype") or res_struct.dtype.names is None:
        return _matlab_struct_to_dict(res_struct)

    # Structured array: convert to dict
    res_dict = {}
    for name in res_struct.dtype.names:
        field_val = res_struct[name]
        if isinstance(field_val, np.ndarray):
            if field_val.ndim == 0:
                res_dict[name] = field_val.item()
            elif field_val.size == 1:
                res_dict[name] = field_val.flat[0]
            else:
                res_dict[name] = field_val
        else:
            res_dict[name] = field_val
    return res_dict


# ---------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------
def parse_subject_stage(filename: str) -> Tuple[str, str]:
    """
    Parse subject ID and stage from filenames like:
      EPCTL01_N1_ssdi_spectral_domain_results.mat
    """
    m = re.match(r"EPCTL(\d+)_(\w+)_ssdi_.*\.mat", filename)
    if not m:
        raise ValueError(f"Cannot parse subject/stage from {filename}")
    subj = f"{int(m.group(1)):02d}"
    stage = m.group(2)
    return subj, stage


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


def extract_subspace(L_data, n: int) -> np.ndarray:
    """
    Extract macroscopic subspace (first n columns) from the L-matrix
    cell corresponding to macro size n.

    L_data can be:
    - A list/array where index n-1 corresponds to macro size n
    - A dict/struct where keys are macro sizes
    - A single 2D matrix
    """
    # Handle different possible structures
    if isinstance(L_data, dict):
        if n in L_data:
            L = np.asarray(L_data[n])
        else:
            raise IndexError(f"Requested macro size n={n} not in L_dict.")
    elif isinstance(L_data, (list, tuple, np.ndarray)):
        idx = lopto_cell_index(n)
        if not (0 <= idx < len(L_data)):
            raise IndexError(f"Requested macro size n={n} not in L_list.")
        L = np.asarray(L_data[idx])
    else:
        L = np.asarray(L_data)

    if L.ndim == 0:
        raise ValueError("L cell contains a scalar, not a matrix.")
    elif L.ndim == 1:
        L = L.reshape(-1, 1)
    elif L.ndim > 2:
        L = L.reshape(L.shape[0], -1)

    if L.shape[1] < n:
        return L
    return L[:, :n]


def subspace_distance(U: np.ndarray, V: np.ndarray) -> float:
    """
    Distance between subspaces spanned by columns of U and V.

    Uses SciPy's principal-angle ('subspace_angles') algorithm, and returns
    RMS(sin(theta)) in [0,1]:

      0  -> identical subspaces
      1  -> orthogonal subspaces
    """
    U = np.asarray(U)
    V = np.asarray(V)
    angles = subspace_angles(U, V)  # radians
    dist = float(np.sqrt(np.mean(np.sin(angles) ** 2)))
    return dist


# ---------------------------------------------------------------------
# Dataclasses for collected records
# ---------------------------------------------------------------------
@dataclass
class AdditivityRecord:
    stage: str
    subject: str
    n: int
    dd_broad: float
    dd_sum_bands: float
    sf_broad: float
    sf_mean_bands: float


@dataclass
class AlignmentRecord:
    stage: str
    subject: str
    n: int
    band_i: str
    band_j: str
    dist: float


@dataclass
class BandDDRecord:
    stage: str
    subject: str
    n: int
    band: str
    dd_band: float


# ---------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------
def gather_records() -> Tuple[List[AdditivityRecord],
                              List[AlignmentRecord],
                              List[BandDDRecord]]:
    add_records: List[AdditivityRecord] = []
    align_records: List[AlignmentRecord] = []
    band_records: List[BandDDRecord] = []

    spectral_files = sorted(SPECTRAL_DIR.glob("EPCTL*_ssdi_spectral_domain_results.mat"))

    for path in spectral_files:
        subj, stage = parse_subject_stage(path.name)
        if stage not in STAGES:
            continue

        print(f"Loading spectral SSDI for subject {subj}, stage {stage} ...")
        res_spec = load_results_mat(path)

        # Optional: time-domain broadband file for this subject/stage
        if TIME_DIR is not None:
            time_file = TIME_DIR / TIME_PATTERN.format(subject=subj, stage=stage)
            if time_file.exists():
                res_time = load_results_mat(time_file)
            else:
                print(f"  [WARN] No time-domain file at {time_file}; "
                      f"broadband comparison skipped for this case.")
                res_time = None
        else:
            res_time = None

        for n in MACRO_SIZES:
            # ---------- Additivity part (broadband vs band sum) ----------
            if res_time is not None:
                if BROADBAND_PREFIX:
                    dopto_broad_key = f"{BROADBAND_PREFIX}_dopto_m{n}"
                    gopto_broad_key = f"{BROADBAND_PREFIX}_gopto_m{n}"
                else:
                    dopto_broad_key = f"dopto_m{n}"
                    gopto_broad_key = f"gopto_m{n}"
                if dopto_broad_key in res_time and gopto_broad_key in res_time:
                    dopto_broad = np.asarray(res_time[dopto_broad_key]).ravel()
                    dd_broad = float(dopto_broad[0])
                    gopto_broad = np.asarray(res_time[gopto_broad_key])
                    sf_broad = structure_function_index(gopto_broad)
                else:
                    dd_broad = np.nan
                    sf_broad = np.nan
            else:
                dd_broad = np.nan
                sf_broad = np.nan

            dd_sum = 0.0
            sf_bands = []

            for band in BANDS:
                dopto_key = f"{band}_dopto_m{n}"
                gopto_key = f"{band}_gopto_m{n}"
                if dopto_key not in res_spec or gopto_key not in res_spec:
                    dd_band = np.nan
                    sf_band = np.nan
                else:
                    dopto_band = np.asarray(res_spec[dopto_key]).ravel()
                    dd_band = float(dopto_band[0])
                    gopto_band = np.asarray(res_spec[gopto_key])
                    sf_band = structure_function_index(gopto_band)

                dd_sum += dd_band
                sf_bands.append(sf_band)

                # store per-band DD for later band-contribution analysis
                band_records.append(
                    BandDDRecord(
                        stage=stage,
                        subject=subj,
                        n=n,
                        band=band,
                        dd_band=dd_band,
                    )
                )

            add_records.append(
                AdditivityRecord(
                    stage=stage,
                    subject=subj,
                    n=n,
                    dd_broad=dd_broad,
                    dd_sum_bands=dd_sum,
                    sf_broad=sf_broad,
                    sf_mean_bands=float(np.nanmean(sf_bands)),
                )
            )

            # ---------- Macro alignment across bands ----------
            subspaces: Dict[str, np.ndarray] = {}
            for band in BANDS:
                lopto_key = f"{band}_Lopto_m"
                if lopto_key not in res_spec:
                    continue
                L_list = res_spec[lopto_key]
                try:
                    U = extract_subspace(L_list, n)
                except Exception as e:
                    print(f"  [WARN] Cannot extract subspace for {band}, n={n}: {e}")
                    continue
                subspaces[band] = U

            bands_present = sorted(subspaces.keys())
            for i in range(len(bands_present)):
                for j in range(i + 1, len(bands_present)):
                    bi = bands_present[i]
                    bj = bands_present[j]
                    U = subspaces[bi]
                    V = subspaces[bj]
                    if U.shape[0] != V.shape[0]:
                        continue
                    dist = subspace_distance(U, V)
                    align_records.append(
                        AlignmentRecord(
                            stage=stage,
                            subject=subj,
                            n=n,
                            band_i=bi,
                            band_j=bj,
                            dist=dist,
                        )
                    )

    return add_records, align_records, band_records


# ---------------------------------------------------------------------
# Scale factor (global)
# ---------------------------------------------------------------------
def compute_scale_factor(records: List[AdditivityRecord]) -> float:
    """Estimate global scale factor c so that DD_broad ≈ c * sum_b DD_b."""
    x = np.array([r.dd_sum_bands for r in records])
    y = np.array([r.dd_broad for r in records])
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0)
    if not np.any(mask):
        raise RuntimeError("No valid records to estimate scale factor.")
    ratios = y[mask] / x[mask]
    c_scale = float(np.median(ratios))
    print(f"[INFO] Estimated global scale factor c ≈ {c_scale:.3e} "
          f"(DD_broad ≈ c * sum_b DD_b)")
    return c_scale


# ---------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------
def plot_additivity(records: List[AdditivityRecord], outdir: Path, c_scale: float) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    recs = [
        r for r in records
        if np.isfinite(r.dd_broad)
        and np.isfinite(r.dd_sum_bands)
        and r.dd_sum_bands > 0
    ]
    if not recs:
        print("[INFO] No valid broadband vs band-sum records; skipping additivity plots.")
        return

    stages = STAGES
    colors = plt.cm.tab10(np.linspace(0, 1, len(stages)))

    # Rescaled x-values
    x_raw = np.array([r.dd_sum_bands for r in recs])
    y_raw = np.array([r.dd_broad for r in recs])
    x_scaled_all = x_raw * c_scale
    y_all = y_raw

    lim_min = float(np.nanmin(np.concatenate([x_scaled_all, y_all])))
    lim_max = float(np.nanmax(np.concatenate([x_scaled_all, y_all])))
    pad = 0.05 * (lim_max - lim_min)
    lim_min -= pad
    lim_max += pad

    # Scatter per macro size
    sizes = sorted(set(r.n for r in recs))
    nrows = 2
    ncols = int(np.ceil(len(sizes) / nrows))

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(14, 6), sharex=True, sharey=True
    )
    axes = np.array(axes).ravel()

    for ax, n in zip(axes, sizes):
        for si, st in enumerate(stages):
            xs = [
                r.dd_sum_bands for r in recs
                if r.n == n and r.stage == st
            ]
            if not xs:
                continue
            xs = np.array(xs) * c_scale
            ys = np.array([
                r.dd_broad for r in recs
                if r.n == n and r.stage == st
            ])
            ax.scatter(xs, ys, label=st, s=20, alpha=0.7, color=colors[si])

        ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=1)
        ax.set_title(f"n = {n}")
        ax.set_xlabel(r"$c \cdot \sum_b DD_b$ (band-limited, rescaled)")
        ax.set_ylabel(r"$DD_{\mathrm{broad}}$")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", frameon=True)
    fig.suptitle("Broadband DD vs scaled sum of band-limited DDs\n"
                 r"(DD$_{\mathrm{broad}} \approx c \sum_b DD_b$)")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(outdir / "DD_broad_vs_scaled_sumBands.png", dpi=300)

    # Difference per stage
    diff = np.array([r.dd_broad - c_scale * r.dd_sum_bands for r in recs])
    stage_labels = np.array([r.stage for r in recs])

    fig, ax = plt.subplots(figsize=(8, 5))
    data_by_stage = [diff[stage_labels == st] for st in stages]
    ax.boxplot(data_by_stage, labels=stages, showfliers=False)
    ax.axhline(0, color="k", linestyle="--", lw=1)
    ax.set_ylabel(r"$DD_{\mathrm{broad}} - c \sum_b DD_b$")
    ax.set_title("Difference between broadband DD and scaled sum of band-limited DDs")
    fig.tight_layout()
    fig.savefig(outdir / "DD_diff_by_stage_scaled.png", dpi=300)


def plot_alignment(records: List[AlignmentRecord], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if not records:
        print("[INFO] No alignment records; skipping alignment plots.")
        return

    sizes = sorted(set(r.n for r in records))
    band_pairs = sorted(
        set((r.band_i, r.band_j) for r in records),
        key=lambda x: (x[0], x[1]),
    )

    for n in sizes:
        rec_n = [r for r in records if r.n == n]
        if not rec_n:
            continue

        fig, ax = plt.subplots(figsize=(12, 5))
        pair_labels = [f"{bi}-{bj}" for (bi, bj) in band_pairs]
        x = np.arange(len(band_pairs))
        width = 0.8 / len(STAGES)
        colors = plt.cm.tab10(np.linspace(0, 1, len(STAGES)))

        for si, st in enumerate(STAGES):
            offsets = x - 0.4 + width / 2 + si * width
            data = []
            for (bi, bj) in band_pairs:
                vals = [
                    r.dist
                    for r in rec_n
                    if r.stage == st and r.band_i == bi and r.band_j == bj
                ]
                data.append(vals if vals else [np.nan])
            clean_data = [np.array(d)[np.isfinite(d)] for d in data]
            bp = ax.boxplot(
                clean_data,
                positions=offsets,
                widths=width * 0.9,
                patch_artist=True,
                showfliers=False,
            )
            for patch in bp["boxes"]:
                patch.set_facecolor(colors[si])
                patch.set_alpha(0.6)

        ax.set_xticks(x)
        ax.set_xticklabels(pair_labels, rotation=45, ha="right")
        ax.set_ylabel("Subspace distance (RMS sin θ)")
        ax.set_title(f"Inter-band macro subspace distances (n={n})")
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        handles = [plt.Line2D([0], [0], color=colors[si], lw=6)
                   for si in range(len(STAGES))]
        ax.legend(handles, STAGES, title="Stage", loc="upper right")
        fig.tight_layout()
        fig.savefig(outdir / f"macro_alignment_n{n}.png", dpi=300)


# ------------------------ NEW PLOTS -----------------------------------
def plot_band_contributions(
    band_records: List[BandDDRecord],
    add_records: List[AdditivityRecord],
    c_scale: float,
    outdir: Path,
    target_n: int = TARGET_N,
) -> None:
    """
    For a chosen macro size n, plot per-band contributions:
        contribution_b = (c * DD_b) / DD_broad

    as barplots (median + IQR) for each stage × band.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # map (stage, subject) -> DD_broad at target_n
    dd_broad_map: Dict[Tuple[str, str], float] = {}
    for r in add_records:
        if r.n == target_n and np.isfinite(r.dd_broad):
            dd_broad_map[(r.stage, r.subject)] = r.dd_broad

    contrib: Dict[Tuple[str, str], List[float]] = {}
    for br in band_records:
        if br.n != target_n:
            continue
        key_ss = (br.stage, br.subject)
        if key_ss not in dd_broad_map:
            continue
        dd_broad = dd_broad_map[key_ss]
        if not np.isfinite(dd_broad) or dd_broad <= 0:
            continue
        if not np.isfinite(br.dd_band):
            continue
        val = (c_scale * br.dd_band) / dd_broad
        contrib.setdefault((br.stage, br.band), []).append(val)

    # organise into matrix [band, stage]
    med = np.full((len(BANDS), len(STAGES)), np.nan)
    q1 = np.full_like(med, np.nan)
    q3 = np.full_like(med, np.nan)

    for bi, band in enumerate(BANDS):
        for si, stage in enumerate(STAGES):
            vals = np.array(contrib.get((stage, band), []), float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            med[bi, si] = np.median(vals)
            q1[bi, si] = np.percentile(vals, 25)
            q3[bi, si] = np.percentile(vals, 75)

    x = np.arange(len(BANDS))
    width = 0.8 / len(STAGES)
    colors = plt.cm.tab10(np.linspace(0, 1, len(STAGES)))

    fig, ax = plt.subplots(figsize=(10, 5))
    for si, stage in enumerate(STAGES):
        offsets = x - 0.4 + width / 2 + si * width
        m = med[:, si]
        lo = m - q1[:, si]
        hi = q3[:, si] - m
        ax.bar(offsets, m, width=width, color=colors[si], alpha=0.8, label=stage)
        ax.errorbar(
            offsets, m, yerr=[lo, hi], fmt="none", ecolor="k", elinewidth=1, capsize=3
        )

    ax.axhline(1.0, color="k", linestyle="--", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(BANDS)
    ax.set_ylabel(r"$(c \, DD_b) / DD_{\mathrm{broad}}$")
    ax.set_title(
        rf"Relative band contributions at $n={target_n}$ "
        r"(median ± IQR across subjects)"
    )
    ax.legend(title="Stage", loc="upper right")
    fig.tight_layout()
    fig.savefig(outdir / f"band_contributions_n{target_n}.png", dpi=300)


def plot_synergy_single_band(
    band_records: List[BandDDRecord],
    add_records: List[AdditivityRecord],
    c_scale: float,
    outdir: Path,
    target_n: int = TARGET_N,
) -> None:
    """
    For a chosen macro size n, compute a "single-band synergy" index:

        S = DD_broad - max_b (c * DD_b)

    and show as stage-wise boxplots.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # map (stage, subject) -> DD_broad at target_n
    dd_broad_map: Dict[Tuple[str, str], float] = {}
    for r in add_records:
        if r.n == target_n and np.isfinite(r.dd_broad):
            dd_broad_map[(r.stage, r.subject)] = r.dd_broad

    # collect c*DD_b for each (stage, subject)
    best_band_scaled: Dict[Tuple[str, str], float] = {}
    for br in band_records:
        if br.n != target_n:
            continue
        key_ss = (br.stage, br.subject)
        if key_ss not in dd_broad_map:
            continue
        if not np.isfinite(br.dd_band):
            continue
        val = c_scale * br.dd_band
        if key_ss not in best_band_scaled:
            best_band_scaled[key_ss] = val
        else:
            best_band_scaled[key_ss] = max(best_band_scaled[key_ss], val)

    # compute synergy S per subject/stage
    synergy_by_stage: Dict[str, List[float]] = {st: [] for st in STAGES}
    for (stage, subj), best_val in best_band_scaled.items():
        dd_broad = dd_broad_map[(stage, subj)]
        if not np.isfinite(dd_broad):
            continue
        S = dd_broad - best_val
        synergy_by_stage[stage].append(S)

    fig, ax = plt.subplots(figsize=(8, 5))
    data = [synergy_by_stage[st] for st in STAGES]
    ax.boxplot(data, labels=STAGES, showfliers=False)
    ax.axhline(0, color="k", linestyle="--", linewidth=1)
    ax.set_ylabel(r"$DD_{\mathrm{broad}} - \max_b (c \, DD_b)$")
    ax.set_title(
        rf"Single-band synergy index at $n={target_n}$ "
        r"(broadband minus best single band)"
    )
    fig.tight_layout()
    fig.savefig(outdir / f"single_band_synergy_n{target_n}.png", dpi=300)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    add_rec, align_rec, band_rec = gather_records()

    c_scale = compute_scale_factor(add_rec)

    plot_additivity(add_rec, FIG_DIR, c_scale)
    plot_alignment(align_rec, FIG_DIR)
    plot_band_contributions(band_rec, add_rec, c_scale, FIG_DIR, target_n=TARGET_N)
    plot_synergy_single_band(band_rec, add_rec, c_scale, FIG_DIR, target_n=TARGET_N)

    print(f"Saved cross-frequency figures to: {FIG_DIR}")


if __name__ == "__main__":
    main()
