"""Find low-dimensional subspaces with low dynamical dependence."""

from dataclasses import dataclass
import numpy as np

from .model import orthonormalise, prepare_var, itransform_subspace
from .objective import cak2ddx, cak2ddxgrad, canonical_band_objective_and_grad


@dataclass
class Optimisation:
    values: np.ndarray
    bases: np.ndarray
    converged: np.ndarray
    iterations: np.ndarray
    step_sizes: np.ndarray
    histories: list


def _factors(gdls):
    if np.isscalar(gdls):
        return float(gdls), 1/float(gdls)
    return float(gdls[0]), float(gdls[1])


def _tolerances(tol, proxy=False):
    if np.isscalar(tol):
        return float(tol), float(tol), float(tol)/(10 if proxy else 1)
    return tuple(float(x) for x in tol)


def _pre_single(cak, initial, niters, gdes, sig0, gdls, tol, history):
    """Direct SSDI-1 opt_gd{1,2}_ddx update order and convergence tests."""
    if gdes not in (1, 2):
        raise ValueError("gdes must be 1 or 2")
    increase, decrease = _factors(gdls)
    stol, dtol, gtol = _tolerances(tol, proxy=True)
    l = initial.copy()
    gmat, g = cak2ddxgrad(l, cak)
    dd = cak2ddx(l, cak)
    sig = sig0
    hist = [[dd, sig, g]] if history else None
    conv, iters = 0, 1
    for iters in range(2, niters+1):
        trial = orthonormalise(l-sig*gmat/g)
        ddnew = cak2ddx(trial, cak)
        if gdes == 2:
            l = trial
            gmat, g = cak2ddxgrad(l, cak)
            if ddnew < dd:
                dd, sig = ddnew, increase*sig
            else:
                sig = decrease*sig
        elif ddnew < dd:
            l, dd = trial, ddnew
            gmat, g = cak2ddxgrad(l, cak)
            sig = increase*sig
        else:
            sig = decrease*sig
        if history:
            hist.append([dd, sig, g])
        if sig < stol:
            conv = 1
        elif dd < dtol:
            conv = 2
        elif g < gtol:
            conv = 3
        if conv:
            break
    return dd, l, conv, iters, sig, np.asarray(hist) if history else None


def _final_single(h, freqs, idx, initial, niters, sig0, gdls, tol, history):
    """Optimise band-limited DD, accepting only steps that improve it."""
    increase, decrease = _factors(gdls)
    stol, dtol, gtol = _tolerances(tol)
    l = orthonormalise(initial)
    sig = sig0
    dd, gmat, g, ok = canonical_band_objective_and_grad(l, h, freqs, idx)
    if not ok:
        dd, gmat, g = float("inf"), np.zeros_like(l), float("inf")
    hist = [[dd, sig, g, 1]] if history else None
    conv, iters = 0, 1
    for iters in range(2, niters+1):
        accepted = 0
        if np.isfinite(g) and g > 0:
            trial = orthonormalise(l-sig*gmat/g)
            candidate, gc, mc, valid = canonical_band_objective_and_grad(trial, h, freqs, idx)
            if valid and np.isfinite(candidate) and candidate < dd:
                l, dd, gmat, g = trial, candidate, gc, mc
                sig *= increase
                accepted = 1
            else:
                sig *= decrease
        else:
            sig *= decrease
        if history:
            hist.append([dd, sig, g, accepted])
        if sig < stol:
            conv = 1
        elif dd < dtol:
            conv = 2
        elif g < gtol:
            conv = 3
        if conv:
            break
    return dd, l, conv, iters, sig, np.asarray(hist) if history else None


def _multirun(initial, run):
    initial = np.asarray(initial, float)
    if initial.ndim == 2:
        initial = initial[:, :, None]
    results = [run(initial[:, :, k]) for k in range(initial.shape[2])]
    order = np.argsort([r[0] for r in results], kind="stable")
    ordered = [results[i] for i in order]
    return Optimisation(np.array([r[0] for r in ordered]),
                        np.stack([r[1] for r in ordered], axis=2),
                        np.array([r[2] for r in ordered]),
                        np.array([r[3] for r in ordered]),
                        np.array([r[4] for r in ordered]),
                        [r[5] for r in ordered])


def opt_gd_ddx_mruns(cak, initial, niters=10000, gdes=2, sig0=1.0,
                      gdls=2.0, tol=1e-10, history=False):
    return _multirun(initial, lambda l: _pre_single(cak, l, niters, gdes, sig0,
                                                     gdls, tol, history))


def canonical_opt_dds_mruns(h, freqs, band_idx, initial, niters=10000,
                              gdes=2, sig0=0.1, gdls=2.0, tol=1e-10,
                              history=False):
    if gdes not in (1, 2):
        raise ValueError("gdes must be 1 or 2")
    return _multirun(initial, lambda l: _final_single(h, freqs, band_idx, l,
                                                       niters, sig0, gdls, tol, history))


def _distance(a, b):
    # Knyazev-Argentati small-angle branch used by SSDI-1 subspacea.m.
    af = a.copy()
    bf = b.copy()
    threshold = np.finfo(float).eps**0.981
    for x in (af, bf):
        for col in range(x.shape[1]):
            scale = np.max(np.abs(x[:, col]))
            if scale > threshold:
                x[:, col] /= scale
    qf, qg = orthonormalise(af), orthonormalise(bf)
    u, singular, vh = np.linalg.svd(qf.T @ qg, full_matrices=False)
    singular = np.minimum(singular, 1)
    angles = np.maximum(np.arccos(singular), 0)
    small = singular > np.sqrt(2)/2
    if np.any(small):
        rg = (qg @ vh.T)[:, small]
        residual = rg-qf @ (qf.T @ rg)
        x = np.linalg.svd(residual, compute_uv=False)
        angles[small] = np.maximum(np.arcsin(np.minimum(x, 1)), 0)[::-1]
    return float(np.max(angles))/(np.pi/2)


def cluster_representatives(bases, tolerance=1e-10):
    """SSDI-1 Lcluster: greedy representatives in ascending objective order."""
    count = bases.shape[2]
    available = np.ones(count, bool)
    representatives, sizes = [], []
    for i in range(count):
        if not available[i]:
            continue
        representatives.append(i)
        available[i] = False
        size = 1
        for j in range(count):
            if available[j] and _distance(bases[:, :, i], bases[:, :, j]) < tolerance:
                available[j] = False
                size += 1
        sizes.append(size)
    return np.asarray(representatives, int), np.asarray(sizes, int)


def run_two_stage(a, v, fs, fres, m, bands, *, initial_bases,
                  pre_iters=10000, final_iters=10000, cluster_tol=1e-10,
                  history=False):
    """Optimise a fast proxy, keep distinct solutions, then refine band DD.

    ``initial_bases`` is an (n,m,runs) array or a mapping from band name to
    such arrays. Use the same starting bases to compare optimisation runs.
    ``bands`` maps names to inclusive (low_hz, high_hz) tuples.
    """
    _, _, _, cak, h, freqs = prepare_var(a, v, fs, fres)
    if not 1 <= m < a.shape[0]:
        raise ValueError("m must be between 1 and n-1")
    out = {}
    for name, (lo, hi) in bands.items():
        idx = np.flatnonzero((freqs >= lo) & (freqs <= hi))
        if idx.size == 0:
            continue
        starts = initial_bases[name] if isinstance(initial_bases, dict) else initial_bases
        if starts.shape[0] != a.shape[0] or starts.shape[1] != m:
            raise ValueError("initial basis dimensions must be (n,m,runs)")
        pre = opt_gd_ddx_mruns(cak, starts, niters=pre_iters, history=history)
        # Cluster after mapping bases back to the original signal coordinates.
        untransformed_pre = itransform_subspace(pre.bases, v)
        reps, sizes = cluster_representatives(untransformed_pre, cluster_tol)
        final = canonical_opt_dds_mruns(h, freqs, idx, pre.bases[:, :, reps],
                                         niters=final_iters, history=history)
        out[name] = {"pre": pre, "pre_original_bases": untransformed_pre,
                     "representatives": reps, "cluster_sizes": sizes,
                     "final": final, "final_original_bases": itransform_subspace(final.bases, v),
                     "frequencies": freqs[idx]}
    return out


def random_initial_bases(n, m, runs, seed=None):
    """SVD-orthonormalised Gaussian restarts for new Python analyses."""
    rng = np.random.default_rng(seed)
    bases = np.empty((n, m, runs))
    for i in range(runs):
        bases[:, :, i] = orthonormalise(rng.standard_normal((n, m)))
    return bases
