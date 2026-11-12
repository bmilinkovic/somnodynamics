"""Direct state-space temporal DD optimisation from SSDI-1's ES path."""

from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from .model import orthonormalise
from .optimise import Optimisation
from .riccati import covariance_balance_iss, iss2dd


_SHARED_ES = None


def _es_worker_init(a, c, k, maxiters, sig, factors, tol, history):
    global _SHARED_ES
    _SHARED_ES = (a, c, k, maxiters, sig, factors, tol, history)


def _es_worker(task):
    index, initial, seed, mutations = task
    a, c, k, maxiters, sig, factors, tol, history = _SHARED_ES
    result = _single(a, c, k, initial, maxiters, sig, factors, tol, history,
                     np.random.default_rng(seed), mutations)
    return index, result


@dataclass
class TemporalESResult:
    value: float
    basis: np.ndarray
    converged: int
    step_size: float
    iterations: int
    history: np.ndarray | None


def es_parms(rule, dimension):
    """SSDI-1 ``es_parms`` success/failure multipliers for a 1+1 strategy."""
    if dimension < 1:
        raise ValueError("dimension must be positive")
    if rule < np.finfo(float).eps:
        return 1.5, 1.5**(-0.25)
    gain = 1/np.sqrt(dimension+1)
    return float(np.exp((1-rule)*gain)), float(np.exp(-rule*gain))


def _single(a, c, k, initial, maxiters, sig, factors, tol, history, rng,
            mutations):
    l = np.asarray(initial, float).copy()
    n, m = l.shape
    if not np.allclose(l.T @ l, np.eye(m), atol=1e-8):
        raise ValueError("initial basis must be orthonormal")
    ifac, nfac = factors
    stol, dtol = (tol, tol) if np.isscalar(tol) else tol
    dd = iss2dd(l, a, c, k)
    hist = [[dd, sig]] if history else None
    converged, iters = 0, 1
    for iters in range(2, maxiters+1):
        noise = (rng.standard_normal((n, m)) if mutations is None
                 else mutations[:, :, iters-2])
        trial = orthonormalise(l+sig*noise)
        candidate = iss2dd(trial, a, c, k)
        if candidate < dd:
            l, dd, sig = trial, candidate, ifac*sig
        else:
            sig = nfac*sig
        if history:
            hist.append([dd, sig])
        if sig < stol:
            converged = 1
        elif dd < dtol:
            converged = 2
        if converged:
            break
    return TemporalESResult(float(dd), l, converged, float(sig), iters,
                            np.asarray(hist) if history else None)


def opt_es_dd(a, c, k, initial, *, maxiters=100, sig=1e-4, rule=0.2,
              esrule=None, tol=1e-8, seed=None, mutations=None,
              history=False, prebalance=True):
    """One direct DARE-DD run using SSDI-1's 1+1 evolution strategy.

    A mutation is accepted only when its temporal DD is lower. The step size
    uses SSDI-1's ``es_parms`` success/failure factors. ``mutations`` can
    supply an ``(n,m,maxiters-1)`` stream for exact trajectory comparisons.
    ``prebalance`` changes state coordinates before DARE evaluation while
    leaving the underlying process and projection unchanged.
    """
    initial = np.asarray(initial, float)
    if initial.ndim != 2:
        raise ValueError("initial basis must have shape (n,m)")
    if maxiters < 1 or sig <= 0:
        raise ValueError("maxiters and sig must be positive")
    if mutations is not None and np.shape(mutations) != (*initial.shape, maxiters-1):
        raise ValueError("mutations must have shape (n,m,maxiters-1)")
    if prebalance:
        a, c, k = covariance_balance_iss(a, c, k)
    factors = es_parms(rule, initial.shape[1]*(initial.shape[0]-initial.shape[1])) \
        if esrule is None else tuple(esrule)
    return _single(a, c, k, initial, maxiters, sig, factors, tol, history,
                   np.random.default_rng(seed), mutations)


def opt_es_dd_mruns(a, c, k, initial, *, maxiters=100, sig=1e-4, rule=0.2,
                    esrule=None, tol=1e-8, seed=None, mutations=None,
                    history=False, prebalance=True, workers=1):
    """Run and sort direct temporal-DD ES optimisations over all starts.

    ``initial`` has shape ``(n,m,runs)``. An optional exact mutation stream
    has shape ``(n,m,maxiters-1,runs)``. This method is usually much slower
    than the transfer-function gradient optimiser.
    """
    initial = np.asarray(initial, float)
    if initial.ndim != 3:
        raise ValueError("initial bases must have shape (n,m,runs)")
    n, m, runs = initial.shape
    if mutations is not None and np.shape(mutations) != (n, m, maxiters-1, runs):
        raise ValueError("mutations must have shape (n,m,maxiters-1,runs)")
    if maxiters < 1 or sig <= 0 or workers < 1:
        raise ValueError("maxiters, sig, and workers must be positive")
    if prebalance:
        a, c, k = covariance_balance_iss(a, c, k)
    factors = es_parms(rule, m*(n-m)) if esrule is None else tuple(esrule)
    seeds = np.random.SeedSequence(seed).spawn(runs)
    tasks = [(run, initial[:, :, run], seeds[run],
              None if mutations is None else mutations[:, :, :, run])
             for run in range(runs)]
    _es_worker_init(a, c, k, maxiters, sig, factors, tol, history)
    if workers == 1:
        indexed = [_es_worker(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers,
                                 initializer=_es_worker_init,
                                 initargs=(a, c, k, maxiters, sig, factors,
                                           tol, history)) as pool:
            indexed = list(pool.map(_es_worker, tasks))
    results = [result for _, result in sorted(indexed)]
    order = np.argsort([r.value for r in results], kind="stable")
    ordered = [results[i] for i in order]
    return Optimisation(np.array([r.value for r in ordered]),
                        np.stack([r.basis for r in ordered], axis=2),
                        np.array([r.converged for r in ordered]),
                        np.array([r.iterations for r in ordered]),
                        np.array([r.step_size for r in ordered]),
                        [r.history for r in ordered])
