"""Parallel SSDI grid search across independent bands and macro dimensions."""

from concurrent.futures import ProcessPoolExecutor, as_completed
import os

import numpy as np

from .model import prepare_var, itransform_subspace
from .optimise import (opt_gd_ddx_mruns, canonical_opt_dds_mruns,
                       cluster_representatives, random_initial_bases)

_SHARED = None


def _initialise(cak, h, freqs, v):
    global _SHARED
    _SHARED = (cak, h, freqs, v)


def _run_block(task):
    m, name, low, high, starts, seed, runs, pre_iters, final_iters, history, ctol = task
    cak, h, freqs, v = _SHARED
    if starts is None:
        starts = random_initial_bases(v.shape[0], m, runs, seed=seed)
    indices = np.flatnonzero((freqs >= low) & (freqs <= high))
    if indices.size == 0:
        raise ValueError(f"band {name!r} has no bins on the frequency grid")
    pre = opt_gd_ddx_mruns(cak, starts, niters=pre_iters, history=history)
    pre_original = itransform_subspace(pre.bases, v)
    representatives, sizes = cluster_representatives(pre_original, ctol)
    final = canonical_opt_dds_mruns(h, freqs, indices,
                                     pre.bases[:, :, representatives],
                                     niters=final_iters, history=history)
    result = dict(pre=pre, pre_original_bases=pre_original,
                  representatives=representatives, cluster_sizes=sizes,
                  final=final,
                  final_original_bases=itransform_subspace(final.bases, v),
                  frequencies=freqs[indices])
    return m, name, result


def run_grid_parallel(a, v, fs, fres, dimensions, bands, *, runs=100, seed=0,
                      initial_bases=None, workers=None, pre_iters=10000,
                      final_iters=10000, cluster_tol=1e-10, history=False,
                      progress=None, pairs=None):
    """Optimise every requested ``(dimension, band)`` block independently.

    Parameters
    ----------
    a, v : VAR coefficients (n,n,p) and innovations covariance (n,n)
    dimensions : iterable of integers in ``1..n-1``
    bands : mapping of name to inclusive ``(low_hz, high_hz)`` bounds
    initial_bases : optional mapping ``(dimension, band_name) -> (n,m,runs)``
        Supply fixed bases here when comparing optimisation runs.
    workers : processes; one runs serially. Defaults to at most eight CPUs.
    progress : optional callback ``(completed, total, dimension, band, result)``
    pairs : optional iterable of ``(dimension, band_name)``
        Run only these blocks. Seeds still follow the complete grid, so a
        resumed subset has the same starts as a full uninterrupted grid.

    Returns a dictionary keyed by ``(dimension, band_name)``. The model
    transform and transfer function are prepared once, then supplied once per
    worker process rather than rebuilt for every block. Random seeds are
    deterministic across worker counts when ``initial_bases`` is omitted.
    """
    a, v = np.asarray(a, float), np.asarray(v, float)
    n = a.shape[0]
    dimensions = tuple(int(m) for m in dimensions)
    if any(not 1 <= m < n for m in dimensions):
        raise ValueError("every dimension must satisfy 1 <= m < n")
    bands = dict(bands)
    if not dimensions or not bands:
        raise ValueError("dimensions and bands must be non-empty")
    if workers is None:
        workers = max(1, min(8, (os.cpu_count() or 1)-1))
    if workers < 1:
        raise ValueError("workers must be positive")
    _, _, _, cak, h, freqs = prepare_var(a, v, fs, fres)
    combinations = [(m, name, bounds) for m in dimensions
                    for name, bounds in bands.items()]
    wanted = None if pairs is None else set(pairs)
    valid = {(m, name) for m, name, _ in combinations}
    if wanted is not None and (not wanted or not wanted <= valid):
        raise ValueError("pairs must be a non-empty subset of the dimension-band grid")
    seeds = np.random.SeedSequence(seed).spawn(len(combinations))
    tasks = []
    for (m, name, bounds), child_seed in zip(combinations, seeds):
        if wanted is not None and (m, name) not in wanted:
            continue
        low, high = bounds
        starts = None if initial_bases is None else np.asarray(initial_bases[(m, name)])
        if starts is not None and (starts.ndim != 3 or starts.shape[:2] != (n, m)):
            raise ValueError(f"initial bases for {(m, name)} must have shape (n,m,runs)")
        tasks.append((m, name, low, high, starts, child_seed, runs,
                      pre_iters, final_iters, history, cluster_tol))
    results = {}
    if workers == 1:
        _initialise(cak, h, freqs, v)
        for completed, task in enumerate(tasks, start=1):
            m, name, result = _run_block(task)
            results[(m, name)] = result
            if progress is not None:
                progress(completed, len(tasks), m, name, result)
    else:
        with ProcessPoolExecutor(max_workers=workers,
                                 initializer=_initialise,
                                 initargs=(cak, h, freqs, v)) as pool:
            futures = [pool.submit(_run_block, task) for task in tasks]
            for completed, future in enumerate(as_completed(futures), start=1):
                m, name, result = future.result()
                results[(m, name)] = result
                if progress is not None:
                    progress(completed, len(tasks), m, name, result)
    return results
