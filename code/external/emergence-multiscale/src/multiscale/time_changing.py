"""Follow dynamical dependence as a recording moves through named windows.

Each window gets its own VAR fit and its own optimisation. A result therefore
describes the dynamics *within* that part of the recording. Window boundaries
are sample indices, with the end index excluded.
"""

from dataclasses import dataclass

import numpy as np

from .fit import fit_var_auto
from .optimise import random_initial_bases
from .parallel import run_grid_parallel
from .pipeline import run_full_dd


@dataclass
class WindowResult:
    """The fitted model and DD searches for one named part of a recording."""

    start: int
    stop: int
    model: object
    blocks: dict


def run_time_changing_dd(signals, windows, *, domain, fs, fres, dimensions,
                         bands=None, pmax=20, criterion="bic", regmode="LWR",
                         runs=100, seed=0, workers=1, pre_iters=10000,
                         final_iters=10000, state_space_iters=100,
                         cluster_tol=1e-10, progress=None):
    """Fit and optimise DD independently in each recording window.

    ``signals`` has shape ``(nodes, samples[, trials])``. ``windows`` maps
    names to ``(start, stop)`` sample indices. Windows may overlap; each uses
    the same sample range in every trial. Keep enough samples in every window
    for model fitting, and compare windows only when their preprocessing and
    sampling rate are consistent.

    Choose ``domain='temporal'`` for the direct state-space DD optimiser, or
    ``domain='spectral'`` for frequency-band canonical DD. Spectral mode needs
    ``bands={name: (low_hz, high_hz)}``. Both modes first optimise the fast DDx
    proxy and then optimise actual DD from representative proxy solutions.
    The model order is chosen separately in each window using AIC, BIC, HQC,
    or LRT. The regression method is LWR or OLS. ``progress`` receives
    ``(window_name, completed_windows, total_windows)`` after each window.

    Returns a dictionary of ``WindowResult`` objects. In spectral mode,
    ``blocks[(dimension, band)]`` is the same result as ``run_grid_parallel``.
    In temporal mode, ``blocks[dimension]`` is the result of ``run_full_dd``;
    its ``state_space`` entry contains the optimised DD values and bases.
    """
    x = np.asarray(signals, float)
    if x.ndim == 2:
        x = x[:, :, None]
    if x.ndim != 3 or not np.isfinite(x).all():
        raise ValueError("signals must be a finite (nodes, samples[, trials]) array")
    n, samples, _ = x.shape
    dimensions = tuple(int(m) for m in dimensions)
    if not dimensions or any(m < 1 or m >= n for m in dimensions):
        raise ValueError("dimensions must contain integers from 1 to nodes-1")
    if domain not in ("temporal", "spectral"):
        raise ValueError("domain must be 'temporal' or 'spectral'")
    if domain == "spectral" and not bands:
        raise ValueError("spectral DD needs at least one frequency band")
    if not isinstance(windows, dict) or not windows:
        raise ValueError("windows must be a non-empty mapping of names to sample ranges")
    if not isinstance(pmax, int) or pmax < 1:
        raise ValueError("pmax must be a positive integer")
    if runs < 1 or workers < 1:
        raise ValueError("runs and workers must be positive")
    for name, bounds in windows.items():
        if (not isinstance(bounds, (tuple, list)) or len(bounds) != 2 or
                not all(isinstance(v, (int, np.integer)) for v in bounds)):
            raise ValueError(f"window {name!r} must be (start, stop) sample indices")
        start, stop = bounds
        if start < 0 or stop > samples or stop <= start:
            raise ValueError(f"window {name!r} is outside the recording")
        if stop-start <= pmax+1:
            raise ValueError(f"window {name!r} is too short for pmax={pmax}")
    seeds = np.random.SeedSequence(seed).spawn(len(windows))
    output = {}
    for index, ((name, (start, stop)), window_seed) in enumerate(zip(windows.items(), seeds), 1):
        model = fit_var_auto(x[:, start:stop, :], pmax, criterion=criterion,
                             regmode=regmode)
        if model.order == 0:
            raise ValueError(f"window {name!r} selected VAR order zero; DD needs dynamics")
        a, v = model.fit.coefficients, model.fit.covariance
        if domain == "spectral":
            blocks = run_grid_parallel(
                a, v, fs, fres, dimensions, bands, runs=runs,
                seed=int(window_seed.generate_state(1)[0]), workers=workers, pre_iters=pre_iters,
                final_iters=final_iters, cluster_tol=cluster_tol)
        else:
            blocks = {}
            for m, dimension_seed in zip(dimensions, window_seed.spawn(len(dimensions))):
                starts = random_initial_bases(n, m, runs, seed=dimension_seed)
                blocks[m] = run_full_dd(
                    a, v, fs, fres, starts, final_method="state_space",
                    pre_iters=pre_iters, state_space_iters=state_space_iters,
                    cluster_tol=cluster_tol,
                    seed=int(dimension_seed.generate_state(1)[0]),
                    state_space_workers=workers)
        output[name] = WindowResult(start, stop, model, blocks)
        if progress is not None:
            progress(name, index, len(windows))
    return output


def minimum_dd_by_window(results, *, domain):
    """Return simple records for plotting the best DD in every window.

    Each record has ``window``, ``start``, ``stop``, ``dimension``, ``band``
    (``None`` for temporal DD), ``dd_min``, ``n_solutions``, and ``model_order``.
    Convert the returned list to a pandas DataFrame if desired.
    """
    if domain not in ("temporal", "spectral"):
        raise ValueError("domain must be 'temporal' or 'spectral'")
    rows = []
    for name, window in results.items():
        for key, block in window.blocks.items():
            if domain == "temporal":
                m, band, values = key, None, block["state_space"].values
            else:
                m, band = key
                values = block["final"].values
            rows.append(dict(window=name, start=window.start, stop=window.stop,
                             dimension=m, band=band,
                             dd_min=float(np.min(values)),
                             n_solutions=len(values), model_order=window.model.order))
    return rows
