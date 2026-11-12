"""Proxy preoptimization followed by either actual DD calculation route."""

import numpy as np

from .integrated import opt_gd_dds_mruns
from .model import itransform_subspace, prepare_var
from .optimise import opt_gd_ddx_mruns, cluster_representatives
from .temporal import opt_es_dd_mruns


def run_full_dd(a, v, fs, fres, initial_bases, *,
                final_method="spectral", pre_iters=10000,
                spectral_iters=1000, state_space_iters=100,
                cluster_tol=1e-10, spectral_sig=0.1,
                state_space_sig=1e-4, es_rule=0.2, seed=None,
                history=False, state_space_workers=1):
    """Run SSDI-1 proxy DDx and select either or both actual DD optimisers.

    ``spectral`` uses ``trfun2dd``/``trfun2ddgrad`` and gradient descent on
    the full transfer-function spectrum. ``state_space`` uses ``iss2dd`` and
    the original 1+1 evolution strategy on the state-space model. ``both``
    starts each from the same proxy representatives; ``spectral_then_state_space``
    follows the original simulation's sequential refinement path.

    This is full-spectrum DD. For frequency-band canonical DD use
    ``run_two_stage`` or ``run_grid_parallel`` instead.
    """
    methods = {"spectral", "state_space", "both", "spectral_then_state_space"}
    if final_method not in methods:
        raise ValueError(f"final_method must be one of {sorted(methods)}")
    a, v, initial_bases = (np.asarray(a, float), np.asarray(v, float),
                           np.asarray(initial_bases, float))
    if initial_bases.ndim != 3 or initial_bases.shape[0] != a.shape[0]:
        raise ValueError("initial bases must have shape (regions,m,runs)")
    m = initial_bases.shape[1]
    if not 1 <= m < a.shape[0]:
        raise ValueError("m must satisfy 1 <= m < regions")
    ass, c, k, cak, h, freqs = prepare_var(a, v, fs, fres)
    pre = opt_gd_ddx_mruns(cak, initial_bases, niters=pre_iters,
                            history=history)
    original_pre = itransform_subspace(pre.bases, v)
    reps, sizes = cluster_representatives(original_pre, cluster_tol)
    starts = pre.bases[:, :, reps]
    out = dict(pre=pre, pre_original_bases=original_pre,
               representatives=reps, cluster_sizes=sizes,
               frequencies=freqs)
    if final_method in ("spectral", "both", "spectral_then_state_space"):
        spectral = opt_gd_dds_mruns(h, starts, niters=spectral_iters,
                                    sig0=spectral_sig, tol=1e-10,
                                    history=history)
        out["spectral"] = spectral
        out["spectral_original_bases"] = itransform_subspace(spectral.bases, v)
    if final_method in ("state_space", "both", "spectral_then_state_space"):
        state_starts = (out["spectral"].bases if final_method ==
                        "spectral_then_state_space" else starts)
        state = opt_es_dd_mruns(ass, c, k, state_starts,
                                maxiters=state_space_iters,
                                sig=state_space_sig, rule=es_rule,
                                seed=seed, history=history,
                                workers=state_space_workers)
        out["state_space"] = state
        out["state_space_original_bases"] = itransform_subspace(state.bases, v)
    return out
