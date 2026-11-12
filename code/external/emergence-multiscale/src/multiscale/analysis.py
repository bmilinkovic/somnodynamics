"""Turn optimised DD values and subspaces into tidy analysis tables."""

import numpy as np
import pandas as pd

from .optimise import _distance
from .statistics import sf_index


def similarity_matrix(bases):
    """Distance between every pair of optimised subspaces.

    ``bases`` has shape ``(nodes, dimension, solutions)``. Distances use the
    largest principal angle, scaled so zero means the same subspace and one
    means orthogonal subspaces.
    """
    bases = np.asarray(bases, float)
    if bases.ndim != 3:
        raise ValueError("bases must have shape (nodes, dimension, solutions)")
    count = bases.shape[2]
    out = np.zeros((count, count))
    for i in range(count):
        for j in range(i+1, count):
            out[i, j] = out[j, i] = _distance(bases[:, :, i], bases[:, :, j])
    return out


def summarise_grid_results(results, *, subject, condition, threshold=1e-4):
    """Summarise a spectral grid search for group analysis and plotting.

    ``results`` is the mapping returned by ``run_grid_parallel``. Returns
    ``(dd_table, sf_table, matrices)``. The tables have one row per band and
    macro dimension; ``matrices[(condition, band, dimension, subject)]`` keeps
    each run-distance matrix for optional visualisation.
    """
    dd_rows, sf_rows, matrices = [], [], {}
    for (m, band), block in results.items():
        values = np.asarray(block["final"].values, float)
        g = similarity_matrix(block["final_original_bases"])
        pair_count = len(values)*(len(values)-1)//2
        shared = dict(subject=subject, condition=condition, band=band, m=m)
        dd_rows.append(shared | dict(dd_min=float(np.min(values)),
                                     dd_median=float(np.median(values)),
                                     dd_mean=float(np.mean(values)),
                                     n_runs=len(values)))
        sf_rows.append(shared | dict(sf_index=sf_index(g, threshold),
                                     n_runs=len(values), n_pairs=pair_count,
                                     near_zero_threshold=threshold))
        matrices[(condition, band, m, subject)] = g
    return pd.DataFrame(dd_rows), pd.DataFrame(sf_rows), matrices


def summarise_window_results(results, *, domain):
    """Build a table of best DD values across named recording windows.

    This accepts ``run_time_changing_dd`` results. It reports the number of
    final solutions and selected VAR order so changes across windows can be
    interpreted alongside model-fit changes.
    """
    from .time_changing import minimum_dd_by_window

    return pd.DataFrame(minimum_dd_by_window(results, domain=domain))


# Keep earlier public names available for existing analysis scripts.
summarize_grid_results = summarise_grid_results
summarize_window_results = summarise_window_results
