"""Original SSDI-1 frequency-integrated temporal DI optimisation.

This is the ``trfun2dd`` / ``trfun2ddgrad`` route used by older SSDI analyses.
This full-spectrum objective differs from the band-limited objective in
``objective.py``.
"""

import numpy as np
from scipy.linalg import cholesky

from .model import orthonormalise
from .optimise import Optimisation, _factors, _multirun, _tolerances


def trfun2dd(l, h):
    """SSDI-1 half-logdet samples and doubled trapezoidal mean."""
    vals = np.empty(h.shape[2])
    for i in range(h.shape[2]):
        q = h[:, :, i].conj().T @ l
        upper = cholesky(q.conj().T @ q, lower=False, check_finite=False)
        vals[i] = float(np.real(np.sum(np.log(np.diag(upper)))))
    return float(np.sum(vals[:-1]+vals[1:])/(len(vals)-1)), vals


def trfun2ddgrad(l, h):
    """SSDI-1 integrated Grassmann gradient with the explicit -2L term."""
    terms = np.empty((h.shape[2], *l.shape))
    for i in range(h.shape[2]):
        hk = h[:, :, i]
        q = hk.conj().T @ l
        terms[i] = np.real(np.linalg.solve((q.conj().T @ q).T,
                                           (hk @ q).T).T)
    grad = np.sum(terms[:-1]+terms[1:], axis=0)/(len(terms)-1)-2*l
    return grad, float(np.sqrt(np.sum(grad*grad)))


def _single(h, initial, niters, gdes, sig0, gdls, tol, history):
    if gdes not in (1, 2):
        raise ValueError("gdes must be 1 or 2")
    increase, decrease = _factors(gdls)
    stol, dtol, gtol = _tolerances(tol)
    l = initial.copy()
    grad, g = trfun2ddgrad(l, h)
    dd, _ = trfun2dd(l, h)
    sig = sig0
    hist = [[dd, sig, g]] if history else None
    conv, iters = 0, 1
    for iters in range(2, niters+1):
        trial = orthonormalise(l-sig*grad/g)
        ddnew, _ = trfun2dd(trial, h)
        if gdes == 2:
            l = trial
            grad, g = trfun2ddgrad(l, h)
            if ddnew < dd:
                dd, sig = ddnew, increase*sig
            else:
                sig = decrease*sig
        elif ddnew < dd:
            l, dd = trial, ddnew
            grad, g = trfun2ddgrad(l, h)
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


def opt_gd_dds_mruns(h, initial, niters=10000, gdes=2, sig0=1e-3,
                     gdls=2.0, tol=1e-9, history=False) -> Optimisation:
    """SSDI-1 multirun integrated spectral optimiser, variants 1 or 2."""
    return _multirun(initial, lambda l: _single(h, l, niters, gdes, sig0,
                                                gdls, tol, history))
