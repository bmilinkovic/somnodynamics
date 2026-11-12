"""Fast proxy and band-limited spectral DD objectives."""

import numpy as np
from scipy.linalg import cholesky, solve_triangular


def cak2ddx(l, cak):
    d = 0.0
    for k in range(cak.shape[2]):
        a = l.T @ cak[:, :, k]
        b = a @ l
        d += np.sum(a*a) - np.sum(b*b)
    return float(d)


def cak2ddxgrad(l, cak):
    p = l @ l.T
    g = np.zeros((l.shape[0], l.shape[0]))
    for k in range(cak.shape[2]):
        q = cak[:, :, k]
        g += q @ q.T - q.T @ p @ q - q @ p @ q.T
    grad = 2*g @ l
    grad -= p @ grad
    return grad, float(np.sqrt(np.sum(grad*grad)))


def _weights(freqs):
    n = len(freqs)
    if n == 1:
        return np.ones(1)
    span = freqs[-1]-freqs[0]
    if span <= 0:
        return np.ones(n)/n
    df = np.diff(freqs)
    w = np.empty(n)
    w[0], w[-1] = df[0]/2, df[-1]/2
    w[1:-1] = (df[:-1]+df[1:])/2
    return w/span


def canonical_band_objective_and_grad(l, h, freqs, band_idx):
    """Band-limited DD objective and its analytic Grassmann gradient.

    Follows canonical_band_objective_and_grad.m: Cholesky of C, congruence
    transform, Hermitian eigenvalues, Cholesky inverse of A, and B solve.
    Invalid terms return (inf, zeros, inf, False).
    """
    l = np.asarray(l, float)
    idx = np.asarray(band_idx, int)
    bad = (float("inf"), np.zeros_like(l), float("inf"), False)
    if idx.size == 0 or np.any(idx < 0) or np.any(idx >= h.shape[2]):
        return bad
    vals = np.empty(len(idx))
    grads = np.empty((len(idx), *l.shape))
    for t, i in enumerate(idx):
        hk = h[:, :, i]
        s = hk @ hk.conj().T
        a = l.T @ s @ l
        a = (a+a.conj().T)/2
        b = l.T @ hk @ l
        c = b @ b.conj().T
        c = (c+c.conj().T)/2
        try:
            lc = cholesky(c, lower=True, check_finite=False)
            temp = solve_triangular(lc, a, lower=True, check_finite=False)
            temp = solve_triangular(lc, temp.conj().T, lower=True, check_finite=False).conj().T
            temp = (temp+temp.conj().T)/2
            e = np.linalg.eigvalsh(temp).real
            if not np.all(np.isfinite(e)) or np.any(e <= 0):
                return bad
            vals[t] = np.sum(np.log(e))
            la = cholesky(a, lower=True, check_finite=False)
            eye = np.eye(a.shape[0], dtype=complex)
            ainv = solve_triangular(la.conj().T, solve_triangular(la, eye, lower=True,
                                      check_finite=False), lower=False, check_finite=False)
            if 1/np.linalg.cond(b, p=1) < 1e-14:
                return bad
            binv = np.linalg.solve(b, eye)
            ga = 2*np.real(s @ l @ ainv)
            gb = 2*np.real(hk @ l @ binv + hk.conj().T @ l @ binv.conj().T)
            grads[t] = ga-gb
        except (np.linalg.LinAlgError, ValueError):
            return bad
    w = _weights(np.asarray(freqs)[idx])
    dd = float(w @ vals)
    ge = np.einsum("k,kij->ij", w, grads)
    g = np.real(ge-l @ (l.T @ ge))
    mg = float(np.linalg.norm(g, "fro"))
    if not np.isfinite(dd) or not np.isfinite(mg):
        return bad
    return dd, g, mg, True


def canonical_band_value(l, h, freqs, band_idx):
    """The identical canonical Cholesky/eigenvalue objective without gradient work."""
    idx = np.asarray(band_idx, int)
    if idx.size == 0 or np.any(idx < 0) or np.any(idx >= h.shape[2]):
        return float("inf")
    vals = np.empty(idx.size)
    for t, i in enumerate(idx):
        hk = h[:, :, i]
        a = l.T @ (hk @ hk.conj().T) @ l
        a = (a+a.conj().T)/2
        b = l.T @ hk @ l
        c = b @ b.conj().T
        c = (c+c.conj().T)/2
        try:
            lc = cholesky(c, lower=True, check_finite=False)
            temp = solve_triangular(lc, a, lower=True, check_finite=False)
            temp = solve_triangular(lc, temp.conj().T, lower=True,
                                    check_finite=False).conj().T
            temp = (temp+temp.conj().T)/2
            e = np.linalg.eigvalsh(temp).real
            if not np.all(np.isfinite(e)) or np.any(e <= 0):
                return float("inf")
            vals[t] = np.sum(np.log(e))
        except (np.linalg.LinAlgError, ValueError):
            return float("inf")
    return float(_weights(np.asarray(freqs)[idx]) @ vals)
