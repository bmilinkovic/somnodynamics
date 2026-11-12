"""MVGC2's extended-pencil QZ discrete Riccati solver.

This follows ``mdare.m``: build H/J, deflate by QR, order the *outside*
unit-circle generalised eigenvalues, then recover the invariant subspace.
It deliberately does not call a generic discrete-ARE wrapper.
"""

import numpy as np
from scipy.linalg import (cholesky, ordqz, qr, solve_discrete_lyapunov,
                          solve_triangular)


def mdare(a, c, q, r=None, s=None):
    a, c, q = np.asarray(a, float), np.asarray(c, float), np.asarray(q, float)
    dim = a.shape[0]
    obs = c.shape[0]
    r = np.eye(obs) if r is None else np.asarray(r, float)
    s = np.zeros((dim, obs)) if s is None else np.asarray(s, float)
    if a.shape != (dim, dim) or c.shape != (obs, dim) or q.shape != a.shape:
        raise ValueError("invalid DARE dimensions")
    h = np.block([[a.T, np.zeros((dim, dim)), c.T],
                  [-q, np.eye(dim), -s],
                  [s.T, np.zeros((obs, dim)), r]])
    j = np.block([[np.eye(dim), np.zeros((dim, dim+obs))],
                  [np.zeros((dim, dim)), a, np.zeros((dim, obs))],
                  [np.zeros((obs, dim)), -c, np.zeros((obs, obs))]])
    # MATLAB qr(...), then q(:,n+1:end)' times the pencil.
    basis, _ = qr(h[:, 2*dim:], mode="full")
    projected_h = basis[:, obs:].T @ h[:, :2*dim]
    projected_j = basis[:, obs:].T @ j[:, :2*dim]
    jj, hh, alpha, beta, _, z = ordqz(projected_j, projected_h,
                                        sort="ouc", output="real")
    with np.errstate(divide="ignore", invalid="ignore"):
        eigen = alpha/beta
    outside = np.abs(eigen) > 1
    if not np.all(outside[:dim]) or np.any(outside[dim:]):
        return None, None, -1.0, None
    p1, p2 = z[:dim, :dim], z[dim:, :dim]
    if np.linalg.cond(p1) > 1/np.finfo(float).eps:
        return None, None, -2.0, None
    p = np.linalg.solve(p1.T, p2.T).T
    p = (p+p.T)/2
    u = a @ p @ c.T+s
    v = c @ p @ c.T+r
    k = np.linalg.solve(v.T, u.T).T
    apa = a @ p @ a.T-p
    uk = u @ k.T
    rep = np.linalg.norm(apa-uk+q, 1)/(1+np.linalg.norm(apa, 1)+
                                       np.linalg.norm(uk, 1)+np.linalg.norm(q, 1))
    return k, v, float(rep), p


def covariance_balance_iss(a, c, k):
    """Similarity transform to unit stationary state covariance.

    This changes only the state coordinates. The input-output process and
    innovations covariance are identical, while the companion realisation's
    very large state covariance no longer degrades the QZ DARE residual.
    """
    a, c, k = np.asarray(a, float), np.asarray(c, float), np.asarray(k, float)
    sigma = solve_discrete_lyapunov(a, k @ k.T)
    sigma = (sigma + sigma.T) / 2
    root = cholesky(sigma, lower=True)
    balanced_a = solve_triangular(root, a @ root, lower=True)
    balanced_c = c @ root
    balanced_k = solve_triangular(root, k, lower=True)
    return balanced_a, balanced_c, balanced_k


def iss2dd(l, a, c, k, *, return_diagnostics=False):
    """Temporal DI from the ordered-QZ DARE, with state-coordinate fallback.

    The original SSDI/MVGC2 call is tried first. If its residual fails the
    SSDI threshold, a stationary-covariance Cholesky similarity transform is
    applied and the same ``mdare`` routine is called again. Set
    ``return_diagnostics=True`` to also get the method and both residuals.
    """
    l = np.asarray(l, float)
    q = k @ k.T
    try:
        _, v, rep, _ = mdare(a, l.T @ c, q, np.eye(l.shape[1]), k @ l)
    except (np.linalg.LinAlgError, ValueError):
        rep, v = -3.0, None
    original_rep = rep
    method = "original"
    if rep < 0 or rep > 1e-8:
        try:
            ba, bc, bk = covariance_balance_iss(a, c, k)
            _, v, rep, _ = mdare(ba, l.T @ bc, bk @ bk.T,
                                 np.eye(l.shape[1]), bk @ l)
            method = "covariance_balanced"
        except (np.linalg.LinAlgError, ValueError):
            rep, v = -3.0, None
            method = "failed"
    value = float("nan")
    if v is not None and 0 <= rep <= 1e-8:
        sign, logdet = np.linalg.slogdet(v)
        if sign > 0:
            value = float(logdet)
    if return_diagnostics:
        return value, {"method": method, "original_residual": float(original_rep),
                       "residual": float(rep)}
    return value
