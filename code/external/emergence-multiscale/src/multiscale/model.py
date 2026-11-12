"""Convert a fitted VAR model to the innovations form used for DD.

Array convention is (output, input, lag), as in MATLAB. All transforms use
the original Cholesky and SVD constructions, not covariance whitening by an
eigendecomposition or QR basis extraction.
"""

import numpy as np
from scipy.linalg import cholesky


def orthonormalise(x):
    return np.linalg.svd(np.asarray(x, dtype=float), full_matrices=False)[0]


def transform_var(a, v):
    a, v = np.asarray(a, float), np.asarray(v, float)
    root = cholesky(v, lower=True)
    inverse = np.linalg.inv(root)
    transformed = np.empty_like(a)
    for k in range(a.shape[2]):
        transformed[:, :, k] = inverse @ a[:, :, k] @ root
    return transformed


def itransform_subspace(bases, v):
    """SSDI-1 pullback through inverse transpose Cholesky, then SVD."""
    inverse = np.linalg.inv(cholesky(np.asarray(v, float), lower=True).T)
    bases = np.asarray(bases, float)
    if bases.ndim == 2:
        return orthonormalise(inverse @ bases)
    out = np.empty_like(bases)
    for i in range(bases.shape[2]):
        out[:, :, i] = orthonormalise(inverse @ bases[:, :, i])
    return out


def var_to_ss(a):
    a = np.asarray(a, float)
    n, n2, p = a.shape
    if n != n2 or p < 1:
        raise ValueError("VAR coefficients must have shape (n,n,p) with p >= 1")
    c = np.concatenate([a[:, :, k] for k in range(p)], axis=1)
    ass = np.zeros((n*p, n*p))
    ass[:n, :] = c
    ass[n:, :n*(p-1)] = np.eye(n*(p-1))
    k = np.zeros((n*p, n))
    k[:n, :] = np.eye(n)
    return ass, c, k


def iss2cak(a, c, k):
    n, r = c.shape
    cak = np.empty((n, n, r))
    ak = np.eye(r)
    for lag in range(r):
        cak[:, :, lag] = c @ ak @ k
        ak = ak @ a
    return cak


def sfreqs(fres, fs):
    return np.arange(fres + 1, dtype=float) * (fs / (2*fres))


def ss2trfun(a, c, k, fres):
    """MVGC2 ss2trfun: I + C(zI-A)^(-1)K at z=exp(i*pi*j/fres)."""
    n, r = c.shape
    h = np.empty((n, n, fres+1), dtype=complex)
    ir, inn = np.eye(r), np.eye(n)
    for j in range(fres+1):
        z = np.exp(1j*np.pi*j/fres)
        h[:, :, j] = inn + c @ np.linalg.solve(z*ir-a, k)
    return h


def var_to_cak(a):
    """Exact CAK impulse coefficients for the VAR companion innovations model.

    For a VAR of order p, the companion state has n*p coordinates. Its first
    n*p CAK slices are the first n*p coefficients of
    (I - sum(A_k z^-k))^-1. The recurrence below is the block elimination of
    C @ A_ss**j @ K, avoiding dense companion matrix powers.
    """
    a = np.asarray(a, float)
    n, n2, p = a.shape
    if n != n2 or p < 1:
        raise ValueError("VAR coefficients must have shape (n,n,p) with p >= 1")
    count = n*p
    impulse = [np.eye(n)]
    cak = np.empty((n, n, count))
    for j in range(1, count+1):
        term = np.zeros((n, n))
        for lag in range(1, min(p, j)+1):
            term += a[:, :, lag-1] @ impulse[j-lag]
        impulse.append(term)
        cak[:, :, j-1] = term
    return cak


def var_to_trfun(a, fres):
    """Exact transfer function of the VAR companion innovations model.

    Solving the companion state-space equations by block elimination gives
    H(z) = (I - sum(A_k z^-k))^-1. This is algebraically identical to
    ss2trfun(var_to_ss(a)), while solving an n-by-n system at each frequency.
    """
    a = np.asarray(a, float)
    n, n2, p = a.shape
    if n != n2 or p < 1:
        raise ValueError("VAR coefficients must have shape (n,n,p) with p >= 1")
    h = np.empty((n, n, fres+1), dtype=complex)
    identity = np.eye(n, dtype=complex)
    for j in range(fres+1):
        zinv = np.exp(-1j*np.pi*j/fres)
        polynomial = identity.copy()
        power = zinv
        for lag in range(p):
            polynomial -= a[:, :, lag]*power
            power *= zinv
        h[:, :, j] = np.linalg.solve(polynomial, identity)
    return h


def prepare_var(a, v, fs, fres):
    """Return transformed model, CAK proxy, transfer function, frequencies."""
    transformed = transform_var(a, v)
    ass, c, k = var_to_ss(transformed)
    return ass, c, k, var_to_cak(transformed), var_to_trfun(transformed, fres), sfreqs(fres, fs)
