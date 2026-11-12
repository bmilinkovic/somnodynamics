"""MVGC2 VAR fitting paths and model order selection for multitrial data."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky, lstsq, solve_triangular
from scipy.stats import f as f_distribution


@dataclass
class VARFit:
    coefficients: np.ndarray
    covariance: np.ndarray
    residuals: np.ndarray


@dataclass
class VARModelOrder:
    aic: int
    bic: int
    hqc: int
    lrt: int
    aic_curve: np.ndarray
    bic_curve: np.ndarray
    hqc_curve: np.ndarray

    def select(self, criterion="aic"):
        """Choose a model order by AIC, BIC, HQC, or sequential LRT."""
        name = criterion.lower()
        if name not in ("aic", "bic", "hqc", "lrt"):
            raise ValueError("criterion must be AIC, BIC, HQC, or LRT")
        return getattr(self, name)


@dataclass
class SelectedVAR:
    criterion: str
    order: int
    candidates: VARModelOrder
    fit: VARFit


def _demean(x):
    return x-np.mean(x, axis=(1, 2), keepdims=True)


def _lwr_prediction_errors(x, af, ab, p1n, k):
    """Form Morf forward and backward errors without a full lag tensor."""
    n, m, trials = x.shape
    kn = k*n
    effective = trials*(m-k)
    ef = np.empty((n, effective))
    eb = np.empty((n, effective))
    for start in range(0, trials, 16):
        stop = min(trials, start+16)
        width = (m-k)*(stop-start)
        forward = np.empty((kn, width))
        backward = np.empty((kn, width))
        for lag in range(k):
            rows = slice(lag*n, (lag+1)*n)
            forward[rows] = x[:, k-lag:m-lag, start:stop].reshape(
                n, width, order="F")
            backward[rows] = x[:, k-1-lag:m-1-lag, start:stop].reshape(
                n, width, order="F")
        columns = slice((m-k)*start, (m-k)*stop)
        ef[:, columns] = af[:, :kn] @ forward
        eb[:, columns] = ab[:, p1n-kn:] @ backward
    return ef, eb


def tsdata_to_var(x, p, regmode="LWR"):
    """Fit MVGC2's Morf LWR or QR-based OLS VAR model.

    Data have shape (nodes,time,trials). LWR retains the original lag-stack
    layout, Cholesky solves, recursion order, and residual covariance divisor.
    The latter uses the final loop variable M, one time point less per trial
    than the residual count; this is kept for result parity. OLS uses the
    original least-squares QR route and divides by N*(m-p)-1.
    """
    x = np.asarray(x, float)
    if x.ndim == 2:
        x = x[:, :, None]
    if x.ndim != 3:
        raise ValueError("data must have shape (nodes,time[,trials])")
    n, m, trials = x.shape
    if not 0 <= p < m:
        raise ValueError("require 0 <= p < time")
    mode = regmode.upper()
    if mode not in ("LWR", "OLS"):
        raise ValueError("regmode must be LWR or OLS")
    if mode == "LWR" and p >= m-1:
        raise ValueError("LWR requires p < time-1 to match MVGC2's covariance divisor")
    if p == 0:
        centred = _demean(x).reshape(n, m*trials, order="F")
        v = centred @ centred.T/(m*trials-1)
        return VARFit(np.zeros((n, n, 0)), v, x.copy())
    x = _demean(x)
    if mode == "OLS":
        obs = np.arange(p, m)
        count = trials*(m-p)
        current = x[:, obs, :].reshape(n, count, order="F")
        lagged = np.empty((n*p, count))
        for lag in range(p):
            lagged[lag*n:(lag+1)*n] = x[:, obs-lag-1, :].reshape(n, count, order="F")
        # MATLAB X0/XL: solve XL.T @ A.T = X0.T with rank-revealing QR.
        flat = lstsq(lagged.T, current.T, lapack_driver="gelsy")[0].T
        coefficients = np.stack([flat[:, lag*n:(lag+1)*n] for lag in range(p)], axis=2)
        e = current-flat @ lagged
        v = e @ e.T/(count-1)
        return VARFit(coefficients, v, e.reshape(n, m-p, trials, order="F"))
    p1 = p+1
    p1n = p1*n
    ee = x.transpose(0, 1, 2).reshape(n, trials*m, order="F")
    c = cholesky(ee @ ee.T, lower=True)
    ic = np.linalg.inv(c)
    af = np.zeros((n, p1n))
    ab = np.zeros((n, p1n))
    af[:, :n] = ic
    ab[:, -n:] = ic
    eye = np.eye(n)
    k = 1
    while k <= p:
        kn = k*n
        # Build the exact MATLAB lag-stack columns in small groups of trials.
        # The former four-dimensional XX array repeated every observation
        # p+1 times, which exceeded RAM for long recordings. Each block has
        # the same node/lag/time/trial ordering as XX's Fortran reshape.
        ef, eb = _lwr_prediction_errors(x, af, ab, p1n, k)
        cf = cholesky(ef @ ef.T, lower=True)
        cb = cholesky(eb @ eb.T, lower=True)
        reflection = solve_triangular(cf, ef, lower=True) @ solve_triangular(cb, eb, lower=True).T
        k += 1
        kn = k*n
        effective = trials*(m-k)
        previous_f = af[:, :kn].copy()
        previous_b = ab[:, p1n-kn:].copy()
        scale_f = cholesky(eye-reflection @ reflection.T, lower=True)
        scale_b = cholesky(eye-reflection.T @ reflection, lower=True)
        af[:, :kn] = solve_triangular(scale_f, previous_f-reflection @ previous_b, lower=True)
        ab[:, p1n-kn:] = solve_triangular(scale_b, previous_b-reflection.T @ previous_f, lower=True)
    a0 = af[:, :n]
    flat = -np.linalg.solve(a0, af[:, n:p1n])
    coefficients = np.stack([flat[:, lag*n:(lag+1)*n] for lag in range(p)], axis=2)
    e = np.linalg.solve(a0, ef)
    v = e @ e.T/(effective-1)
    residuals = e.reshape(n, m-p, trials, order="F")
    return VARFit(coefficients, v, residuals)


def _lwr_loglikelihoods(x, pmax):
    """Compute every LWR order in one recursion, as MVGC2 does."""
    x = _demean(np.asarray(x, float))
    n, m, trials = x.shape
    p1n = (pmax+1)*n
    ee = x.transpose(0, 1, 2).reshape(n, trials*m, order="F")
    c = cholesky(ee @ ee.T, lower=True)
    ic = np.linalg.inv(c)
    af = np.zeros((n, p1n))
    ab = np.zeros((n, p1n))
    af[:, :n] = ic
    ab[:, -n:] = ic
    eye = np.eye(n)
    likelihood = np.empty(pmax)
    for order in range(1, pmax+1):
        kn = order*n
        ef, eb = _lwr_prediction_errors(x, af, ab, p1n, order)
        cf = cholesky(ef @ ef.T, lower=True)
        cb = cholesky(eb @ eb.T, lower=True)
        reflection = solve_triangular(cf, ef, lower=True) @ solve_triangular(
            cb, eb, lower=True).T
        next_kn = (order+1)*n
        previous_f = af[:, :next_kn].copy()
        previous_b = ab[:, p1n-next_kn:].copy()
        scale_f = cholesky(eye-reflection @ reflection.T, lower=True)
        scale_b = cholesky(eye-reflection.T @ reflection, lower=True)
        af[:, :next_kn] = solve_triangular(
            scale_f, previous_f-reflection @ previous_b, lower=True)
        ab[:, p1n-next_kn:] = solve_triangular(
            scale_b, previous_b-reflection.T @ previous_f, lower=True)
        e = np.linalg.solve(af[:, :n], ef)
        effective = trials*(m-order-1)
        covariance_ml = e @ e.T/effective
        likelihood[order-1] = -np.linalg.slogdet(covariance_ml)[1]/2
    return likelihood


def tsdata_to_varmo(x, pmax, regmode="LWR", alpha=0.01):
    """MVGC2 VAR order selection by AIC, BIC, HQC, and sequential LRT.

    Mirrors ``tsdata_to_varll`` and ``tsdata_to_varmo`` including order zero,
    each order's loop-variable effective sample size, and Bonferroni LRT.
    """
    x = np.asarray(x, float)
    if x.ndim == 2:
        x = x[:, :, None]
    n, m, trials = x.shape
    mode = regmode.upper()
    if mode not in ("LWR", "OLS"):
        raise ValueError("regmode must be LWR or OLS")
    max_order = m-2 if mode == "LWR" else m-1
    if not 1 <= pmax <= max_order:
        raise ValueError("pmax is too large for the data and regression mode")
    likelihood = np.empty(pmax+1)
    parameters = np.arange(pmax+1)*n*n
    sizes = np.empty(pmax+1)
    centred = _demean(x).reshape(n, m*trials, order="F")
    cov0 = centred @ centred.T/(m*trials)
    likelihood[0] = -np.linalg.slogdet(cov0)[1]/2
    sizes[0] = m*trials
    if mode == "LWR":
        # MATLAB tsdata_to_varll advances one Morf recursion through pmax.
        # Each order's covariance uses the loop's post-increment sample size.
        sizes[1:] = trials*(m-np.arange(1, pmax+1)-1)
        likelihood[1:] = _lwr_loglikelihoods(x, pmax)
    else:
        for p in range(1, pmax+1):
            fit = tsdata_to_var(x, p, regmode)
            size = trials*(m-p)
            sizes[p] = size
            covariance_ml = fit.residuals.reshape(n, -1, order="F")
            covariance_ml = covariance_ml @ covariance_ml.T/size
            likelihood[p] = -np.linalg.slogdet(covariance_ml)[1]/2
    ratio = parameters/sizes
    aic = -2*likelihood+2*ratio
    bic = -2*likelihood+ratio*np.log(sizes)
    hqc = -2*likelihood+2*ratio*np.log(np.log(sizes))
    statistic = 2*np.diff(sizes*likelihood)
    df = n*n
    df2 = sizes[1:]-n*np.arange(1, pmax+1)-1
    pvalues = np.full(pmax, np.inf)
    good = df2 > 0
    pvalues[good] = f_distribution.sf(statistic[good]/df, df, df2[good])
    hits = np.flatnonzero(pvalues < alpha/pmax)
    lrt = int(hits[-1]+1) if len(hits) else 0
    return VARModelOrder(int(np.nanargmin(aic)), int(np.nanargmin(bic)),
                         int(np.nanargmin(hqc)), lrt, aic, bic, hqc)


def fit_var_auto(x, pmax, criterion="aic", *, regmode="LWR", alpha=0.01):
    """Select one of the four MVGC2 orders and fit an LWR or OLS VAR model.

    Returns the selected order, fitted coefficients/covariance/residuals, and
    all four candidate orders for inspection. LRT may select order zero.
    """
    candidates = tsdata_to_varmo(x, pmax, regmode=regmode, alpha=alpha)
    order = candidates.select(criterion)
    fit = tsdata_to_var(x, order, regmode=regmode)
    return SelectedVAR(criterion.upper(), order, candidates, fit)
