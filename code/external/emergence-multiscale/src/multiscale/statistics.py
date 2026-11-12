"""Study-independent paired DD and structure-fragmentation statistics."""

import numpy as np
import pandas as pd
from scipy.stats import norm, ranksums, wilcoxon
from itertools import combinations


def bh_fdr(pvalues, alpha=0.05):
    """Benjamini-Hochberg correction, preserving NaN positions."""
    p = np.asarray(pvalues, float)
    q = np.full(p.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    if valid.size:
        order = valid[np.argsort(p[valid], kind="stable")]
        ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
        q[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0, 1)
    return q, np.isfinite(q) & (q <= alpha)


def sign_balance(delta):
    """Fraction of positive changes minus fraction of negative changes."""
    d = np.asarray(delta, float)
    d = d[np.isfinite(d) & (d != 0)]
    return float(np.sign(d).mean()) if d.size else np.nan


def paired_tests(table, value, condition_a, condition_b, groups,
                 correction_groups, alpha=0.05):
    """Paired two-sided Wilcoxon test using an approximate p value.

    ``condition_b - condition_a`` defines every reported effect. BH is
    applied separately within ``correction_groups``. One row per subject,
    condition and requested group is required.
    """
    cols = ["subject", "condition", value, *groups]
    frame = table[cols].dropna(subset=[value])
    if frame.duplicated(["subject", "condition", *groups]).any():
        raise ValueError("duplicate subject/condition/group rows")
    rows = []
    for key, sub in frame.groupby(list(groups), sort=False, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        left = sub[sub.condition == condition_a][["subject", value]]
        right = sub[sub.condition == condition_b][["subject", value]]
        pair = left.merge(right, on="subject", suffixes=("_a", "_b"))
        if len(pair) < 3:
            continue
        a, b = pair[f"{value}_a"].to_numpy(), pair[f"{value}_b"].to_numpy()
        delta = b-a
        try:
            stat, p = wilcoxon(b, a, alternative="two-sided",
                               zero_method="wilcox", method="approx")
        except ValueError:
            stat, p = np.nan, np.nan
        rows.append(dict(zip(groups, key)) | dict(
            n_subjects=len(pair), delta_mean=float(np.mean(delta)),
            delta_median=float(np.median(delta)),
            effect_rbc=sign_balance(delta), stat=float(stat), p=float(p)))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q"] = np.nan
    out["significant"] = False
    for _, indexes in out.groupby(list(correction_groups), sort=False).groups.items():
        q, significant = bh_fdr(out.loc[indexes, "p"], alpha)
        out.loc[indexes, "q"] = q
        out.loc[indexes, "significant"] = significant
    return out


def grand_average_dd(table, condition_a, condition_b, *, paired=False):
    """Full-run pooled rank-sum or time-window subject-mean paired test."""
    rows = []
    group_cols = ["window", "band"] if "window" in table else ["band"]
    for key, sub in table.groupby(group_cols, sort=False):
        key = key if isinstance(key, tuple) else (key,)
        if paired:
            mean = sub.groupby(["subject", "condition"], as_index=False)["dd_min"].mean()
            a = mean[mean.condition == condition_a][["subject", "dd_min"]]
            b = mean[mean.condition == condition_b][["subject", "dd_min"]]
            pair = a.merge(b, on="subject", suffixes=("_a", "_b"))
            if len(pair) < 3:
                continue
            av, bv = pair.dd_min_a.to_numpy(), pair.dd_min_b.to_numpy()
            stat, p = wilcoxon(bv, av, alternative="two-sided",
                               zero_method="wilcox", method="approx")
            n_a = n_b = len(pair)
            effect = sign_balance(bv-av)
        else:
            av = sub[sub.condition == condition_a].dd_min.to_numpy()
            bv = sub[sub.condition == condition_b].dd_min.to_numpy()
            if min(len(av), len(bv)) < 3:
                continue
            stat, p = ranksums(av, bv)
            n_a, n_b, effect = len(av), len(bv), np.nan
        rows.append(dict(zip(group_cols, key)) | dict(
            n_a=n_a, n_b=n_b, mean_a=float(np.mean(av)), mean_b=float(np.mean(bv)),
            delta_mean=float(np.mean(bv)-np.mean(av)), effect_rbc=effect,
            stat=float(stat), p=float(p)))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q"] = np.nan
    out["significant"] = False
    families = out.groupby("window").groups if "window" in out else {None: out.index}
    for _, indexes in families.items():
        q, significant = bh_fdr(out.loc[indexes, "p"])
        out.loc[indexes, "q"] = q
        out.loc[indexes, "significant"] = significant
    return out


def sf_index(distance_matrix, threshold=1e-4):
    """Fraction of off-diagonal Grassmann distances below the threshold."""
    g = np.asarray(distance_matrix, float)
    if g.ndim != 2 or g.shape[0] != g.shape[1]:
        raise ValueError("similarity matrix must be square")
    tri = g[np.triu_indices(len(g), 1)]
    return float(np.mean(tri < threshold)) if tri.size else np.nan


def pooled_sf_stats(sf_table, condition_a, condition_b):
    """Compare pooled near-zero run-pair proportions with a z test.

    Run pairs from one subject are dependent, so this is descriptive. Use a
    subject-level paired test for the main group inference.
    """
    rows = []
    for (band, m), sub in sf_table.groupby(["band", "m"], sort=False):
        counts = {}
        for condition in (condition_a, condition_b):
            s = sub[sub.condition == condition]
            total = int(s.n_pairs.sum())
            near = int(np.rint((s.sf_index*s.n_pairs).sum()))
            counts[condition] = (near, total)
        (x1, n1), (x2, n2) = counts[condition_a], counts[condition_b]
        if min(n1, n2) == 0:
            continue
        p1, p2 = x1/n1, x2/n2
        pooled = (x1+x2)/(n1+n2)
        se = np.sqrt(pooled*(1-pooled)*(1/n1+1/n2))
        z = (p2-p1)/se if se > 0 else np.nan
        p = 2*norm.sf(abs(z)) if np.isfinite(z) else np.nan
        rows.append(dict(band=band, m=m, n_pairs_a=n1, n_pairs_b=n2,
                         sf_a=p1, sf_b=p2, delta=p2-p1, z=z, p=p))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q"] = np.nan
    out["significant"] = False
    for _, indexes in out.groupby("band").groups.items():
        q, significant = bh_fdr(out.loc[indexes, "p"])
        out.loc[indexes, "q"] = q
        out.loc[indexes, "significant"] = significant
    return out


def window_state_table(table, value, *, baseline="rest", baseline_label="Rest"):
    """Deduplicate a shared baseline and name each nonbaseline window as a state."""
    base = table[table.condition == baseline].copy()
    base["state"] = baseline_label
    base = base.drop_duplicates(["subject", "band", "m"])
    active = table[table.condition != baseline].copy()
    active["state"] = active["window"]
    return pd.concat([base, active], ignore_index=True)[
        ["subject", "state", "band", "m", value]]


def state_pairwise_tests(table, value, states, *, detail_correction_groups=("m",)):
    """All paired state contrasts, with caller-selected BH families.

    Use ``("m",)`` for DD comparisons across bands, or ``("band",)``
    for fragmentation comparisons across macro dimensions.
    """
    detailed, grand = [], []
    for state_a, state_b in combinations(states, 2):
        work = table.rename(columns={"state": "condition"})
        detail = paired_tests(work, value, state_a, state_b,
                              ["m", "band"], list(detail_correction_groups))
        detail.insert(0, "state_a", state_a)
        detail.insert(1, "state_b", state_b)
        detailed.append(detail)
        averaged = work.groupby(["subject", "condition", "band"],
                                as_index=False)[value].mean()
        summary = paired_tests(averaged, value, state_a, state_b,
                               ["band"], ["band"])
        # Paired_tests corrects within each band; redo BH across all bands.
        if len(summary):
            q, significant = bh_fdr(summary.p)
            summary["q"] = q
            summary["significant"] = significant
        summary.insert(0, "state_a", state_a)
        summary.insert(1, "state_b", state_b)
        grand.append(summary)
    return pd.concat(detailed, ignore_index=True), pd.concat(grand, ignore_index=True)
