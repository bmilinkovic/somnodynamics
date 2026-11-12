# Multiscale · dynamical dependence

Find a small set of variables that captures a system's dynamics. This
standalone Python package fits multivariate autoregressive (VAR) models and
searches for macroscopic variables with low dynamical dependence (DD). It can
measure DD directly in the **time domain**, by integration in the **spectral
domain**, or within named recording windows to follow **time-changing DD**.

The mathematical routines follow Lionel Barnett's SSDI and MVGC2 MATLAB
implementations. The Python package needs no MATLAB installation.

## Install

```bash
python -m pip install ./multiscale
```

For group tables and figures, install `./multiscale[plots]`. The package has no
study data, atlases, or fixed cohort names. The
[quickstart notebook](notebooks/quickstart.ipynb) runs the main workflows on a
small synthetic recording.

## Prepare a model

Pass signals as `(nodes, samples, trials)`. A single trial can be a two-dimensional
`(nodes, samples)` array. Decide how to clean and segment your recordings
**before** fitting. EEG and fMRI differ in sampling rate, filtering, artefacts,
stationarity, and the meaning of a sample. Keep preprocessing and node order
consistent across subjects and conditions.

### fMRI: use particular caution

**For now, treat DI/DD computed from fMRI BOLD recordings as exploratory, and
do not interpret it as a direct estimate of neural-level dynamical
independence.** The BOLD signal samples neural activity indirectly and much
more slowly than many neural interactions. Severe downsampling can change
which predictive relationships are detectable or even produce misleading
Granger-causality results. The haemodynamic response can also affect those
results: invariance to its filtering depends on a *minimum-phase* condition.
Recent modelling finds that many plausible response functions satisfy this
condition, but some parameter choices, including ones producing an initial
dip, do not. Even a minimum-phase response does not remove the sampling
problem. These papers study Granger causality rather than this package's DI
estimator directly; the warning for DI follows from its use of predictive
dependence. The package does not test the response function or correct for
these measurement effects. See [Seth, Chorley & Barnett (2013)](https://doi.org/10.1016/j.neuroimage.2012.09.049),
[Barnett & Seth (2017)](https://doi.org/10.1016/j.jneumeth.2016.10.016), and
[Novelli, Barnett, Seth & Razi (2025)](https://doi.org/10.1002/hbm.70285).

```python
from multiscale import fit_var_auto

selected = fit_var_auto(signals, pmax=20, criterion="bic", regmode="LWR")
A = selected.fit.coefficients  # nodes × nodes × lags
V = selected.fit.covariance    # innovation covariance
print(selected.order, selected.candidates)
```

Choose AIC, BIC, HQC, or LRT for the model order; choose LWR or OLS for fitting.
LWR follows Morf's multichannel recursion. OLS solves the lagged least-squares
problem by rank-revealing QR. An order-zero model has no dynamics to optimise.

## Spectral DD

Name the frequency bands that make sense for your signals and sampling rate.
The bands below are examples, not prescribed ranges.

```python
from multiscale import run_grid_parallel

results = run_grid_parallel(
    A, V, fs=128, fres=128,
    dimensions=(1, 2),
    bands={"slow": (1, 4), "fast": (8, 20)},
    runs=100, workers=4, seed=42,
)
print(results[(2, "slow")]["final"].values[0])
```

Each dimension × band search first optimises a fast proxy, groups near-identical
solutions, and then optimises band-limited DD. `final.values` are sorted from
lowest DD upward. The associated bases are in decorrelated coordinates;
`final_original_bases` maps them back to the original signals. Increase the
number of restarts and iteration limits for a serious analysis. Dimensions
must satisfy `1 <= dimension < number of nodes`.

## Time-domain DD

Direct temporal DD uses a state-space Riccati calculation. It is more expensive
than the spectral gradient route, since each candidate subspace needs another
Riccati solve.

```python
from multiscale import random_initial_bases, run_full_dd

starts = random_initial_bases(A.shape[0], m=2, runs=100, seed=42)
result = run_full_dd(A, V, fs=128, fres=128, initial_bases=starts,
                     final_method="state_space", state_space_workers=4)
print(result["state_space"].values[0])
```

`run_full_dd` also supports `final_method="spectral"` for full-spectrum
transfer-function integration, `"both"` to compare the two methods from the
same proxy solutions, and `"spectral_then_state_space"` for sequential
refinement. Band-limited spectral DD uses `run_grid_parallel` instead.

## Time-changing DD

A recording can change state. Give each period a name and a `(start, stop)`
range in **sample indices**; `stop` is excluded. The package fits a fresh VAR
model in every window and performs a fresh DD optimisation there. This is a
series of local, stationary-model estimates, rather than a single model whose
parameters vary continuously. Windows may overlap and can have different
lengths, but each must contain enough samples for its VAR order.

```python
from multiscale import run_time_changing_dd, minimum_dd_by_window

windows = {"before": (0, 1000), "during": (1000, 2000),
           "after": (2000, 3000)}
trajectory = run_time_changing_dd(
    signals, windows, domain="spectral", fs=128, fres=128,
    dimensions=(1, 2), bands={"slow": (1, 4)},
    pmax=20, criterion="bic", regmode="LWR", runs=100, workers=4,
    progress=lambda name, done, total: print(done, total, name),
)
rows = minimum_dd_by_window(trajectory, domain="spectral")
```

Set `domain="temporal"` and omit `bands` for direct state-space DD in each
window. `trajectory["during"].model` contains that window's selected fit;
`trajectory["during"].blocks` contains its optimisation results. The summary
rows include best DD, model order, and number of final solutions. Treat a
change in DD alongside a change in model order or fit quality. Use the same
preprocessing, node definitions, and frequency grid in every window.

For a quick figure with `multiscale[plots]`:

```python
from multiscale.analysis import summarise_window_results
from multiscale.visualisation import plot_time_changing_dd

window_table = summarise_window_results(trajectory, domain="spectral")
plot_time_changing_dd(window_table, "figures/dd_over_time")
```

## Group analysis and figures

`multiscale.analysis` converts optimised solutions to tidy subject tables and
run-similarity matrices. It uses British spellings such as
`summarise_grid_results` and `summarise_window_results`; the earlier
`summarize_*` names remain available for existing scripts.
`multiscale.statistics` provides paired Wilcoxon tests,
Benjamini–Hochberg correction, and descriptive pooled comparisons.
`multiscale.visualisation` saves publication figures as PDF and PNG.

```python
from multiscale.analysis import summarise_grid_results
from multiscale.statistics import paired_tests
from multiscale.visualisation import plot_dd_violin_by_macro

# Repeat this for each participant and concatenate the resulting tables.
dd, sf, matrices = summarise_grid_results(results,
                                          subject="participant_01",
                                          condition="control")
# After combining all participants and both conditions into all_dd:
stats = paired_tests(all_dd, "dd_min", "control", "practice",
                     groups=["m", "band"], correction_groups=["m"])
plot_dd_violin_by_macro(all_dd, stats, "figures/dd_by_dimension",
                        conditions=("control", "practice"))
```

The structure-fragmentation (SF) index is the fraction of pairs of final
solutions whose subspace distance is below a chosen threshold. Compare SF only
with the same threshold and restart settings. The pooled pair-count test is
descriptive because pairs from one subject are dependent; use subject-level
paired tests for the main inference. Specify the family of comparisons for
multiple-testing correction before interpreting adjusted p values.

## Licence and credit

Copyright (c) 2026 Borjan Milinkovic. Multiscale is distributed under
[GPL-3.0-only](LICENSE). Core routines were translated and adapted from
Lionel Barnett's [SSDI](https://github.com/lcbarnett/ssdi) and
[MVGC2](https://github.com/lcbarnett/MVGC2) toolboxes; see
[NOTICE.md](NOTICE.md) for attribution. The separate Synergy package has its
own MIT licence.

## Source papers

Our applications of this framework:

- Milinkovic, B., Barnett, L., Carter, O., Seth, A. K., & Andrillon, T.
  (2025). “Capturing the emergent dynamical structure in biophysical neural
  models.” *PLOS Computational Biology, 21*(5), e1012572.
  [doi:10.1371/journal.pcbi.1012572](https://doi.org/10.1371/journal.pcbi.1012572).
- Milinkovic, B., Seth, A. K., Barnett, L., Carter, O., & Andrillon, T.
  (2026). “Emergent multiscale organisation of neural dynamics fragments in
  anaesthesia.” *Imaging Neuroscience, 4*, IMAG.a.1364.
  [doi:10.1162/IMAG.a.1364](https://doi.org/10.1162/IMAG.a.1364).

Foundational methods and source toolboxes:

- Barnett, L., & Seth, A. K. (2023). “Dynamical independence: Discovering
  emergent macroscopic processes in complex dynamical systems.” *Physical
  Review E, 108*, 014304. [doi:10.1103/PhysRevE.108.014304](https://doi.org/10.1103/PhysRevE.108.014304).
  The paper underlying the SSDI approach.
- Barnett, L., & Seth, A. K. (2014). “The MVGC multivariate Granger causality
  toolbox: A new approach to Granger-causal inference.” *Journal of
  Neuroscience Methods, 223*, 50–68.
  [doi:10.1016/j.jneumeth.2013.10.018](https://doi.org/10.1016/j.jneumeth.2013.10.018).
  This describes the MVGC toolbox lineage from which MVGC2 develops; it is
  not a separate publication of the MVGC2 repository.
