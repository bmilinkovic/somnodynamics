# SomnoDynamics Figure Build

SomnoDynamics creates stage-specific figures during analysis and then consolidates
publication-facing assets under `results/figures`.

## GC DMT-style figures

```bash
MPLCONFIGDIR=/tmp/mpl python3 code/gc/plot_gc_stage_figures.py
```

Primary outputs:
- `results/gc/stage_mvgc2_anatomical7/figures/time_*_wake_contrasts.(png|svg|pdf)`
- `results/gc/stage_mvgc2_anatomical7/figures/spectral_*_wake_contrast_*.{png,svg,pdf}`

## SSDI DMT-style figures

```bash
MPLCONFIGDIR=/tmp/mpl python3 code/ssdi/plot_ssdi_stage_figures.py
```

Primary outputs:
- `results/ssdi/figures/ssdi_dd_median_violin_by_macro.(png|svg|pdf)`
- `results/ssdi/figures/ssdi_dd_median_grand_average.(png|svg|pdf)`
- `results/ssdi/figures/ssdi_dd_median_wake_contrasts.(png|svg|pdf)`

## Build consolidated figure set

```bash
python3 code/figures/build_publication_figures.py
```

Outputs:
- `results/figures/<stage>/*`
- `results/figures/figure_manifest.csv`
- `results/figures/figure_summary.txt`

Supported source stages:
- sensor power
- GC statistics
- SSDI
- GC-SSDI cross-analysis
