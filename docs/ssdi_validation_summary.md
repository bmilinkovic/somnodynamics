# SSDI Validation Summary (DMT-Faithful Stage Pipeline)

Generated: 2026-03-17T14:26:51

## Coverage
- Input subject-stage cases found: **140**
- Output subject-stage cases found: **1**
- Stage input counts:
  - W_: 28
  - N1: 28
  - N2: 28
  - N3: 28
  - R_: 28
- Stage output counts:
  - W_: 0
  - N1: 1
  - N2: 0
  - N3: 0
  - R_: 0

## Optimizer Validation
- Macro sizes per case: [2, 3, 4, 5, 6] (count=5)
- Bands per case: ['delta', 'theta', 'alpha', 'sigma', 'beta', 'gamma'] (count=6)
- Expected optimization blocks per subject-stage case:
  - preopt calls: 30
  - second-stage calls: 30
- Representative config from output meta: nrunsp=100, nitersp=10000, niterso=10000
- Final objective: `canonical_band_dd_eq25_26` (optimiser `canonical_opt_dds_mruns`)
- Pre objective: `proxy_ddx_cak` (optimiser `canonical_opt_ddx_mruns`)

## Runtime
- Status rows: 1
- Runtime seconds (min/median/max): 112.07 / 112.07 / 112.07

## Structural Checks
- Files with missing required DMT-style key families: 0
- All inspected outputs contained required key families for all band/m pairs.

## Representative Output
- File: `results/ssdi/stage_ssdi_anatomical7/N1/EPCTL01_ssdi_spectral_domain_results.mat`
- Subject: `EPCTL01`
- Stage: `N1`
- Shape summary:
  - delta_dopto_m2: {"shape": [92], "dtype": "float64"}
  - delta_doptp_m2: {"shape": [100], "dtype": "float64"}
  - delta_Lopto_m2: {"shape": [7, 2, 92], "dtype": "float64"}
  - delta_Loptp_m2: {"shape": [7, 2, 100], "dtype": "float64"}
  - delta_node_weights_m2: {"shape": [7, 92], "dtype": "float64"}
  - delta_beta_m2: {"shape": [7, 92], "dtype": "float64"}
  - gamma_dopto_m6: {"shape": [28], "dtype": "float64"}
  - gamma_Lopto_m6: {"shape": [7, 6, 28], "dtype": "float64"}

## DMT Structure Comparison
- Reference file: `/Users/borjan/code/dmt_project/results/ssdi_analysis/rest/P2M63_ssdi_spectral_domain_results.mat`
- Somno key count: 638
- DMT key count: 125
- Shared keys: 125
- Somno-only keys: 513
- DMT-only keys: 0
- Note: available DMT reference output appears from an older result set; key-count differences are expected.
