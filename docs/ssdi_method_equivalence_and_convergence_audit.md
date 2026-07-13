# SSDI Method Equivalence and Convergence Audit (SomnoDynamics vs DMT)

Date: 2026-03-19

## Scope
Audit requested for:
- spectral SSDI equation-level equivalence to DMT canonical implementation
- preoptimisation + final optimisation structure equivalence
- convergence behavior diagnostics from current `momax4 / 400epochs` run logs

Run context audited:
- Somno run root:  
  `/Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_anatomical7_artifact_removed_256hz_complete5_center400_momax4_aic`
- Main run log:  
  `/Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_anatomical7_artifact_removed_256hz_complete5_center400_momax4_aic/logs/stage_ssdi_pipeline_20260319_163006.log`
- Run was manually stopped on request; 13 subject-stage files completed before stop.

## 1) Spectral SSDI structure (Eq. 25/26 canonical form)

### Result
Spectral canonical objective implementation matches DMT.

### Evidence
- File-level identity (`md5`) between Somno and DMT for:
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/canonical_band_objective_and_grad.m`
  - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/canonical_band_objective_and_grad.m`
- Same for:
  - `canonical_opt_dds_mruns.m`
  - `canonical_opt_ddx_mruns.m`

### Equation-level details (Somno)
- Canonical per-frequency objective and gradient:
  - `f(L) = log|A| - log|C|` with
    - `A = L'(H H*)L`
    - `B = L'HL`
    - `C = B B*`
  - implemented in:
    - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/canonical_band_objective_and_grad.m:85`
- Band objective is weighted frequency average of per-frequency canonical DD terms:
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/canonical_band_objective_and_grad.m:56`
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/canonical_band_objective_and_grad.m:138`

## 2) Preoptimisation and optimisation structure

### Result
Somno pipeline structure matches DMT canonical structure (sleep-stage adaptation only in I/O axis and band definitions).

### Stage pipeline structure (Somno)
- Preoptimisation:
  - random orthonormal multistarts (`nrunsp=100`)
  - proxy DDX objective via `canonical_opt_ddx_mruns`
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/sd_ssdi01_compute_stage_spectral.m:179`
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/sd_ssdi01_compute_stage_spectral.m:182`
- Restart reduction:
  - similarity clustering via `Lcluster` from pre-stage outputs
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/sd_ssdi01_compute_stage_spectral.m:188`
- Final canonical optimisation:
  - canonical Eq.25/26 objective via `canonical_opt_dds_mruns`
  - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/sd_ssdi01_compute_stage_spectral.m:245`

### DMT canonical reference structure
- Same pre → cluster → final structure:
  - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/run_spectral_domain_optimisation_dmt_anatomical7_canonical.m:267`
  - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/run_spectral_domain_optimisation_dmt_anatomical7_canonical.m:272`
  - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/run_spectral_domain_optimisation_dmt_anatomical7_canonical.m:329`

### Hyperparameter parity
- Somno and DMT use matching canonical defaults:
  - `nrunsp=100`, `nitersp=10000`, `niterso=10000`, `gdtol*=1e-10`, `sig0p=1`, `sig0o=0.1`
  - Somno config:
    - `/Users/borjan/code/somnodynamics/code/ssdi/stage_ssdi_anatomical7/sd_ssdi_config_anatomical7.m:32`
  - DMT config:
    - `/Users/borjan/code/dmt_project/code/ssdi_anatomical7/ssdi_anatomical7_config.m:39`

## 3) Convergence behavior diagnostics

Diagnostics were extracted from the active run log before stop.

### Log-level convergence summary
- Pre-stage runs parsed: `41,700`
  - iterations: min `35`, median `163`, p90 `481`, p95 `763`, max `8409`
  - fast runs (`<=50` iterations): `2,005` (`~4.8%`)
  - all labeled `converged(1)` in log output
- Final-stage runs parsed: `36,070`
  - iterations: min `61`, median `503`, p90 `1873`, p95 `3236`, max `10000`
  - fast runs (`<=50` iterations): `0`
  - convergence labels:
    - `converged(1)`: `35,813`
    - `converged(3)`: `11`
    - `unconverged`: `246` (`~0.68%`)

### Restart diversity and collapse checks (from completed files)
- Completed files before stop: `13` (all currently in stage `N1`)
- Across 390 `(subject,band,m)` entries in these files:
  - pre-stage cluster count (`Lcluster`) median: `93` (min `28`, max `100`)
  - largest cluster share median: `0.05`, p90: `0.28`, max: `0.51` (1 case)
  - interpretation: most runs do **not** collapse into a single near-identical basin; diversity is generally high.

### Representative objective trajectories (EPCTL10, N1)
- File:
  - `/Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_anatomical7_artifact_removed_256hz_complete5_center400_momax4_aic/N1/EPCTL10_ssdi_spectral_domain_results.mat`
- `gamma, m=2` final-stage trajectories:
  - best run: `dd 8.92e-3 -> 5.77e-4` in `159` iterations
  - median run: `dd 7.29e-3 -> 5.77e-4` in `487` iterations
  - worst run: `dd 5.69e-3 -> 1.04e-3` in `331` iterations
- `gamma, m=2` pre-stage trajectories:
  - best run: `dd 43.11 -> 6.23` in `105` iterations
  - worst run: `dd 43.93 -> 16.57` in `39` iterations
  - fast-pre examples exist, but these are filtered by clustering and followed by full canonical final-stage optimisation.

## 4) Suspicious behavior flags

### Methodological mismatch
- Not detected for canonical spectral SSDI equations or optimisation structure.

### Operational issues detected
- Parallel worker instability during run:
  - warning count (`worker aborted ... parfor`): `11`
  - serial fallback warnings in final optimiser: `8`
- One explicit subject-stage failure due pool shutdown in status CSV:
  - `/Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_anatomical7_artifact_removed_256hz_complete5_center400_momax4_aic/logs/ssdi_anatomical7_stage_status.csv`

## 5) Conclusion

1. **Spectral SSDI method**: matches DMT canonical implementation (Eq.25/26 structure and decomposition).
2. **Optimisation structure**: matches DMT canonical (preoptimisation + clustering + final multistart canonical optimisation).
3. **Convergence health**:
   - final-stage behavior is generally healthy (median iterations high; zero ultra-fast runs; low unconverged fraction),
   - but runtime robustness is limited by MATLAB parallel-worker instability in this environment.

## 6) Optional improvements (only operational, not methodological)

If desired, next safe steps without changing method:
- force serial (`parallel_enable=false`) or lower workers to reduce worker-abort failures,
- keep canonical objective and restart structure unchanged.

These are execution-stability adjustments, not algorithm redesign.
