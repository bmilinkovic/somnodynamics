# Pipeline Profile: Artifact-Removed Preprocessed Sensor Dataset (Setup Only)

This profile is preconfigured to run **from**:
- `/Users/borjan/code/somnodynamics/results/preprocessing/preprocessed_sensor_space_artifact_removed_256hz`

and therefore skips preprocessing.

## Fixed assumptions
- Sensor epoching is unchanged: 10 s epochs at 256 Hz.
- Source windows are unchanged: 2 s non-overlapping windows.
- Downstream branch: DK5/DK68 + GC (LRT, momax=30) + SSDI.

## 1) Source reconstruction (DK5 + DK68)
```bash
NUMBA_DISABLE_JIT=1 conda run -n mne python /Users/borjan/code/somnodynamics/code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/src_recon/run_source_reconstruction_dk5_dk68.py \
  --sensor-dir /Users/borjan/code/somnodynamics/results/preprocessing/preprocessed_sensor_space_artifact_removed_256hz \
  --raw-root /Volumes/ex_data/data_somnoemergence/raw \
  --target-sfreq 256 \
  --epoch-seconds 2 \
  --min-epochs 15
```

## 2) GC (DK5, LRT)
```matlab
run('/Users/borjan/code/somnodynamics/code/run_all/setup_matlab_toolboxes.m');
addpath('/Users/borjan/code/somnodynamics/code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc');
run_gc_sleep_stage_pipeline_dk5('mode','stage_mvgc2','atlas','dk5','momax',30,'mosel','LRT','fs_original',256,'fs_target',256,'fres',256);
```

## 3) GC (DK68, LRT) — only if you want SSDI-B
```matlab
run_gc_sleep_stage_pipeline_dk5('mode','stage_mvgc2','atlas','dk68','momax',30,'mosel','LRT','fs_original',256,'fs_target',256,'fres',256);
```

## 4) SSDI (DK5)
```matlab
addpath('/Users/borjan/code/somnodynamics/code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi');
run_ssdi_sleep_stage_pipeline_dk('mode','stage_ssdi','atlas','dk5','nrunsp',100,'nitersp',10000,'niterso',10000);
```

## 5) SSDI (DK68)
```matlab
run_ssdi_sleep_stage_pipeline_dk('mode','stage_ssdi','atlas','dk68','nrunsp',100,'nitersp',10000,'niterso',10000);
```

## 6) Figures and summaries
```bash
python3 /Users/borjan/code/somnodynamics/code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc/plot_gc_stage_figures_dk5.py

python3 /Users/borjan/code/somnodynamics/code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/summarize_ssdi_results_dk.py \
  --spectral-dir /Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz \
  --out-dir /Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz

python3 /Users/borjan/code/somnodynamics/code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/plot_ssdi_stage_figures_dk.py \
  --long-csv /Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/ssdi_subject_long.csv \
  --spectral-dir /Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz \
  --out-dir /Users/borjan/code/somnodynamics/results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/figures
```
