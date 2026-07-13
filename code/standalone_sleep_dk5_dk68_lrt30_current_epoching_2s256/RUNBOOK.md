# Standalone Runbook (Current Epoching Preserved)

## 1) Source reconstruction (DK5 + DK68 outputs)
```bash
python3 code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/src_recon/run_source_reconstruction_dk5_dk68.py \
  --sensor-dir /Users/borjan/code/somnodynamics/results/preprocessing/preprocessed_sensor_space_artifact_removed_256hz \
  --target-sfreq 256 \
  --epoch-seconds 2 \
  --min-epochs 15
```

## 2) GC (DK5, LRT momax=30)
```matlab
run('code/run_all/setup_matlab_toolboxes.m');
addpath('code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc');
run_gc_sleep_stage_pipeline_dk5('mode','stage_mvgc2','atlas','dk5','momax',30,'mosel','LRT','fs_original',256,'fs_target',256,'fres',256);
```

## 3) GC figures (DK5)
```bash
python3 code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/gc/plot_gc_stage_figures_dk5.py
```

## 4) SSDI analysis A (DK5)
```matlab
addpath('code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi');
run_ssdi_sleep_stage_pipeline_dk('mode','stage_ssdi','atlas','dk5','nrunsp',100,'nitersp',10000,'niterso',10000);
```

## 5) SSDI analysis B (DK68)
```matlab
% Requires GC subject-measures from DK68 branch first:
run_gc_sleep_stage_pipeline_dk5('mode','stage_mvgc2','atlas','dk68','momax',30,'mosel','LRT','fs_original',256,'fs_target',256,'fres',256);
run_ssdi_sleep_stage_pipeline_dk('mode','stage_ssdi','atlas','dk68','nrunsp',100,'nitersp',10000,'niterso',10000);
```

## 6) SSDI summaries/figures
```bash
python3 code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/summarize_ssdi_results_dk.py \
  --spectral-dir results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz \
  --out-dir results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz

python3 code/standalone_sleep_dk5_dk68_lrt30_current_epoching_2s256/ssdi/plot_ssdi_stage_figures_dk.py \
  --long-csv results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/ssdi_subject_long.csv \
  --spectral-dir results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz \
  --out-dir results/ssdi/stage_ssdi_dk5_lrt30_current_epoching_2s_256hz/figures
```
