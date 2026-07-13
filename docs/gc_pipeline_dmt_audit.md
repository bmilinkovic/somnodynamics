# SomnoDynamics GC Pipeline Audit and Port (DMT -> Sleep Stages)

Date: 2026-03-17

## Scope
This document audits the DMT GC implementation and documents the SomnoDynamics port for stage-based sleep analysis.

Reference implementation audited:
- `/Users/borjan/code/dmt_project/code/paper_style_mvgc_anatomical7`

SomnoDynamics target:
- stage-based within-participant analysis across `W_`, `N1`, `N2`, `N3`, `R_`
- no time-windowed-within-state contrast logic

## DMT files audited
Primary DMT GC files inspected:
- `run_paper_style_mvgc_anatomical7_pipeline.m`
- `ps_config_anatomical7.m`
- `ps01_compute_subject_measures_anatomical7.m`
- `ps02_group_stats_paired_anatomical7.m`
- `psw_config_anatomical7.m`
- `psw01_compute_subject_measures_timewindows_anatomical7.m`
- `psw02_group_stats_timewindows_paired_anatomical7.m`
- `README.md`

## Upstream dependencies identified
DMT dependencies used by the GC pipeline:
- MVGC2 toolbox (`startup.m`, `tsdata_to_varmo`, `tsdata_to_var`, `moselect`, GC/MI transforms)
- Source-export step writing anatomical7 macro time-series MAT files
- MATLAB toolboxes for stats/signal/control

Confirmed MVGC2 load pattern in DMT:
- `run_paper_style_mvgc_anatomical7_pipeline.m` initializes `MVGC2` from project-local path and checks `tsdata_to_varmo`.

## DMT execution flow (audited)
1. Source export to anatomical7 MAT files
2. Subject-level measure extraction (`ps01`)
   - VAR order selection
   - VAR fit
   - directed GC (time + spectral)
   - undirected MI-family measures (time + spectral)
3. Group stats (`ps02`)
   - paired Rest vs DMT signed-rank tests
   - BH/BY correction by family
4. Optional time-window branch (`psw*`)
   - paired Rest vs W1/W2/W3/W4 windows

## Key assumptions in original DMT code
- Two-condition paired design (rest vs dmt)
- Time-window branch still fundamentally paired against rest
- Condition-centric folder structure (`rest/`, `dmt/`)
- Group stats are paired two-condition tests, not multi-level repeated measures

## Required adaptation for SomnoDynamics
SomnoDynamics scientific design differs:
- Analysis compares sleep stages within participants (5 conditions) rather than windows within one state.

Critical logic changes implemented:
- Replaced rest/dmt condition assumptions with stage IDs: `W_`, `N1`, `N2`, `N3`, `R_`
- Replaced paired two-condition group stats with repeated-measures multi-stage stats:
  - Friedman omnibus per feature
  - pairwise signed-rank between stage pairs
  - BH/BY FDR within correction families
- Updated naming/metadata/output organization to stage-aware labels and standardized display labels (`Wake`, `N1`, `N2`, `N3`, `REM`).

## SomnoDynamics implementation added
New module:
- `code/gc/stage_mvgc2_anatomical7/`

Files:
- `sd_gc_config_anatomical7.m`
- `run_stage_mvgc2_anatomical7_pipeline.m`
- `sd_gc01_compute_subject_measures_by_stage.m`
- `sd_gc02_group_stats_repeated_measures.m`

Entrypoint wiring:
- `code/gc/run_gc_sleep_stage_pipeline.m`
  - default mode: new stage-based MVGC2 pipeline
  - optional `mode='legacy'` compatibility wrapper

Pipeline runner fix:
- `code/run_all/run_pipeline.py`
  - GC stage now calls MATLAB with `addpath('code/gc'); run_gc_sleep_stage_pipeline();`

## MVGC2 local vendoring
MVGC2 is now package-local:
- `code/external/MVGC2`

Path setup updated:
- `code/run_all/setup_matlab_toolboxes.m`
  - first MVGC2 candidate now `code/external/MVGC2`

This removes hard dependency on an external machine path for MVGC2.

## Input/output contract (SomnoDynamics GC)
Input:
- source MAT files from:
  - `results/source_recon/anatomical7/EPCTL*/mat/*_source_time_series.mat`
- expected fields: `source_ts` and `labels`
- expected shape: `(7, 256, n_epochs)`

Subject-level outputs:
- `results/gc/stage_mvgc2_anatomical7/data_intermediate/subject_measures/<stage>/<subject>_subject_measures.mat`

Group-level outputs:
- `results/gc/stage_mvgc2_anatomical7/stats/repeated_measures_stats.mat`
- `results/gc/stage_mvgc2_anatomical7/stats/omnibus_long.csv`
- `results/gc/stage_mvgc2_anatomical7/stats/pairwise_long.csv`
- `results/gc/stage_mvgc2_anatomical7/reports/repeated_measures_summary.txt`
- logs in `results/gc/stage_mvgc2_anatomical7/logs/`

## Run evidence (current repository)
Executed full new pipeline run:
- command: `run_gc_sleep_stage_pipeline('mode','stage_mvgc2')`
- input files detected: 140 (28 per stage)
- stage1 subject-measure status: 140/140 success
- repeated-measures tables generated:
  - omnibus rows: 735
  - pairwise rows: 7350

## Ambiguities / risks / notes
- MVGC2 C MEX binaries (`mvfilter`, `findin`) were not built in this environment, so MATLAB used slower scripted fallbacks.
- The current group-level repeated-measures framework is feature-wise nonparametric (`friedman` + `signrank`) with FDR families; if you want a strict one-model repeated-measures formulation (e.g., mixed-effects/RMANOVA variants), that should be added as a separate stats mode.
- Legacy GC wrappers remain available (`mode='legacy'`) for compatibility comparisons but are not the default scientific path.

## Recommendation
Use the new default stage-based MVGC2 pipeline in SomnoDynamics for GC on sleep data.
It preserves DMT measure computation philosophy while correctly adapting inferential logic from two-condition/window comparisons to within-participant across-stage repeated measures.
