# SomnoDynamics Audit: Sleep Preprocessing, Source Reconstruction, and GC Preparation

Date: 2026-03-16  
Repository audited: `/Users/borjan/code/somnodynamics`

## Scope and constraints
This audit answers Questions A–G from the request, using code-level evidence and direct inspection of raw and preprocessed files.

- Raw data inspected in-place (no copying):
  - `/Volumes/ex_data/data_somnoemergence/raw`
- Legacy outputs inspected in-place:
  - `/Volumes/dataSets/healthy-sleep-project-datasets/data/preprocessed_source_reconstructed`
  - `/Volumes/dataSets/healthy-sleep-project-datasets/data/preprocessed_sensor_space`
- Current outputs inspected:
  - `results/preprocessing/preprocessed_sensor_space`
  - `results/source_recon/anatomical46`
  - `results/source_recon/anatomical7`

## A) How the raw sleep data actually looks

### A.1 Storage and stage representation
- Sleep staging used by this project is in per-subject TXT files (`<subject>.txt`), with columns:
  - `label`, `start_time`, `duration`
- Internal labels present: `W`, `N1`, `N2`, `N3`, `R`, plus non-stage `L` rows.
- Across all 29 subjects, TXT `duration` values are only `30.0` seconds.
- Therefore the stage files encode 30-second scored rows (standard sleep-scoring granularity).

### A.2 EDF metadata and annotations
For 2 example subjects directly inspected (EPCTL01, EPCTL03):
- EDF sampling rate: 1000 Hz
- EDF channels: 93
- EDF internal annotations exist, but are sparse/non-stage descriptors (examples: `Segment: REC START...`, `A1+A2 OFF`, posture labels like `supine`, `left side`).
- Conclusion: stage labels for analysis come from TXT scoring files, not EDF annotations.

### A.3 Stage durations and expected non-overlapping 2-second windows
Computed from raw TXT durations:

| Subject | Stage | Total duration (s) | Expected 2s windows |
|---|---|---:|---:|
| EPCTL01 | W_ | 1380 | 690 |
| EPCTL01 | N1 | 2220 | 1110 |
| EPCTL01 | N2 | 14130 | 7065 |
| EPCTL01 | N3 | 5550 | 2775 |
| EPCTL01 | R_ | 5280 | 2640 |
| EPCTL03 | W_ | 1680 | 840 |
| EPCTL03 | N1 | 1980 | 990 |
| EPCTL03 | N2 | 9150 | 4575 |
| EPCTL03 | N3 | 6210 | 3105 |
| EPCTL03 | R_ | 3810 | 1905 |

---

## B) What the legacy pipeline did

### B.1 Legacy preprocessing script (sensor-space)
`code/legacy/source/01_preprocessing.py` creates 10-second events by expanding each annotation row by duration:
- event tiling logic: lines 162–173
- 10-second epochs: lines 177–179

This explains old `preprocessed_sensor_space` counts (e.g., 30s row -> 3 ten-second epochs).

### B.2 Legacy source-reconstruction scripts
Two legacy source scripts create source epochs from row onsets only:
- `code/legacy/source/load-edf-data.py`
  - one event per row onset: lines 145–149
  - 2-second epoch extraction: lines 152–153
- `code/legacy/source/source_reconstruction_yeo7.py`
  - one event per row onset: lines 266–270
  - epoch creation with `tmin=0`, `tmax=epoch_seconds`: lines 275–283

Interpretation: source reconstruction used one short snippet per scored row, not full tiling of each scored-row duration.

### B.3 Direct validation against saved old outputs
For EPCTL01 and EPCTL03, old source epoch count equals TXT row count per stage, and current new output equals full expected 2s windows:

| Subject | Stage | Raw duration (s) | Expected full 2s | TXT row count | Old epochs | Current epochs |
|---|---|---:|---:|---:|---:|---:|
| EPCTL01 | W_ | 1380 | 690 | 46 | 46 | 690 |
| EPCTL01 | N1 | 2220 | 1110 | 74 | 74 | 1110 |
| EPCTL01 | N2 | 14130 | 7065 | 471 | 471 | 7065 |
| EPCTL01 | N3 | 5550 | 2775 | 185 | 185 | 2775 |
| EPCTL01 | R_ | 5280 | 2640 | 176 | 176 | 2640 |
| EPCTL03 | W_ | 1680 | 840 | 56 | 56 | 840 |
| EPCTL03 | N1 | 1980 | 990 | 66 | 66 | 990 |
| EPCTL03 | N2 | 9150 | 4575 | 305 | 305 | 4575 |
| EPCTL03 | N3 | 6210 | 3105 | 207 | 207 | 3105 |
| EPCTL03 | R_ | 3810 | 1905 | 127 | 127 | 1905 |

Why they differ: with 30s scoring rows, full non-overlap 2s tiling gives 15 windows per row. New:Old = 15x for EPCTL03 across all stages.

### B.4 Current preprocessing/source/GC-prep behavior (code-level)
Current sensor preprocessing (`code/preprocessing/preprocess_sensor_space.py`) applies:
- notch filter 50/100 Hz: line 150
- band-pass 1–45 Hz: line 151
- resample to 256 Hz: line 152
- channel cleaning + average reference projection: lines 164–175
- duration-aware event tiling from TXT rows: lines 119–138
- stage FIF output (`*_epochs.fif`): lines 183–201

Current source reconstruction (`code/src_recon/run_source_reconstruction_anatomical46.py`) then:
- reads preprocessed sensor epochs (does not reread raw EDF): lines 4–6, 57–61, 322–334
- decimates/segments to 2-second source windows at 128 Hz: lines 160–184
- computes LCMV with `reg=0.05`, `pick_ori='max-power'`, `weight_norm='nai'`, `rank='info'`: lines 300–309
- writes anatomical46 then anatomical7: lines 343–376

Current GC preparation expects source arrays shaped `(regions, time, epochs)` and treats epoch/trial axis as multitrial realizations:
- `code/legacy/run_mvgc_atlas.m` lines 82–116

---

## C) Was legacy behavior correct for the intended GC goal?
Intended GC goal: use full stage signal, then produce non-overlapping 2-second source trials as multitrial realizations.

Finding:
- Legacy source behavior does **not** do this. It samples one 2-second snippet per scored row onset.
- It therefore under-samples each stage by ~15x relative to full 30s-row coverage.

So for the stated GC/VAR/state-space objective, legacy source epoch definition is methodologically mismatched.

---

## D) Correct order of operations for source reconstruction and GC preparation

### Strategy comparison
1. **Strategy 1 (2s sensor epochs first, then LCMV from those short epochs):**
   - feasible, but covariance and filter estimation are based on short windows and may be less stable.
2. **Strategy 2 (stage-continuous source reconstruction first, then 2s source segmentation):**
   - methodologically preferable for this use-case.
   - estimate beamformer from more complete stage data, then segment source outputs into GC-ready trials.

### Evidence from tooling/docs
- MNE LCMV is covariance-driven (`make_lcmv`) and can be applied to `Raw` or `Epochs` outputs (`apply_lcmv_raw` / `apply_lcmv_epochs`), supporting “reconstruct first, segment later” workflows.
- MNE tutorial guidance indicates covariance quality depends on enough samples relative to channels.

### Practical recommendation
Use Strategy 2 as primary:
- preprocess once,
- build stage-specific continuous data objects,
- estimate filters/covariance with adequate samples,
- reconstruct source time series,
- then tile to non-overlap 2s source trials for GC.

---

## E) Common filter vs stage-specific filters

### What authoritative docs state
- MNE condition-comparison beamformer tutorial recommends a **common filter** for condition contrasts (to avoid confounding differences caused by separate filters).
- FieldTrip beamformer tutorials similarly recommend computing one common spatial filter from concatenated/all-condition data and applying it per condition.

### Interpretation for sleep stages
- If the main inferential target is **between-stage comparability**, common filter is usually preferable.
- If the target is **best stage-specific source capture** and covariance structure likely differs strongly across stages, stage-specific filters can be informative but are less directly comparable.

### Principled compromise
- Primary analysis: common filter across all included sleep stages per subject.
- Sensitivity analysis: stage-specific filters; report whether conclusions (directionality/significant contrasts) are robust.

---

## F) Recommended SomnoDynamics pipeline (precise)

1. Read stage annotations from raw TXT (`label`, `start_time`, `duration`) and canonicalize labels to `W_`, `N1`, `N2`, `N3`, `R_`.
2. Build stage-specific continuous sensor data streams from raw (or equivalently all stage samples after continuous preprocessing), not one-snippet-per-row sampling.
3. Source reconstruction:
   - use anatomical46 first,
   - then coarse-grain to canonical anatomical7 mapping from DMT.
4. Filter strategy:
   - default common filter across stages per subject,
   - optional stage-specific sensitivity run.
5. Segmentation timing:
   - segment into non-overlapping 2-second trials **after** source reconstruction.
6. Output shapes:
   - anatomical46: `(46, 256, n_trials_stage)` at 128 Hz and 2 s
   - anatomical7: `(7, 256, n_trials_stage)`
   - where `n_trials_stage = floor(total_stage_duration_seconds / 2)`.
7. Non-overlap enforcement:
   - deterministic tiling by fixed 2s window length, no overlap, no random sampling.
8. Guardrails against hidden errors:
   - save per-stage audit metadata (`duration_s`, expected window count, actual trial count),
   - fail if actual != expected,
   - keep one saved output per subject-stage-atlas target path.

Current code status relevant to this:
- Current source script already reads preprocessed sensor epochs (no second preprocessing pass):
  - `code/src_recon/run_source_reconstruction_anatomical46.py` lines 4–6, 57–61, 322–334
- Current source script transforms to 2-second, 128 Hz windows and writes anatomical46 then anatomical7:
  - transformation: lines 160–184
  - fine/coarse saves: lines 343–376

---

## G) Deliverables created in this audit

1. Audit report:
- `docs/source_recon_gc_audit.md`

2. Validation scripts:
- `code/validation/inspect_raw_sleep_annotations.py`
- `code/validation/compare_old_vs_new_source_epoch_counts.py`

3. Concise behavior summary:

### Legacy behavior
- Sensor preprocessing: 10-second tiling from 30-second rows.
- Source reconstruction: one onset snippet per row (2-second epoch), under-sampling for full-stage GC purpose.

### Current behavior
- Preprocessed sensor-space matches old 10-second counts for checked subject.
- Source outputs are 2-second windows covering full stage durations (for checked subjects), with expected trial counts.

### Recommended behavior
- Keep current full-coverage 2-second source trial strategy.
- Formalize stage-continuous reconstruction logic and common-filter primary analysis with stage-specific sensitivity checks.
- Enforce count checks in validation for every run.

---

## External references (authoritative)
- MNE LCMV beamformer tutorial (condition comparison and common filter context):  
  https://mne.tools/stable/auto_tutorials/inverse/50_beamformer_lcmv.html
- MNE API (`make_lcmv`, covariance-driven filter construction):  
  https://mne.tools/stable/generated/mne.beamformer.make_lcmv.html
- FieldTrip beamformer tutorial (common filter for condition contrasts):  
  https://www.fieldtriptoolbox.org/tutorial/source/beamformer/

## Uncertainties and limits
- This is an audit pass, not a full end-to-end rerun of all downstream GC/SSDI outputs.
- Recommendations on common vs stage-specific filtering are methodologically grounded, but should be finalized against the exact inferential question (contrast-focused vs stage-optimized reconstruction).
- Current source output filenames for `W_`/`R_` include a double underscore artifact (e.g., `..._W__source_time_series.mat`), which is a naming hygiene issue but not a count/shape error.
