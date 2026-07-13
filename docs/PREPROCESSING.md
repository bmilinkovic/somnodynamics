# SomnoDynamics Preprocessing

## Legacy audit outcome

The existing preprocessing pipeline was located and preserved at:
- `code/legacy/source/01_preprocessing.py` (sensor-space preprocessing)
- `code/legacy/source/load-edf-data.py` (source-oriented preprocessing)

These scripts were audited for:
- filtering and resampling operations
- stage segmentation logic
- channel exclusion and referencing
- file output behavior

## Integrated preprocessing entrypoint

SomnoDynamics uses:
- `code/preprocessing/preprocess_sensor_space.py`

Default raw input:
- `/Volumes/ex_data/data_somnoemergence/raw`

Default output:
- `results/preprocessing/preprocessed_sensor_space/`

Outputs per subject/stage:
- `EPCTLxx_<stage>_epochs.fif`

where `<stage>` uses internal standardized labels (`W_`, `N1`, `N2`, `N3`, `R_`).

## Run

```bash
python3 code/preprocessing/preprocess_sensor_space.py
```

Optional overrides:

```bash
python3 code/preprocessing/preprocess_sensor_space.py \
  --raw-root /Volumes/ex_data/data_somnoemergence/raw \
  --output-dir results/preprocessing/preprocessed_sensor_space \
  --subject EPCTL01 --subject EPCTL02
```

## Reproducibility notes

- Raw data is read in-place and never moved.
- Derived outputs are written only inside the SomnoDynamics repository.
- Stage labels are canonicalized through `code/run_all/somnodynamics_config.toml`.
