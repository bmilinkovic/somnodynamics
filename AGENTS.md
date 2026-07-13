# AGENTS.md

## Project purpose

This repository is being migrated into a clean, reproducible analysis toolbox for the healthy sleep EEG dataset.

The goal is to rebuild the healthy sleep analysis pipeline so that it follows the same architecture, methodology, and output structure as the DMT project.

Reference DMT project:

/Users/borjan/code/dmt_project

Key DMT reference directories:

SSDI / Dynamical Dependence pipeline:
/Users/borjan/code/dmt_project/code/ssdi_anatomical7

Granger causality pipeline:
/Users/borjan/code/dmt_project/code/paper_style_mvgc_anatomical7

SSDI figure outputs:
/Users/borjan/code/dmt_project/results/ssdi_anatomical7_canonical/figures

GC figure outputs:
/Users/borjan/code/dmt_project/results/paper_style_mvgc_yeo7_rest_vs_5meodmt_hierarchical_stats/time_resolved_diagnostics/figures_python

---

## Raw dataset location

The healthy sleep dataset is located outside this repository:

/Volumes/ex_data/data_somnoemergence/raw

Rules:

- Raw data must remain in this location
- Do not move, rename, duplicate, or overwrite raw data
- All pipelines must read raw data from this path
- All derived outputs must be saved inside this repository

---

## Sleep stage naming

Internal dataset labels:

W_
N1
N2
N3
R_

Standardized labels used in all downstream outputs:

W_ → Wake  
N1 → N1  
N2 → N2  
N3 → N3  
R_ → REM  

These standardized labels must be used in:

- figures
- statistics tables
- reports
- figure annotations
- saved outputs

---

## Preprocessing requirements

The healthy sleep project already contains a preprocessing pipeline.

This preprocessing generates:

preprocessed_sensor_space

This preprocessing must:

- be located in the existing code
- be audited before any changes
- be preserved and integrated into the new toolbox

The preprocessing pipeline should not be rewritten unless necessary.

---

## Sensor-space analysis branch

The toolbox must include sensor-level spectral analyses.

These operate on:

preprocessed_sensor_space

Required analyses:

- absolute power
- relative power

Statistics:

- repeated-measures within-subject across sleep stages

Outputs must include:

- statistical tables
- publication-ready figures
- saved intermediate outputs

---

## Source-space analysis branch

The toolbox must reproduce the DMT-style source-space analysis pipeline.

Required components:

1. source reconstruction
2. canonical anatomical coarse-graining
3. Granger causality analysis
4. Dynamical Dependence / SSDI analysis
5. statistics
6. figure generation
7. GC–SSDI cross-analysis

---

## Canonical anatomical mapping

The canonical macro-regions used in the DMT project must be extracted directly from the DMT code.

Do not guess.

The mapping must be:

- reproduced exactly
- documented explicitly
- stored as a configuration or mapping file

---

## Repository structure expectations

The refactored repository should resemble:

healthy-sleep-project/

AGENTS.md  
README.md  

docs/

code/
    legacy/
    preprocessing/
    sensor_power/
    src_recon/
    coarse_graining/
    gc/
    ssdi/
    stats/
    figures/
    cross_analysis/
    run_all/

data/
    interim/
    processed/

results/
    preprocessing/
    sensor_power/
    source_recon/
    coarse_grained/
    gc/
    ssdi/
    stats/
    figures/
    cross_analysis/

tests/

Notes:

- legacy code must be preserved
- raw dataset must not be copied into the repository

---

## Coding guidelines

Prefer:

- modular code
- small functions
- clear docstrings
- explicit configuration files

Avoid:

- monolithic scripts
- hidden global state
- hard-coded paths

---

## Documentation requirements

The repository must include:

- README explaining the full pipeline
- migration plan describing repository changes
- canonical region mapping documentation
- preprocessing documentation
- documentation for each analysis stage

---

## Statistical design

The sleep dataset uses a repeated-measures within-subject design across sleep stages.

All statistical analyses must respect this design.

---

## Output requirements

All analyses must save:

- machine-readable result tables
- publication-ready figures
- human-readable summary reports

Figures should be saved in:

- PNG
- vector formats where appropriate (SVG or PDF)

---

## Workflow required for the agent

Before modifying the repository the agent must:

1. Audit the healthy-sleep-project repository
2. Audit the DMT project pipelines
3. Write a migration plan

Then perform the following steps:

4. Restructure the repository
5. Integrate preprocessing
6. Implement sensor power analyses
7. Implement source reconstruction
8. Implement canonical coarse-graining
9. Implement the GC pipeline
10. Implement the SSDI pipeline
11. Implement cross-analysis
12. Implement statistics
13. Implement figure generation
14. Add documentation and validation

---

## Deliverables

The repository must ultimately contain:

- a clean modular analysis toolbox
- preserved legacy code
- migration plan markdown file
- implementation summary
- canonical region mapping documentation
- sensor power analyses
- GC analyses
- SSDI analyses
- cross-analysis results
- publication-ready figures
- statistical result tables

---

## Definition of done

The migration is complete only when:

- preprocessing pipeline has been located and integrated
- sensor-space analyses are implemented
- source-space analyses are implemented
- GC pipeline matches the DMT methodology
- SSDI pipeline matches the DMT methodology
- cross-analysis between GC and SSDI is implemented
- figures and statistics are reproducible
- documentation explains how to run the entire pipeline
