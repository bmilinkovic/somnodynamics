
# Healthy Sleep Project - Code Repository

This repository contains MATLAB code and scripts for the **Healthy Sleep Project**, which focuses on the analysis of brain dynamics during different sleep stages using Granger Causality (GC), Dynamical Independence, and other computational techniques. The analyses are performed across time and spectral domains for different regions of interest (ROIs) using coarse-grained participant data.

## Dual-Atlas Note

The repository now includes a side-by-side Yeo7 pipeline (without replacing legacy Desikan scripts). See:

- `code/YEO7_PIPELINE.md`

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Data Preparation](#data-preparation)
- [Running the Code](#running-the-code)
  - [Time-Domain Optimisation](#time-domain-optimisation)
  - [Spectral-Domain Optimisation](#spectral-domain-optimisation)
  - [Running on SPARTAN](#running-on-spartan)
- [Results](#results)
- [Contributors](#contributors)

## Overview

This project uses the Desikan-Killiany atlas with 68 regions of interest (ROIs) to perform inter-ROI pairwise conditional analysis between different sleep stages. The analysis is conducted both in the time domain and spectral domain, and key metrics such as node weights, Granger causality (GC) statistics, and beta statistics are computed and saved.

### Key Analyses:

1. **Time-Domain Optimisation**: Optimises for macroscopic variables using Dynamical Independence across different sleep stages.
2. **Spectral-Domain Optimisation**: Optimises for macroscopic variables in specific frequency bands (Delta, Theta, Alpha, Beta, Gamma) using spectral Granger causality analysis.

## Repository Structure

```
healthy-sleep-project/
│
├── code/
│   ├── time_domain_optimisation.m           # Time-domain optimisation script
│   ├── spectral_domain_optimisation.m       # Spectral-domain optimisation script
│   └── results/
│       ├── pwcgc/                           # Directory for participant data
│       ├── ssdi_results/                    # Directory for storing results
│       └── jobarrays/                       # Directory for job-specific scratch data
├── README.md                                # Project documentation
└── slurm_scripts/
    ├── bm_time_domain_optimisation.sh       # SLURM script for time-domain optimisation
    └── bm_spectral_domain_optimisation.sh   # SLURM script for spectral-domain optimisation
```

## Requirements

### Software:
- MATLAB (including Parallel Computing Toolbox)
- SPARTAN access at The University of Melbourne
- SLURM for job submission on SPARTAN

### MATLAB Toolboxes:
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox

### Data Files:
Ensure that the following directories and data files are in place:
- `/code/results/pwcgc/`: Contains the `.mat` files for each participant's brain data (e.g., `pwcgc_params_*.mat`).
- `/code/results/ssdi_results/`: The directory where the optimisation results are saved.

## Data Preparation

Before running any scripts, make sure your participant data files are stored in the `pwcgc/` directory. The `.mat` files should contain variables such as `V`, `CAK`, and `Vss`, which are used in the analysis.

## Running the Code

### Time-Domain Optimisation

To run the time-domain optimisation, use the script `time_domain_optimisation.m`. This script performs optimisation for macroscopic variables across different sleep stages for each participant and saves the results.

**Command:**
```matlab
run_time_domain_optimisation('participant_id');
```
- Replace `'participant_id'` with the participant's ID (e.g., `'EPCTL01_N1'`).

**Results:**
- The results are stored in `/code/results/ssdi_results/time/` as `.mat` files.

### Spectral-Domain Optimisation

To run the spectral-domain optimisation, use the script `spectral_domain_optimisation.m`. This script performs spectral-domain analysis in specific frequency bands (Delta, Theta, Alpha, Beta, Gamma).

**Command:**
```matlab
run_spectral_domain_optimisation('participant_id');
```
- Replace `'participant_id'` with the participant's ID.

**Results:**
- The results are stored in `/code/results/ssdi_results/spectral/` as `.mat` files.

### Running on SPARTAN

The project is set up to run large-scale analyses on the SPARTAN cluster at The University of Melbourne. SLURM scripts are provided to facilitate job submission.

#### SLURM Scripts:
- `bm_time_domain_optimisation.sh`: SLURM script to submit time-domain jobs.
- `bm_spectral_domain_optimisation.sh`: SLURM script to submit spectral-domain jobs.

**To run on SPARTAN:**

1. **Copy the SLURM scripts and data:**
   Ensure your participant data is copied from your local machine to SPARTAN.
   ```bash
   rsync -avz /Users/borjan/code/healthy-sleep-project/code/results/pwcgc_lowres/ bmilinkovic@spartan.unimelb.edu.au:/data/gpfs/projects/punim1761/healthy-sleep-project/code/results/pwcgc
   ```

2. **Submit jobs:**
   ```bash
   sbatch bm_time_domain_optimisation.sh
   sbatch bm_spectral_domain_optimisation.sh
   ```

The SLURM scripts are designed to run analyses in parallel using 48 CPUs per job.

## Results

Results for both time-domain and spectral-domain optimisations are saved in the following directories:
- `/code/results/ssdi_results/time/`: Time-domain results.
- `/code/results/ssdi_results/spectral/`: Spectral-domain results.

Each `.mat` file contains:
- `Lopto_m`: Optimised subspace matrices (stored in cell arrays due to varying dimensions).
- `dopto_m`, `gopto_m`: Optimisation metrics.
- `node_weights_m`: Node weights across regions.
- `beta_m`: Beta statistics.
- `beta_pvals_m`: P-values for statistical significance.
- `beta_sigs_m`: Statistical significance flags.

## Contributors

- **Borjan Milinkovic** - Study design, data analysis, optimisation code development.

For more details or questions, contact **bmilinkovic@student.unimelb.edu.au**.
