
## Required Duration Estimation Framework for HFO Characterisation

### Overview

This repository accompanies the under-review manuscript:

Chen Z*, Yu W*, et al. (2026)
The influence of recording duration and vigilance state on high-frequency oscillation characterization in epilepsy

This repository implements a three-stage analytical framework to estimate the required recording duration necessary to capture stable high-frequency oscillation (HFO) spatial distributions across vigilance states.

The framework quantifies when truncated intracranial EEG recordings reliably approximate full-recording HFO spatial patterns.

Please cite the final published version when available.

Scientific Rationale

Reliable characterization of HFO spatial distributions is essential for epilepsy biomarker development and surgical outcome prediction. However, the minimum recording duration required to obtain stable spatial estimates remains uncertain.

The reqDur framework systematically evaluates temporal stability of HFO spatial structure using matrix factorisation and time-alignment similarity metrics.

### Analytical Pipeline

The reqDur framework comprises three sequential stages:

#### Stage 1 – Extraction of Noise-Reduced HFO Spatial Distributions

Method: Non-negative matrix factorization (NMF)
Purpose: Derive stable, low-dimensional spatial representations of HFO rate distributions across vigilance states

Main script:

analyse_nmfByEpochByState.m
#### Stage 2 – Similarity Quantification Using Best-Match Similarity Padding (BSP)

Method: Best-match Similarity Padding (BSP)

Purpose:
Quantify similarity between HFO spatial distributions derived from truncated recordings and the full recording.

This stage assesses how closely shorter recordings reproduce the spatial structure of the complete dataset.

Main script:

analyse_bsp_nmf.m
#### Stage 3 – Stable Plateau Detection (reqDur Identification)

Method: Stability plateau detection

Purpose:
Identify the earliest time point at which similarity:

Reaches a high level

Remains consistently stable

This time point is defined as the required duration (reqDur).

Main script:

analyse_bsp_nmf_findStabPoint.m

Repository Structure
.
├── analyse_nmfByEpochByState.m
├── analyse_bsp_nmf.m
├── analyse_bsp_nmf_findStabPoint.m
├── functions/                  % Supporting utility functions
├── data/               % demonstration data
├── results/               % results
├── figure/               % figure
└── README.md
Requirements

MATLAB (tested on R2023)

Statistics and Machine Learning Toolbox

Signal Processing Toolbox


### Citation

If you use this framework in your research, please cite:

Chen Z*, Yu W*, et al. (2026).
The influence of recording duration and vigilance state on high-frequency oscillation characterization in epilepsy.

DOI: (to be updated upon publication)

### Contact

For methodological questions or collaboration inquiries, please contact:

Dr. Zhuying Chen
email: zhuychen@unimelb.edu.au

Dr. William Stacey
email: wstacey@umich.edu
