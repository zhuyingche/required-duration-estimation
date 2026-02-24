%% ======================================================================
% Code accompanying the manuscript:
% "The influence of recording duration and vigilance state on 
% high-frequency oscillation characterization in epilepsy"
% Chen Z, et al. (2026)
%
% This repository implements a three-stage analytical framework to estimate
% the required recording duration (reqDur) necessary to capture stable
% high-frequency oscillation (HFO) spatial distributions across vigilance
% states.
%
% Please cite the final published version when available.
% ----------------------------------------------------------------------
% Required Duration (reqDur) Estimation Pipeline
%
% The reqDur framework comprises three sequential stages:
%
% Stage 1: Extraction of noise-reduced HFO spatial distributions using
%          non-negative matrix factorization (NMF)
%          -> analyse_nmfByEpochByState.m
%
% Stage 2: Quantification of similarity between distributions derived from
%          truncated recordings and the full recording using dynamic time
%          warping with best-match similarity padding (BSP)
%          -> analyse_bsp_nmf.m
%
% Stage 3: Identification of the earliest time point at which similarity
%          reaches and maintains a high, stable plateau (defined as reqDur)
%          -> analyse_bsp_nmf_findStabPoint.m
%
% ----------------------------------------------------------------------
% CURRENT SCRIPT: Stage 2
%
% PURPOSE
% For each patient, this script loads Stage 1 NMF outputs and computes
% similarity trajectories across increasing recording duration relative
% to the full recording. These trajectories are subsequently used in
% Stage 3 to determine reqDur.
%% ======================================================================

%% ------------------------ Patient Selection ----------------------------

patientIds1 = [18:85,93];

% Exclusions (documented for reproducibility)
patientIds2 = [57,58,61];              % files not generated / metadata issues
patientIds3 = [33,77,50,54,55,93];     % anomalous/insufficient data, missing sleep scoring

patientIds = setdiff(patientIds1, patientIds2);
patientIds = setdiff(patientIds, patientIds3);

%% ------------------------ Paths ----------------------------------------

filePathData   = 'C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data\';
filePathResult = 'C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data analysis\result\';
filePathFig    = 'C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data analysis\figure';

% HPC paths (uncomment if needed)
% filePathData   = '/data/gpfs/projects/punim1181/Postdoc/data/';
% filePathResult = '/data/gpfs/projects/punim1181/Postdoc/data analysis/result/';
% filePathFig    = '/data/gpfs/projects/punim1181/Postdoc/data analysis/figure/';

%% ------------------------ Core Parameters ------------------------------

hfoversion = "qHFO_v2.3_running";
winMin     = 5;  % HFO rate window (minutes)

%% ------------------------ Stage 1 Consistency Options ------------------
% These parameters must match those used in Stage 1.

RmAbnEpOptions = {'','RmHighEpoch','RmHighLowEpoch','RmHighLowEpoch_bold','RmLowHighEpoch'};
RmAbnEpOption  = RmAbnEpOptions{5};

RmFewValidSampleEpOptions = {'','RmFewValidSampleEpoch'};
RmFewValidSampleEpOption  = RmFewValidSampleEpOptions{2};

normalizeMethods = {'','zscore','minMax','L1norm'};
normalizeMethod  = normalizeMethods{1};

InterictalDataOnlyOptions = {'','InterictalDataOnly'};
InterictalDataOnly = InterictalDataOnlyOptions{2};

%% ------------------------ Stage 2 Execution ----------------------------

for iPt = 1 % use 1 for testing; change to 1:length(patientIds) for batch
    patientId = patientIds(iPt);

    cosine_dist_nmf = bestMatchCosine_mainFunction( ...
        patientId, filePathResult, filePathFig, ...
        hfoversion, winMin, ...
        RmAbnEpOption, RmFewValidSampleEpOption, ...
        normalizeMethod, InterictalDataOnly);

end