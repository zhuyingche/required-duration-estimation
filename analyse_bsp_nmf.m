%% ======================================================================
% Required Duration (reqDur) Estimation Framework
%
% OVERVIEW
% The reqDur framework quantifies the minimum recording duration required
% to reliably capture stable HFO spatial organization. It consists of three
% sequential stages:
%
%   Stage 1: NMF-based extraction of noise-reduced HFO spatial and temporal
%            patterns from windowed HFO-rate matrices.
%            -> analyse_nmfByEpochByState_mainFunction.m
%
%   Stage 2: Similarity quantification between partial-recording patterns
%            (up to epoch n) and the full-recording reference using
%            best-match cosine similarity with padding to accommodate
%            unequal recording lengths.
%            -> bestMatchCosine_mainFunction.m
%
%   Stage 3: Stability-point detection to identify the earliest recording
%            time at which similarity reaches and maintains a high plateau
%            (defined as reqDur).
%            -> analyse_dtw_nmf_findStabPoint.m
%
% ----------------------------------------------------------------------
% CURRENT SCRIPT: Stage 2
%
% PURPOSE
% For each patient, this script loads Stage 1 NMF outputs and computes
% similarity trajectories across increasing recording duration relative
% to the full recording. These trajectories are subsequently used in
% Stage 3 to determine reqDur.
%
% AUTHOR
% Zhuying Chen, 20 Aug 2024
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