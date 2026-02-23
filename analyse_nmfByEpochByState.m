%% ======================================================================
% Code accompanying the manuscript:
% "The influence of recording duration and vigilance state on 
% high-frequency oscillation characterization in epilepsy"
% Chen Z*, Yu W*, et al. (2026)
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
%          -> analyse_nmfByEpochByState_mainFunction.m
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
% The present script executes Stage 1. Outputs are stored for subsequent
% similarity analysis (Stage 2) and stability detection (Stage 3).
%% ======================================================================

%% Add library path
addpath('C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data analysis\code\lib');
% addpath('/data/gpfs/projects/punim1181/Postdoc/data analysis/code/lib'); % HPC

%% Patient selection
patientIds = [26];   % modify for batch processing

%% Paths and core parameters
filePathData   = 'C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data\';
filePathResult = 'C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data analysis\result\';

hfoversion = "qHFO_v2.3_running";
winMin     = 5;      % HFO rate window (minutes)
normalizeMethod = ''; % define if required

%% Epoch quality control
RmAbnEpOptions = {'','RmHighLowEpoch','RmLowHighEpoch'};
RmAbnEpOption  = RmAbnEpOptions{2};

RmFewValidSampleEpOptions = {'','RmFewValidSampleEpoch'}; % remove epochs <90% valid samples
RmFewValidSampleEpOption  = RmFewValidSampleEpOptions{2};

InterictalDataOnlyOptions = {'','InterictalDataOnly'};
InterictalDataOnly = InterictalDataOnlyOptions{2};

%% Vigilance states
varNames = {'All','SWS','Awake','REM'};

%% ---------------------- Stage 1: NMF Extraction of Noise-reduced HFO distribution ------------------------
for iPt = 1:length(patientIds)

    curPt = sprintf('UMHS-00%d', patientIds(iPt));

    for iState = 1:numel(varNames)

        analyse_nmfByEpochByState_mainFunction( ...
            varNames{iState}, curPt, ...
            filePathData, filePathResult, ...
            hfoversion, winMin, ...
            RmAbnEpOption, RmFewValidSampleEpOption, InterictalDataOnly);

    end

    % Combine state-specific outputs (used in Stage 2 BSP analysis)
    combineDiffStateIntoOne( ...
        curPt, varNames, ...
        filePathResult, hfoversion, winMin, ...
        RmAbnEpOption, RmFewValidSampleEpOption, InterictalDataOnly);

end