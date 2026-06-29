%% ======================================================================
% Code accompanying the manuscript:
% Chen Z, Gliske SV, Alsammani A, et al. (2026).
% The influence of recording duration and vigilance state on
% high-frequency oscillation characterization in epilepsy.
% Neurology, 107(2), e218225.
%
% This repository implements a three-stage analytical framework to estimate
% the required recording duration (reqDur) necessary to adequately capture
% each vigilance state's high-frequency oscillation (HFO) distribution.
%
% Please cite the final published version when available.
% ----------------------------------------------------------------------
% Required Duration (reqDur) Estimation Pipeline
%
% The reqDur framework comprises three sequential stages:
%
% Stage 1: Extraction of noise-reduced HFO spatial distributions using
%          non-negative matrix factorization (NMF)
%          -> stage1_nmf.m
%
% Stage 2: Quantification of similarity between distributions derived from
%          truncated recordings and the full recording using best-match similarity padding (BSP)
%          -> stage2_bsp.m
%
% Stage 3: Identification of the earliest time point at which similarity
%          reaches and maintains a high, stable plateau (defined as reqDur)
%          -> stage3_findStabPoint.m
%
%  Created in 2026 by Zhuying Chen (zhuychen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
%% ======================================================================

%% Predefine file paths & core parameters
addpath(genpath('lib'));% Add library path

filePathData   = 'data\'; % data file path
filePathResult = 'result\'; % result file path
filePathFig = 'figure\'; % figure file path

winMin = 5; % HFO rate calculation window
patientId = 26; % Patient ID (e.g., 26)
varNames = {'All','NREM','Awake','REM'}; % Vigilance states: All: all state; NREM: NREM2+NREM3; Awake: Awake; REM: REM state

%% ---------- Stage 1: Extraction of Noise-reduced HFO distribution using NMF ----------
for iState = 2%1:numel(varNames) 
    varName = varNames{iState};
    stage1_nmf(varName, curPt, filePathData, filePathResult);
end

%% ---------- Stage 2: Quantification of similarity between truncated and full-recording distributions using BSP ----------
for iState = 2%1:numel(varNames) 
    varName = varNames{iState};
    stage2_bsp(varName, patientId, filePathResult, filePathFig, winMin);
end

%% ---------- Stage 3: Identification of the earliest time point with stable and high similarity value ---------- 
for iState = 2%1:numel(varNames) 
    varName = varNames{iState};
    stage3_findStabPoint(varName, patientId, filePathResult, filePathFig, winMin);
end

