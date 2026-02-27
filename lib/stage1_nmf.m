function stage1_nmf(varName, patientId, filePathData, filePathResult)
%===============================================================================
% stage1_nmf
%
% PURPOSE
% This function implements Stage 1 of the reqDur estimation framework:
%
% Stage 1: Extraction of noise-reduced HFO spatial distributions using
%          non-negative matrix factorization (NMF)
%
% For the specified vigilance state (varName), the function:
%   1) Loads precomputed HFO rate matrices (window length = winMin).
%   2) Executes NMF to obtain noise-reduced spatial HFO distributions.
%   3) Saves results
%
% INPUTS
%   varName                  - Vigilance state ('All','NREM','Awake','REM')
%   patientId                - Patient Id
%   filePathData             - Root directory containing HFO rate data
%   filePathResult           - Root directory for saving outputs
%
% OUTPUTS
%   Saves nmfEpochResult (e.g., W/H factors and metadata) to filePathResult.
%
%  Created in 2026 by Zhuying Chen (zhuychen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
% ==========================================================================================

    %% Construct input file path
    curPt = sprintf('UMHS-00%d', patientId);  % Patient identifier (e.g., 'UMHS-0026')
    dataFile = fullfile(filePathData, ....
            [curPt, '.HFOrate.', varName, '.mat']);

    %% Proceed only if HFO rate file is available
    if ~isfile(dataFile)
        warning('Data file %s does not exist.', dataFile);
        return;
    end

    %% Load precomputed HFO rate data
    load(dataFile);

    %% Execute NMF to obtain noise-reduced HFO distributions
    nmfEpochResult = analyseNMF( ...
        Rates, nValidSamples, fs, winMin, ...
        timeDayRaw, timeEndDayRaw,filePathData, curPt);

    %% Define structured output directory
    resultPath = fullfile(filePathResult, 'nmf');

    if ~exist(resultPath, 'dir')
        mkdir(resultPath);
    end

    %% Save NMF results and available timing metadata
    resultFile = fullfile(resultPath, strcat(curPt, '.nmf.', varName, '.mat'));
    save(resultFile, 'nmfEpochResult', 'fileStartTimeMin', 'fileStopTimeMin','SOZ');%fileStartTimeMin：the clock time that the recording start since the midnight of electrode implantation

end