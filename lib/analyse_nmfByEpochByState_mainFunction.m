function analyse_nmfByEpochByState_mainFunction( ...
    varName, curPt, filePathData, filePathResult, ...
    hfoversion, winMin, RmAbnEpOption, ...
    RmFewValidSampleEpOption, InterictalDataOnly)
%===============================================================================
% analyse_nmfByEpochByState_mainFunction
%
% PURPOSE
% Performs state-specific non-negative matrix factorization (NMF) on windowed
% HFO rate matrices for a given patient. This function implements Stage 1 of
% the reqDur estimation framework.
%
% For the specified vigilance state (varName), the function:
%   1) Loads precomputed HFO rate matrices (window length = winMin).
%   2) Reconstructs epoch-level timing information when not explicitly stored.
%   3) Applies quality-control filters (abnormal epochs, insufficient
%      valid samples, and applies interictal-only restriction).
%   4) Executes NMF to obtain noise-reduced spatial HFO distributions.
%   5) Saves results using a version- and option-specific directory structure
%      to ensure analytical traceability and reproducibility.
%
% INPUTS
%   varName                  - Vigilance state ('All','SWS','Awake','REM')
%   curPt                    - Patient identifier (e.g., 'UMHS-0026')
%   filePathData             - Root directory containing HFO rate data
%   filePathResult           - Root directory for saving outputs
%   hfoversion               - HFO detector/version label
%   winMin                   - Window size (minutes) used for HFO rate estimation
%   RmAbnEpOption            - Abnormal-epoch removal strategy
%   RmFewValidSampleEpOption - Remove epochs with <90% valid samples (optional)
%   InterictalDataOnly       - Restrict analysis to interictal data (optional)
%
% OUTPUTS
%   Saves nmfEpochResult (e.g., W/H factors and metadata) to filePathResult.
%
% NOTES
%   - If time vectors are missing in the input file, they are reconstructed
%     from available indexing variables (pts/timeIdx).
% ==========================================================================================

    %% Construct input file path according to vigilance state and detector version
    if strcmp(varName, 'All') && strcmp(hfoversion, "qHFO_v2.3_running")
        dataFile = fullfile(filePathData, 'HFOrate', hfoversion, ...
            [num2str(winMin) 'min'], [curPt, '.HFOrate.mat']);
    else
        dataFile = fullfile(filePathData, 'HFOrate', hfoversion, ...
            [num2str(winMin) 'min'], 'SleepStage', ...
            [curPt, '.HFOrate.', varName, '.mat']);
    end

    %% Proceed only if HFO rate file is available
    if ~isfile(dataFile)
        warning('Data file %s does not exist.', dataFile);
        return;
    end

    %% Load precomputed HFO rate data
    load(dataFile);

    % Initialise timing arrays
    timeDayRaw    = [];
    timeEndDayRaw = [];

    %% Reconstruct epoch timing information when not explicitly stored
    % Some qHFO_v2.3_running exports may not include timeDay/timeEndDay.
    if ~exist('timeDay', 'var') && strcmp(hfoversion, "qHFO_v2.3_running")

        if strcmp(varName, 'All')
            % 'All' typically uses continuous indexing (pts) plus file start time.
            % Remove terminal element if it is an artefact of segmentation.
            if exist('pts', 'var') && ~isempty(pts)
                pts(end) = [];
            end

            timeDayRaw    = pts./fs./(60*60*24) + fileStartTimeMin./(60*24);
            timeEndDayRaw = timeDayRaw + winMin/60/24;

        else
            % Sleep-stage files may store epoch timing via timeIdx (start/stop minutes).
            uniqueIdx = unique(timeIdx.idx(:,1));
            timeDayRaw    = zeros(length(uniqueIdx), 1);
            timeEndDayRaw = zeros(length(uniqueIdx), 1);

            for j = 1:length(uniqueIdx)
                matchingRows = (timeIdx.idx == uniqueIdx(j));

                firstStartMin     = timeIdx.startMin(find(matchingRows, 1));
                timeDayRaw(j)     = firstStartMin./(60*24);

                lastStopMin       = timeIdx.stopMin(find(matchingRows, 1, 'last'));
                timeEndDayRaw(j)  = lastStopMin./(60*24);
            end
        end

    else
        % Use stored timing if available
        timeDayRaw = timeDay;

        % If timeEndDay is not stored, approximate using winMin
        if exist('timeEndDay', 'var')
            timeEndDayRaw = timeEndDay;
        else
            timeEndDayRaw = timeDayRaw + winMin/60/24;
        end
    end

    %% Execute NMF to obtain noise-reduced spatial HFO distributions
    nmfEpochResult = analyseNMF( ...
        Rates, nValidSamples, fs, winMin, ...
        RmAbnEpOption, RmFewValidSampleEpOption, ...
        timeDayRaw, hfoversion, InterictalDataOnly, ...
        filePathData, curPt, timeEndDayRaw);

    %% Define structured output directory (version- and option-specific)
    resultPath = fullfile(filePathResult, hfoversion, ...
        [num2str(winMin) 'min'], 'nmf', ...
        RmAbnEpOption, RmFewValidSampleEpOption, InterictalDataOnly);

    if ~exist(resultPath, 'dir')
        mkdir(resultPath);
    end

    %% Save NMF results and available timing metadata
    resultFile = fullfile(resultPath, strcat(curPt, '.nmfEpochResult.', varName, '.mat'));

    if exist('fileStartTimeMin', 'var') && exist('fileStopTimeMin', 'var')
        save(resultFile, 'nmfEpochResult', 'fileStartTimeMin', 'fileStopTimeMin');
    else
        save(resultFile, 'nmfEpochResult');
    end

end