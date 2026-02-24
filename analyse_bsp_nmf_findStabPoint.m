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
% CURRENT SCRIPT: Stage 3 (single-patient demonstration)
%
% PURPOSE
% Demonstrates Stage 3 by computing reqDur for one patient (optionally
% across vigilance states). The script loads Stage 2 similarity trajectories
% and identifies the earliest epoch at which similarity stabilizes at a
% high level for a minimum sustained duration. Each subplot additionally
% shows a DayID strip aligned to epoch index to aid interpretation.
%
% NOTES
% - Minimal release version: single-patient computation only.
% - Cohort-level summaries and statistical tests are intentionally omitted.
%
%% ======================================================================

%% ------------------------ Paths ----------------------------------------
filePathResult = 'C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data analysis\result\';
filePathFig    = 'C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data analysis\figure';
addpath(genpath('C:\Users\ZHUYCHEN\OneDrive - The University of Melbourne\2024 work\data analysis\code\lib'));

%% ------------------------ Parameters -----------------------------------
patientId   = 26;                    % demonstration patient
hfoversion  = "qHFO_v2.3_running";
winMin      = 5;                     % epoch duration (minutes)

% Stage 1/2 consistency options (must match those used upstream)
RmAbnEpOptions = {'','RmHighEpoch','RmHighLowEpoch','RmHighLowEpoch_bold','RmLowHighEpoch'};
RmAbnEpOption  = RmAbnEpOptions{5};

RmFewValidSampleEpOptions = {'','RmFewValidSampleEpoch'};
RmFewValidSampleEpOption  = RmFewValidSampleEpOptions{2};

normalizeMethods = {'','zscore','minMax','L1norm'};
normalizeMethod  = normalizeMethods{1};

InterictalDataOnlyOptions = {'','InterictalDataOnly'};
InterictalDataOnly = InterictalDataOnlyOptions{2};

% Vigilance-state list (must match Stage 2 outputs)
varNames = {'All','SWS','Awake','REM'};

%% ------------------------ Stage 3 Execution ----------------------------
[stabPoint_nEpoch, stabPoint_ratioEpoch, stabPoint_totalEpoch, ...
 stabPoint_timeDayClock, stabPoint_timeDayRelRecStart, ...
 stabPoint_ratioDurDay, stabPoint_totalDurDay] = ...
    FindStabPoint_bestMatchCosine_cosine( ...
        varNames, patientId, filePathResult, filePathFig, ...
        hfoversion, winMin, RmAbnEpOption, RmFewValidSampleEpOption, ...
        normalizeMethod, InterictalDataOnly);

disp(stabPoint_timeDayRelRecStart)


%% ======================================================================
% Core routine: compute reqDur (stability point) from similarity trajectory
% + plot DayID strip aligned to epoch index in each subplot
%% ======================================================================
function [stabPoint_nEpoch, stabPoint_ratioEpoch, stabPoint_totalEpoch, ...
          stabPoint_timeDayClock, stabPoint_timeDayRelRecStart, ...
          stabPoint_ratioDurDay, stabPoint_totalDurDay] = ...
    FindStabPoint_bestMatchCosine_cosine( ...
        varNames, patientId, filePathResult, filePathFig, hfoversion, winMin, ...
        RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly)

    % Outputs (tables indexed by vigilance state)
    varNamesTbl = [{'PatientID'}, varNames];
    stabPoint_nEpoch            = array2table(nan(1, numel(varNamesTbl)), 'VariableNames', varNamesTbl);
    stabPoint_ratioEpoch        = array2table(nan(1, numel(varNamesTbl)), 'VariableNames', varNamesTbl);
    stabPoint_totalEpoch        = array2table(nan(1, numel(varNamesTbl)), 'VariableNames', varNamesTbl);
    stabPoint_timeDayClock      = array2table(nan(1, numel(varNamesTbl)), 'VariableNames', varNamesTbl);
    stabPoint_timeDayRelRecStart= array2table(nan(1, numel(varNamesTbl)), 'VariableNames', varNamesTbl);
    stabPoint_ratioDurDay       = array2table(nan(1, numel(varNamesTbl)), 'VariableNames', varNamesTbl);
    stabPoint_totalDurDay       = array2table(nan(1, numel(varNamesTbl)), 'VariableNames', varNamesTbl);

    stabPoint_nEpoch.PatientID             = patientId;
    stabPoint_ratioEpoch.PatientID         = patientId;
    stabPoint_totalEpoch.PatientID         = patientId;
    stabPoint_timeDayClock.PatientID       = patientId;
    stabPoint_timeDayRelRecStart.PatientID = patientId;
    stabPoint_ratioDurDay.PatientID        = patientId;
    stabPoint_totalDurDay.PatientID        = patientId;

    % Load Stage 2 output (cosine_dist_nmf struct)
    curPt = sprintf('UMHS-00%d', patientId);
    inDir = fullfile(filePathResult, hfoversion, [num2str(winMin) 'min'], ...
        'bestMatchCosine', RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly);
    inFile = fullfile(inDir, strcat(curPt, '.mat'));

    if ~isfile(inFile)
        warning('Stage 2 output not found: %s', inFile);
        return;
    end

    load(inFile); % expects cosine_dist_nmf and fileStartTimeMin/fileStopTimeMin

    if ~exist('cosine_dist_nmf', 'var')
        warning('cosine_dist_nmf not found in %s', inFile);
        return;
    end

    if ~exist('fileStartTimeMin','var') || ~exist('fileStopTimeMin','var')
        warning('fileStartTimeMin/fileStopTimeMin not found in %s (DayID/time conversions will be limited).', inFile);
    end

    % Global y-limit for consistent visualization
    max_value = 0;
    for i = 1:numel(varNames)
        if isfield(cosine_dist_nmf, varNames{i})
            max_value = max(max_value, max(cosine_dist_nmf.(varNames{i}).overallSimilarity));
        end
    end

    figure
    set(gcf,"Position",[266.3333 142.3333 710 474.0000])
    fontSize   = 12;
    markerSize = 3;
    sgtitle(sprintf('Patient %d: HFO similarity stabilization (reqDur)', patientId), 'FontSize', fontSize)

    for i = 1:numel(varNames)

        state = varNames{i};
        axMain = subplot(ceil(numel(varNames)/2), 2, i);

        if ~isfield(cosine_dist_nmf, state)
            title(axMain, sprintf('%s (not available)', state));
            ylim(axMain, [0 max_value]); box(axMain, 'off');
            continue;
        end

        % Similarity trajectory and timing
        sim     = cosine_dist_nmf.(state).overallSimilarity(:);
        timeDay = cosine_dist_nmf.(state).timeDay(:);

        plot(axMain, 1:numel(sim), sim, 'LineWidth', 1);
        title(axMain, state);
        ylim(axMain, [0 max_value]);
        box(axMain, 'off');
        hold(axMain, 'on');

        % Total epochs/duration
        stabPoint_totalEpoch.(state) = sum(~isnan(sim));

        if exist('fileStartTimeMin','var') && exist('fileStopTimeMin','var')
            recStartDay = fileStartTimeMin/60/24;
            totalDurDay = (fileStopTimeMin - fileStartTimeMin)/60/24;
            stabPoint_totalDurDay.(state) = totalDurDay;

            % DayID relative to recording start (Day 1, Day 2, ...)
            dayID = floor(timeDay - recStartDay) + 1;
        else
            recStartDay = nan;
            totalDurDay = nan;
            dayID = ones(numel(timeDay),1); % fallback
        end

        % reqDur as stability point on similarity curve
        stable_point = FindStabPoint_cosine(sim, RmAbnEpOption);

        if ~isempty(stable_point)
            stabPoint_nEpoch.(state)     = stable_point;
            stabPoint_ratioEpoch.(state) = stable_point / sum(~isnan(sim));

            if ~isnan(recStartDay)
                stabPoint_timeDayClock.(state)       = timeDay(stable_point);
                stabPoint_timeDayRelRecStart.(state) = timeDay(stable_point) - recStartDay;
                if ~isnan(totalDurDay) && totalDurDay > 0
                    stabPoint_ratioDurDay.(state) = stabPoint_timeDayRelRecStart.(state) / totalDurDay;
                end
            end

            plot(axMain, stable_point, sim(stable_point), 'ro', ...
                'MarkerSize', markerSize, 'MarkerFaceColor', 'r');
        end

        % ---------------- DayID strip (aligned to epoch index) -------------
        addDayIDStrip(axMain, dayID);
        % ------------------------------------------------------------------

    end

    % Common axis labels
    han = axes(gcf, 'visible', 'off');
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabel(han, 'Cosine similarity relative to full recording', 'FontSize', fontSize-2);
    xlabel(han, sprintf('Number of recorded %d-min epochs', winMin), 'FontSize', fontSize-2);

    % Optional figure export (disabled by default)
    saveFig = 0;
    if saveFig
        outDir = fullfile(filePathFig, hfoversion, [num2str(winMin) 'min'], ...
            'bestMatchCosine', RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly, ...
            'FindStabPoint_cosine');
        if ~exist(outDir, 'dir'); mkdir(outDir); end
        savefig(gcf, fullfile(outDir, [curPt, '.fig']));
        saveas(gcf,  fullfile(outDir, [curPt, '.tif']));
    end
end


%% ======================================================================
% Helper: add a DayID strip below a main axis, aligned to epoch index
%% ======================================================================
function addDayIDStrip(axMain, dayID)

    if isempty(dayID) || numel(dayID) < 2
        return;
    end

    % Reserve space under main axis for DayID strip
    p = axMain.Position;
    stripH = 0.035;
    gapH   = 0.008;

    axMain.Position = [p(1), p(2) + stripH + gapH, p(3), p(4) - (stripH + gapH)];
    p2 = axMain.Position;

    % Create strip axis directly beneath main axis (same width)
    axStrip = axes('Position', [p2(1), p2(2) - (stripH + gapH), p2(3), stripH]);
    plotDayMarkers(dayID);

    % Match x-limits to the main axis
    xlim(axStrip, xlim(axMain));

    % Clean strip axis formatting
    set(axStrip, 'XTick', [], 'Box', 'off', 'Color', 'none');
end


%% ======================================================================
% Helper: stability point detection on cosine similarity trajectory
%% ======================================================================
function stable_point = FindStabPoint_cosine(similarity, RmAbnEpOption)
    d = similarity(:);
    if numel(d) < 2 || all(isnan(d))
        stable_point = [];
        return;
    end

    slopeThreshold = 1e-2;
    minDuration    = 5;

    if strcmp(RmAbnEpOption,'RmLowHighEpoch')
        similarityThreshold = 1 - 0.075;
    else
        similarityThreshold = 1 - 0.05;
    end

    x  = (1:numel(d))';
    dy = diff(d) ./ diff(x);

    belowSlope = abs(dy) < slopeThreshold;
    belowSlope = [belowSlope; true]; % align lengths

    eligible = belowSlope & (d > similarityThreshold);

    convRes = conv(double(eligible), ones(minDuration,1), 'valid');
    convRes = [convRes; zeros(minDuration-1,1)];

    stable_point = find(convRes == minDuration, 1, 'first');
end


%% ======================================================================
% DayID strip plotting: marks day boundaries and labels day segments
%% ======================================================================
function plotDayMarkers(dayID)
    nEps = numel(dayID);

    y1 = 0.02;
    y2 = 0.20;

    hold on

    lastMark = 1;

    for j = 2:nEps
        if dayID(j-1) ~= dayID(j)
            line([j j], [y1 y2], 'Color', 'k', 'LineWidth', 1);

            text((lastMark + (j-1))/2, y2, sprintf('%d', dayID(j-1)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize', 8, 'Color','k');

            lastMark = j;
        end
    end

    % last segment label
    text((lastMark + nEps)/2, y2, sprintf('%d', dayID(end)), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize', 8, 'Color','k');

    % boundaries
    line([1 1],       [y1 y2], 'Color','k', 'LineWidth', 1);
    line([nEps nEps], [y1 y2], 'Color','k', 'LineWidth', 1);

    xlim([1 nEps]);
    ylim([0 0.25]);
    set(gca, 'YTick', [], 'YColor','none');
    hold off
end