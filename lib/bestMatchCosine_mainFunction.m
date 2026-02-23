%% ======================================================================
% bestMatchCosine_mainFunction
%
% PURPOSE
% Implements Stage 2 of the reqDur framework by quantifying how closely
% NMF-reconstructed HFO spatiotemporal structure derived from partial
% recordings (up to epoch n) matches that derived from the full recording.
% Similarity is computed using a best-match cosine similarity strategy with
% padding to accommodate unequal observation lengths.
%
% INPUTS
%   patientId                - Patient index/identifier (dataset-dependent)
%   filePathResult           - Root directory containing Stage 1 NMF outputs
%   filePathFig              - Output directory for diagnostic figures
%   hfoversion               - HFO detector/version label ('qHFO_v2.3_running' or 'NV')
%   winMin                   - Epoch duration (minutes)
%   RmAbnEpOption            - Abnormal-epoch removal strategy (Stage 1 QC)
%   RmFewValidSampleEpOption - Remove epochs with insufficient valid samples
%   normalizeMethod          - Normalization option used in Stage 1
%   InterictalDataOnly       - Optional restriction to interictal epochs
%
% OUTPUT
%   cosine_dist_nmf          - Structure indexed by vigilance state (e.g., All,
%                             SWS, Awake, REM) containing epoch-wise similarity
%                             trajectories relative to the full recording:
%                               .overallSimilarity
%                               .nonPaddedSimilarity
%                               .paddedSimilarity
%                             and associated timing vectors (.timeDay, .timeDayRaw).
%
% METHOD SUMMARY
%   For each vigilance state:
%     1) Load NMF results computed at increasing recording lengths.
%     2) Reconstruct V = W*H for each partial recording and V_all for the full recording.
%     3) Compute best-match cosine similarity between V and V_all, returning:
%         - nonPaddedSimilarity: similarity without padding
%         - paddedSimilarity: similarity attributable to padded segments
%         - overallSimilarity: combined similarity
%     4) Save similarity trajectories and generate a diagnostic figure.
%
% NOTES
%   - If NMF is not performed for a given epoch (e.g., insufficient HFO activity),
%     similarity values are recorded as NaN for that epoch.
%% ======================================================================

function cosine_dist_nmf = bestMatchCosine_mainFunction( ...
    patientId, filePathResult, filePathFig, hfoversion, winMin, ...
    RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly)

    cosine_dist_nmf = [];

    %% Resolve patient identifier by dataset
    if strcmp(hfoversion,'qHFO_v2.3_running')
        curPt = sprintf('UMHS-00%d', patientId);
    elseif strcmp(hfoversion,'NV')
        Patients = {'23_002', '23_003', '23_004', '23_005', '23_006', '23_007', ...
                    '24_001', '24_002', '24_004', '24_005', '25_001', '25_002', ...
                    '25_003', '25_004', '25_005'};
        curPt = Patients{patientId};
    else
        error('Unsupported hfoversion: %s', hfoversion);
    end

    %% Load Stage 1 outputs (nmfEpochResults)
    nmfPath = fullfile(filePathResult, hfoversion, [num2str(winMin) 'min'], 'nmf', ...
        RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly);

    filename = fullfile(nmfPath, strcat(curPt, '.nmfEpochResult.mat'));

    if ~isfile(filename)
        warning('NMF result file not found: %s', filename);
        return;
    end

    load(filename); % expects nmfEpochResults (and optionally fileStartTimeMin/fileStopTimeMin)

    if ~exist('nmfEpochResults', 'var')
        warning('Variable nmfEpochResults not found in %s', filename);
        return;
    end

    varNames = fieldnames(nmfEpochResults); % vigilance states
    cosine_dist_nmf = struct();

    %% Compute similarity trajectories per vigilance state
    for i = 1:numel(varNames)

        stateName = varNames{i};
        data = nmfEpochResults.(stateName);

        if isempty(data)
            continue;
        end

        nEpochs = length(data);

        % Full-recording reference V_all (only if NMF available at final epoch)
        if ~isfield(data(nEpochs).nmf, 'W')
            continue;
        end

        V_all = data(nEpochs).nmf.W * data(nEpochs).nmf.H;

        % Preallocate vectors for clarity
        cosine_dist_nmf.(stateName).overallSimilarity   = nan(1, nEpochs);
        cosine_dist_nmf.(stateName).nonPaddedSimilarity = nan(1, nEpochs);
        cosine_dist_nmf.(stateName).paddedSimilarity    = nan(1, nEpochs);

        for n = 1:nEpochs

            if isfield(data(n).nmf, 'W')
                V = data(n).nmf.W * data(n).nmf.H;

                [overallSimilarity, nonPaddedSimilarity, paddedSimilarity] = ...
                    best_match_cosine_similarity(V, V_all);

                cosine_dist_nmf.(stateName).overallSimilarity(n)   = overallSimilarity;
                cosine_dist_nmf.(stateName).nonPaddedSimilarity(n) = nonPaddedSimilarity;
                cosine_dist_nmf.(stateName).paddedSimilarity(n)    = paddedSimilarity;
            end

        end

        % Attach timing vectors
        cosine_dist_nmf.(stateName).timeDay    = [data.timeDay];
        cosine_dist_nmf.(stateName).timeDayRaw = data(nEpochs).timeDayRaw;

    end

    %% Diagnostic figure: similarity vs recording length (per vigilance state)
    figure
    set(gcf, "Position", [266.3333 142.3333 710 474.0000])
    fontSize = 12;

    sgtitle(sprintf('Patient %d: HFO spatiotemporal similarity across recording length', patientId), ...
        'FontSize', fontSize)

    for i = 1:numel(varNames)

        stateName = varNames{i};
        subplot(4, 2, i)

        if isfield(cosine_dist_nmf, stateName)

            n = length(cosine_dist_nmf.(stateName).nonPaddedSimilarity);

            plot(1:n, cosine_dist_nmf.(stateName).nonPaddedSimilarity, 'g', 'DisplayName', 'NonPadded Similarity');
            hold on
            plot(1:n, cosine_dist_nmf.(stateName).overallSimilarity,   'b', 'DisplayName', 'Overall Similarity');
            plot(1:n, cosine_dist_nmf.(stateName).paddedSimilarity,    'r', 'DisplayName', 'Padded Similarity');

            title(stateName)
            ylim([0 1])
            legend('NonPadded Similarity','Overall Similarity','Padded Similarity', 'Location','best')
            box off

        else
            title(sprintf('%s (not available)', stateName))
            ylim([0 1])
            box off
        end

    end

    % Common axis labels
    han = axes(gcf, 'visible', 'off');
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabel(han, 'Cosine similarity relative to full recording', 'FontSize', fontSize-2);
    xlabel(han, sprintf('Number of recorded %d-min epochs', winMin), 'FontSize', fontSize-2);

    %% Save similarity trajectories
    outDataDir = fullfile(filePathResult, hfoversion, [num2str(winMin) 'min'], 'bestMatchCosine', ...
        RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly);

    if ~exist(outDataDir, 'dir')
        mkdir(outDataDir);
    end

    outFile = fullfile(outDataDir, strcat(curPt, '.mat'));

    if exist('fileStartTimeMin', 'var') && exist('fileStopTimeMin', 'var')
        save(outFile, 'cosine_dist_nmf', 'fileStartTimeMin', 'fileStopTimeMin');
    else
        save(outFile, 'cosine_dist_nmf');
    end

    %% Save diagnostic figure
    saveFig = 1;
    if saveFig

        outFigDir = fullfile(filePathFig, hfoversion, [num2str(winMin) 'min'], 'bestMatchCosine', ...
            RmAbnEpOption, RmFewValidSampleEpOption, normalizeMethod, InterictalDataOnly);

        if ~exist(outFigDir, 'dir')
            mkdir(outFigDir);
        end

        savefig(gcf, fullfile(outFigDir, [curPt, '.fig']));
        saveas(gcf, fullfile(outFigDir, [curPt, '.tif']));
    end

end