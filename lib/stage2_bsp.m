function stage2_bsp(varName, patientId, filePathResult, filePathFig, winMin)
%% ======================================================================
% stage2_bsp
%
% PURPOSE
% Implements Stage 2 of the reqDur framework by quantifying how closely
% NMF-reconstructed noise-reduced HFO distribution derived from truncated
% recordings (up to epoch n) matches that derived from the full recording.
% Similarity is computed using a best-match cosine similarity strategy with
% padding to accommodate unequal observation lengths.
%
% INPUTS
%   varName                  - Vigilance state ('All','NREM','Awake','REM')
%   patientId                - Patient Id
%   filePathResult           - Root directory containing Stage 1 NMF outputs
%   filePathFig              - Output directory for figures
%   winMin                   - Epoch duration (minutes)
%
% OUTPUT
%   cosine_nmf          - cosine similarity trajectories relative to the full recording
%
% METHOD SUMMARY
%     1) Load NMF results computed at increasing recording lengths.
%     2) Reconstruct V = W*H for each shorter recording and V_all for the full recording.
%     3) Compute best-match cosine similarity between V and V_all, returning:
%         - nonPaddedSimilarity: similarity without padding
%         - paddedSimilarity: similarity attributable to padded segments
%         - overallSimilarity: combined similarity
%     4) Save similarity trajectories and generate a figure.
%
%
%  Created in 2026 by Zhuying Chen (zhuying.chen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
%% ======================================================================

    % Format patient ID
    curPt = sprintf('UMHS-00%d', patientId);

    % Initialize output structure
    cosine_nmf = struct();

    %% -------------------- Load Stage 1 NMF results --------------------
    % Construct file path for NMF results
    filename = fullfile(filePathResult, 'nmf', ...
        strcat(curPt, '.nmf.', varName, '.mat'));

    % Exit if file does not exist
    if ~isfile(filename)
        warning('NMF result file not found: %s', filename);
        return;
    end

    % Load nmfEpochResult (contains W (spatial pattern), H (temporal pattern) information)
    load(filename);

    % Exit if expected variable missing
    if ~exist('nmfEpochResult', 'var')
        warning('Variable nmfEpochResult not found in %s', filename);
        return;
    end   

    %% -------------------- Compute similarity trajectories --------------------
    % Extract Stage 1 results
    data = nmfEpochResult;
    nEpochs = size(data,1);
    nCh = size(data,1);

    % Use final epoch (full recording) as reference
    if ~isfield(data(nEpochs).nmf, 'W')
       return; % Skip if no NMF at full recording
    end

    % Reconstruct full-recording matrix
    V_all = data(nEpochs).nmf.W * data(nEpochs).nmf.H;

    % Preallocate similarity vectors
    cosine_nmf.overallSimilarity   = nan(1, nEpochs);
    cosine_nmf.nonPaddedSimilarity = nan(1, nEpochs);
    cosine_nmf.paddedSimilarity    = nan(1, nEpochs);

    % Loop over increasing recording lengths
    for n = 1:nEpochs

        % Only compute similarity if NMF exists at this epoch
        if isfield(data(n).nmf, 'W')

            % Reconstruct partial-recording matrix
            V = data(n).nmf.W * data(n).nmf.H;

            % Compute best-match cosine similarity with padding
            [overallSimilarity, nonPaddedSimilarity, paddedSimilarity] = ...
                best_match_similarity_padding(V, V_all);

            % Store similarity values
            cosine_nmf.overallSimilarity(n)   = overallSimilarity;
            cosine_nmf.nonPaddedSimilarity(n) = nonPaddedSimilarity;
            cosine_nmf.paddedSimilarity(n)    = paddedSimilarity;
        end
    end

    % Get timing and rate information
    timeDay    = [data.timeDay];
    timeEndDay = [data.timeEndDay];
    rate = data(nEpochs).rate;

    %% -------------------- Save similarity results --------------------
    % Define output directory
    resultPath = fullfile(filePathResult, 'bsp');

    if ~exist(resultPath, 'dir')
        mkdir(resultPath);
    end

    % Save similarity structure
    resultFile = fullfile(resultPath, ...
        strcat(curPt, '.bsp.', varName, '.mat'));

    save(resultFile, 'cosine_nmf', 'fileStartTimeMin', 'fileStopTimeMin','SOZ','rate','timeDay','timeEndDay');

    %% -------------------- Plot HFO heatmap --------------------
    figure

    subplot(2,1,1)

    % Adjust axis width
    h0 = gca;
    p0 = h0.Position;
    set(h0, 'Position', [p0(1), p0(2), p0(3) - 0.08, p0(4)]);
    hMain = gca;
    data = rate.*60;

    % Plot HFO rate matrix (channels × epochs)
    imagesc(hMain, 1:nEpochs, 1:nCh, data)

    % Robust color scaling
    clim_upper = prctile(data(:), 99.9);
    set(hMain, 'CLim', [0 clim_upper]);

    box off
    ylabel('Channel')

    %% Add SOZ side bar
    p0 = hMain.Position;
    hSOZ  = axes('Position', [p0(1) + p0(3), p0(2), 0.02, p0(4)]);
    imagesc(1, 1:nCh, SOZ);
    colormap(hSOZ, [1, 1, 1; 0.5, 0, 0.5]);
    set(hSOZ, 'YTick', [], 'XTick', [], ...
        'YColor', 'none', 'XColor', 'none');
    title(hSOZ, 'SOZ','Color',[0.5, 0, 0.5]);
    box off

    % Add colorbar
    c = colorbar(hMain, 'Position', ...
        [p0(1) + p0(3) + 0.043, p0(2), 0.015, p0(4)]);
    c.Label.String = 'HFO rate/min';

    %% -------------------- Plot similarity curve --------------------
    subplot(2,1,2)

    % Adjust subplot position
    h0 = gca;
    p0 = h0.Position;
    set(h0, 'Position', [p0(1), p0(2)-0.01, p0(3) - 0.08, p0(4)]);

    % Plot only if valid similarity exists
    if ~all(isnan(cosine_nmf.overallSimilarity(:)))

        n = length(cosine_nmf.overallSimilarity);

        % Plot similarity vs number of epochs
        plot(1:n, cosine_nmf.overallSimilarity);
        ylim([0 1])
        xlim([0 n])
        box off

        title(sprintf('%s; %s', curPt, varName))

        % Add day markers relative to recording start
        recStartDay = fileStartTimeMin/60/24;
        dayID = floor(timeDay - recStartDay) + 1;
        plotDayMarkers(dayID);

    else
        % If similarity unavailable
        title(sprintf('%s; %s (not available)', curPt, varName))
        ylim([0 1])
        xlim([0 n])
        box off
    end

    ylabel('Cosine similarity relative to full recording');
    xlabel(sprintf('Number of recorded %d-min epochs', winMin));

    %% -------------------- Save figure --------------------
    saveFig = 1;

    if saveFig
        outFigDir = fullfile(filePathFig, 'bsp');

        if ~exist(outFigDir, 'dir')
            mkdir(outFigDir);
        end

        savefig(gcf, fullfile(outFigDir, ...
            strcat(curPt, '.bsp.', varName, '.fig')));

        saveas(gcf, fullfile(outFigDir, ...
            strcat(curPt, '.bsp.', varName, '.tif')));
    end

end