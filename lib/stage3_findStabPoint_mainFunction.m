function  [stabPoint_timeDaySinceRecStart] = stage3_findStabPoint_mainFunction (varName, patientId, filePathResult, filePathFig, winMin)
%% ======================================================================
% stage3_findStabPoint_mainFunction
%
% PURPOSE
%   Implements Stage 3 of the reqDur framework:
%   Identification of the earliest time point at which the similarity
%   between truncated and full-recording distributions reaches and
%   maintains a high, stable plateau (defined as reqDur).
%
% DESCRIPTION
%   This function loads Stage 2 (BSP-based similarity) outputs for a
%   given patient and analyzes the cosine similarity trajectory derived
%   from progressively truncated recordings.
%
%   The earliest sustained high-similarity plateau is identified using
%   findStabPoint_cosine(), which determines the stabilization point
%   in the similarity curve.
%
%   The function:
%       1. Loads Stage 2 BSP similarity results.
%       2. Extracts the overall cosine similarity trajectory.
%       3. Identifies the earliest stable plateau (reqDur).
%       4. Converts the stabilization time to days since recording start.
%       5. Generates diagnostic visualization of:
%            - HFO rate distribution across channels and epochs
%            - Cosine similarity trajectory
%            - Stabilization point annotation
%
% INPUTS
%   varName        - Vigilance state ('All','NREM','Awake','REM')
%   patientId      - Patient Id
%   filePathResult - Directory containing Stage 2 BSP results
%   filePathFig    - Directory to save diagnostic figures
%   winMin         - Epoch duration (minutes)
%
% OUTPUT
%   stabPoint_timeDaySinceRecStart
%       - Required recording duration (reqDur) in days since
%         recording start at which similarity stabilizes.
%
% DEPENDENCIES
%   - cosine_nmf structure from Stage 2 BSP output
%   - findStabPoint_cosine()
%   - plotDayMarkers()
%
% NOTES
%   If BSP outputs or cosine similarity are unavailable, the function
%   exits gracefully with a warning.
%
%% ======================================================================
   
    % Format patient ID string
    curPt = sprintf('UMHS-00%d', patientId);

    %% -------------------- Load Stage 2 BSP results --------------------
    % Construct file path for BSP output
    filename = fullfile(filePathResult, 'bsp', ...
        strcat(curPt, '.bsp.', varName, '.mat'));

    % Exit if file does not exist
    if ~isfile(filename)
        warning('BSP result file not found: %s', filename);
        return;
    end

    % Load BSP result file (expects cosine_nmf structure)
    load(filename);

    % Exit if similarity structure is missing
    if ~exist('cosine_nmf', 'var')
        warning('Variable cosine_nmf not found in %s', filename);
        return;
    end   

    %% -------------------- Identify stabilization point --------------------
    % Extract similarity trajectory
    d = cosine_nmf.overallSimilarity';

    % Find earliest stable high-similarity plateau
    stable_point = findStabPoint_cosine(d);

    % Convert stabilization index to time (absolute and relative)
    stabPoint_timeDayClock = timeDay(stable_point);
    stabPoint_timeDaySinceRecStart = timeDay(stable_point) - fileStartTimeMin/60/24;

    %% -------------------- Save results --------------------
    % Define output directory
    resultPath = fullfile(filePathResult, 'findStablePoint');

    if ~exist(resultPath, 'dir')
        mkdir(resultPath);
    end

    % Save similarity structure
    resultFile = fullfile(resultPath, ...
        strcat(curPt, '.findStablePoint.', varName, '.mat'));

    save(resultFile, 'cosine_nmf','timeDay','timeEndDay', 'fileStartTimeMin', 'fileStopTimeMin','SOZ','rate','stable_point','stabPoint_timeDaySinceRecStart','stabPoint_timeDayClock');

    %% -------------------- Plot HFO rate heatmap --------------------
    figure
    data = rate.*60;
    nEpochs = size(data,1);
    nCh = size(data,1);


    subplot(2,1,1)

    % Adjust axis width to leave space for SOZ and colorbar
    h0 = gca;
    p0 = h0.Position;
    set(h0, 'Position', [p0(1), p0(2), p0(3) - 0.08, p0(4)]);
    hMain = gca;

    % Plot channel × epoch HFO rate matrix
    imagesc(hMain, 1:nEpochs, 1:nCh, data)

    % Limit color scale to 99.9th percentile for robustness
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
    colorbar(hMain, 'Position', ...
        [p0(1) + p0(3) + 0.043, p0(2), 0.015, p0(4)]);

    %% -------------------- Plot similarity curve --------------------
    subplot(2,1,2)

    % Adjust subplot position
    h0 = gca;
    p0 = h0.Position;
    set(h0, 'Position', [p0(1), p0(2)-0.01, p0(3) - 0.08, p0(4)]);

    % Plot only if valid similarity exists
    if ~all(isnan(cosine_nmf.overallSimilarity(:)))

        n = length(cosine_nmf.nonPaddedSimilarity);

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

        % Mark stabilization point
        hold on;
        plot(stable_point, d(stable_point), ...
            'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

        hold on
        % Annotate required duration
        text(stable_point*1.01, d(stable_point), ...
            sprintf('%.2f days', stabPoint_timeDaySinceRecStart), ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 8, 'Color', 'r');

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
        outFigDir = fullfile(filePathFig, 'findStablePoint');

        if ~exist(outFigDir, 'dir')
            mkdir(outFigDir);
        end

        savefig(gcf, fullfile(outFigDir, ...
            strcat(curPt, '.findStablePoint.', varName, '.fig')));

        saveas(gcf, fullfile(outFigDir, ...
            strcat(curPt, '.findStablePoint.', varName, '.tif')));
    end


end