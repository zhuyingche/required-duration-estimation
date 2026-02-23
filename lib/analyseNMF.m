function nmfEpochResult = analyseNMF(Rates, nValidSamples, fs, winMin, ...
    RmAbnEpOption, RmFewValidSampleEpOption, timeDayRaw, hfoversion, ratioVS)
%======================================================================
% analyseNMF
%
% PURPOSE
% Implements Stage 1 of the reqDur framework by extracting noise-reduced
% HFO spatial distributions from windowed HFO-rate matrices using
% non-negative matrix factorization (NMF).
%
% OVERVIEW
% Given an epoch-by-channel HFO rate matrix, this function:
%   1) Applies epoch-level quality control (QC), including removal of epochs
%      with insufficient valid samples and/or globally abnormal rate profiles.
%   2) Removes any remaining NaN-contaminated epochs to ensure numerical stability.
%   3) Performs iterative NMF on progressively longer prefixes of the recording
%      (from the first epoch up to the full recording), saving decomposition
%      results for each time point. The full-recording decomposition is used
%      as the reference in downstream similarity/stability analyses (Stages 2–3).
%
% INPUTS
%   Rates                  - HFO rate matrix [nEpoch x nChannel]
%   nValidSamples          - Valid sample counts per epoch [nEpoch x nChannel] (qHFO)
%   fs                     - Sampling frequency (Hz)
%   winMin                 - Epoch duration (minutes)
%   RmAbnEpOption          - Abnormal-epoch removal option (e.g., 'RmHighLowEpoch')
%   RmFewValidSampleEpOption - Remove epochs with <90% valid samples (optional)
%   timeDayRaw             - Epoch start times (days; relative to recording start)
%   hfoversion             - HFO detector/version label ('qHFO_v2.3_running' or 'NV')
%   ratioVS                - Valid-sample ratio matrix (NV only), same size as Rates
%
% OUTPUT
%   nmfEpochResult         - Struct array (length = nEpoch) where each element
%                            contains NMF decomposition results for data up to
%                            that epoch, plus timing metadata.
%
% REPRODUCIBILITY NOTES
%   - QC decisions are explicitly parameterized (thresholds, percentiles).
%   - NMF is repeated (nReps = 10) and the most frequently observed model order
%     (K) is selected to reduce sensitivity to random initialization.
%======================================================================


    % Use a filtered copy of timing vector aligned to post-QC Rates
    timeDay = timeDayRaw;

    %% -------------------- Epoch QC: insufficient valid samples --------------------
    % Exclude epochs with valid-sample ratio < 0.9 (i.e., <90% of expected samples).
    if strcmp(RmFewValidSampleEpOption,'RmFewValidSampleEpoch')
        threshold = 0.9;

        if strcmp(hfoversion,"qHFO_v2.3_running")
            validRatio = nanmean(nValidSamples./(fs*60*winMin), 2);
            I = find(validRatio < threshold & ~isnan(nanmean(Rates,2)));

        elseif strcmp(hfoversion,'NV')
            validRatio = nanmean(ratioVS, 2);
            I = find(validRatio < threshold & ~isnan(nanmean(Rates,2)));
        end

        Rates(I,:) = [];
        timeDay(I) = [];
    end

    %% -------------------- Epoch QC: globally abnormal rate profiles --------------------
    % Remove epochs where the majority of channels show unusually high rates.
    if strcmp(RmAbnEpOption,'RmHighEpoch') || strcmp(RmAbnEpOption,'RmHighLowEpoch') || ...
       strcmp(RmAbnEpOption,'RmHighLowEpoch_bold')

        thrHigh = prctile(Rates, 65, "all");
        I = find(nansum(Rates > thrHigh, 2) > 0.80 * sum(~isnan(Rates), 2));
        Rates(I,:) = [];
        timeDay(I) = [];

        % Optional: additionally remove globally low-rate epochs.
        if strcmp(RmAbnEpOption,'RmHighLowEpoch')
            thrLow = 0.001; % empirically defined (approx. 1st percentile scale)
            I = find(nanmean(Rates,2) < thrLow | prctile(Rates,75,2) == 0);
            Rates(I,:) = [];
            timeDay(I) = [];
        end
    end

    %% -------------------- Remove NaN-contaminated epochs --------------------
    I = find(any(isnan(Rates), 2));
    Rates(I,:) = [];
    timeDay(I) = [];

    %% -------------------- Iterative NMF across increasing recording duration --------------------
    if ~isempty(Rates) && size(Rates,1) > 1

        data = Rates'; % [nChannel x nEpoch]

        % Store one result per truncation time point (epoch 1...nEpoch)
        nmfEpochResult = repmat(struct('perK',[],'K',[],'nmf',[],'rate',[]), size(data,2), 1);

        for nEpoch = 1:size(data,2)

            rate = data(:,1:nEpoch);  % progressively increasing prefix

            [perK, results] = run_nmf(rate);

            % Store full-recording rate and original timing only at final epoch
            if nEpoch == size(data,2)
                nmfEpochResult(nEpoch).rate       = rate;
                nmfEpochResult(nEpoch).timeDayRaw = timeDayRaw;
            end

            nmfEpochResult(nEpoch).nmf    = results;
            nmfEpochResult(nEpoch).perK   = perK;
            nmfEpochResult(nEpoch).K      = results.K;
            nmfEpochResult(nEpoch).timeDay = timeDay(nEpoch);
        end

    else
        nmfEpochResult = [];
    end

    %% ======================== Nested helper: repeated NMF =========================
    function [perK, results] = run_nmf(rate)
        % Repeat NMF to mitigate initialization sensitivity; select the
        % decomposition whose model order (K) occurs most frequently.
        nReps = 10;

        results_all  = cell(nReps,1);
        K_all        = nan(nReps,1);
        category_all = cell(nReps,1);

        for nRep = 1:nReps
            tmp = blind_source_separation(rate);
            results_all{nRep} = tmp;
            K_all(nRep) = tmp.K;

            if tmp.K > 0
                category_all{nRep} = tmp.label;
            end
        end

        idx_mostFreq = find(K_all == mode(K_all), 1, 'first');
        perK = sum(K_all == mode(K_all)) / sum(~isnan(K_all)) * 100;
        results = results_all{idx_mostFreq};
    end

end