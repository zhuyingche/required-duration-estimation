function nmfEpochResult = analyseNMF(Rates, nValidSamples, fs, winMin, timeDayRaw)
%======================================================================
% analyseNMF
%
% PURPOSE
% Implements Stage 1 of the reqDur framework by extracting noise-reduced
% HFO distributions from windowed HFO-rate matrices using
% non-negative matrix factorization (NMF).
%
% OVERVIEW
% Given an epoch-by-channel HFO rate matrix, this function:
%   1) Applies epoch-level quality control (QC)
%   2) Performs iterative NMF on progressively longer prefixes of the recording
%      (from the first epoch up to the full recording), saving decomposition
%      results for each time point. The full-recording decomposition is used
%      as the reference in downstream similarity/stability analyses (Stages 2–3).
%
% INPUTS
%   Rates                  - HFO rate matrix [nEpoch x nChannel]
%   nValidSamples          - Valid sample counts per epoch [nEpoch x nChannel] (qHFO)
%   fs                     - Sampling frequency (Hz)
%   winMin                 - Epoch duration (minutes)
%   timeDayRaw             - Epoch start times (days; relative to recording start)
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
%
%
%  Created in 2026 by Zhuying Chen (zhuychen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
%======================================================================

    %% -------------------- Epoch QC --------------------
    % Use a filtered copy of timing vector aligned to post-QC Rates
    timeDay = timeDayRaw;
    timeEndDay = timeEndDayRaw;

    % perform epoch QC
    [Rates, timeDay, timeEndDay] = applyEpochQC( ...
    Rates, timeDay, timeEndDay, nValidSamples, fs, winMin, ...
    filePathData, curPt);

    %% -------------------- Iterative NMF across increasing recording duration --------------------
    if ~isempty(Rates) && size(Rates,1) > 1
        data = Rates'; % [nChannel x nEpoch]

        % Store one result per truncation time point (epoch 1...nEpoch)
        nmfEpochResult = repmat(struct('perK',[],'K',[],'nmf',[],'rate',[],'timeDay',[],'timeEndDay',[]), size(data,2), 1);

        for nEpoch = 1:size(data,2)
            rate = data(:,1:nEpoch);  % progressively increasing prefix
            [perK, results] = run_nmf(rate);

            % Store full-recording rate at final epoch
            if nEpoch == size(data,2)
                nmfEpochResult(nEpoch).rate       = rate;
            end

            nmfEpochResult(nEpoch).nmf    = results;
            nmfEpochResult(nEpoch).perK   = perK;
            nmfEpochResult(nEpoch).K      = results.K;
            nmfEpochResult(nEpoch).timeDay = timeDay(nEpoch);
            nmfEpochResult(nEpoch).timeEndDay = timeEndDay(nEpoch);
        end

    else
        nmfEpochResult = [];
    end

  
end