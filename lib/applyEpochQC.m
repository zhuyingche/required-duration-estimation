function [Rates, timeDay, timeEndDay] = applyEpochQC( ...
    Rates, timeDay, timeEndDay, nValidSamples, fs, winMin, ...
    filePathData, curPt)
%==========================================================================
% applyEpochQC
%
% PURPOSE:
%   Perform epoch-level quality control (QC) for HFO rate analysis.
%
% INPUTS:
%   Rates           - [nEpoch x nChannel] HFO rate matrix
%   timeDay         - [nEpoch x 1] epoch start time (in days)
%   timeEndDay      - [nEpoch x 1] epoch end time (in days)
%   nValidSamples   - [nEpoch x nChannel] number of valid samples
%   fs              - sampling frequency (Hz)
%   winMin          - epoch window length (minutes)
%   filePathData    - base data directory
%   curPt           - patient ID
%
% OUTPUTS:
%   Rates           - cleaned rate matrix
%   timeDay         - cleaned start times
%   timeEndDay      - cleaned end times
%
% QC STEPS:
%   1. Remove epochs with <90% valid samples
%   2. Remove globally high-rate epochs
%   3. Remove globally low-rate epochs
%   4. Remove NaN-contaminated epochs
%   5. Remove epochs within 5 min before and 30 min after seizures
%
%  Created in 2026 by Zhuying Chen (zhuying.chen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
%==========================================================================

%% -------------------- Epoch QC (single-pass, all criteria) --------------------
nEpoch = size(Rates,1);

%% --- 1. Insufficient valid samples ---
threshold = 0.9;
validRatio = nanmean(nValidSamples ./ (fs * 60 * winMin), 2);
mask_validSample = validRatio < threshold & ~isnan(nanmean(Rates,2));

%% --- 2. Globally high-rate epochs ---
thrHigh = prctile(Rates, 65, "all");
mask_highRate = nansum(Rates > thrHigh, 2) > 0.80 * sum(~isnan(Rates), 2);

%% --- 3. Globally low-rate epochs ---
thrLow = 0.001; % count/sec (~0.06 count/min)
mask_lowRate = nanmean(Rates,2) < thrLow | prctile(Rates,75,2) == 0;

%% --- 4. NaN-contaminated epochs ---
mask_nan = any(isnan(Rates), 2);

%% --- 5. Peri-seizure epochs (5 min before, 30 min after) ---

% Load seizure data
dataFile = fullfile(filePathData,[curPt, '.HFOrate.All.mat']);
d = load(dataFile);

szStartTimesDay = d.szTimes(:,1) ./ 60 ./ 24;
szEndTimesDay   = d.szTimes(:,2) ./ 60 ./ 24;

preSzWinDay  = 5  / 60 / 24;
postSzWinDay = 30 / 60 / 24;

% Construct expanded seizure windows
exclusionStart = szStartTimesDay - preSzWinDay;
exclusionEnd   = szEndTimesDay   + postSzWinDay;

% Vectorized interval overlap test
mask_periSz = false(nEpoch,1);

for i = 1:length(exclusionStart)
    mask_periSz = mask_periSz | ...
        (timeEndDay >= exclusionStart(i) & timeDay <= exclusionEnd(i));
end

%% --- Combine ALL exclusion criteria ---
excludeMask = mask_validSample | mask_highRate | ...
              mask_lowRate | mask_nan | mask_periSz;

%% --- Apply filtering ONCE ---
Rates(excludeMask,:)   = [];
timeDay(excludeMask)   = [];
timeEndDay(excludeMask)= [];

end