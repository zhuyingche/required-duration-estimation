function [perK, results] = run_nmf(rate)
%==========================================================================
% run_nmf
%
% PURPOSE:
%   Perform repeated non-negative matrix factorization (NMF) to reduce
%   sensitivity to random initialization and select a stable model order (K).
%
% DESCRIPTION:
%   NMF solutions can vary depending on random initialization. To improve
%   robustness, this function repeats the blind source separation procedure
%   multiple times (default = 10 repetitions). For each repetition, the
%   estimated model order (K) is recorded.
%
%   The final solution is selected based on the most frequently occurring
%   K across repetitions (majority rule). The decomposition corresponding
%   to the first occurrence of that modal K is returned.
%
%   Additionally, the function computes the percentage of repetitions that
%   yielded the selected K (perK), which serves as a stability metric for
%   model order selection.
%
% INPUT:
%   rate    - [nChannel × nEpoch] non-negative matrix (e.g., HFO rates)
%
% OUTPUTS:
%   perK    - Percentage of repetitions supporting the selected K
%   results - Structure containing the selected NMF decomposition
%             (as returned by blind_source_separation), including:
%             • W  (spatial components)
%             • H  (temporal activations)
%             • K  (model order)
%
%  Created in 2026 by Zhuying Chen (zhuychen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
%==========================================================================

    nReps = 10;

    results_all  = cell(nReps,1);
    K_all        = nan(nReps,1);

    for nRep = 1:nReps
        tmp = blind_source_separation(rate);
        results_all{nRep} = tmp;
        K_all(nRep) = tmp.K;
    end

    idx_mostFreq = find(K_all == mode(K_all), 1, 'first');
    perK = sum(K_all == mode(K_all)) / sum(~isnan(K_all)) * 100;
    results = results_all{idx_mostFreq};

end