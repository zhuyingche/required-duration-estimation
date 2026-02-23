function results = blind_source_separation( rate )
%% ======================================================================
% blind_source_separation
%
% PURPOSE
% Decomposes a channel-by-epoch HFO rate matrix into spatial and temporal
% components using robust non-negative matrix factorization (NNMF). This
% procedure identifies latent spatial HFO patterns (W) and their temporal
% activation profiles (H), and estimates the optimal number of channel
% groups/patterns (K).
%
% INPUT
%   rate : [nChannel × nEpoch] HFO rate matrix.
%
% OUTPUT (structure "results")
%   W0, H0  - Raw NNMF factors (unnormalized)
%   W,  H   - Normalized spatial (W) and temporal (H) components
%   K       - Estimated number of spatial patterns/channel groups
%
% METHOD SUMMARY
%   1) Validity check: requires detectable HFO activity (max(rate) ≥ 0.5/min).
%   2) Initialize model order: K0 = min(nChannel, nEpoch, 12).
%   3) Robust NNMF with repeated attempts to avoid degenerate solutions.
%   4) Model-order refinement: K is reduced iteratively until inter-component
%      redundancy is low (Spearman ρ ≤ 0.3 for both columns of W and rows of H).
%   5) Normalization: W and H are rescaled to preserve W*H while improving
%      interpretability (H approximates an effective HFO rate).
%
% INTERPRETATION
%   Columns of W represent spatial HFO patterns across channels.
%   Rows of H represent temporal modulation of each spatial pattern.
%   The decomposition approximates:
%       rate ≈ W × H
%
% ORIGIN
%   Created by S. Gliske (2016; CC-BY-NC-4.0).
%   Notes and refinements by Z. Chen (2024) for the reqDur framework.
%% ======================================================================

%% Validity check: require at least minimal detectable HFO activity
if ( numel(rate) < 1 || ~sum(rate(:),'omitmissing') || max(rate(:)) < 0.5/60 )
    results.K = 0;
    return;
end

%% Preparation
[nChan, nEp] = size(rate);

%% Initial model order
K = min(size(rate));
K = min(K, 12);   % K0 = min(nChan, nEpoch, 12)

%% Iterative factorization with redundancy control
maxC = 1;

while ( maxC > 0.3 )

    [W,H] = robust_nnmf(rate, K);

    % Ensure W has K components; repeat if necessary
    k = 0;
    while ( size(W,2) < K && k < 1000 )
        [W,H] = robust_nnmf(rate, K);
        k = k + 1;
    end

    % Handle non-finite solutions
    if ( ~all(isfinite(W(:))) || ~all(isfinite(H(:))) )
        if ( K > 1 ), K = K - 1; end
        continue;
    end

    % Handle degenerate (all-zero) factors
    if any(all(W == 0, 1)) || any(all(H == 0, 2))
        if ( K > 1 ), K = K - 1; end
        continue;
    end

    % If fewer components than expected, update K
    if ( size(W,2) < K )
        K = size(W,2);
        continue;
    end

    % Redundancy assessment (Spearman correlation)
    C = corr(H', 'type', 'spearman');
    C(eye(size(C)) == 1) = 0;
    maxC1 = max(C(:));

    C = corr(W, 'type', 'spearman');
    C(eye(size(C)) == 1) = 0;
    maxC2 = max(C(:));

    maxC = max(maxC1, maxC2);
    K = K - 1;

end

K = size(W,2);

%% Normalization: rescale W and H while preserving W*H
W0 = W;
H0 = H;

mu = sum(W,1) / nChan * K;
W  = W * diag(1./mu);
H  = diag(mu) * H;

%% Output structure
results.W0 = W0;
results.H0 = H0;
results.W  = W;
results.H  = H;
results.K  = K;

end