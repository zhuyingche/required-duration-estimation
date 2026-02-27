function [overallSimilarity, nonPaddedSimilarity, paddedSimilarity, paddedX, paddedY, pathX, pathY] = best_match_similarity_padding(x, y)
%==========================================================================
% best_match_similarity_padding
%
% PURPOSE:
%   Computes similarity between two multichannel signals of unequal length
%   by padding the shorter signal using a best-match segment selection
%   strategy. 
%
% INPUT:
%   x, y   - Multichannel signals of size [nChannels × nObservations].
%            Both signals must have the same number of channels.
%
% OUTPUT:
%   overallSimilarity      - Mean cosine similarity across the fully padded signals
%   nonPaddedSimilarity    - Mean cosine similarity over the originally aligned portion
%   paddedSimilarity       - Mean cosine similarity over the padded portion only
%   paddedX, paddedY       - Length-matched signals after padding
%   pathX, pathY           - Index mappings indicating segment selection during padding
%
% METHOD:
%   1. Identify the shorter signal.
%   2. For each missing observation, search within the shorter signal for
%      candidate segments that best match the corresponding segment in the
%      longer signal.
%   3. Rank candidate segments by cosine similarity.
%   4. From the top-ranking subset (e.g., top 10%), select the segment with
%      minimal absolute dissimilarity.
%   5. Append the selected segment to the shorter signal.
%   6. Compute mean cosine similarity for:
%        (i)  the entire padded signals,
%        (ii) the originally aligned portion,
%        (iii) the padded portion.
%
% RATIONALE:
%   This hybrid similarity strategy prioritizes angular alignment (cosine
%   similarity) while refining matches using magnitude-based dissimilarity.
%   The approach is particularly suitable when structural pattern similarity
%   is of primary interest.
%
% NOTE:
%   If both signals have equal length, no padding is performed and the
%   padded similarity is returned as NaN.
%
%  Created in 2026 by Zhuying Chen (zhuychen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
%==========================================================================

    % Ensure x and y have the same number of channels (rows)
    [mCh1, nObs1] = size(x);
    [mCh2, nObs2] = size(y);
    
    if mCh1 ~= mCh2
        error('x and y must have the same number of rows (channels)');
    end
    
    % Initialize paddedX, paddedY, and paths
    paddedX = x;
    paddedY = y;
    pathX = 1:nObs1;
    pathY = 1:nObs2;
    
    if nObs1 < nObs2
        % Pad x to match the length of y
        paddingLength = nObs2 - nObs1;
        
        for i = 1:paddingLength
            segmentToMatch = y(:, nObs1 + i);
            cosineSimilarities = nan(1, nObs1);
            absDissimilarities = nan(1, nObs1);
            
            for j = 1:nObs1
                currentSegment = x(:, j);
                cosineSimilarities(j) = dot(currentSegment, segmentToMatch) / ...
                    (norm(currentSegment) * norm(segmentToMatch));
                absDissimilarities(j) = sum(abs(currentSegment - segmentToMatch), 'all');
            end
            
            % Find the indices of the best matches with the highest cosine similarity
            [~, sortedIndices] = sort(cosineSimilarities, 'descend');
            nBest =  ceil(min(nObs1, nObs2) * 0.1);
            topIndices = sortedIndices(1:nBest);
            
            % Select the best match among the top matches based on the shortest absolute dissimilarity
            [~, bestIndex] = min(absDissimilarities(topIndices));
            bestSegmentIndex = topIndices(bestIndex);
            
            % Pad with the best matching segment and update pathX
            paddedX = [paddedX, x(:, bestSegmentIndex)];
            pathX = [pathX, bestSegmentIndex];
        end
        pathY = 1:nObs2;
        
    elseif nObs1 > nObs2
        % Pad y to match the length of x
        paddingLength = nObs1 - nObs2;
        
        for i = 1:paddingLength
            segmentToMatch = x(:, nObs2 + i);
            cosineSimilarities = nan(1, nObs2);
            absDissimilarities = nan(1, nObs2);
            
            for j = 1:nObs2
                currentSegment = y(:, j);
                cosineSimilarities(j) = dot(currentSegment, segmentToMatch) / ...
                    (norm(currentSegment) * norm(segmentToMatch));
                absDissimilarities(j) = sum(abs(currentSegment - segmentToMatch), 'all');
            end
            
            % Find the indices of the best matches with the highest cosine similarity
            [~, sortedIndices] = sort(cosineSimilarities, 'descend');
            nBest = ceil(min(nObs1, nObs2) * 0.1);
            topIndices = sortedIndices(1:nBest);
            
            % Select the best match among the top matches based on the shortest absolute dissimilarity
            [~, bestIndex] = min(absDissimilarities(topIndices));
            bestSegmentIndex = topIndices(bestIndex);
            
            % Pad with the best matching segment and update pathY
            paddedY = [paddedY, y(:, bestSegmentIndex)];
            pathY = [pathY, bestSegmentIndex];
        end
        pathX = 1:nObs1;
        
    else
        % No padding needed, pathX and pathY are simple mappings
        paddedX = x;
        paddedY = y;
    end
    
    % Calculate the average cosine similarity over the entire padded signals
    overallSimilarity = mean(sum(paddedX .* paddedY, 1) ./ ...
        (sqrt(sum(paddedX.^2, 1)) .* sqrt(sum(paddedY.^2, 1))));
    
    % Calculate the cosine similarity only for the non-padded part (no padding)
    nAligned = min(nObs1, nObs2);
    nonPaddedSimilarity = mean(sum(x(:, 1:nAligned) .* y(:, 1:nAligned), 1) ./ ...
        (sqrt(sum(x(:, 1:nAligned).^2, 1)) .* sqrt(sum(y(:, 1:nAligned).^2, 1))));
    
    % Calculate the cosine similarity only for the padded part
    if nObs1 < nObs2
        paddedPartX = paddedX(:, nObs1 + 1:end);
        paddedPartY = paddedY(:, nObs1 + 1:end);
    elseif nObs1 > nObs2
        paddedPartX = paddedX(:, nObs2 + 1:end);
        paddedPartY = paddedY(:, nObs2 + 1:end);
    else
        paddedPartX = [];
        paddedPartY = [];
    end
    
    if isempty(paddedPartX) || isempty(paddedPartY)
        paddedSimilarity = NaN; % No padding was applied
    else
        paddedSimilarity = mean(sum(paddedPartX .* paddedPartY, 1) ./ ...
            (sqrt(sum(paddedPartX.^2, 1)) .* sqrt(sum(paddedPartY.^2, 1))));
    end

end
