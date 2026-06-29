function stable_point = findStabPoint_cosine(d)
%==========================================================================
% findStabPoint_cosine
%
% PURPOSE:
% Identifies the earliest stabilization point in a cosine similarity
% trajectory. A stable point is defined as the first index at which the
% signal becomes both sufficiently high (above a cosine threshold) and
% exhibits minimal change (derivative below a predefined threshold) for a
% specified number of consecutive points.
%
% INPUT:
%   d   - Vector of cosine similarity values (e.g., similarity trajectory
%         relative to full-recording reference).
%
% OUTPUT:
%   stable_point   - Index of the first detected stable point satisfying
%                    all criteria. Returns empty if no stable segment is found.
%
% METHOD:
%   1. Compute the first derivative of the trajectory.
%   2. Identify indices where the derivative magnitude is below a small
%      threshold (indicating flattening).
%   3. Apply an additional constraint requiring cosine similarity to exceed
%      a predefined threshold.
%   4. Use convolution to detect the first segment where the criteria are
%      satisfied for a minimum number of consecutive points.
%
% NOTE:
%   Thresholds and minimum duration parameters can be adjusted depending
%   on the desired sensitivity for stabilization detection.
%
%  Created in 2026 by Zhuying Chen (zhuying.chen@unimelb.edu.au)
%  Released under the CC-BY-NC-4.0 License
%  http://creativecommons.org/licenses/by-nc/4.0/
%==========================================================================

% --- Parameters controlling stabilization criteria ---
threshold = 1e-2;          % Derivative threshold for considering the curve "flat"
min_duration = 5;         % Required number of consecutive stable points
cosine_threshold = 1-0.05; % Minimum cosine similarity required (e.g., ≥0.95)

% --- Define x-axis as index positions ---
x = 1:length(d);  % Index vector
y = d;            % Cosine similarity trajectory

% --- Compute first derivative (rate of change) ---
dy_dx = diff(y) ./ diff(x');  % Numerical derivative

% --- Identify points where derivative magnitude is small ---
below_threshold = abs(dy_dx) < threshold;

% Append last index as true (diff reduces length by 1)
below_threshold = [below_threshold; true];

% --- Apply additional condition: cosine similarity must be high ---
valid_indices = below_threshold & y(1:end) > cosine_threshold;

% --- Detect consecutive stable segments using convolution ---
% Count how many consecutive points satisfy the condition
convolution_result = conv(double(valid_indices), ...
                          ones(min_duration,1), 'valid');

% Pad the end (last few indices cannot be fully assessed by convolution)
convolution_result = [convolution_result; ...
                      zeros(min_duration-1,1)];

% --- Find first index satisfying stability for required duration ---
stable_point = find(convolution_result == min_duration, 1);

% --- Optional visualization ---
plotFig = 1;
if plotFig
    
    figure;
    plot(x, y, '-o');      % Plot original trajectory
    hold on;
    
    if ~isempty(stable_point)
        % Mark detected stable point
        stable_point_x = x(stable_point);
        stable_point_y = y(stable_point_x);
        
        plot(stable_point_x, stable_point_y, 'ro', ...
             'MarkerSize', 10, 'MarkerFaceColor', 'r');
        
        legend('Original Data', 'Stable Point', ...
               'Location', 'southeast');
        
        title(['Stable Point at Index: ', ...
               num2str(stable_point_x)]);
    else
        disp('No Stable point found within the given criteria.');
    end
    
    hold off;
    xlabel('Index');
    ylabel('Cosine Similarity');
end

end