function  stable_point = findStabPoint_cosine(d)
% Define thresholds and duration for detecting continuous low derivatives
threshold = 1e-2; % Threshold for derivative to consider it "flat"
min_duration = 5; % Minimum number of consecutive points below the threshold
cosine_threshold = 1-0.05; % cosine dissimilarity threshold for a point to be eligible for considering

x = 1:length(d); % Assuming x represents the indices
y = d; % Replace 'd' with your actual data

% Calculate the first derivative (rate of change)
dy_dx = diff(y) ./ diff(x');

% Find where the first derivative is below the threshold
below_threshold = abs(dy_dx) < threshold;
below_threshold =[below_threshold;true];%append the last index 

% Apply the additional criterion that y must be below 0.1
valid_indices = below_threshold & y(1:end) > cosine_threshold;

% Perform convolution to find segments where the condition is met for at least min_duration points
convolution_result = conv(double(valid_indices), ones(min_duration, 1), 'valid');
convolution_result = [convolution_result;zeros(min_duration-1,1)];%append the last 4 indexes as 0 as last 4 element can not be assessed

% Identify the first index where the condition is satisfied for min_duration consecutive points
stable_point = find(convolution_result == min_duration, 1);


plotFig = 1;
if(plotFig)
% Plot the original data and mark the Stable point
figure;
plot(x, y, '-o');
hold on;
% Mark the Stable point on the plot if found
if ~isempty(stable_point)
    stable_point_x = x(stable_point ); % Adjust for smoothing and derivative
    stable_point_y = y(stable_point_x);
    plot(stable_point_x, stable_point_y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    legend('Original Data', 'Stable Point','location','southeast');
    title(['Stable Point at Index: ', num2str(stable_point_x)]);
else
    disp('No Stable point found within the given criteria.');
end

hold off;
xlabel('x');
ylabel('y');
end


end
