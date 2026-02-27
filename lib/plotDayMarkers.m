function plotDayMarkers(dayID)
%% ======================================================================
% DayID strip plotting: marks day boundaries and labels day segments
%
% Created in 2026 by Zhuying Chen (zhuychen@unimelb.edu.au)
% Released under the CC-BY-NC-4.0 License
% http://creativecommons.org/licenses/by-nc/4.0/
%% ======================================================================
    nEps = numel(dayID);

    y1 = 0.002;
    y2 = 0.095;

    hold on

    lastMark = 1;

    for j = 2:nEps
        if dayID(j-1) ~= dayID(j)
            line([j j], [y1 y2], 'Color', 'k', 'LineWidth', 1);

            text((lastMark + (j-1))/2, y2, sprintf('%d', dayID(j-1)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                 'Color','k');

            lastMark = j;
        end
    end

    % last segment label
    text((lastMark + nEps)/2, y2, sprintf('%d', dayID(end)), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'Color','k');

    % boundaries
    line([0.5 0.5],       [y1 y2], 'Color','k', 'LineWidth', 1);
    line([nEps nEps], [y1 y2], 'Color','k', 'LineWidth', 1);

    text( nEps/2, 0.15, sprintf('Day'), ...
                'horizontalalign', 'right', 'verticalalign', 'middle', 'color','k')


    hold off
end