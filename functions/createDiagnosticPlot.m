function createDiagnosticPlot(s, originalkVal, kvalMin, indexStart, indexEnd, index)
% createDiagnosticPlot Generates diagnostic plots for k-value comparison and model fitting
%
% Syntax:
%   createDiagnosticPlot(s, originalkVal, kvalMin, indexStart, indexEnd, index)
%
% Description:
%   This function creates diagnostic plots to visualize the changes in k-values and
%   the fitting results of the model. It plots the difference between the original 
%   and the current best k-values, the current best k-values, and the fit plot 
%   over a specified window of time.
%
% Input Arguments:
%   s           - The simulation object containing data and configuration.
%   originalkVal - The original k-values before fitting.
%   kvalMin     - The optimized k-values after fitting.
%   indexStart  - The starting index of the fit window.
%   indexEnd    - The ending index of the fit window.
%   index       - The index of the current fit window for file naming.

    % Initialize the figure and layout
    figure('Visible', 'off');
    %figure;
    tl = tiledlayout('flow');

    % Plot the change in k-value from original to current best
    nexttile
    bar(originalkVal - kvalMin);
    title('Delta k-value');

    % Plot current best k-value
    nexttile
    hold on;
    bar(kvalMin);
    hold off;
    title('Current best k-value');

    % Generate the fit plot based on specified indices
    s.fitPlot(true, [indexStart, indexEnd], false);

    % Set the overall title for the tiled layout
    title(tl, sprintf('Last window from: %s to: %s', datestr(s.dateNumList(indexStart, :)), datestr(s.dateNumList(indexEnd, :))));

    % Save the plot if saveDiagnosticFlag is true
    if s.saveDiagnosticFlag
        saveas(gcf, fullfile(s.folder, s.diagnosticFolder, sprintf('fitPlotIndex_%d.png', index)));
    end
end
