% ****************************************************************************************************
% The purpose of this function is to display an error message in the center of the plot when
% calibration data is not available for the sensor.
% ****************************************************************************************************
function displayNoDataMessage(strMessage)
    grid off;

    set(gca, 'XLim', [1 100], 'YLim', [0 100]);
    set(gca, 'ytick', []);
    set(gca, 'xtick', []);

    v = axis;
    xtxt = (v(2) - v(1)) / 2.0;
    ytxt = (v(4) - v(3)) / 2.0;
    text(xtxt, ytxt, strMessage, 'HorizontalAlignment', 'center', 'FontSize', 12);
end
