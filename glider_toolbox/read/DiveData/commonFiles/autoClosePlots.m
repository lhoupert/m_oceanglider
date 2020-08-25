% ********************************************************************************
% GUI callback
% ********************************************************************************
function autoClosePlots()
    [plotFigureHandles] = getPlotFigureHandles();

    if (isempty(plotFigureHandles) == 0)
        close(plotFigureHandles);
    end
end
