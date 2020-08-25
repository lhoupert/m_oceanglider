% ********************************************************************************
% Retreive handles for all plots
% ********************************************************************************
function [figureHandleList, figureNameList] = getPlotFigureHandles()
    figureHandleList = [];
    
    figure_handles = findobj('Type','figure');
    figureCount = 0;

    for i=1:length(figure_handles)
        figInfo=get(figure_handles(i));
        %if isempty(strfind(figInfo.Name, 'Plot_')) == 0
        if isempty(strfind(figInfo.UserData, 'Plot_')) == 0
            figureHandleList = [figureHandleList figure_handles(i)];
            figureCount = figureCount + 1;
            figureNameList{figureCount} = figInfo.Name;
        end
    end
    
    if (figureCount == 0)
        figureNameList = [];
    else
        figureNameList = char(figureNameList');
    end
