% ********************************************************************************
%
% ********************************************************************************
function [selectionState] = getPlotSelectionState(handles, plotNumber)
    strEval = sprintf( 'isfield(handles, ''checkboxSelectPlot_%d'')', plotNumber);
    if (eval(strEval) == 1)
        strEval = sprintf( 'get(handles.checkboxSelectPlot_%d, ''Value'');', plotNumber);
        selectionState = eval(strEval);
    else
        selectionState = 1;
    end
end