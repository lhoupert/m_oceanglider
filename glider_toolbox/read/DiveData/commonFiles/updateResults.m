% ****************************************************************************************************
% Output the results messages to the log view in the GUI
% ****************************************************************************************************
function updateResults(handles, strResults)
    strResultsData = strtrim(get(handles.textResults, 'String'));
    
    strNewResultsData = [];
    [m, n] = size(strResultsData);
    for i=1:m
        strNewResultsData = sprintf('%s\n%s', strNewResultsData, char(strResultsData(i,:)));
    end
    
    strNewResultsData = strtrim(sprintf('%s\n%s\n', strNewResultsData, strResults));
    
    set(handles.textResults, 'String', strNewResultsData);
end
