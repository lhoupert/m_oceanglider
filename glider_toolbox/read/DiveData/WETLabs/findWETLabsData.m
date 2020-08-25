% ********************************************************************************
%
% ********************************************************************************
function [list] = findWETLabsData(eng)
    % SEARCH FOR WET Labs DATA IN THE COLUMN HEADER. ALL WET Labs COLUMN
    % LABELS ARE PREFACED WITH "WL"
    list = find(not(cellfun('isempty', (regexpi(eng.columns, 'WL*')))));
end
