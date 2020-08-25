% *************************************************************************
%
% *************************************************************************
function [WETLabsData] = getWETLabsData(eng, WETLabsCalData)
    WETLabsData = [];
  
    [list] = findWETLabsData(eng);
    if ~isempty(list)
        [WETLabsData] = extractWETLabsData(eng, length(eng.elaps_t), eng.columns(list), WETLabsCalData);
    end
end
