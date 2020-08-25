% ********************************************************************************
% ********************************************************************************
function [scaledScatteringData]=scaleScatteringData(scatteringData, scaleFactor, darkCounts)
    scaledScatteringData = scaleFactor * (scatteringData - darkCounts);
end
