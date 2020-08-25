% ****************************************************************************************************
%
% ****************************************************************************************************
function [newWaitbarIndex]=incrementWaitbar(hWaitbar, waitbarIndex, numWaitStages, strMsg)

    waitbar(waitbarIndex/numWaitStages, hWaitbar, strMsg);
    
    newWaitbarIndex = waitbarIndex + 1;
end
