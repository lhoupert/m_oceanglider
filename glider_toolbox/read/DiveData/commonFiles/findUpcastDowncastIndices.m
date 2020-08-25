% **********************************************************************************************************************************
% Function: findUpcastDowncastIndices
%
% Purpose: the purpose of this function is to find the indices into the data which correspond to downcast and upcast data
%
% Description: the function filters out bad data that can cause data from one cast to be confused for data from the other cast.
%
% Author           Date             Remarks
% ----------------------------------------------------------------------------------------------------------------------------------
% John E. Stranzl	04-Oct-2011     INITIAL INCORPORATION
% 
% **********************************************************************************************************************************
function [iwn, iwp] = findUpcastDowncastIndices(depth, pitchCtl)
    % FIND INDEXES FOR DOWNCAST/DESECENT DATA POINTS
    iwn = 1:find(depth==max(depth));

    % FIND INDEXES FOR UPCAST/ASCENT DATA POINTS
    iwp = find(depth==max(depth)):length(depth);

    % 1. SEARCH FOR NEGATIVE PITCH CONTROL VALUES AFTER MAX DEPTH IS IDENTIFIED
    % 2. REMOVE ANY DATA ON THE UPCAST WHERE THE PITCHCTL IS NEGATIVE. THIS IS A CARTE BLANCHE APPROACH FOR NOW. IT WILL REMOVE ANY
    %    DATA ON THE UPCAST WITH A NEGATIVE PITCH CONTROL VALUE ASSOCIATED WITH IT.
    %
    % NOTE: THIS IS ONLY AN ISSUE WITH ALI SENSORS SINCE THEY LOG DATA INDEPENDENT OF THE DATA LOGGED IN THE ENGINEERING FILE.
    if ~isempty(iwp)
        pitchPos = find(pitchCtl(iwp)<0);
        if ~isempty(pitchPos)
            iwp(pitchPos) = [];
        end
    end
end
