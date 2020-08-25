% *************************************************************************
%
% *************************************************************************
function [PARCalData] = getPARCalData(results)

    %if isfield(results.sg_calib_constants, 'PARCalData')
    %    PARCalData = results.sg_calib_constants.PARCalData;
    %    return;
    %end
    
    PARCalData = [];
    list = fieldnames(results.sg_calib_constants);
    for i = 1:length(list)
        if (strncmp(list(i), 'PARCalData_', 11) == 1)
            [prefix, remain] = strtok(list(i), '_');
            fieldStr = char(remain);
            fieldStr = fieldStr(2:end);
            if (strncmp(fieldStr, 'manufacturer', 12) == 1)
                PARCalData.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'serialNumber', 12) == 1)
                PARCalData.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'calDate', 7) == 1)
                PARCalData.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'darkOffset', 10) == 1)
                PARCalData.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'scaleFactor', 11) == 1)
                PARCalData.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            else
                display(sprintf('Unexpected field in PAR cal data (%s%s)', char(prefix), char(remain)));
            end
        end
    end
    
end
