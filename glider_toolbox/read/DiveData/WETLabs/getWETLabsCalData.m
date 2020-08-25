% *************************************************************************
%
% *************************************************************************
function [WETLabsCalData] = getWETLabsCalData(results)

    if isfield(results.sg_calib_constants, 'WETLabsCalData')
        WETLabsCalData = results.sg_calib_constants.WETLabsCalData;
        return;
    end
    
    WETLabsCalData = [];
    list = fieldnames(results.sg_calib_constants);
    for i = 1:length(list)
        if (strncmp(list(i), 'WETLabsCalData_', 15) == 1)
            % New approach - sensor type included in cal fields, e.g.,
            % WETLabsCalData_wlbbfl2_Chlorophyll_<XXX>.
            %
            % Old approach - sensor type not included, e.g.,
            % WETLabsCalData_Chlorophyll_<XXX>
            line = char(list(i));

            if (strncmp(list(i), 'WETLabsCalData_wlbbfl2_', 23) == 1)
                type = 'wlbbfl2';
                fieldStr = line(24:end);
            elseif (strncmp(list(i), 'WETLabsCalData_wlbb2fl_', 23) == 1)
                type = 'wlbb2fl';
                fieldStr = line(24:end);
            elseif (strncmp(list(i), 'WETLabsCalData_wlbb3_', 21) == 1)
                type = 'wlbb3';
                fieldStr = line(22:end);
            elseif (strncmp(list(i), 'WETLabsCalData_wlfl3_', 21) == 1)
                type = 'wlfl3';
                fieldStr = line(22:end);
            else
                % Old approach - puck type not included in calibration
                % fields.  Field is everything after 'WETLabsCalData_'.
                type = 'wl';
                fieldStr = line(16:end);
            end
            
            if (strncmp(fieldStr, 'Chlorophyll_', 12) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Chlorophyll.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'CDOM_', 5) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).CDOM.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Uranine_', 8) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Uranine.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Phycoerythrin_', 14) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Phycoerythrin.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Rhodamine_', 10) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Rhodamine.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Scatter470_', 11) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Scatter470.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Scatter532_', 11) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Scatter532.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Scatter650_', 11) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Scatter650.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Scatter700_', 11) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Scatter700.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'Scatter880_', 11) == 1)
                [prefix, remain] = strtok(fieldStr, '_');
                fieldStr = char(remain);
                fieldStr = fieldStr(2:end);
                WETLabsCalData.(type).Scatter880.(fieldStr) = results.sg_calib_constants.(char(list(i)));
            elseif (strncmp(fieldStr, 'calinfo', 7) == 1)
                WETLabsCalData.(type).calinfo = results.sg_calib_constants.(char(list(i)));
            else
                %display(sprintf('Unexpected field in WETLabs cal data (%s%s)', char(prefix), char(remain)));
                display(sprintf('Unexpected field in WETLabs cal data: %s', char(fieldStr)));
            end
        end
    end
    
end
