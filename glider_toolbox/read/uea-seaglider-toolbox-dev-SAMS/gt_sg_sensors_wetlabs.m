function data = gt_sg_sensors_wetlabs(data)
% -=-=-=-=-=-=-=-=-=-=-=-=-=-
% COMPUTE
%  WETLABS SENSOR
% -=-=-=-=-=-=-=-=-=-=-=-=-=-
%
% calibs = sg_calib_constants struct
% eng    = eng struct
% wps    = struct to which optics data will be appended (can be eng, or a new struct)
%

%% Check for calibration data
if ~isfield(data.gt_sg_settings,'WETLabsCalData')
    gt_sg_sub_echo({'No Wetlabs calibration data. Cannot continue with processing'});
    return;
end

%% Identify sensor(s) that need to be processed
sensorPosition = find(strcmp(data.gt_sg_settings.sensorList(:,3),'gt_sg_sensors_wetlabs.m'));

%% List dives to process, no need for log files or currents for processing,
%% so all eng files should work
dives = [data.eng.dive];

%% For each sensor need to be processed...
for sPos = [sensorPosition]
    % List all subfields containing the "ref" substring
    subFields = cellfun(@(x) strsplit(x,'.'), data.gt_sg_settings.sensorList{sPos,4},'Uni',0);
    ref_ind = find(cellfun(@(x) ~isempty(strfind(x{2},'ref')),subFields));
    
    for rstep = [ref_ind]
        %Check subfields with "ref" substring for reference wavelength and set
        %names
        switch mode(data.eng(dives(1)).(subFields{rstep}{1}).(subFields{rstep}{2}))
            case 695
                Name = 'Chlorophyll';
            case 460
                Name = 'CDOM';
            case 530
                Name = 'Uranine';
            case 680
                Name = 'Phycocyanin';
            case 570
                Name = 'Phycoerythrin';
            otherwise
                Name = ['Scatter_',num2str(mode(data.eng(dives(1)).(subFields{rstep}{1}).(subFields{rstep}{2})))];
        end
        % Identify name of data column
        Sig = strrep(subFields{rstep}{2}, 'ref', 'sig');

        % Process data if calibration data exists
        if isfield(data.gt_sg_settings.WETLabsCalData,Name)
            if any(~isfield(data.gt_sg_settings.WETLabsCalData.(Name),{'scaleFactor','darkCounts'}))
                gt_sg_sub_echo({['WARNING: No calibration info for the ',Name,' Wetlabs data.']});
                data.gt_sg_settings.WETLabsCalData.(Name).scaleFactor = NaN;
                data.gt_sg_settings.WETLabsCalData.(Name).darkCounts = NaN;
            else
                gt_sg_sub_echo({['Calculating Wetlabs ',Name,' data.']});
                for dstep = [dives]
                    data.hydrography(dstep).(Name) = data.gt_sg_settings.WETLabsCalData.(Name).scaleFactor .* (data.eng(dstep).(subFields{rstep}{1}).(Sig) - data.gt_sg_settings.WETLabsCalData.(Name).darkCounts);
                end
            end
        else
            gt_sg_sub_echo({['WARNING: No calibration info for the ',Name,' Wetlabs data.']});
        end %isfield(wetlabs cal data)
        
    end % for rstep (ie. for each reference subfield)
end % for sPos (ie. for each wetlabs sensor)