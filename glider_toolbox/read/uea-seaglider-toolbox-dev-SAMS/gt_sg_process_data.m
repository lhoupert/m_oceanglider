function data = gt_sg_process_data(data)
%% Check inputs
if isempty(data) || ~isfield(data,'log') || ~isfield(data,'eng')
    gt_sg_sub_echo({'Input data does not contain parsed .eng and .log data.','ABORTING.'});
    return;
elseif ~isstruct(data.log) || ~isstruct(data.eng)
    gt_sg_sub_echo({'Input data does not contain parsed .eng and .log data.','ABORTING.'});
    return;
end

%% Set processing settings from gloder toolbox defaults and sg_calib_constants
% Load processing constants from gt_sg_toolbox folder or current directory
if ~isfield(data,'gt_sg_settings')
    data.gt_sg_settings = gt_sg_import_constants('gt_sg_settings.m');
else
    gt_sg_sub_echo({'gt_sg_settings already loaded in structure.','Not updating constants and keeping pre-loaded data.'});
end

% Merge sg_calib_constants into gt_sg_settings. sg_calib_constants
% overwrite to allow setting specific settings for each glider individually
% if desired.
% Remove from gt_sg_settings fields that will be overwritten by sg_calib_const
data.gt_sg_settings = rmfield(data.gt_sg_settings, intersect(fieldnames(data.gt_sg_settings), fieldnames(data.sg_calib_const)));
% Merge sg_calib_const into gt_sg_settings
data.gt_sg_settings = cell2struct([struct2cell(data.gt_sg_settings); struct2cell(data.sg_calib_const)], [fieldnames(data.gt_sg_settings); fieldnames(data.sg_calib_const)], 1);

%% General context calculations (date, time, gps)
% Fix density units
if data.gt_sg_settings.rho0 >= 1 && data.gt_sg_settings.rho0 <= 2
    data.gt_sg_settings.rho0 =  data.gt_sg_settings.rho0 .* 1000;
elseif data.gt_sg_settings.rho0 <= 1000 || data.gt_sg_settings.rho0 >= 1045
    gt_sg_sub_echo({'Odd rho0 in sg_calib_constants.m. Units should be kg.m-3.'});
end


%--------------------------------------------------------------------------------------	

% Identify dives to be processed
processedDives = intersect([data.log.dive],[data.eng.dive]);

% Create date array
data.date = nan(size(data.log));
data.date(processedDives) = arrayfun(@(x) datenum(x.start(3)+1900,x.start(1),x.start(2),x.start(4),x.start(5),x.start(6)), data.log(processedDives));

% Load GPS data
[data.gps_postpreviousdive,data.gps_predive,data.gps_postdive] = gt_sg_sub_parse_GPS(data.log);

%% Process sensor names
    engNames = fieldnames(data.eng)';
if ~isfield(data.gt_sg_settings,'sensorList')
    % List registered sensors
    % Start GUI to get user input of matching devices, eng headers and
    % scripts
    gt_sg_sub_queryColumnsHeaders;
end
while ~any(strcmp(data.gt_sg_settings.sensorList(:,3),'gt_sg_sensors_sbect.m'));
    h = errordlg('You must process Seabird CT data.');
    uiwait(h);
    gt_sg_sub_queryColumnsHeaders;
end
% Get all subheaders in eng file for each sensor (ie. cond and temp for the
% seabird sensors).
for istep=1:numel(data.gt_sg_settings.sensorList(:,1))
    if ~any(strcmp(data.gt_sg_settings.sensorList{istep,2},{'Not in .eng file','Select .eng header name.'}))
        data.gt_sg_settings.sensorList{istep,4} = cellfun(@(x) [data.gt_sg_settings.sensorList{istep,2},x],...
            gt_sg_sub_recurse_fieldnames(data.eng(processedDives(1)).(data.gt_sg_settings.sensorList{istep,2})),'Uni',0);
    end
end

% Write output
sensorMessage = {'Matching log file sensor timings to eng file columns:'};
for istep = 1:numel(data.gt_sg_settings.sensorList(:,1))
    if ~isempty(data.gt_sg_settings.sensorList{istep,1}) && ~any(strcmp(data.gt_sg_settings.sensorList{istep,2},{'Not in .eng file','Select .eng header name.'}))
        sensorMessage{numel(sensorMessage)+1} = ['Log $SENSORS: "',data.gt_sg_settings.sensorList{istep,1},'" to eng header: "',data.gt_sg_settings.sensorList{istep,2},'"'];
    end
end

processVars = find(cellfun(@(x) ~any(strcmp(x,{'Do not process','Select processing script.'})),data.gt_sg_settings.sensorList(:,3)));
for istep = 1:numel(processVars)
    sensorMessage{numel(sensorMessage)+1} = ['Processing ',data.gt_sg_settings.sensorList{processVars(istep),1},' through ',data.gt_sg_settings.sensorList{processVars(istep),3}];
    if any(strcmp(data.gt_sg_settings.sensorList{processVars(istep),2},{'Not in .eng file','Select .eng header name.'}))
       gt_sg_sub_echo({['WARNING: No .eng header selected for processing ',data.gt_sg_settings.sensorList{processVars(istep),1},' data.'],'Assuming this is a logger device with external data storage.'});
    end
end
gt_sg_sub_echo(sensorMessage);

%% Correct pressure and create time array
% The glider converts pressure readings to PSI onboard using Z = PSI * 0.685.
% So we backcalculate to PSI, then convert to dbar.
sgdpth2psi = 1/0.685;
psi2bar = 1.0/14.5037738007;
data.hydrography = cell2struct(arrayfun(@(x) x.depth ./100 .*sgdpth2psi .*psi2bar .*10, data.eng,'Uni',0), 'pressure',2); % dbar
[data.hydrography.time] = deal([]);
tmp = deal(arrayfun(@(x) (datenum(x.start(3)+1900,x.start(1),x.start(2)) + x.elaps_t_0000/60/60/24),data.eng([data.eng.dive]),'Uni',0)); %TODO: verify why there is a difference in duration between hydro.time and eng.elaps_t
[data.hydrography([data.eng.dive]).time] = tmp{:}; % matlab datenum
tmp = arrayfun(@(x,y) -gsw_z_from_p(x.pressure,y),data.hydrography([data.eng.dive]),data.gps_postdive([data.eng.dive],1),'Uni',0);
[data.hydrography([data.eng.dive]).depth] = tmp{:}; % m
% TODO : the two-liner with tmp then deal(tmp) can probably be reduced to a
% one-line with subsref or something, can't figure it out though.
eng_dives = [data.eng.dive];

[max_pressure,max_pressure_index] = arrayfun(@(x) max(x.depth),data.hydrography(eng_dives));%,'un',0); % add the option 'Uni',0 to arrayfun 
% by L. Houpert 25/01/2016 due to an error when processing sg602 

max_pressure = num2cell(max_pressure); max_pressure_index = num2cell(max_pressure_index);
[data.hydrography(eng_dives).max_pressure] = max_pressure{:};
[data.hydrography(eng_dives).max_pressure_index] = max_pressure_index{:};
clear max_pressure max_pressure_index

%% Create offset time and pressure arrays for each sensor
%% TODO: Calculate the time and pressure offsets for each sensor
% keyboard
%gt_sg_sub_echo('Beginning sensor time offset calculations.');
%
% % Get total sensor run time from all dives
% data.gt_sg_settings.sensorList{4,:} = nansum(vertcat(data.log.SENSOR_SECS));
%
% for istep = 1:length(data.gt_sg_settings.sensorList)
%     if ~isstruct(data.eng(processedDives(1)).(data.gt_sg_settings.sensorList{2,istep}))
%         data.gt_sg_settings.sensorList{5,istep} = sum(~isnan([data.eng.(data.gt_sg_settings.sensorList{2,istep})]));
%     end
% end
% %sensorNumSamples =
% TODO: Check on deployments it t_elaps_0 is always 5 ish or if there are
% irregular delays before between sampling...

%% Begin scientific data processing
[data.hydrography.flagged_samples] = deal([]);
[data.hydrography.used_samples] = deal([]);

%% Regress flight model
% First naive approximation of flight base on observed vertical
% velocity (pressure), vehicle pitch and hydrodynamic constants.
data = gt_sg_sub_slope(data);

% First rough pass of T&S data,
data = gt_sg_sensors_sbect(data);

% Verify volmax and hydrodynamic parameters
% Calculate flight speed
data = gt_sg_sub_flightmodelregression(data);

%% Clean pass of T&S data
data = gt_sg_sensors_sbect(data);

%% Calculate vertical velocities, depth-averaged currents and interpolated GPS positions
if isfield(data.flight,'model_spd')
    gt_sg_sub_echo({'Estimating up/downwelling, depth-average currents and underwater coordinates.'});
    data = gt_sg_sub_currents(data);
else
    gt_sg_sub_echo({'No flight model output available.',...
        'Cannot estimate up/downwelling, depth-averaged currents or',...
        'coordinates of the glider underwater.',...
        'WARNING: Flight model data should be present, something has gone wrong'});
end

%% TODO: Estimate piloting parameters (C_VBD, C_ROLL, C_PITCH)

%% Process other sensor data
for istep = 1:length(data.gt_sg_settings.sensorList(:,1))
    % Identify routines that need to be run
    switch data.gt_sg_settings.sensorList{istep,3}
        case 'Do not process'
            % Do nothing
        case 'gt_sg_sensors_sbect.m'
            % Do nothing as it should already have been processed, if they
            % have two seabird CT sails, I can't help them, they're just
            % being difficult.
        case 'gt_sg_sensors_wetlabs.m'
            gt_sg_sub_echo({['Beginning processing on ',data.gt_sg_settings.sensorList{istep,1},' through the ',data.gt_sg_settings.sensorList{istep,3},' routine.']});
            data = gt_sg_sensors_wetlabs(data);
        case 'gt_sg_sensors_aa4330.m'
            gt_sg_sub_echo({['Beginning processing on ',data.gt_sg_settings.sensorList{istep,1},' through the ',data.gt_sg_settings.sensorList{istep,3},' routine.']});
            data = gt_sg_sensors_aa4330(data);
        case 'gt_sg_sensors_qsppar.m'
            gt_sg_sub_echo({['Beginning processing on ',data.gt_sg_settings.sensorList{istep,1},' through the ',data.gt_sg_settings.sensorList{istep,3},' routine.']});
            data = gt_sg_sensors_qsppar(data);
    end
end

%% Finish up
data = orderfields(data);
data.gt_sg_log = gt_sg_sub_echo('FINISHED PROCESSING OF SEAGLIDER DATA');


%%%%%%%%%%%%%%% END OF MAIN ROUTINES  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% START OF SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%

    function gt_sg_sub_queryColumnsHeaders
        f = figure('Position',[100 100 650 300],'Visible','off');
        set(f,'MenuBar','none')
        
        deviceNames = strsplit(data.log(processedDives(1)).SENSORS{:},',');
        [~,deviceNames_ind] = setxor(deviceNames,{'nil','empty'});
        deviceNames = deviceNames(deviceNames_ind);        
        sensorFunctions = dir([fileparts(which('gt_sg_process_data.m')),filesep,'gt_sg_sensors_*.m']);
        sensorFunctions = {'Do not process',sensorFunctions.name};
        
        % Column names and column format
        columnname = {'Device Name','eng Header','Process?'};
        columnformat = {'char',{'Not in .eng file',engNames{:}},sensorFunctions};
        
        % Define the data
        UItable = deviceNames';
%         for sens_step=1:numel(UItable)
%             UItable{sens_step,2} = gt_sg_sub_find_nearest_string(UItable{sens_step,1},engNames);
%         end
        [UItable{:,2}] = deal('Select .eng header name.');
        [UItable{:,3}] = deal('Select processing script.');
        
        % Create the uitable
        t = uitable('Data', UItable,...
            'ColumnName', columnname,...
            'ColumnFormat', columnformat,...
            'ColumnEditable', [false true true],...
            'ColumnWidth',{90 200 350},...
            'RowName',[],...
            'Units','normalized','Position',[0 0.2 1 0.8]);
        
        pb = uicontrol(f,'Style','pushbutton','String','Save & Close',...
            'Units','normalized','Position',[0 0 1 0.2],...
            'CallBack',@button_clickcallback);
        
        function button_clickcallback(button_data,event_data)
            data.gt_sg_settings.sensorList = get(t,'Data');
            close(gcf);
            delete(gcf);
        end
        
        movegui(f,'center')
        set(f,'Visible','on');
        
        uiwait(f)
    end
end
