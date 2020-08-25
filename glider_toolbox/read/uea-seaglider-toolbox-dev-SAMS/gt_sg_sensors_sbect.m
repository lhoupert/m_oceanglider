function data = gt_sg_sensors_sbect(data)

%% Output messages
if ~isfield(data.hydrography,'temp')
    % List of necessary constants
    sbectCalibrationConstants = {...
        't_g',...
        't_h',...
        't_i',...
        't_j',...
        'c_g',...
        'c_h',...
        'c_i',...
        'c_j',...
        'cpcor',...
        'ctcor'};
    % Verify necessary calibration coeffcients are present
    if any(~isfield(data.gt_sg_settings,sbectCalibrationConstants))
        gt_sg_sub_echo({'Missing the following SBE CT calibration constants in sg_calib_constants.m:',...
            sbectCalibrationConstants{~isfield(data.gt_sg_settings,sbectCalibrationConstants)}});
        return
    end
    gt_sg_sub_echo({'Beginning first-pass calculation of temperature and salinity before the flight regression.','Using naive velocity estimate from glide_spd for data correction.'});
else
    if isfield(data.flight,'model_spd')
        gt_sg_sub_echo({'Beginning calculations based on improved flight model data to correct temperature and salinity profiles.'});
    else
        gt_sg_sub_echo({'No flight model output available, using glide model instead.',...
            'Doing naive calculation of temperature and salinity using velocity estimate from glide_spd.',...
            'WARNING: Flight model data should be present, something has gone wrong'});
    end
end

%% Load general coefficients.
% Temperature lag
if isfield(data.gt_sg_settings,'sbe_temp_tau')
    tempTau =  data.gt_sg_settings.sbe_temp_tau;
else
    tempTau =  0.6; % From UW/iRobot code back in the day (circa. 2010)
end

% Parameters for conductivity cell thermal lag
%% IF EXIST FLIGHT
%% TODO: ADD BLOCK FOR ESTIMATE OF FLOW PARAMETERS
if isfield(data.gt_sg_settings,'sbe_cond_coefficients')
    a_offset = data.gt_sg_settings.sbe_cond_coefficients(1);
    a_slope  = data.gt_sg_settings.sbe_cond_coefficients(2);
    t_offset = data.gt_sg_settings.sbe_cond_coefficients(3);
    t_slope  = data.gt_sg_settings.sbe_cond_coefficients(4);
else
    % These values are proposed in Morison 1994
    a_offset = 0.0135;
    a_slope  = 0.0264;
    t_offset = 7.1499;
    t_slope  = 2.7858;
    % % TODO: These values are proposed in Garau 2011, but the SOCIB toolbox is
    % % still using the Garau parameters as a default. I suspect because they
    % % haven't updated the default parameters because they recalculate
    % % everytime. Would be worth asking Tomeu at some point.
    % a_offset = 0.1328;
    % a_slope  = 0.0208;
    % t_offset = 9.7492;
    % t_slope  = 4.6128;
end

% Check for offsets in sg_calib_constants.m
if isfield(data.gt_sg_settings,'sbe_temp_offset')
    tempOffset =  data.gt_sg_settings.sbe_temp_offset;
    gt_sg_sub_echo(['Applying a temperature offset of ',num2str(tempOffset),' ?C']);
else
    tempOffset = 0;
end
if isfield(data.gt_sg_settings,'sbe_cond_offset')
    condOffset =  data.gt_sg_settings.sbe_cond_offset;
    gt_sg_sub_echo(['Applying a conductivity offset of ',num2str(condOffset),' mS/cm']);
else
    condOffset = 0;
end
if isfield(data.gt_sg_settings,'sbe_sal_offset')
    salOffset =  data.gt_sg_settings.sbe_sal_offset;
    gt_sg_sub_echo(['Applying a salinity offset of ',num2str(salOffset),' (PSS-78)']);
else
    salOffset = 0;
end

% Identify correct eng field to use
sbectInd = find(strcmp(data.gt_sg_settings.sensorList(:,3),'gt_sg_sensors_sbect.m'));
[~, condInd] = gt_sg_sub_find_nearest_string('condFreq',data.gt_sg_settings.sensorList{sbectInd,4});
[~, tempInd] = gt_sg_sub_find_nearest_string('tempFreq',data.gt_sg_settings.sensorList{sbectInd,4});

% Initialise temp and salinity arrays
[data.hydrography.temp] = deal([]);
[data.hydrography.cons_temp] = deal([]);
[data.hydrography.salinity] = deal([]);
[data.hydrography.abs_salinity] = deal([]);
[data.hydrography.sigma0] = deal([]);

%% Start calculations
% The data.log.dive only exists for processed dives.
for istep = [data.eng.dive]
    %% TODO: Extract corrected time and pressure stamps, standard one. And verify it improves things.
    [data.hydrography(istep).temp, data.hydrography(istep).conductivity] = ...
        gt_sg_sensors_sub_sbect(...
        eval(['data.eng(istep).',data.gt_sg_settings.sensorList{sbectInd,4}{condInd}]),...
        eval(['data.eng(istep).',data.gt_sg_settings.sensorList{sbectInd,4}{tempInd}]),...
        data.hydrography(istep).time,data.hydrography(istep).pressure);
end

%% END OF MAIN FUNCTION %%
%% BEGINNIGN OF SUBFUNCTIONS %%

    function [tempCorr, condCorr] = gt_sg_sensors_sub_sbect(condFreq,tempFreq,time,pressure)
        
        % Initialise flagging mask. (1 = keep, 0 = bad)
        flagBool = ones(size(condFreq+tempFreq));
        
        % Flag values outside of set range from sg_calib_constants.m
        if  isfield(data.gt_sg_settings,'sbe_cond_freq_min')...
                && isfield(data.gt_sg_settings,'sbe_cond_freq_max')
            flagBool = flagBool & condFreq > data.gt_sg_settings.sbe_cond_freq_min*1000 ...
                & condFreq < data.gt_sg_settings.sbe_cond_freq_max*1000;
        end
        if  isfield(data.gt_sg_settings,'sbe_temp_freq_min')...
                && isfield(data.gt_sg_settings,'sbe_temp_freq_max')
            flagBool = flagBool & tempFreq > data.gt_sg_settings.sbe_temp_freq_min ...
                & tempFreq < data.gt_sg_settings.sbe_temp_freq_max;
        end
        
        % Get glider flight speed for cell flushing
        if isfield(data.flight,'model_spd')
            gliderSpeed = data.flight(istep).model_spd;
            if isempty(gliderSpeed)
                gliderSpeed = data.flight(istep).glide_spd; %% TODO: cm.s to m.s? Check Garau for proper units.
                flagBool(:) = 0;
                gt_sg_sub_echo({'Flight model speed is empty, using glide model instead.',...
                    'Naive calculation of temperature and salinity using velocity estimate from the slope model.',...
                    'WARNING: Flight model data should be present, something has gone wrong'});
            end
        else
            gliderSpeed = data.flight(istep).glide_spd; %% TODO: cm.s to m.s? Check Garau for proper units.
            %flagBool(:) = 0;
        end
        
        flagBool(gliderSpeed < 0.001) = 0;
        gliderSpeed(gliderSpeed < 0.001) = 0.001;
        
        % Flag glider speeds outside of set range from gt_sg_settings.m
        if  isfield(data.gt_sg_settings,'sbe_flush_speed_min')...
                && isfield(data.gt_sg_settings,'sbe_flush_speed_max')
            flagBool = flagBool & gliderSpeed > data.gt_sg_settings.sbe_flush_speed_min ...
                & condFreq < data.gt_sg_settings.sbe_flush_speed_max;
        end
        
        % Check for finite and imaginary numbers (occasionally occurs)
        flagBool = flagBool & isfinite(tempFreq+condFreq) & isreal(tempFreq+condFreq);
        
        % First estimate of temperature as per Seabird manual and UW
        % original code.
        tempPrelim = log(1000./tempFreq) ;
        tempPrelim = (1 ./ ( data.gt_sg_settings.t_g + tempPrelim.*(data.gt_sg_settings.t_h + tempPrelim.*(data.gt_sg_settings.t_i + tempPrelim.*data.gt_sg_settings.t_j)))) -273.15;
        
        % Correct for first order lag
        tempPrelim = tempPrelim + tempTau*gt_sg_sub_diff({tempPrelim,time*86400});
        
        % Apply temperature offset
        tempPrelim = tempPrelim + tempOffset;
        
        % Interpolate over missing values (?) and accept as final temperature
        %% TODO: implement interpolation nans or leaving as is as an option 
        tempCorr = tempPrelim;
        
        % First estimate of conductivity
        condPrelim = condFreq/1000;
        condPrelim = (data.gt_sg_settings.c_g + ...
            condPrelim.*condPrelim.*(data.gt_sg_settings.c_h + condPrelim.*(data.gt_sg_settings.c_i + condPrelim.*data.gt_sg_settings.c_j))) ./ ...
            ( 10.*( 1.0 + data.gt_sg_settings.ctcor*tempCorr + data.gt_sg_settings.cpcor*pressure./100 ));
        
        % Correct for thermal inertia of the CT cell
        % The relative coefficient between the flow speed inside and outside
        % of the conductivity cell. Here are three choices for first three orders polynomial.
        speedFactorPols = [0.00, 0.00, 0.40;  % 0th order degree
            0.00, 0.03, 1.15;  % 1st order degree - 3rd element now 1.15 from 0.45 as corrected by Gerd Krahmann.
            1.58, 1.15, 0.70]; % 2nd order degree
        
        selectedDegree = 1; % First order approximation, second row of the matrix
        speedFactor = polyval(speedFactorPols(selectedDegree+1, :), gliderSpeed);
        flowSpeed = abs(speedFactor .* gliderSpeed) + eps; % Avoid division by zero
        
        % The alpha and tau parameters, as suggested in the reference paper,
        % depend on the flow with tne next expressions
        alpha = a_offset + a_slope ./      flowSpeed ;
        tau   =   t_offset +   t_slope ./ sqrt(flowSpeed);
        
        % Relation between a and b coefficients with respect to alpha and tau
        freq = 1 ./ gt_sg_sub_diff(time*86400);
        coefa = 4 .* freq .* alpha .* tau ./ (1 + 4 .* freq .* tau);
        coefb = 1 - 2 .* coefa ./ alpha;
        
        % Sensitivity of conductivity with respect to temperature,
        % approximation suggested by SeaBird: SBE Data Processing Users Manual
        % at Section 6: Data Processing Modules, Cell Thermal Mass
        % Software Release 7.16a and later. Date: 01/18/08
        % dCdT = 0.1 .* (1 + 0.006 .* (temp - 20));
        dCdT = 0.088 + 0.0006 * tempCorr;
        dT = gt_sg_sub_diff(tempCorr);
        
        % Recursive processing of the corrections
        condCorrection = zeros(size(condPrelim));
        tempCorrection = zeros(size(condPrelim));
        
        for depthLevel = 1:numel(flowSpeed)-1,
            % Compute corrections for next depth level
            condCorrection(depthLevel+1) = ...
                - coefb(depthLevel) .* condCorrection(depthLevel) + ...
                coefa(depthLevel) .* dCdT(depthLevel) .* dT(depthLevel);
            tempCorrection(depthLevel+1) = ...
                - coefb(depthLevel) .* tempCorrection(depthLevel) + ...
                coefa(depthLevel) .* dT(depthLevel);
        end
            
            % Apply corrections and save them as fields in the output struct
            condPrelim = condPrelim + condCorrection;
            % tempSBEcell  = tempCorr - tempCorrection;
            
            % Apply conductivity offset
            condPrelim = condPrelim + condOffset;
            
            % Interpolate over missing values (?) and accept as final
            % conductivity
            condCorr = condPrelim.*10;
            %% TODO: figure out this *10 and see if it needs to be applied way
            % before instead
            
            %% TODO: implement interpolation nans or leaving as is as an option
            
            % Calculate practical salinity, and apply any offsets that need to
            % be applied.
            data.hydrography(istep).salinity = gsw_SP_from_C(condCorr,tempCorr,pressure) + salOffset;
            
            % Calculate TEOS-10 variables.
            data.hydrography(istep).abs_salinity = gsw_SA_from_SP(data.hydrography(istep).salinity,pressure,data.gps_postdive(istep,2),data.gps_postdive(istep,1));
            data.hydrography(istep).cons_temp = gsw_CT_from_t(data.hydrography(istep).abs_salinity,tempCorr,pressure);
            data.hydrography(istep).sigma0 = gsw_rho(data.hydrography(istep).abs_salinity,data.hydrography(istep).cons_temp,0);
            data.hydrography(istep).rho = gsw_rho(data.hydrography(istep).abs_salinity,data.hydrography(istep).cons_temp,pressure);
            
            % Export flagging, I don't like these two arrays, but keep them for
            % backwards compatability?
            data.hydrography(istep).used_samples = flagBool;
            data.hydrography(istep).flagged_samples = ~data.hydrography(istep).used_samples;
        end
    end