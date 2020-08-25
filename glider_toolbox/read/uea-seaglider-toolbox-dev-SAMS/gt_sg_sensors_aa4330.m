function data = gt_sg_sensors_aa4330(data)

%% Identify sensor(s) that need to be processed
sensorPosition = find(strcmp(data.gt_sg_settings.sensorList(:,3),'gt_sg_sensors_aa4330.m'));

%% Identify lag correction coefficients
if isfield(data.gt_sg_settings,'aa4330_tau') && isfield(data.gt_sg_settings,'aa4330_tau_degC')
    tau_DO = [data.gt_sg_settings.aa4330_tau data.gt_sg_settings.aa4330_tau_degC];
    gt_sg_sub_echo({['Using a AA4330 lag correction of ',num2str(tau_DO(1)),' s ',num2str(tau_DO(2)),' s.K^-1 (ref 20degC) as per Johannes Hahn''s thesis (2013).']});
else
    gt_sg_sub_echo({'No AA4330 lag correction coefficients defined in sg_calib_const.m','Not performing any lag correction on the oxygen data','Suggest using aa4330_tau = 14.8 seconds and aa4330_tau_degC = -0.4 s.K^-1 as per Johannes Hahn''s thesis (2013).'})
    tau_DO = [14.8 -0.4];
end

%% Check for missing coefficients
% List of necessary constants
aa4330CalibrationConstants = {...
    'aa4330_TPhaseCoef0',...
    'aa4330_TPhaseCoef1',...
    'aa4330_TPhaseCoef2',...
    'aa4330_TPhaseCoef3',...
    'aa4330_FoilCoefA0',...
    'aa4330_FoilCoefA1',...
    'aa4330_FoilCoefA2',...
    'aa4330_FoilCoefA3',...
    'aa4330_FoilCoefA4',...
    'aa4330_FoilCoefA5',...
    'aa4330_FoilCoefA6',...
    'aa4330_FoilCoefA7',...
    'aa4330_FoilCoefA8',...
    'aa4330_FoilCoefA9',...
    'aa4330_FoilCoefA10',...
    'aa4330_FoilCoefA11',...
    'aa4330_FoilCoefA12',...
    'aa4330_FoilCoefA13',...
    'aa4330_FoilCoefB0',...
    'aa4330_FoilCoefB1',...
    'aa4330_FoilCoefB2',...
    'aa4330_FoilCoefB3',...
    'aa4330_FoilCoefB4',...
    'aa4330_FoilCoefB5',...
    'aa4330_FoilCoefB6',...
    'aa4330_FoilCoefB7',...
    'aa4330_FoilCoefB8',...
    'aa4330_FoilCoefB9',...
    'aa4330_FoilCoefB10',...
    'aa4330_FoilCoefB11',...
    'aa4330_FoilCoefB12',...
    'aa4330_FoilCoefB13'};
% Verify necessary calibration coeffcients are present
if any(~isfield(data.gt_sg_settings,aa4330CalibrationConstants))
    gt_sg_sub_echo({'Missing the following AA4330 calibration constants in sg_calib_constants.m:',...
        aa4330CalibrationConstants{~isfield(data.gt_sg_settings,aa4330CalibrationConstants)}});
    return
end

%% For each sensor needing to be processed...
for sPos = [sensorPosition]
    % Identify fields
    Fields = strsplit(gt_sg_sub_find_nearest_string('TCPhase',data.gt_sg_settings.sensorList{sPos,4}),'.');
    optfieldname = Fields{1};
    tcphasefield = Fields{2};
    gt_sg_sub_echo({['Proceeding with aa4330 oxygen calculations using TPhase from the ',optfieldname,'.',tcphasefield,' column.']})
    
    for istep=[data.eng.dive]
        data.hydrography(istep).oxygen = gt_sg_sub_calc_oxygen;
    end
end

    function oxygen = gt_sg_sub_calc_oxygen
        
        if tau_DO(1) > 0
            tau = tau_DO(1) + tau_DO(2)*(data.hydrography(istep).temp - 20); % Sensor lag as per Johannes Hahn.
            ts = (data.hydrography(istep).time - data.hydrography(istep).time(1)) * 24 * 60 * 60; % Time array in seconds
            ts_hr = ts(1):1:ts(end); % Making high-res 1Hz time array
            
            % Interpolate TPhase to 1 Hz, apply filtering/smoothing in case of spikes
            tphase_hr = interp1(ts,gt_sg_sub_filter(data.eng(istep).(optfieldname).(tcphasefield),2,0),ts_hr);
            tphase_hr = gt_sg_sub_filter(tphase_hr,0,5);
            
            % Interpolate tau to 1 Hz, apply filtering/smoothing in case of spikes
            tau_hr = interp1(ts,gt_sg_sub_filter(tau,2,0),ts_hr);
            tau_hr = gt_sg_sub_filter(tau_hr,0,5);
            
            % Pre-allocate for speed
            tphase_new = tphase_hr;
            
            % Apply step-wise lag correction
            for jstep = 2:numel(tphase_hr)
                tphase_new(jstep) = ...
                    (tphase_hr(jstep) - (tphase_hr(jstep-1) .* ...
                    (exp(-(ts_hr(jstep)-ts_hr(jstep-1)) ./ tau_hr(jstep-1))))) ./ ...
                    (1-sum(exp(-(ts_hr(jstep)-ts_hr(jstep-1))./tau_hr(jstep-1))));
            end
            
            % Downsample back to original time array
            tphase_new = interp1(ts_hr,tphase_new,ts);
            
            calPhase = ...
                data.gt_sg_settings.aa4330_TPhaseCoef0 + tphase_new .*...
                (data.gt_sg_settings.aa4330_TPhaseCoef1 + tphase_new .*...
                (data.gt_sg_settings.aa4330_TPhaseCoef2 + tphase_new .*...
                data.gt_sg_settings.aa4330_TPhaseCoef3));
        else
            calPhase = ...
                data.gt_sg_settings.aa4330_TPhaseCoef0 + data.eng(istep).(optfieldname).(tcphasefield) .*...
                (data.gt_sg_settings.aa4330_TPhaseCoef1 + data.eng(istep).(optfieldname).(tcphasefield) .*...
                (data.gt_sg_settings.aa4330_TPhaseCoef2 + data.eng(istep).(optfieldname).(tcphasefield) .*...
                data.gt_sg_settings.aa4330_TPhaseCoef3));
        end
        
        % Calculate DELTA-p
        dp = data.gt_sg_settings.aa4330_FoilCoefB1 +...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefB2 + ...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefB3 + ...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefB4 + ...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefB5 + ...
            data.hydrography(istep).temp .* data.gt_sg_settings.aa4330_FoilCoefB6)))) + ...
            calPhase.* (...
            (data.gt_sg_settings.aa4330_FoilCoefA10 +...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefA11 +...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefA12 +...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefA13 +...
            data.hydrography(istep).temp .* data.gt_sg_settings.aa4330_FoilCoefB0)))) +...
            calPhase.* (...
            (data.gt_sg_settings.aa4330_FoilCoefA6 +...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefA7 +...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefA8 +...
            data.hydrography(istep).temp .* data.gt_sg_settings.aa4330_FoilCoefA9))) +...
            calPhase .* (...
            (data.gt_sg_settings.aa4330_FoilCoefA3 +...
            data.hydrography(istep).temp .* (data.gt_sg_settings.aa4330_FoilCoefA4 +...
            data.hydrography(istep).temp .* data.gt_sg_settings.aa4330_FoilCoefA5)) +...
            calPhase .* (...
            data.gt_sg_settings.aa4330_FoilCoefA2 +...
            data.hydrography(istep).temp .* data.gt_sg_settings.aa4330_FoilCoefA0 + ...
            calPhase .* data.gt_sg_settings.aa4330_FoilCoefA1))));
        
        % Benson and Krause coefficients (1984) for umol.kg output (recommended for oceanographic work)
        A0 = 5.80871;
        A1 = 3.20291;
        A2 = 4.17887;
        A3 = 5.10006;
        A4 = -9.86643-2;
        A5 = 3.80369;
        B0 = -7.01577e-3;
        B1 = -7.70028e-3;
        B2 = -1.13864e-2;
        B3 = -9.51519e-3;
        C0 = -2.75915e-7;
        
        % Calcate vapour pressure
        pvapour = exp(52.57 - (6690.9./(data.hydrography(istep).temp + 273.15))-(4.681.*log(data.hydrography(istep).temp + 273.15)));
        
        % Calculate Air Saturation
        airsat = (dp .* 100)./((1013.25 - pvapour) .* 0.20946);
        
        % Calculate scaled temperature
        Ts = log((298.15 - data.hydrography(istep).temp)./(273.15 + data.hydrography(istep).temp));
        
        % Calculate oxygen solubility
        lnC = A0 + Ts.*(A1 + Ts.*(A2 + Ts.*(A3 + Ts.*(A4 + A5.*Ts)))) +...
            data.hydrography(istep).salinity .*(B0 + Ts.*(B1 + Ts.*(B2 + B3.*Ts))) + ...
            C0 .* data.hydrography(istep).salinity .* data.hydrography(istep).salinity;
        solub = exp(lnC);
        
        % Oxygen concentration calculated as:
        oxygen = (solub .* airsat) ./ 100;
        
        % Pressure correction:
        oxygen = oxygen .* (1.0 + ((0.032 .* data.hydrography(istep).pressure) ./ 1000.0));
        
        if isfield(data.gt_sg_settings,'aa4330_offset')
            oxygen = oxygen + data.gt_sg_settings.aa4330_offset;
        end
        
    end % gt_sg_sub_calc_oxygen.m

end % gt_sg_sensors_aa4330.m