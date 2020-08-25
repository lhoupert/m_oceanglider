function data = gt_sg_sensors_qsppar(data)
% -=-=-=-=-=-=-=-=-=-=-=-=-=-
% COMPUTE
%  QSP PAR DATA
% -=-=-=-=-=-=-=-=-=-=-=-=-=-
%
% calibs = sg_calib_constants struct
% eng    = eng struct
%

%% Check for calibration data
if ~isfield(data.gt_sg_settings,'PAR_CalData')
    gt_sg_sub_echo({'No PAR calibration data. Cannot continue with processing'});
    return;
end

%% Identify sensor(s) that need to be processed
sensorPosition = strcmp(data.gt_sg_settings.sensorList(:,3),'gt_sg_sensors_qsppar.m');

%% List dives to process, no need for log files or currents for processing,
%% so all eng files should work
dives = [data.eng.dive];


Fields = strsplit(gt_sg_sub_find_nearest_string('PARuV',data.gt_sg_settings.sensorList{sensorPosition,4}),'.');
sensorfieldname = Fields{1};
datafieldname = Fields{2};

gt_sg_sub_echo({['Calculating QSP PAR data.']});

for dstep = [dives]
    data.hydrography(dstep).PAR = data.gt_sg_settings.PAR_CalData.scalingFactor .* (data.eng(dstep).(sensorfieldname).(datafieldname) - data.gt_sg_settings.PAR_CalData.sensorDark/1000);
end

end