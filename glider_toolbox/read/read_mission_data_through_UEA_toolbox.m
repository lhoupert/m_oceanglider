function [glider,data] = read_mission_data_through_UEA_toolbox(missioninfo,proc_param)
%
%=========================================================================%
% Function to read seaglider data from eng and log files using UEA toolbox 
% routines (include in the toolbox). The version of the UEA included in 
% this toolbox is the branch devSAMS which integrate additional support for
% the Kongsberg Seaglider
%
% L. Houpert, SAMS 25/01/2016
%=========================================================================%
%
% [glider,data] = read_mission_data_through_UEA_toolbox(missioninfo,proc_param)
% Function that read seaglider data (TS, flightmodel) from the basestation
% netcdf and concatenate all the dive in a structure. 
%   Inputs: o missioninfo: [1x1] structure with mission details
%                where - missioninfo.glname      : name of the glider (e.g. sg605);
%                      - missioninfo.glmission   : name of the mission (e.g. OSNAP1);
%                      - missioninfo.rawdata.dataDir : path of the eng and log files;
%                      - missioninfo.rawdata.ueadivetoproces  : the divenumbers to process 
%                      - missioninfo.rawdata.flight_regression_passes :   number of 
%                            regression passes to optimize the flight 
%                            model parameters (hd_a, hd_b, hd_c, volmax, abs_compress, 
%                            therm_expan).
%           o proc_param [1x1] structure containing:
%                      - (optional: proc_param.length)     : the length of the hull
%                      - (optional: proc_param.sensorList) : the list of the sensors mounted on the glider 
%                            with the associated matlab script (e.g. 
%                            proc.sensorList={'SBE_CT' 'sbect' 'gt_sg_sensors_sbect.m' })
%                            By default, a gui will ask the user. 

%  Output:  - glider : structure with the concatenate and process glider
%                      data convert to a uniform format (same that the data
%                      structure from the netcdf conversion)
%           - data   : original structure from the UEA processing

% created by L. Houpert (houpertloic@gmail.com), 25/01/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 
% If one wants to add others field to the output structure (biogeochemical sensor data, flight model output,...), just modify the
% script. Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote
% repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)
%


disp(['Flight regression passes: ' ...
                            num2str(missioninfo.rawdata.flight_regression_passes)])

glname    = missioninfo.glname;
glmission = missioninfo.glmission;
dataDir   = missioninfo.rawdata.dataDir;
gnumber   = str2double(glname(3:end));

if exist('./matlab_diary','dir')~=7
    mkdir('matlab_diary')
end

diaryfile1 = ['matlab_diary/diary_gt_sg_load_merge_data_' glmission '_' glname '.log'];
diaryfile2 = ['matlab_diary/diary_gt_sg_process_data_' glmission '_' glname '.log'];
logfile1   = ['matlab_diary/' glmission '_' glname '_dives_with_negative_depth_values.log'];
logfile2   = ['matlab_diary/' glmission '_' glname '_dives_with_nan_CT.log '];


% Load and merge the eng and log files through the UEA toolbox:
data = gt_sg_load_merge_data(gnumber,dataDir ,missioninfo.rawdata.ueadivetoprocess,[dataDir '/..']);
eval(['!/bin/mv dives_with_negative_depth_values.log ' logfile1])
eval(['!/bin/mv dives_with_nan_CT.log ' logfile2])

eval(['!/bin/mv gt_sg_processing.log ' diaryfile1])
%save temploadmerge.mat data

% Process the merge data through the UEA toolbox:
data.sg_calib_const.flight_regression_passes = missioninfo.rawdata.flight_regression_passes;
if isfield(proc_param,'length')==1
    data.sg_calib_const.length                   = proc_param.length;
end
if isfield(proc_param,'sensorList')==1
    data.sg_calib_const.sensorList               = proc_param.sensorList;
end
data=gt_sg_process_data(data);
eval(['!/bin/mv gt_sg_processing.log ' diaryfile2])

%save tempprocess.mat data
%load tempprocess

%==========================================================================
% Conversion to generic matlab structure for the glider mission (similar to
% the one loaded with read_seaglider_nc_basestation.m:
glider = convert_UEA_struc_to_common_glider_struc(data,missioninfo,proc_param);


end