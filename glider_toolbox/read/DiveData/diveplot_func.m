% ************************************************************************
% Plot sg depth, pitch, roll, pitch control,  buoyancy, etc. vs time
%	includes corrections for sensor lags and offsets
%	Includes estimate of glider speed from pitch and buoyancy alone
%
% dataFolder = 'Z:\Workspace\Data\sg525\PortSusan_5-6Oct2010'
% diveNumber = 20
% bias = 0
% dive_out_dir = 'Z:\Workspace\Data\sg525\PortSusan_5-6Oct2010'
% ************************************************************************

function [logp, eng, loginfo, err] = diveplot_func(dataFolder, dive_out_dir, diveNumber, vbdbias, settings, filenames, handles)

    global seagliderID;
    global NUM_PLOTS;
    global COUNTS_PER_CC;
    global gSpeedOfSound;

    COUNTS_PER_CC = 4.0767;
    numWaitStages = 25;
    xl = 1.8;

    closeOrphanedWaitbars();

    % ------------------------------------------------------------------------------------------------------------------------
    %
    % ------------------------------------------------------------------------------------------------------------------------
    if exist('settings', 'var') == 0
        close all;
        settings = [];
        settings.plotSelectionState(1:NUM_PLOTS) = 1;
        % DISABLE SPEED-OF-SOUND
        %settings.plotSelectionState(13) = 0;
    else
        if (length(settings.plotSelectionState) < NUM_PLOTS)
            for i=length(settings.plotSelectionState):NUM_PLOTS
                %settings.plotSelectionState(1:NUM_PLOTS) = 1;
            end
        end
    end

    % ------------------------------------------------------------------------------------------------------------------------
    % CONFIGURE PATH TO POINT TO FOLDERS CONTAINING REQUIRED FUNCTIONS
    % ------------------------------------------------------------------------------------------------------------------------
    % jf configurePaths();

    % ------------------------------------------------------------------------------------------------------------------------
    % INITIALIZE VARIABLES
    % ------------------------------------------------------------------------------------------------------------------------
    waitbarIndex = 0;
    hWaitbar = waitbar(waitbarIndex/numWaitStages, 'Initializing variables...');
    mp = [];
    sigma_apogee = [];
    sigma_1m = [];
    vbd_max_cnts = [];
    vbd_1m = [];
    c_vbd_true = [];
    mass_apogee_implied = [];
    cc_surf_min = [];
    vbd_min_cnts = [];
    istdy = [];
    logp = [];
    eng  = [];
    loginfo = [];

    % ------------------------------------------------------------------------------------------------------------------------
    % READ THE SEAGLIDER DATA FILES
    % ------------------------------------------------------------------------------------------------------------------------
    err = 0;
    
    [waitbarIndex] = incrementWaitbar(hWaitbar, waitbarIndex, numWaitStages, 'Reading Seaglider data files...');

    seagliderID = settings.seagliderID;
    
    [filenames] = buildFilenames(dataFolder, seagliderID, diveNumber);

    try
        id_dive = sprintf('%03d%04d',seagliderID,diveNumber);
        
        % Extract the data from the netCDF file.
        [loginfo,eng,results] = dd_get_dive_data(id_dive,1,dataFolder);
        if (isfield(results,'processing_error'))
            vizMsgBox('Inconsistent dive data - basestation marked .nc file as "bad". May not be able to generate the plots for this dive.', 'ERROR: dive data', 'modal');
            closeOrphanedWaitbars();
            err = 1;
            %return;
        end
        sg_cal = results.sg_calib_constants;
    catch
        vizMsgBox('Unable to load Seaglider data file', 'ERROR: Data read', 'modal');
        closeOrphanedWaitbars();
        err = 1;
        return;
    end
    
    base_file = sprintf('p%s',id_dive);
    
    serial_date_time = loginfo.log_dn;
    timestr = strcat(datestr(serial_date_time),'Z'); % UTC
    startstr = strcat(' ', timestr);

    titlestring = sprintf('SG%d - Dive %d - %s\nMission Start Time: %s\n', settings.seagliderID, diveNumber, mission_title, startstr);

    diveplot_body(loginfo, eng, results, sg_cal, base_file, hWaitbar, waitbarIndex, numWaitStages, settings, titlestring, handles);

    %err = 0;
    return;
    
    
% ************************************************************************************************************************
% ************************************************************************************************************************
% ************************************************************************************************************************
% utility functions
% ************************************************************************************************************************
% ************************************************************************************************************************
% ************************************************************************************************************************

% ************************************************************************************************************************
%
% ************************************************************************************************************************
function ReportState(fid, states, energy, statename, varname)
   j = find(strcmp(states, statename) == 1);
   if ~isempty(j)
      fprintf(fid, '%s = %5.0f\n', varname, energy(j));
   else
      fprintf(fid, '%s = 0\n', varname);
   end

% ********************************************************************************
%
% ********************************************************************************
function [zscl] = get_zscl(dataStore)
    [~, vmdepth] = dataStore.queryDouble('depth');
        
    if max(vmdepth) > 500
        zscl = 10;
    elseif max(vmdepth) > 300
        zscl = 5;
    elseif max(vmdepth) > 200
        zscl = 3;
    elseif max(vmdepth) > 100
        zscl = 2;
    else
        zscl = 1;
    end

% ********************************************************************************
%
% ********************************************************************************
function [hv_voltage, hv_used, lv_voltage, lv_used] = getVoltages(logp)
    if isfield(logp, 'a24V_AH')
        hv_voltage = logp.a24V_AH(1);
        hv_used = logp.a24V_AH(2) / logp.AH0_24V;
    else
        hv_voltage = 0;
        hv_used = 0;
    end

    if isfield(logp, 'a10V_AH')
        lv_voltage = logp.a10V_AH(1);
        lv_used = logp.a10V_AH(2) / logp.AH0_10V;
    else
        lv_voltage = 0;
        lv_used = 0;
    end

% ********************************************************************************
%
% ********************************************************************************
function [sm_amp_sec] = getAmps(logp)
    if isfield(logp, 'SM_CCo')
        sm_amp_sec = logp.SM_CCo(2) * logp.SM_CCo(3);
    else
        sm_amp_sec = 0;
    end

% ********************************************************************************
%
% ********************************************************************************
function [climb_right_rolls, climb_left_rolls,...
          dive_right_rolls, dive_left_rolls,...
          climb_vbd_pumps, dive_vbd_pumps,...
          climb_pitch, dive_pitch] = getManeuveringStatistics(roll_control, pitch_control, pitch, vbddiff)

    j1 = find(roll_control > 10);
    j1(end + 1) = 100000;
    jp = j1(diff(j1) > 1);

    j1 = find(roll_control < -10);
    j1(end + 1) = 100000;
    jn = j1(diff(j1) > 1);

    iwn = find( pitch_control < 0 );
    iwp = find( pitch_control > 0 );

    ivb = find( vbddiff > 0.5 );

    climb_right_rolls = length(intersect(iwp, jp));
    climb_left_rolls = length(intersect(iwp, jn));
    dive_right_rolls = length(intersect(iwn, jp));
    dive_left_rolls = length(intersect(iwn, jn));
    climb_vbd_pumps = length(intersect(iwp, ivb));
    dive_vbd_pumps = length(intersect(iwn, ivb));
    climb_pitch = mean(pitch(iwp));
    dive_pitch = mean(pitch(iwn));

% ********************************************************************************
%
% convert to cartesian vectors to determine current direction
%
% ********************************************************************************
function [Udir, Umag] = calculateCurrentDirection(cog_diving, dog_diving, ave_brg, displacement, vmtime, mp)
    [dogN, dogE] = pol2cart(cog_diving*pi/180, dog_diving);
    % [dtwN, dtwE] = pol2cart((ave_brg + decl)*pi/180, displacement);
    [dtwN, dtwE] = pol2cart((ave_brg)*pi/180, displacement);

    % Ue = dogE/gps_dt - dtwE/(vmtime(mp) - vmtime(1));
    Ue = (dogE - dtwE)/(vmtime(mp) - vmtime(1));
    Un = (dogN - dtwN)/(vmtime(mp) - vmtime(1));
    [Udir, Umag] = cart2pol(Un, Ue);
    Udir = Udir * 180/pi;
    if(Udir < 0);
        Udir = Udir + 360;
    end;
    
