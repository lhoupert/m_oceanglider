% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% spread the arguments...

% see get_dive_data()
% assumes results are in variables [loginfo,eng,results] = get_dive_data('GGGYYYY')

%UNUSED but available globals = results.GLOBALS;
% dive_scalars
avg_lat = results.avg_latitude; % average latitude between GPS2 (dive) and GPS (surfacing)

% depth-average currents
if (str2num(results.depth_avg_curr_qc) == QC_GOOD)
  u_h2o = results.depth_avg_curr_north; % m/s northward
  v_h2o = results.depth_avg_curr_east; % m/s eastward
  u_av = results.flight_avg_speed_east;
  v_av = results.flight_avg_speed_north;
else
  u_h2o = 0;
  v_h2o = 0;
  u_av = 0;
  v_av = 0;
end

iter =  1; % MADE UP ?? dive_scalars(9); % number of TSV iterations
converged = (str2num(results.hdm_qc) == QC_GOOD); % did TSV converge?
CTD_qc = str2num(results.CTD_qc);

% dive_results, per data point
dive_start_time = eng.start_ts; % unix epoch seconds (see note below)
sg_epoch_time = results.time;
sg_time = sg_epoch_time - dive_start_time; % elapsed
time = eng.elaps_t; % elapsed time in secs since start of dive
% vmtime = time;
% vmdepth = eng.depth; % meters
ctd_epoch_time = results.ctd_time;
ctd_time = ctd_epoch_time - dive_start_time; % elapsed

ctd_np = length(ctd_time); % total number of CTD points
sg_np = length(time); % total number of SG data points
mp = sg_np;  % BACKCOMPAT
unknown_sg_v = nan(sg_np,1);
qc_missing_sg_v = QC_MISSING*ones(sg_np,1);

% various indices for different events.
%?? apogee_pump_start_i = directives.pump_start;

sg_press    = results.pressure; % [dbar]
sg_depth_m  = results.depth; % [m] the measured depth
ctd_depth_m = results.ctd_depth;
ct_depth_m  = ctd_depth_m; % BACKCOMPAT

% w = m2cm.*ctr1stdiffderiv(-sg_depth_m,sg_time);

DAC_qc = str2num(results.depth_avg_curr_qc); % trust the depth-average current info?

serial_date_time = unix_to_datenum(dive_start_time);
secs_per_day  = 86400.0;
day = serial_date_time + ctd_time/secs_per_day; % the fraction of days since start of dive
lat = results.latitude; % latitude
lon = results.longitude; % longitude
                         
% if GPS2_ok:
dtgps2E = loginfo.GPS_t - loginfo.GPS2_t; % elapsed time between $GPS and $GPS2
midpoint_i = find(ctd_time > dtgps2E/2,1,'first');
if (isempty(midpoint_i))
  midpoint_i = find(ctd_time > ctd_time(end)/2,1,'first');
end
midpoint_day = day(midpoint_i);
midpoint_lat = lat(midpoint_i);
midpoint_lon = lon(midpoint_i);

if (isfield(results,'start_of_climb_time'))
  start_of_climb_time = results.start_of_climb_time; % elapsed
  start_of_climb_i = find(ctd_time >= start_of_climb_time,1,'first');
else
  fprintf(1,'Unable to determine start of climb index. Old NC file?\n');
  start_of_climb_i = find(sg_depth_m == max(sg_depth_m),1,'first');
  start_of_climb_time = sg_time(start_of_climb_i); % elapsed secs
end

spd_model = results.speed; % [m/s] predicted speed through the water
glideangle_model = results.glide_angle;  % [deg] glide angle
w_model = results.vert_speed; % [m/s] vertical velocity

temp_raw = results.temperature_raw; % [deg]
temp_raw_qc = str2num(results.temperature_raw_qc);
cond_raw = results.conductivity_raw; % [...] lag-corrected temperature
cond_raw_qc = str2num(results.conductivity_raw_qc);
salin_raw = results.salinity_raw; % [PSU] lag-corrected temperature
salin_raw_qc = str2num(results.salinity_raw_qc); 

temp = results.temperature; % [deg] lag-corrected temperature
temp_qc = str2num(results.temperature_qc);
cond = results.conductivity; % [...] lag-corrected conductivity
cond_qc = str2num(results.conductivity_qc);
salin = results.salinity; % [PSU] lag-corrected salinity
salin_qc = str2num(results.salinity_qc); 
salinity_error = zeros(ctd_np,1); % TODO replace this if we return this

% createa a fake directives structure
directives = ([]); 
directives.start_of_dive = 1;
directives.start_of_climb = start_of_climb_i;
directives.end_of_climb = sg_np;
directives.interpolated_salinity = find(salin_qc == QC_INTERPOLATED);
directives.stable_ts = []; % sigh..

% flags: do we have any of these optional sensors?
solubility_micromoles_per_kg = results.dissolved_oxygen_sat; % ctd 
sbe43_present = isfield(results,'SBE43_qc');
if (sbe43_present)
  sbe43_qc = str2num(results.SBE43_qc);
  sbe43_results_qc = str2num(results.sbe43_dissolved_oxygen_qc);
  sbe43_time = results.sbe43_results_time - results.sbe43_results_time(1); % elapsed
  oxygen = results.sbe43_dissolved_oxygen;
end

% TODO look for optode instrument O2 rather than qc, ensure at least that is reported
% REQUIRES basestation to echo those (and optode_time) even if calib constants are not present
optode_present = isfield(results,'aanderaa3830_qc');
if (optode_present)
  optode_qc = str2num(results.aanderaa3830_qc);
  optode_time = results.aanderaa3830_results_time - results.aanderaa3830_results_time(1); % elapsed
  optode_instrument_oxygen = results.aanderaa3830_instrument_dissolved_oxygen;
  optode_oxygen_qc = results.aanderaa3830_dissolved_oxygen_qc;
  optode_oxygen = results.aanderaa3830_dissolved_oxygen;
else
  optode_present = isfield(results,'aanderaa4330_qc');
  if (optode_present)
    optode_qc = str2num(results.aanderaa4330_qc);
    optode_time = results.aanderaa4330_results_time - results.aanderaa4330_results_time(1); % elapsed
    optode_instrument_oxygen = results.aanderaa4330_instrument_dissolved_oxygen;
    optode_oxygen_qc = results.aanderaa4330_dissolved_oxygen_qc;
    optode_oxygen = results.aanderaa4330_dissolved_oxygen;
  end
end

% We only look for red/blue/fluor counts and ignore ref signals
% NOTE: red and blue have been inverted for forever!!
% The tags are in red/blue/fluor/temp order always
% At the moment we always expect a redCount?  What about BB3 or FL3? (the later loses, we drop green on the former)
% Should we trigger on fluor counts since that is present in BB2FL and BBFL2?
% Unfortunately we drop the CDOM signal on the later instrument
wetlabs_data = ([]); % empty struct
if isfield(eng,'redCount')
	wetlabs_data = eng;
	wetlabs_time = sg_time;
	wetlabs_tags = {'redCount','blueCount','fluorCount','VFtemp'};
elseif isfield(eng,'wlbb2f_redCount')  % alternative names on eng files
	wetlabs_data = eng;
	wetlabs_time = sg_time;
	wetlabs_tags = {'wlbb2f_redCount','wlbb2f_blueCount','wlbb2f_fluorCount','wlbb2f_VFtemp'};
elseif isfield(eng,'wlbb2flvmt_wl470sig')  % alternative names on eng files
	wetlabs_data = eng;
	wetlabs_time = sg_time;
	wetlabs_tags = {'wlbb2flvmt_wl470sig','wlbb2flvmt_wl700sig','wlbb2flvmt_Chlsig','wlbb2flvmt_VFtemp'};
elseif isfield(results,'wlbb2f_redCount') % try scicon
	wetlabs_data = results;
	wetlabs_time = results.wlbb2f_time - dive_start_time;
	wetlabs_tags = {'wlbb2f_redCount','wlbb2f_blueCount','wlbb2f_fluorCount','wlbb2f_VFtemp'};
elseif isfield(results,'wlbb2flvmt_wl470sig')  % alternative name on scicon
	wetlabs_data = results;
	wetlabs_time = results.wlbb2flvmt_time - dive_start_time;
	wetlabs_tags = {'wlbb2flvmt_wl470sig','wlbb2flvmt_wl700sig','wlbb2flvmt_Chlsig','wlbb2flvmt_VFtemp'};
elseif isfield(results,'wlbb2flvmt_470sig')  % yet another alternative name on scicon
	wetlabs_data = results;
	wetlabs_time = results.wlbb2flvmt_time - dive_start_time;
	wetlabs_tags = {'wlbb2flvmt_470sig','wlbb2flvmt_700sig','wlbb2flvmt_Chlsig','wlbb2flvmt_temp'};
end

if ~isempty(wetlabs_data)
	clear red_scttr blue_scttr fluor fluor_temp
	wetlabs_present = 1;
	wetlabs_qc = ones(length(wetlabs_time),1)*QC_GOOD;
	if isfield(wetlabs_data,wetlabs_tags{1})
		red_scttr = wetlabs_data.(wetlabs_tags{1});
	end
	if isfield(wetlabs_data,wetlabs_tags{2})
		blue_scttr = wetlabs_data.(wetlabs_tags{2});
	end
	if isfield(wetlabs_data,wetlabs_tags{3})
		fluor = wetlabs_data.(wetlabs_tags{3});
	end
	if isfield(wetlabs_data,wetlabs_tags{4})
		fluor_temp = wetlabs_data.(wetlabs_tags{4});
	end
else
	wetlabs_present = 0; % Jim Bennet fix 
end

%
%  KUTI WET Labs detection, based on iRobot's approach.
%
wetlabs_bb2fl_present = 0;
wetlabs_bbfl2_present = 0;
wetlabs_fl3_present = 0;
wetlabs_bb3_present = 0;

if (isfield(eng,'wlbb2fl_BB1sig') || ...
        isfield(eng,'wlbb2flvmt_BB1sig'))           % old - used on SG558 refurb
    % KUTI - the BB1sig field name is generated by wlbb2fl.cnf and
    % passed thru to us by the basestation. Use iRobot approach to plotting
    % the data.
    wetlabs_bb2fl_present = 1;
end
if (isfield(eng,'wlbbfl2_BB1sig') || ...
        isfield(eng,'wlbbfl2_BBsig') || ...         % old
        isfield(eng,'wlbbfl2vmt_BBsig'))            % old - used on SG602-SG606
    % KUTI - the BBSig field name is generated by wlbbfl2.cnf and
    % passed thru to us by the basestation. Use iRobot approach to plotting
    % the data.
    wetlabs_bbfl2_present = 1;
end
if isfield(eng,'wlfl3_FL1sig')
    % KUTI - the FL1Sig field name is generated by wlfl3.cnf and
    % passed thru to us by the basestation. Use iRobot approach to plotting
    % the data.
    wetlabs_fl3_present = 1;
end
if isfield(eng,'wlbb3_BB1sig')
    % KUTI - the BB1Sig field name is generated by wlbb3.cnf and
    % passed thru to us by the basestation. Use iRobot approach to plotting
    % the data.
    wetlabs_bb3_present = 1;
end
if (isfield(eng,'wlbbfl2_650sig') && ~isfield(eng, 'wlbbfl2_650ref') && ...
        isfield(eng,'wlbbfl2_Chlsig') && ~isfield(eng, 'wlbbfl2_Chlref') && ...
        isfield(eng,'wlbbfl2_Cdomsig') && ~isfield(eng, 'wlbbfl2_Cdomref') && ...
        isfield(eng,'wlbbfl2_temp'))
    %
    % Fixup for Pedro (SG576, June 2014).  Add the three reference value columns.
    %
    eng.wlbbfl2_650ref(1:length(eng.wlbbfl2_650sig)) = 650;
    eng.wlbbfl2_Chlref(1:length(eng.wlbbfl2_Chlsig)) = 695;
    eng.wlbbfl2_Cdomref(1:length(eng.wlbbfl2_Cdomsig)) = 460;
    wetlabs_bbfl2_present = 1;
end


compute_dive_climb_indices
