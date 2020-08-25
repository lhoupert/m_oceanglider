% ********************************************************************************
% ********************************************************************************
% ********************************************************************************
% ********************************************************************************
function diveplot_body(loginfo, eng, results, sg_cal, base_file, hWaitbar, waitbarIndex, numWaitStages, settings, titlestring, handles)

% from here down could be called after dive_process
cm2m = 0.01;
radpd = pi/180;
degpr = 180./pi;
secs_per_hour = 3600.0;
dbar_per_psi = 1.0/1.4503774;
psi_per_meter = 1.0/0.685; % this needs to match the value in glider/constants.h

if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
  fprintf(1,'Processing error or skipped profile; skipping\n');
  return;
end

qc_declarations;
clear_sg_config_constants;

%
% Do not read in sg_calib_constants - use sg_cal instead.
%
%jf sg_calib_constants;
%evalin('caller','sg_config_constants'); % update these constants based on sg_calib_constant values
hd_a = sg_cal.hd_a;
hd_b = sg_cal.hd_b;
hd_c = sg_cal.hd_c;
volmax = sg_cal.volmax;
rho0 = sg_cal.rho0;
mass = sg_cal.mass;
mission_title = sg_cal.mission_title;
calibcomm = sg_cal.calibcomm;

sg_config_constants;

% assumes loginfo, eng, and results are available from caller
unpack_data;
% reduce nearest ctd data to sg grid
% we need these indices to plot data against different scales in dd_diveplot_graphs
[ctd_depth_ignored,ctd2sg_i,sg_ctd_valid_i,ctd2sg_valid_i,ctd_dive_i,ctd_climb_i] = compute_sensor_sg_indices(ctd_time,salin_qc,sg_time,sg_depth_m,sg_dive_i,sg_climb_i);
if (sbe43_present)
  [oxy_depth_m,oxy2sg_i,sg_oxy_valid_i,oxy2sg_valid_i,oxy_dive_i,oxy_climb_i] = compute_sensor_sg_indices(sbe43_time,sbe43_results_qc,sg_time,sg_depth_m,sg_dive_i,sg_climb_i);
end
if (optode_present)
  [opt_depth_m,opt2sg_i,sg_opt_valid_i,opt2sg_valid_i,opt_dive_i,opt_climb_i] = compute_sensor_sg_indices(optode_time,optode_oxygen_qc,sg_time,sg_depth_m,sg_dive_i,sg_climb_i);
end
if (wetlabs_present)
  [wet_depth_m,wet2sg_i,sg_wet_valid_i,wet2sg_valid_i,wet_dive_i,wet_climb_i] = compute_sensor_sg_indices(wetlabs_time,wetlabs_qc,sg_time,sg_depth_m,sg_dive_i,sg_climb_i);
end

if (~exist('calibcomm','var'))
  calibcomm = 'Unknown calibration';
end

cc_surf = loginfo.SM_CC;
c_pitch = loginfo.C_PITCH;
% always use what the glider used, not what the pilot thinks!
pitch_max_cnts = loginfo.PITCH_MAX;
pitch_min_cnts = loginfo.PITCH_MIN;
c_roll_dive = loginfo.C_ROLL_DIVE;
c_roll_climb = loginfo.C_ROLL_CLIMB;
c_roll = 0.5*(c_roll_dive + c_roll_climb);
% always use what the glider used, not what the pilot thinks!
roll_max_cnts = loginfo.ROLL_MAX;
roll_min_cnts = loginfo.ROLL_MIN;
c_vbd = loginfo.C_VBD;
speed_factor = loginfo.SPEED_FACTOR;
head_errband = loginfo.HEAD_ERRBAND;
target_depth = loginfo.D_TGT;
dive_time = loginfo.T_DIVE;
if( isfield(loginfo, 'D_GRID'))
  grid_depth = loginfo.D_GRID;

  mhead = loginfo.MHEAD_RNG_PITCHd_Wd; % get array
  % What the glider thought it should do...
  mag_heading = mhead(1); % best heading (magnetic)
  intended_pitch = mhead(3); % best overall pitch for distance and speed
  w_d = mhead(4); % vertical velocity desired

  magvar = results.magnetic_variation;
  pitch_control_intended = loginfo.gc_pitch_ctl;
  apogee_control_index = find(loginfo.gc_vbd_ctl == 0);
  
  pitch_d = -intended_pitch;
  kalman_control_u = loginfo.KALMAN_CONTROL(1);
  kalman_control_v = loginfo.KALMAN_CONTROL(2);

  mc = length(pitch_control_intended);
  pitch_intended = zeros(1,mc);
  pitch_intended(1:apogee_control_index-1) = intended_pitch; 
  pitch_intended(apogee_control_index) = loginfo.APOGEE_PITCH; 
  pitch_intended(apogee_control_index+1:mc) = -intended_pitch;
else
  grid_depth = 1e4;
  mhead = zeros(1,4);
  mag_heading = 0;
  magvar = 0.;
  kalman_control_u = 0;
  kalman_control_v = 0;
  pitch_d = 0;
  pitch_control_intended = 0;
  apogee_control_index = 1;
  mc = length(pitch_control_intended);
  pitch_intended = loginfo.APOGEE_PITCH*ones(1,mc);
end

max_buoy = loginfo.MAX_BUOY;
course_bias = loginfo.COURSE_BIAS;

dive_depth = min(target_depth, grid_depth);
% w_target should be -w_d*cm2m
w_target = 2*target_depth/(dive_time*60);
% compute desired vertical and horizontal speed
% we use pitch_d/w_target  (scalars) in a version of the glide_slope model to
% estimate th_d.  code assumes things never go complex on us...
a_onbd = hd_a;
b_onbd = hd_b;
c_onbd = hd_c;
s_onbd = hd_s;

[th_d,spdg_d] = glide_slope(w_target/cm2m,pitch_d,a_onbd,b_onbd,c_onbd,s_onbd,rho0);
spdg_d = spdg_d * cm2m; % convert to m/s
hspdg_d = spdg_d*cos(th_d);

timestr = strcat(datestr(serial_date_time),'Z'); % UTC
startstr = strcat(' ', timestr);

% when we made these plots
timelabel = datestr(clock);

% vmtime   = eng.elaps_t;
vmtime = sg_time;
sgtime_m = vmtime/60; % minutes

mp = length(vmtime);
total_dive_time_s = vmtime(mp) - vmtime(1);
dive_duration = total_dive_time_s/secs_per_hour; % hrs

vmdepth  = eng.depth;	% depth in meters
hdg = eng.head;
pitch = eng.pitchAng;
roll = eng.rollAng;
pitch_control = eng.pitchCtl;
roll_control = eng.rollCtl;
vbd = eng.vbdCC;

if (sbe43_present)
  if (exist('calibcomm_oxygen','var'))
    ox_comment = calibcomm_oxygen;
  else
    ox_comment = 'Unknown SBE43 calibration';
  end
end

if (optode_present)
  if (exist('calibcomm_optode','var'))
    ox_comment = calibcomm_optode;
  else
      if (isfield(results.sg_calib_constants, 'comm_oxy_type') == 1 && isfield(results.sg_calib_constants, 'calibcomm_optode') == 1)
          ox_comment = sprintf('%s %s', results.sg_calib_constants.comm_oxy_type, results.sg_calib_constants.calibcomm_optode);
      else
          if (isfield(results.sg_calib_constants, 'calibcomm_optode') == 1)
              ox_comment = results.sg_calib_constants.calibcomm_optode;
          else
	          ox_comment = 'Unknown optode calibration';
          end
      end
  end
end

if (wetlabs_present || wetlabs_bb2fl_present || wetlabs_bbfl2_present || wetlabs_bb3_present || wetlabs_fl3_present )
      if (isfield(results.sg_calib_constants, 'calibcomm_wetlabs') == 1)
          wetlabs_comment = results.sg_calib_constants.calibcomm_wetlabs;
      else
          wetlabs_comment = 'Unknown WETLabs serial num and calibration date';
      end
end
    


delta_t =  zeros(mp,1);
delta_t(2:end) = diff(vmtime);

%
% find changes in control signals
%
pitchdiff = zeros(mp,1);
pitchdiff(2:mp) = diff(pitch_control);
rolldiff = zeros(mp,1);
rolldiff(2:mp)  = diff(roll_control);

hdgdiff = zeros(mp,1);
hdgdiff(2:mp) = diff(hdg);
hdgdiff = mod(hdgdiff,360);
in_hdg = find(hdgdiff > 180);
hdgdiff(in_hdg) = hdgdiff(in_hdg) - 360;
hdg_wrapped = hdg(1) + cumsum(hdgdiff);
turn_rate = ctr1stdiffderiv(hdg_wrapped, vmtime);

vbddiff = zeros(mp,1);
vbddiff(2:mp)  = diff(vbd)./delta_t(2:mp);

zpos = -vmdepth;
w  = 100.*ctr1stdiffderiv(zpos, vmtime); % observed w

%
% get derived/corrected quantities
%

% valid_i, dive_i, climb_i is set in unpack_data
start_valid_i = min(valid_i); % was 1
end_valid_i = max(valid_i); % was mp

% temp_raw = results.temperature_raw;
% salin_raw = results.salinity_raw;
% temp = results.temperature;
% salin = results.salinity;

% always use what the glider used, not what the pilot thinks!
vbd_min_cnts = loginfo.VBD_MIN;
vbd_max_cnts = loginfo.VBD_MAX;
vbd_cnv = loginfo.VBD_CNV;
vbd_cnts_per_cc = 1.0/loginfo.VBD_CNV;

vbd_min = (vbd_max_cnts - c_vbd)*vbd_cnv;
vbd_max = (vbd_min_cnts - c_vbd)*vbd_cnv;
vol0 = volmax + (c_vbd - vbd_min_cnts)/vbd_cnts_per_cc;
density = results.density;
density_raw = sw_dens0(salin_raw,temp_raw);
sigma_t = results.sigma_t;

% recompute speeds, displacements based on possibly changed volmax, mass, etc.
% pilot often does this to see the impact on a past dive of new parameters in sg_calib_constants
ctd_w = interp1d(sg_epoch_time,w,ctd_epoch_time,'linear');
ctd_pitch = interp1d(sg_epoch_time,pitch,ctd_epoch_time,'linear');
ctd_roll = interp1d(sg_epoch_time,roll,ctd_epoch_time,'linear');
ctd_hdg = interp1d(sg_epoch_time,hdg,ctd_epoch_time,'linear');
ctd_vbd = interp1d(sg_epoch_time,vbd,ctd_epoch_time,'linear');
ctd_sg_press = interp1d(sg_epoch_time,sg_press,ctd_epoch_time,'linear');

ctd_depth = results.ctd_depth;
ctd_delta_t = zeros(ctd_np, 1);
ctd_delta_t(2:end) =  diff(ctd_time);

[glideangle_gsm,speed_gsm] = glide_slope(ctd_w,ctd_pitch,hd_a,hd_b,hd_c,hd_s,rho0);
% glide_slope() tags infinities as stalls
stalled_i = find_stalled(speed_gsm,ctd_pitch,min_stall_speed,max_stall_speed,min_stall_angle*radpd);
speed_gsm(stalled_i) = 0; % mark as stall
glideangle_gsm(stalled_i) = 0; % going nowhere
hspdg = speed_gsm.*cos(glideangle_gsm);
glideangle_stdy = glideangle_gsm*degpr; % [deg]

displaced_vol = (vol0 + ctd_vbd);
% The term below is for the volume compression with temp and presure
% TODO CCE reports that (all types of?) oil in the VBD unit chnages volume by temperature
vol = displaced_vol.*exp(-abs_compress*ctd_sg_press + therm_expan*(temp - temp_ref));
if (mass_comp)
  % DG: adjust volume because of compressee at depth...
  kg2g = 1000;
  vol_comp_ref = kg2g*mass_comp/cml_dens(temp_ref, 0);
  vol_comp = kg2g*mass_comp./cml_dens(temp, ctd_sg_press);
  vol = vol + vol_comp - vol_comp_ref;
end
density_insitu = sw_dens(salin_raw, temp_raw, ctd_sg_press);
buoy = 1000*( -mass + density_insitu.*vol*1.e-6 );  % vehicle buoyancy [g]

% ignore interstital, tail wake effects on buoyancy
[glidespeed_stdy_hdm,glideangle_stdy_hdm] = flightvec(buoy, ctd_pitch, glider_length, hd_a, hd_b, hd_c, rho0, hd_s);
stalled_i = find_stalled(glidespeed_stdy_hdm,ctd_pitch,min_stall_speed,max_stall_speed,min_stall_angle*radpd);
glidespeed_stdy_hdm(stalled_i) = 0; % mark as stall
glideangle_stdy_hdm(stalled_i) = 0; % going nowhere
hspd_stdy = glidespeed_stdy_hdm.*cos(glideangle_stdy_hdm*radpd);
w_stdy    = glidespeed_stdy_hdm.*sin(glideangle_stdy_hdm*radpd);

% recompute displacements and avg flight speed and direction
[vv_s, uu_s] = pol2cart( (ctd_hdg + results.magnetic_variation)*radpd, hspd_stdy );
dx = uu_s.*ctd_delta_t*cm2m; % east displacement increment [m]
dy = vv_s.*ctd_delta_t*cm2m; % north displacement increment [m]
inan = find(isnan(hspd_stdy));
dx(inan) = 0;
dy(inan) = 0;
x_s = cumtrapz(dx); % rather than cumsum
y_s = cumtrapz(dy); % rather than cumsum
dx = x_s; % i.e., results.east_displacement
dy = y_s; % i.e., results.north_displacement
flight_time = sum(ctd_delta_t);
u_av = (x_s(end) - x_s(1))/flight_time; % why not dx(end)/flight_time
v_av = (y_s(end) - y_s(1))/flight_time; % why not dy(end)/flight_time

[ave_brg, ave_spd] = cart2pol(v_av, u_av);
ave_brg = ave_brg*degpr;
if(ave_brg < 0) 
  ave_brg = ave_brg + 360;
end
displacement = ave_spd*flight_time;

%fprintf(1, 'Avg speed = %.4f cm/s\n', ave_spd/cm2m);
%fprintf(1, 'Avg bearing %.2f deg True\n', ave_brg);
%fprintf(1, 'Displacement %.2f meters\n', displacement);
%set(handles.textResults, 'String', '');
strResults = sprintf('Avg speed = %.4f cm/s\nAvg bearing %.2f deg True\nDisplacement %.2f meters\n', ave_spd/cm2m, ave_brg, displacement);
updateResults(handles, strResults);

net_glide_slope = 2*dive_depth/displacement;

% do the plots 

% Use slightly-modified version of UW plotting code.
% Renamed to dd_diveplot_graphs to ensure we don't mistakenly
% pickup the unmodified stock UW version if path is messed up
dd_diveplot_graphs

close(hWaitbar);


%============================================================
% get_zscl()
%============================================================
function [zscl] = get_zscl(results)

    maxd = max(results.depth);
    % use 1100 meters for SG so if we go a little over nomimal depth
    % we don't jump to 2km scale
    if maxd > 1100
      zscl = fix(maxd/100); % 1000 - 6000m
      % round up to the nearest 10
      zscl = (zscl - rem(zscl,10)) + 10;
    elseif maxd > 500
      zscl = 10; % 500 - 1000m
    elseif maxd > 300
      zscl = 5; % 300 - 500m
    elseif maxd > 200
      zscl = 3; % 200 - 300m
    elseif maxd > 100
      zscl = 2; % 100 - 200m
    else
      zscl = 1; % 0 - 100m
end



