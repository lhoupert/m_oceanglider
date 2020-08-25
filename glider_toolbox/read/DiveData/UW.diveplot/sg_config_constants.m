% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% the set of default configuration constants for different versions of
% the Seaglider and Deepglider

% the calling sequence should be
% clear sg_configuration % no hang overs from previous runs
% sg_calib_constants % load constants and overrides
% sg_config_constants % load default configuration constants

% First, an overall configuration parameter from which a cascade of
% definitions can flow by default.  These definitions are broken into
% various sections, each with their own grouping parameter.  This permits
% easy mix-and-match and targetted override to the individiaul parameter
% basis.  These overrides should be defined in the sg_calib_constants
% file.

% Setting this variable should provide a set of defaults for all others
if (~exist('sg_configuration','var')) 
  % Assume your stock, old style, Seaglider
  sg_configuration = 0;
  % there are no 'octal' problems parsing '019' as 19 under matlab...
  if (exist('id_str','var') & str2num(id_str) >= 105)
    % when SFC started at 100, all gliders post 105 have been gun-style
    % the hull numbers for SFC can run to 499 (500+ are iRobot numbers?)
    % CCE's original SG's ran to 23, then started to purchase from SFC starting at 100
    % DG's started at 30.  30-32 were SG fairing style DGs w/ aft compresseee
    % 33+ are cylinderical DGs with fwd and aft compressee
    % all DGs use gun-style CTs
    % Assume stock Seaglider with new gun-style CT
    sg_configuration = 1;
  end
end
% assignments that are vehicle-independent

% data biases and offsets, which sg_calib_constants can override

% setting to non-zero often helps catch the vehicle slowing down into a stall and accelerating to flight
% it also helps with sitting on the bottom bouncing around...
if (~exist('min_stall_speed','var')) min_stall_speed = 1.0; end % cm/s (stalled if speed below this)
% faroes/jun08/sg005 dive 135 attained 1.2m/s speeds in upwelling
if (~exist('max_stall_speed','var')) max_stall_speed = 100.0; end % cm/s (stalled if speed above this...salinity spikes? huge downwelling)
if (~exist('min_stall_angle','var')) min_stall_angle = 5.0; end % [degreees] stalled if abs(pitch) is less than stall angle (0.0 defeats this check)

if (~exist('vbdbias','var')) vbdbias = 0.0; end  % [cc]
if (~exist('sbe_temp_freq_offset','var')) sbe_temp_freq_offset = 0; end % temp frequency offset
if (~exist('temp_bias','var')) temp_bias = 0; end  % [deg]
if (~exist('sbe_cond_freq_offset','var')) sbe_cond_freq_offset = 0; end % cond frequency offset
if (~exist('cond_bias','var')) cond_bias = 0; end  % [S/m] % conductivity bias
if (~exist('pitchbias','var')) pitchbias = 0; end % [deg] % pitch sensor bias
if (~exist('depth_bias','var')) depth_bias = 0; end % [m] depth bias (in meters) because of a flakey or mis-tared pressure sensor

switch (sg_configuration)
 case 0 % stock SG with original CT
  if (~exist('sg_ct_type','var')) sg_ct_type = 0; end
  if (~exist('sg_vehicle_geometry','var')) sg_vehicle_geometry = 0; end
  if (~exist('sg_sensor_geometry','var')) sg_sensor_geometry = 0; end
  if (~exist('sg_ct_geometry','var')) sg_ct_geometry = 0; end
 
 case 1 % stock SG new gun style CT
  if (~exist('sg_ct_type','var')) sg_ct_type = 1; end
  if (~exist('sg_vehicle_geometry','var')) sg_vehicle_geometry = 0; end
  if (~exist('sg_sensor_geometry','var')) sg_sensor_geometry = 1; end
  if (~exist('sg_ct_geometry','var')) sg_ct_geometry = 1; end
 
 case 2 % DG
  if (~exist('sg_ct_type','var')) sg_ct_type = 1; end
  if (~exist('sg_vehicle_geometry','var')) sg_vehicle_geometry = 1; end
  if (~exist('sg_sensor_geometry','var')) sg_sensor_geometry = 2; end
  if (~exist('sg_ct_geometry','var')) sg_ct_geometry = 1; end
  
 case 3 % stock SG with pumped GPCTD
  if (~exist('sg_ct_type','var')) sg_ct_type = 2; end
  if (~exist('sg_vehicle_geometry','var')) sg_vehicle_geometry = 0; end
  if (~exist('sg_sensor_geometry','var')) sg_sensor_geometry = 1; end
  if (~exist('sg_ct_geometry','var')) sg_ct_geometry = 0; end % irrelevant
 
 otherwise % unknown value specified
  fprintf(1,'Unknown Seaglider configuration %d; assuming stock SG with original CT mount!\n',sg_configuration);
  sg_configuration = 0; % assume old
  if (~exist('sg_ct_type','var')) sg_ct_type = 0; end
  if (~exist('sg_vehicle_geometry','var')) sg_vehicle_geometry = 0; end
  if (~exist('sg_sensor_geometry','var')) sg_sensor_geometry = 0; end
  if (~exist('sg_ct_geometry','var')) sg_ct_geometry = 0; end
end

% We would love to be able to report in a function
% like overridden() that a variable was overridden
% but the function call masks all the variables in the
% outer scope (workspace).  Thanks MATLAB.

% overall vehicle geometery (sg_vehicle_geometry)
switch (sg_vehicle_geometry)
 case 0 % stock SG
  if (~exist('glider_length','var')) glider_length = 1.8; end
  if (~exist('hd_s','var')) hd_s =  -0.25; end % how the drag scales by shape (-1/4 for SG per Eriksen, et al.)

  % Define parameters that influence buoyancy
  % set this to 0 to disable interstitial calculations
  if (~exist('glider_interstitial_length','var')) glider_interstitial_length = 0.2; end %[m]
  if (~exist('glider_interstitial_volume','var')) glider_interstitial_volume = 12e-3; end %[m3]
  if (~exist('glider_r_en','var')) glider_r_en =  0.00635; end % entry radius [m]
  if (~exist('mass_comp','var')) mass_comp = 0; end % mass of compresseee, if any [kg]
  
  if (~exist('glider_wake_entry_thickness','var')) glider_wake_entry_thickness = 0.; end % wake entry region thickness [m]
  if (~exist('glider_vol_wake','var')) glider_vol_wake =  18e-3; end % attached wake volume [m^3]
  % fairing radius
  if (~exist('glider_r_fair','var')) glider_r_fair =  0.3; end
  % thermal expansion of the hull
  if (~exist('therm_expan','var')) therm_expan =  70.5e-6; end % m^3/degree
  % This is the reference temperature against which the therm_expand coefficient was measured
  if (~exist('temp_ref','var')) temp_ref =  15; end
  % compressibility of the hull per dbar
  if (~exist('abs_compress','var')) abs_compress =  4.4e-6; end % m^3/dbar (the compressibility of SW)
 
 case 1 % 'std' DG
  % assume ogive nose and tail, right regular cylinder DG33++
  % TODO verify these with CCE
  if (~exist('glider_length','var')) glider_length = 1.74; end
  if (~exist('hd_s','var')) hd_s =  -0.25; end % how the drag scales by shape (-1/4 for SG per Eriksen, et al.)
  
  % DG has no flow-through volume to speak of...
  if (~exist('glider_interstitial_length','var')) glider_interstitial_length = 0.0; end %[m]
  if (~exist('glider_interstitial_volume','var')) glider_interstitial_volume = 0; end %[m3]
  if (~exist('glider_r_en','var')) glider_r_en =  0.00635; end % entry radius [m]
  % can't float or fly w/o it
  if (~exist('mass_comp','var')) mass_comp = 10; end % mass of compresseee, if any [kg]
  
  if (~exist('glider_wake_entry_thickness','var')) glider_wake_entry_thickness = 0.; end % wake entry region thickness [m]
  if (~exist('glider_vol_wake','var')) glider_vol_wake =  18e-3; end % attached wake volume [m^3]
  % fairing radius
  if (~exist('glider_r_fair','var')) glider_r_fair =  0.3; end
  % thermal expansion of the hull
  if (~exist('therm_expan','var')) therm_expan =  10e-6; end % m^3/degree
  % This is the reference temperature against which the therm_expand coefficient was measured
  if (~exist('temp_ref','var')) temp_ref =  15; end
  % compressibility of the hull per dbar
  if (~exist('abs_compress','var')) abs_compress =  2.45e-6; end % m^3/dbar Boeing hull
  % tstar = 4420;      pstar = 35800;      vstar = 0.9945;    % compressee parameters
  % tstar = 4516.3746; pstar = 35364.4304; vstar = 0.99909; % 25-94C fit
  % tstar = 4587.3064; pstar = 31674.2091; vstar = 1.0067;  % 25-37C fit
end

% TODO these values need to be set fo the GPCTD in final configuration
% vehicle sensor geometry (sg_sensor_geometry)
switch (sg_sensor_geometry)
 case 0 % original SG pressure location and CT sail
  if (~exist('glider_xT','var')) glider_xT =  -1.1800; end % glider x coord of thermistor tip [m]
  if (~exist('glider_zT','var')) glider_zT =   0.1700; end % glider z coord of thermistor tip [m]
  if (~exist('glider_xP','var')) glider_xP =  -0.6870; end % glider x coord of pressure gauge [m]
  % JSB 9/9/9 -0.125 on SG(!);
  if (~exist('glider_zP','var')) glider_zP =  -0.0254; end % glider z coord of pressure gauge [m]
 
 case 1 % original SG pressure location with new gun CT sail
  % Gun measurements assume the stalk location for the old and new CT are the same
  % JSB 9/9/9 Gun: thermistor 0.07m further aft than SG
  if (~exist('glider_xT','var')) glider_xT =  -1.2500; end % glider x coord of thermistor tip [m]
  if (~exist('glider_zT','var')) glider_zT =   0.1650; end % glider z coord of thermistor tip [m]
  if (~exist('glider_xP','var')) glider_xP =  -0.6870; end % glider x coord of pressure gauge [m]
  % JSB 9/9/9 -0.125 on SG(!);
  if (~exist('glider_zP','var')) glider_zP =  -0.0254; end % glider z coord of pressure gauge [m]
 
 case 2 % DG pressure location with new gun CT sail
  % The gun style sail is mounted forward and the pressure sensor is in
  % the nose of the end cap behind the ogive fairing
  % VERIFY THE THERMISTOR LOCATION OFFSETS using ogive
  if (~exist('glider_xT','var')) glider_xT =  -1.2500; end % glider x coord of thermistor tip [m]
  if (~exist('glider_zT','var')) glider_zT =   0.1650; end % glider z coord of thermistor tip [m]
  if (~exist('glider_xP','var')) glider_xP =  -0.3900; end % glider x coord of pressure gauge [m]
  if (~exist('glider_zP','var')) glider_zP =  -0.0000; end % glider z coord of pressure gauge [m]
end

switch (sg_ct_type)
 case 0 % std Seabord unpumped CTD
  if (~exist('sbect_unpumped','var')) sbect_unpumped =  1; end   
 case 1 % gun-style  Seabord unpumped CTD
  if (~exist('sbect_unpumped','var')) sbect_unpumped =  1; end 
 case 2 % pumped Seabird GPCTD
  if (~exist('sbect_unpumped','var')) sbect_unpumped =  0; end 
  % Note: SBE9 u_f is 1.75m/s (see Morison d'Asaro)
  % CCE has 30ml/s in our tube as 2.38m/s flow rate
  % Creed and Janzen report GPCTD runs at 10ml/s
  % the thermistor sting sits proud in the intake which is ~30mm vertical and ~10mm forward from the start of the tube, so tube is ~4cm removed from thermistor in all orientations
  if (~exist('sbect_gpctd_x_d','var')) sbect_gpctd_x_d = 0.03; end % length of duct [m[ before the tube starts
  if (~exist('sbect_gpctd_u_f','var')) sbect_gpctd_u_f =  2.38*(10/30); end % tube flow speed [m/s] for continuous pumped CTD (personal communication w/ SBE) (was 0.9151)
  if (~exist('sbect_gpctd_tau_1','var')) sbect_gpctd_tau_1 = sbect_gpctd_x_d/sbect_gpctd_u_f; end % time delay [s] between thermistor and mouth of conductivity tube assume 4cm distance
 case 3 % unpumped Seabird SAILCT
  if (~exist('sbect_unpumped','var')) sbect_unpumped =  1; end % assume we have an unpumped SBE41
end


% CT type (sg_ct_type) and construction (sg_ct_geometry)
% TODO assume but verify that this applies to ALL CTDs from SBE, including the GPCTD
if (~exist('sbect_tau_T','var')) sbect_tau_T =  0.6; end % thermistor response[s]
if (sbect_unpumped)
  % all these constants apply to the unpumped SBE41 for its corrections
  if (~exist('sbect_x_m,','var')) sbect_x_m =  0.0087; end % length of mouth portion of cell [m]
  if (~exist('sbect_r_m','var')) sbect_r_m =  0.0081; end % radius of mouth portion of cell [m]
  if (~exist('sbect_cell_length','var')) sbect_cell_length = 0.09; end % combined length of 2 narrow (sample) portions of cell [m]
  if (~exist('sbect_x_w','var')) sbect_x_w = 0.0386; end	% length of wide portion of cell [m]
  if (~exist('sbect_r_w','var')) sbect_r_w = 0.0035; end	% radius of wide portion of cell [m]
  if (~exist('sbect_r_n','var')) sbect_r_n = 0.002; end 	% radius of narrow portion of cell [m]

  % encode information about the glass/expoxy response
  if (~exist('sbect_u_r','var')) sbect_u_r =  10^0; end % thermal inertia response rolloff speed [m/s]

  % assume interior potting is the same between models and mounts (bounded between 7.5 and 9.1 secs).  This is 10^.9
  % if (~exist('sbect_tau_w_min','var')) sbect_tau_w_min = 7.9433; end % thermal response time of glass/epoxy tube [s]
  % after changing the eqn in tsv_iter to use the proper form from Schlicting we need to scale tau_w_min by 1.44
  if (~exist('sbect_tau_w_min','var')) sbect_tau_w_min = 7.9433/1.44; end % thermal response time of glass/epoxy tube [s]

  % Define the (variable) geometry of the CT sail components that affect the T/S corrections
  switch (sg_ct_geometry)
   case 0 % original CT sail
    if (~exist('sbect_x_T','var')) sbect_x_T =  -0.014; end % cell mouth to thermistor x offset[m]
    if (~exist('sbect_z_T','var')) sbect_z_T =  -0.015; end % cell mouth to thermistor z offset[m]
    if (~exist('sbect_C_d0','var')) sbect_C_d0 =  1.2; end % cell mouth drag coefficient
    
   case 1 % gun CT sail
    if (~exist('sbect_x_T','var')) sbect_x_T =  -0.011; end % cell mouth to thermistor x offset[m]
    if (~exist('sbect_z_T','var')) sbect_z_T =  -0.030; end % cell mouth to thermistor z offset[m]
    if (~exist('sbect_C_d0','var')) sbect_C_d0 =  2.4; end % cell drag mouth coefficient
  end
else
  % GPCTD numbers
  % TODO verify these numbers, measurements, which are copies of SBE41 above and likely incorrect
  if (~exist('sbect_x_m,','var')) sbect_x_m =  0.0087; end % length of mouth portion of cell [m]
  if (~exist('sbect_r_m','var')) sbect_r_m =  0.0081; end % radius of mouth portion of cell [m]
  if (~exist('sbect_cell_length','var')) sbect_cell_length = 0.09; end % combined length of 2 narrow (sample) portions of cell [m]
  if (~exist('sbect_x_w','var')) sbect_x_w = 0.0386 + sbect_gpctd_x_d; end	% length of wide portion of cell [m] with intake duct
  if (~exist('sbect_r_w','var')) sbect_r_w = 0.0035; end	% radius of wide portion of cell [m]
  if (~exist('sbect_r_n','var')) sbect_r_n = 0.002; end 	% radius of narrow portion of cell [m]

  % encode information about the glass/expoxy response
  if (~exist('sbect_u_r','var')) sbect_u_r =  10^0; end % thermal inertia response rolloff speed [m/s]

  % assume interior potting is the same between models and mounts (bounded between 7.5 and 9.1 secs).  This is 10^.9
  % if (~exist('sbect_tau_w_min','var')) sbect_tau_w_min = 7.9433; end % thermal response time of glass/epoxy tube [s]
  % after changing the eqn in tsv_iter to use the proper form from Schlicting we need to scale tau_w_min by 1.44
  if (~exist('sbect_tau_w_min','var')) sbect_tau_w_min = 7.9433/1.44; end % thermal response time of glass/epoxy tube [s]

  if (~exist('sbect_x_T','var')) sbect_x_T =  -0.014; end % cell mouth to thermistor x offset[m]
  if (~exist('sbect_z_T','var')) sbect_z_T =  -0.015; end % cell mouth to thermistor z offset[m]
  if (~exist('sbect_C_d0','var')) sbect_C_d0 =  1.2; end % cell mouth drag coefficient
end
