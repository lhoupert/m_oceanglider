% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%

% setup QC, configuration and calibration constants for dive_process use

if (0)
  % TODO reference times and locations; these should be in sg_calib_constants.m
  % defaults: 0 everywhere...
  ref_serial_date_time = datenum(2000,1,1);   % 0000Z 1 Jan 2000 reference for WA coast deployments
  ref_lat = 50;
  ref_lon = -145;
end

qc_declarations
clear_sg_config_constants % clear the top workspace
sg_calib_constants
% if you are running converted dives in a v66 subdir, expose this 
% v66_sg_calib_constants
sg_config_constants

fprintf(1,'SG%s: ', id_str);
switch sg_vehicle_geometry
 case 0
  vehicle_type = 'Seaglider';
 case 1
  vehicle_type = 'Ogive Deepglider';
end
switch sg_ct_type
 case 0
  ctd_type = 'original unpumped Seabird CTD';
 case 1
  ctd_type = 'gun-style unpumped Seabird CTD';
 case 2
  ctd_type = 'pumped GPCTD';
 case 3
  ctd_type = 'unpumped SAILCT';
end
switch sg_sensor_geometry
 case 0
  ctd_mount = 'aft';
 case 1
  ctd_mount = 'aft';
 case 2
  ctd_mount = 'forward';
end
fprintf(1, '%s with %s mounted %s.\n',vehicle_type,ctd_type,ctd_mount);


if (rho0 < 1.1) % max in g/cc is 1.041 (Red Sea)
  % sg_calib_constants.m provides rho0, expressed in g/m^3
  % The glider parameter (and log file) provide $RHO, in g/cc
  % It can sometimes happen (SG524 under GPCTD tests) that the sg_calib_constant value
  % is expressed in g/cc, as though copied from $RHO.
  % Warn and correct, otherwise this will produce speeds in the 1-10m/s range, among other problems
  % In any case, we ignore $RHO in this script since it used used by the glider to compute operational
  % buoyancy-based results our speed guesses could be different below
  % NOTE comment in tsv_iter about redoing this change, sigh.
  fprintf(1,'NOTE: Correcting rho0 in sg_calib_constants, which is %g and less than 1000 g/m^3!\n',rho0);
  rho0 = rho0*1000; % assume they gave us g/cc; convert to g/m^3
end

delta_T0 = 1; % thermal boundary layer scale (dimensionless?)
C_d0 = sbect_C_d0; % cell mouth drag coefficient (dimensionsless?)
% tau_w_min, the minimum time to stablize a thermal pulse in the
% conductivity tube, varies according to the material diffusivity of the
% tube (glass) and the surrounding epoxy tau_w_min also varies with how the
% drag of the cell mouth increases wrt to the 'standard' CT arrangement.
% delta_T0 and u_r (typical flow) are order 1 and so contribute little.
% the form of this expression is constructed so that, in tsv_iter, tau_w*alpha will be a constant in case u_r or delta_T0 vary.
tau_w_min = sbect_tau_w_min/(sqrt(sbect_u_r*1.2/C_d0)*delta_T0);  % minimum tau_w [s]

radpd = pi/180; % radians per degree
degpr = 180./pi; % degress per radian
cmpm = 100; % cm per meter
m_per_deg = 1852*60; % meters per degree of lat/lon (1852 meters per nautical mile; 60 nm per degree)
dbar_per_psi = 1.0/1.4503774;
psi_per_meter = 1.0/0.685; % this needs to match the value in glider/constants.h
c3515 = sw_c3515*0.1; % conductivity for 35PSU at 15C
Kelvin_offset = 273.15; % for 0 deg C
