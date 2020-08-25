% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

% the protocol:

% get_dive_data() returns loginfo, eng, and results structures with 'raw'
% data and results consumers should NOT modify these structures but copy
% into variables to manipulate.  Much of this copying into variables happens
% in unpack_data, which see.  Then, optionally, this script,
% save_dive_results, takes those variables and centrally updates the results
% structure and writes it to a mat file.  It updates the variables
% appropriately so the results look as if they were written to an nc file
% ala the basestation.

% This file is deliberately a script in order to have access to results of
% dive_process_funcs.  Assumes mat_file is also available...

% DEAD results = []; % do NOT reset in case it came from an nc file; save its results we don't update

% emitted to nc file
results.pressure = sg_press;
results.depth = sg_depth_m;

results.avg_latitude = avg_lat; % average latitude between GPS2 (dive) and GPS (surfacing)
results.magnetic_variation = mag_variation; % assumed magnetic variation
results.start_of_climb_time = start_of_climb_time;
% directives has other vectors and scalars of interest
% clean directives of redundant information
directives = rmfield(directives,{'eval_debug','make_suggestions','data_points','time','depth','glider_data_points'});
results.directives = directives; % is this what we want?

% We don't save raw data here but they might be ctd_np or mp long so separate them if so...

% TODO ensure these are all ctd_np long
% NOTE: just after dive_process, time and ctd_time can differ!!
% graphing should always use ctd_time...
results.latitude  = lat; % latitude
results.longitude = lon; % longitude

results.time = sg_epoch_time; % epoch time in secs since start of dive
results.ctd_time = ctd_epoch_time; % epoch time in secs since start of dive
results.ctd_depth = ctd_depth_m; % [m] the depth of the CT orifice

%DEAD results.ctd_sg_depth_m = ctd_sg_depth_m; % [m] the measured depth
results.hdm_qc = num2str(hdm_qc);
results.speed_gsm = speed_gsm;
results.glide_angle_gsm = glideangle_gsm; % degrees
results.horz_speed_gsm = speed_gsm.*cos( glideangle_gsm*radpd);
results.vert_speed_gsm = speed_gsm.*sin( glideangle_gsm*radpd);

results.CTD_qc = num2str(CTD_qc);
results.temperature_raw = temp_raw;
results.temperature_raw_qc = num2str(temp_raw_qc');
results.conductivity_raw = cond_raw;
results.conductivity_raw_qc = num2str(cond_raw_qc');
results.salinity_raw = salin_raw;
results.salinity_raw_qc = num2str(salin_raw_qc');

results.temperature = temp;
results.temperature_qc = num2str(temp_qc');
results.conductivity = cond;
results.conductivity_qc = num2str(cond_qc');
results.salinity = salin;
results.salinity_qc = num2str(salin_qc');
% results.salinity_error = salinity_error;

results.speed = spd_model;
results.glide_angle = glideangle_model; % degrees
results.horz_speed = spd_model.*cos( glideangle_model*radpd);
results.vert_speed = w_model;
results.speed_qc = num2str(speed_qc');

results.east_displacement = east_displacement; 
results.north_displacement = north_displacement;
% update u_av, v_av


results.density = density;
results.buoyancy = buoy;
results.sigma_t = sigma_t;
results.theta = theta;
results.sigma_theta = sigma_theta;
results.sound_velocity = sound_velocity;
results.dissolved_oxygen_sat = solubility_micromoles_per_kg;

results.depth_avg_curr_qc = num2str(DAC_qc);
if (DAC_qc == QC_GOOD)
  results.depth_avg_curr_east  = u_h2o; % m/s east
  results.depth_avg_curr_north = v_h2o; % m/s north
  results.flight_avg_speed_east  = u_av;
  results.flight_avg_speed_north = v_av; 
end
%TODO results.surface_curr_qc = num2str(surface_curr_qc);

% Sensor results
if (sbe43_present)
  results.SBE43_qc = num2str(sbe43_qc);
  results.sbe43_results_time = sbe43_time;
  results.sbe43_dissolved_oxygen = oxygen;
  results.sbe43_dissolved_oxygen_qc = num2str(oxygen_qc');
end

if (aa3830_present)
  results.aanderaa3830_qc = num2str(optode_qc');
  results.aanderaa3830_results_time = optode_time;
  results.aanderaa3830_instrument_dissolved_oxygen = optode_oxygen;
  results.aanderaa3830_dissolved_oxygen = optode_dphase_oxygen;
  results.aanderaa3830_dissolved_oxygen_qc = num2str(optode_oxygen_qc');
end

if (aa4330_present)
  results.aanderaa4330_qc = num2str(optode_qc');
  results.aanderaa4330_results_time = optode_time;
  results.aanderaa4330_instrument_dissolved_oxygen = optode_oxygen;
  results.aanderaa4330_dissolved_oxygen = optode_tcphase_oxygen;
  results.aanderaa4330_dissolved_oxygen_qc = num2str(optode_oxygen_qc');
end

if (wetlabs_present)
  %DEAD results.wetlabs_qc = wetlabs_qc;
end


% These are not emitted in nc file but some matlab scripts want
results.midpoint_day = midpoint_day;
results.midpoint_lat = midpoint_lat;
results.midpoint_lon = midpoint_lon;
results.iter = iter;
results.converged = converged;
results.day = day; % the fraction of days since 1 Jan 0000 (sdn)

% if you use -ascii then you need to use a different load command to assign to proper vars
% also directives, as a structure, is not dumped
%DEAD save_options = '-ascii';
save_options = '-mat'; % default; use load directly

save(mat_file,save_options,'results');
