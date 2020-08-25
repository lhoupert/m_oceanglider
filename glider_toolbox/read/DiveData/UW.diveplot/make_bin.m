% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% map over result files and create a bin'd structure
% this version uses nanmean rather than linear interpolation, which can generate extreme values (see interpolate_data.m for explanation)

% define bin depths
% bin_depth = [0:2:150 155:5:300 310:10:10000]'; % PARAMETER
% This matches MakeDiveProfiles (actually MDP forces to 2m bins?)
bin_depth = [0:5:6000]'; % PARAMETER

dbin = diff(bin_depth);
nbins = length(bin_depth);
bin_top = [0; bin_depth(1:nbins-1)];
bin_bot = [bin_depth(2:nbins); bin_depth(nbins) + 0.5*dbin(end)];

Bin = [];
DAC = [];

qc_declarations
clear_sg_config_constants % clear the top workspace
sg_calib_constants
sg_config_constants

fprintf(1,'SG%s: ', id_str);
available_runs = available_profiles();
if (exist('dive_process_runs','var'))
  % A script called us and supplied which runs to work on...
  runs = intersect(dive_process_runs,available_runs); % we could have skipped dives
  clear dive_process_runs;
else
  runs = ask_which_runs(available_runs);
end
nruns = length(runs);
% loop on runs
tic;
for ir = 1:nruns;
  irun = runs(ir);
  id_dive = sprintf('%s%04d',id_str,irun);
  irun = double(irun); % ensure arrays are all double, not int32!
  base_file = sprintf('p%s',id_dive);
  [loginfo,eng,results] = get_dive_data(id_dive,1); % require results since we look at sigma and DAC, etc.
  if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
    fprintf(1,'Processing error or skipped profile for %d; skipping\n',irun);
    continue;
  end
  unpack_data;
  fprintf(1,'%s\n',base_file); % DBEUGGING
  
  high_salinity_error = zeros(sg_np,1);
  high_i = find(abs(salinity_error) > 0.01); % PARAMETER
  high_salinity_error(high_i) = 1;
  suspect_thermal_points = zeros(sg_np,1);
  if (isfield(directives,'suspect_thermal_inertia')) % DEAD -- in new versions this is always computed
    suspect_thermal_points(directives.suspect_thermal_inertia) = 1;
  end

  % recompute...
  density = sw_dens0(salin, temp);
  sigma_t = density - 1000;
  max_depth = max(ct_depth_m);
  
  DAC = [DAC;[irun;avg_lat;u_h2o;v_h2o;u_av;v_av;midpoint_day;midpoint_lat;midpoint_lon]'];
  % Fit profiles within each bin 
  % Linear fit within bin if 2 or more pts, otherwise a constant.
  % tic
  nfit = 1;   % linear fit with depth
  Bin_d = [];
  Bin_c = [];
  xbins = find(bin_top <= max_depth); % just the bins that matter to this dive
  for jb = xbins'
    % all points on dive and climb in the depth bin
    inbin = find(ct_depth_m >= bin_top(jb) & ct_depth_m < bin_bot(jb));
    
    % fit dive bins
    id = intersect(inbin,dive_i); % the dive points
    if (~isempty(id))
      % we will add info to this bin...
      pitch_sign = -1;
      
      bin_d_salin_err = length(find(high_salinity_error(id)));
      bin_d_suspect_ti = length(find(suspect_thermal_points(id)));
      
      bin_d_temp = nanmean(temp(id));
      bin_d_salin = nanmean(salin(id));
      bin_d_sigma_t = nanmean(sigma_t(id));
      bin_d_lon = nanmean(lon(id));
      bin_d_lat = nanmean(lat(id));
      bin_d_day = nanmean(day(id));
      bin_d = [irun; pitch_sign; bin_d_day; bin_d_lon; ...
               bin_d_lat; bin_depth(jb); bin_d_temp; bin_d_salin; ...
               bin_d_sigma_t; bin_d_salin_err; bin_d_suspect_ti];

      if (sbe43_present)
        bin_d_oxygen = nanmean(oxygen(id));
        % solubility_micromoles_per_kg as well?
        bin_d = [bin_d; bin_d_oxygen];
      end
      if(optode_present) % start fit optode dive bins
          bin_d_optode_oxygen = nanmean(optode_oxygen(id));
          bin_d_optode_dphase_oxygen = nanmean(optode_dphase_oxygen(id));
          bin_d = [bin_d; bin_d_optode_oxygen; bin_d_optode_dphase_oxygen];
      end % end fit optode dive bins
      if (wetlabs_present)
        bin_d_red_scttr = nanmean(red_scttr(id));
        bin_d_blue_scttr = nanmean(blue_scttr(id));
        bin_d_fluor = nanmean(fluor(id));
        bin_d = [bin_d; bin_d_red_scttr; bin_d_blue_scttr; bin_d_fluor];
      end
      
      Bin_d = [Bin_d; bin_d'];
    end
    
    % fit climb bins
    ic = intersect(inbin,climb_i); % the climb points
    if (~isempty(ic))
      % we will add info to this bin...
      pitch_sign = 1.;
      
      bin_c_salin_err = length(find(high_salinity_error(ic)));
      bin_c_suspect_ti = length(find(suspect_thermal_points(ic)));
      
      bin_c_temp = nanmean(temp(ic));
      bin_c_salin = nanmean(salin(ic));
      bin_c_sigma_t = nanmean(sigma_t(ic));
      bin_c_lon = nanmean(lon(ic));
      bin_c_lat = nanmean(lat(ic));
      bin_c_day = nanmean(day(ic));
      bin_c = [irun; pitch_sign; bin_c_day; bin_c_lon; ...
               bin_c_lat; bin_depth(jb); bin_c_temp; bin_c_salin; ...
               bin_c_sigma_t; bin_c_salin_err; bin_c_suspect_ti];
      
      if (sbe43_present)
        bin_c_oxygen = nanmean(oxygen(ic));
        % solubility_micromoles_per_kg as well?
        bin_c = [bin_c; bin_c_oxygen];
      end
      if(optode_present) % start fit optode dive bins
          bin_c_optode_oxygen = nanmean(optode_oxygen(ic));
          bin_c_optode_dphase_oxygen = nanmean(optode_dphase_oxygen(ic));
          bin_c = [bin_c; bin_c_optode_oxygen; bin_c_optode_dphase_oxygen];
      end % end fit optode dive bins
      if (wetlabs_present)
        bin_c_red_scttr = nanmean(red_scttr(ic));
        bin_c_blue_scttr = nanmean(blue_scttr(ic));
        bin_c_fluor = nanmean(fluor(ic));
        bin_c = [bin_c; bin_c_red_scttr; bin_c_blue_scttr; bin_c_fluor];
      end
      
      Bin_c = [Bin_c; bin_c'];
    end
  end % end loop on bins jb
  
  % Append overall data vectors
  if(~isempty(Bin_d))
    Bin = [Bin; Bin_d];
  end
  if(~isempty(Bin_c))
    Bin = [Bin; flipud(Bin_c)]; % preserve time order
  end
  % toc
end % end ir loop on runs

Bin_vars = {'dive_num','pitch_sign','day','lon','lat','depth','temp','salin','density','salin_error','suspect_ti'};
if (sbe43_present)
  Bin_vars{end+1} = 'oxygen';
end
if (optode_present)
  Bin_vars{end+1} = 'optode_oxygen';
  Bin_vars{end+1} = 'dphase_oxygen';
end
if (wetlabs_present)
  Bin_vars{end+1} = 'red_scttr';
  Bin_vars{end+1} = 'blue_scttr';
  Bin_vars{end+1} = 'fluor';
end
save Bin sbe43_present optode_present wetlabs_present Bin_vars Bin;
save DAC DAC;
fprintf(1,'Wrote Bin.mat and DAC.mat\n');
toc % report Elapsed time
