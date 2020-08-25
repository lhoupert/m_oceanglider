function [regress_dives] = dd_regress_vbd_dives(handles, pathname, doPlot, settings)
% Parameters:
%   handles  - Application handles (GUI)
%   pathname - Path where (.nc) files are located
%   doPlot   - If set to 1 then plot
%
% Returns:
%   regress_dives - Array of dives that are recommended for running 
%      regression.
%

existingPath = pwd;
cd(pathname);

if exist('doPlot', 'var') == 0
    doPlot = 0;
end


% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%
% History:
%   6/2014 RTK Incorporated script into SDDA, renamed from
%   regress_vbd_dives.m
%
hWaitBar = waitbar(0, 'Determining best dives for VBD regression...');

msg = sprintf('%s%sDetermining best dives for running regression (dd_regress_vbd_dives)...', char(10), char(10));
updateResults(handles, msg);  

% Find a good set of dives for regress_vbd.m
% For vbdbias we want a set of deep dives
% For abc we want dives that span the buoyancy and pitch desired space to fit the power function properly
% For fairing/hull compression and thermal response we want max temp and pressure ranges
most_dives_weight = 5; % PARAMETER the minimum number of 'most frequently performed dives' to weight
clear_sg_config_constants % clear the top workspace
%sg_calib_constants;  % billr: this blows up the compiled code
sg_config_constants;
qc_declarations;
m2cm = 100;

% Ask for a range to limit search
runs = available_profiles({'nc'}, settings.seagliderID); % billr: added settings.seagliderID to replace global id_str formerly created by sg_calib_constants
nruns = length(runs);

msg = sprintf('Searching %d dives from SG%d %s for good VBD regression candidates...\n',...
		nruns,settings.seagliderID,underscore(settings.mission))
updateResults(handles, msg);    

ir = 1;
data = [];
bad_runs = [];
for ix = 1:nruns;
waitbar(ix / nruns, hWaitBar, sprintf('%.0f%% complete', (ix / nruns) * 100.0));
  irun = runs(ix);
  id_dive = sprintf('%03d%04d',settings.seagliderID,irun); % billr
  base_file = sprintf('p%s',id_dive);
  
  [loginfo,eng,results] = get_dive_data(id_dive,1); % esnure there are results since regress_vbd requires speed, etc.
  
  if (isfield(results,'processing_error') || isfield(results,'skipped_profile') ...
	  || ~ isfield(results,'hdm_qc'))
    msg = sprintf('Processing error or skipped profile for %d; skipping\n',irun);
    updateResults(handles, msg);  
    continue;
  end
  
  if (str2num(results.hdm_qc) == QC_BAD)
	  % Oh course, this is likely to happen because you haven't gotten the vbd regression right!
      msg = sprintf('NOTE: Bad hydrodynamic speed QC for dive %d\n',irun);
      updateResults(handles, msg);  
  end
  
  unpack_data; % get indices
  if (isempty(dive_i))
    msg = sprintf('No points on dive %d; skipping\n',irun);
    updateResults(handles, msg);  
    continue;
  end
  
  if (isempty(climb_i))
    msg = sprintf('No points on climb %d; skipping\n',irun);
    updateResults(handles, msg);  
    continue;
  end
  
  % TODO -- skip dives used for in-situ compass cal (constantly rolled)
  if (0)
	  w = m2cm.*ctr1stdiffderiv(-ctd_depth_m,ctd_time);
	  if (nanmean(abs(w - results.vert_speed)) > 5)
		  % see if the basestation predicted an anomalously large
		  % w difference from observed to the model.  if so
		  % this is a possible indication that volmax was set improperly
		  % and hence the vbdbias regression may fail (because
		  % the model fails to compute the buoyancy properly and
		  % that causes the regression to spin out...)
		  % See OKMC Mar13 sg177 23:25
		  % There were few to no pumps on the climb and initial volamx
		  % was incorrect by ~80cc which lead it to compute very slow speeds
		  
		  % what is interesting is that if it is off that much with those problems
		  % the regression does NOT find that vbdbias...why not?
	  	  
		  % another idea would be make a better, coarser estimate of vbdbias
		  % and provide that as vbdbias_0
		  bad_runs = [bad_runs; irun];
		  continue;
	  end
  end
  max_buoy = loginfo.MAX_BUOY;
  pitch_desired = loginfo.MHEAD_RNG_PITCHd_Wd(3);
  data = [data; [irun max_buoy pitch_desired]];
end
close(hWaitBar);

if doPlot
    hFigPlot = figure('Name','Pitch');
    set(hFigPlot, 'UserData', 'Plot_Regress_VBD_Dives'); 
    hold on;
    plot(data(:,1), data(:,2)/100.);
    plot(data(:,1), data(:,3), 'g');
    grid on;
    lg = legend('Max Buoy', 'Pitch Desired');
    set(lg,'color','none'); % make transparent
end

if (~isempty(bad_runs))
	msg = sprintf('Observed vs. computed w for dives %s substantially disagree; skipping those.\n', succinct_elts(bad_runs));
    updateResults(handles, msg);  
end

unique_max_buoy = unique(data(:,2));
unique_pitch_desired = unique(data(:,3));

% avoid too shallow dives where Kalman probably messed up
unique_pitch_desired = unique_pitch_desired(find(abs(unique_pitch_desired) > 3)); % PARAMETERS [degrees]
msg = sprintf('Unique MAX_BOUY values: %s\n',succinct_elts(unique_max_buoy));
updateResults(handles, msg);  

% unique sorts
trust_abc_fit = 1; % assume the best
if (abs(unique_max_buoy(1) - unique_max_buoy(end)) < 30) % PARAMETER [cc]
	msg = sprintf('Consider a dive with increased or decreased thrust ($MAX_BUOY) for better abc fit.\n')
    updateResults(handles, msg);  
	trust_abc_fit = 0;
end

msg = sprintf('Unique PITCH desired values: %s\n',succinct_elts(unique_pitch_desired));
updateResults(handles, msg);  

if (abs(unique_pitch_desired(1) - unique_pitch_desired(end)) < 5) % PARAMETER [degrees]
	msg = sprintf('Consider doing a steeper dive for better abc fit.\n')
	trust_abc_fit = 0;
end

% CCE says we only need one representative dive for each combination to get a good abc fit
regress_dives = [];
for max_buoy = unique_max_buoy'
	dives_i = find(data(:,2) == max_buoy);
	for pitch_desired = unique_pitch_desired'
		combined_dives_i = dives_i(find(data(dives_i,3) == pitch_desired));
		if (~isempty(combined_dives_i))
			% add most recent dive
			regress_dives = [regress_dives ; data(combined_dives_i(end),1)];
		end
	end
end

regress_dives = unique(regress_dives);
% finally see what the bulk of the dives did and add a few of those to get a better vbdbias
% this weights the regression for vbdbias and the hull parameters
n_max_buoy = histc(data(:,2),unique_max_buoy);
most_max_buoy = unique_max_buoy(find(n_max_buoy == max(n_max_buoy),1,'first'));
n_pitch_desired = histc(data(:,3),unique_pitch_desired);
most_pitch_desired = unique_pitch_desired(find(n_pitch_desired == max(n_pitch_desired),1,'first'));
most_dives = data(find(data(:,2) == most_max_buoy & data(:,3) == most_pitch_desired),1);
most_dives = setdiff(most_dives, regress_dives); % ensure they are different from our unique dives
% Find a set of typical dives to balance the number of extremal dives we've collected
% and get at least weight number of them...
regress_dives = [regress_dives; most_dives(end - min(max(length(regress_dives),most_dives_weight),length(most_dives)-1):end)];
regress_dives = unique(regress_dives)';

msg = sprintf('Run VBD Regression using dives: %s\n',succinct_elts(regress_dives))
updateResults(handles, msg);  

if (length(regress_dives) < 4)
	msg = sprintf('NOTE: There are insufficiently different dives to provide a good abc fit.\n');
    updateResults(handles, msg);  
end


if doPlot
    yL = get(gca,'YLim');
    for i=1:length(regress_dives)
        line([regress_dives(i) regress_dives(i)], yL, 'Color', 'r');
    end
end

% Restore original path
cd(existingPath);


