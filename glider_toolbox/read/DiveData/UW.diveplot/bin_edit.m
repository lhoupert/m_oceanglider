% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%
%   bin_edit.m
%   removes outliers in S, makes scatter plots

% We could be called directly or from bin_edit_dir
if (~exist('bin_fit_runs','var'))
  clear
  % They could have renamed the bin file...
  bin_file = input('Which bin file? ','s');
  if (isempty(bin_file))
    bin_file = 'Bin.mat';
  end
  load(bin_file);
end
if (~exist('Bin_vars','var'))
  fprintf(1,'Old-style bin file.  Please rebuild!\n');
  return;
end
var_cnt = 1;
for var = Bin_vars
  eval(sprintf('%s_index = %d;',var{1},var_cnt));
  var_cnt = var_cnt + 1;
end
% get dive restrictions?
% dive, climb, profile?

dive_climb = 0; % CONTROL  0 - both, 1 - dives 2 - climbs
ts_only = 0; % CONTROL 0 - all possible graphs 1 - TS diagram only
[m,n] = size(Bin);

sg_calib_constants
if(m == 0)
    fprintf(1,'No data in Bin!\n');
    return;
end
runs = sort(unique(Bin(:,dive_num_index))); % all available in Bin
if (exist('bin_fit_runs','var'))
  runs = bin_fit_runs;
  clear bin_fit_runs;
else
  runs = ask_which_runs(runs);
end
runs = runs';
% Bin_i are the indices of entries for the selected dives
All_i = [];
Dive_i = [];
Climb_i = [];
for run = runs
  All_i = [All_i; find(Bin(:,dive_num_index) == run)];
  Dive_i = [Dive_i; find(Bin(:,dive_num_index) == run & Bin(:,pitch_sign_index) == -1)];
  Climb_i = [Climb_i; find(Bin(:,dive_num_index) == run & Bin(:,pitch_sign_index) ==  1)];
end
switch(dive_climb)
 case 0
  tag = 'Profiles';
  Bin_i = All_i;
 case 1 % dives
  tag = 'Dives';
  Bin_i = Dive_i;
 case 2 % climbs
  tag = 'Climbs';
  Bin_i = Climb_i;
end
xruns = succinct_elts(runs);
fprintf(1,'%s: %s\n',tag,xruns); % before truncation...
if (length(xruns) > 20)
  xruns = strcat(xruns(1:20),'...');
end

titlestring = sprintf('SG%s %s %s: %s', id_str, underscore(mission_title),tag,xruns);
% when was this plot made
timelabel = datestr(clock);

% prepare structure for datacursormode function
Bin_data.Bin = Bin;
Bin_data.Bin_i = Bin_i;
user_data.Bin_data = Bin_data; % initialize this common bit

% unique(Bin(Bin_i(find(Bin(Bin_i,salin_index) < 32.2)),dive_num_index))
[depmin,depmax,depdel] = plot_bounds(Bin(Bin_i,depth_index),0,6000);
[tmin,tmax,tdel] = plot_bounds(Bin(Bin_i,temp_index),-3,40,[12,2; 6,1; 0,0.5],0.1);
[smin,smax,sdel] = plot_bounds(Bin(Bin_i,salin_index),19,40,[12,2; 6,1; 0,0.5],0.1);
[sigmin,sigmax,sigdel] = plot_bounds(Bin(Bin_i,density_index),19,40);
% locs_i = salin_error_i(find(Bin(salin_error_i,salin_index) > 33.5))
% unique(Bin(locs_i,dive_num_index))

%   T vs. S
fig = figure;
plot(Bin(Bin_i,salin_index), Bin(Bin_i,temp_index), '.b', 'MarkerSize',1)
grid on
hold on
if (exist('salin_error_index','var'))
  % plot where there is large numeric error
  salin_error_i = Bin_i(find(Bin(Bin_i,salin_error_index)));
  plot(Bin(salin_error_i,salin_index), Bin(salin_error_i,temp_index), '.m', 'MarkerSize',1)
end

if (exist('suspect_ti_index','var'))
  suspect_ti_i= Bin_i(find(Bin(Bin_i,suspect_ti_index)));
  plot(Bin(suspect_ti_i,salin_index), Bin(suspect_ti_i,temp_index), '.r', 'MarkerSize',1)
end


axis([smin smax tmin tmax]);
title(titlestring)
xlabel('Salinity [PSU]')
ylabel('Temperature [\circC]')
v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);
user_data.X = {salin_index,'Salinity',0.01};
user_data.Y = {temp_index,'Temp',0.01};
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@bin_edit_dc);
if (ts_only)
  return; % bail early after just TS diagram
end


%   Temperature vs depth
fig = figure;
plot(Bin(Bin_i,temp_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
set(gca, 'YDir', 'reverse')
grid on
hold on
axis([tmin tmax depmin depmax]);
title(titlestring)
ylabel('Depth [m]')
xlabel('Temperature [\circC]')

v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel); 
user_data.X = {temp_index,'Temp',0.01};
user_data.Y = {depth_index,'Depth',0.1};
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@bin_edit_dc);

%   Salinity vs depth
fig = figure;
plot(Bin(Bin_i,salin_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
set(gca, 'YDir', 'reverse')
grid on
hold on
if (exist('salin_error_index','var'))
  % plot where there is large numeric error
  salin_error_i = Bin_i(find(Bin(Bin_i,salin_error_index)));
  plot(Bin(salin_error_i,salin_index), Bin(salin_error_i,depth_index), '.m', 'MarkerSize',1)
end

if (exist('suspect_ti_index','var'))
  suspect_ti_i= Bin_i(find(Bin(Bin_i,suspect_ti_index)));
  plot(Bin(suspect_ti_i,salin_index), Bin(suspect_ti_i,depth_index), '.r', 'MarkerSize',1)
end
axis([smin smax depmin depmax]);
title(titlestring)
ylabel('Depth [m]')
xlabel('Salinity [PSU]')
v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel); 
user_data.X = {salin_index,'Salinity',0.01};
user_data.Y = {depth_index,'Depth',0.1};
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@bin_edit_dc);

%   Potential density vs depth
fig = figure;
plot(Bin(Bin_i,density_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
set(gca, 'YDir', 'reverse')
grid on
hold on
axis([sigmin sigmax depmin depmax]);
title(titlestring)
ylabel('Depth [m]')
xlabel('\sigma_t [kg/m^3]')

v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel); 

if (sbe43_present)
  [oxymin,oxymax,oxydel] = plot_bounds(Bin(Bin_i,oxygen_index),0,400);
  %   Oxygen vs depth
  figure
  plot(Bin(Bin_i,oxygen_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
  grid on
  hold on
  axis([oxymin oxymax depmin depmax]);
  title(titlestring)
  ylabel('Depth [m]')
  xlabel('Dissolved Oxygen SBE-43 [\mumol/kg]')
  
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel); 
end

if (optode_present)
  [oxymin,oxymax,oxydel] = plot_bounds(Bin(Bin_i,optode_oxygen_index),0,400);
  %   optode Oxygen vs depth
  fig = figure;
  plot(Bin(Bin_i,optode_oxygen_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
  grid on
  hold on
  axis([oxymin oxymax depmin depmax]);
  title(titlestring)
  ylabel('Depth [m]')
  xlabel('Dissolved Oxygen Aanderaa Optode [\mumol/kg]')
  
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel); 
  
  %   optode Oxygen dphase vs depth
  figure
  plot(Bin(Bin_i,dphase_oxygen_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
  grid on
  hold on
  axis([oxymin oxymax depmin depmax]);
  title(titlestring)
  ylabel('Depth [m]')
  xlabel('Dissolved Oxygen Aanderaa Optode Dphase [\mumol/kg]')
  
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel); 
end

if (wetlabs_present)
  % prep for wetlabs
  [redmin,redmax,reddel] = plot_bounds(Bin(Bin_i,red_scttr_index),50,550);
  [bluemin,bluemax,bluedel] = plot_bounds(Bin(Bin_i,blue_scttr_index),50,300);
  [fluormin,fluormax,fluordel] = plot_bounds(Bin(Bin_i,fluor_index),0,500);
  
  %   Red Scatter vs depth
  fig = figure;
  plot(Bin(Bin_i,red_scttr_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
  grid on
  hold on
  axis([redmin redmax depmin depmax]);
  title(titlestring)
  ylabel('Depth [m]')
  xlabel('Red Backscatter [counts]')
  
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel); 
  
  %   Blue Scatter vs depth
  fig = figure;
  plot(Bin(Bin_i,blue_scttr_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
  grid on
  hold on
  axis([bluemin bluemax depmin depmax]);
  title(titlestring)
  ylabel('Depth [m]')
  xlabel('Blue Backscatter [counts]')
  
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel); 
  
  %   Fluorescence vs depth
  fig = figure;
  plot(Bin(Bin_i,fluor_index), Bin(Bin_i,depth_index), '.b', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
  grid on
  hold on
  axis([fluormin fluormax depmin depmax]);
  title(titlestring)
  ylabel('Depth [m]')
  xlabel('Fluorescence [counts]')
  
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel); 
end
