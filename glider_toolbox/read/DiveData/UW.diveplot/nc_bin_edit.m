% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%
%   nc_bin_edit.m
ts_only = 0;
sg_calib_constants
% form default name from MDP:make_mission_profile?
files = dir('*profile.nc');
if (isempty(files))
	bin_file = input('Which bin file? ','s');
else
	bin_file = files(1).name
end
results = get_nc_data(bin_file);
globals = results.GLOBALS;
if (~isfield(globals,'binwidth'))
	fprintf(1,'%s does not appear to be the result of make_mission_profile()\n',bin_file);
	return;
end
% depending on globals.file_data_type we could have only dives, only climbs, or dive and climb data
legend_text = {};
if (strcmp(globals.file_data_type,'Up and Down profile'))
	file_type = 3;
	legend_text = {'Dive','Climb'};
elseif (strcmp(globals.file_data_type,'Down profile only'))
	file_type = 1;
	% legend_text = {'Dive'};
else % Climb profile only
	file_type = 2;
	% legend_text = {'Climb'};
end

dive_numbers = results.dive_number;
runs = ask_which_runs(unique(dive_numbers));
runs = runs';
xruns = succinct_elts(runs);
if (length(xruns) > 20)
  xruns = strcat(xruns(1:20),'...');
end
titlestring = sprintf('%s %s %dm bins %s',globals.platform_id, globals.project, globals.binwidth, xruns);
% when was this plot made
timelabel = datestr(clock);


Dive_i = [];
Climb_i = [];
for run = runs
	dive_i = find(dive_numbers == run,1,'first');
	switch (file_type)
	 case 1
	  Dive_i = [Dive_i; dive_i];
	 case 2
	  Climb_i = [Climb_i; dive_i];
	 case 3
	  Dive_i = [Dive_i; dive_i];
	  Climb_i = [Climb_i; dive_i+1];
	end
end
% prepare structure for datacursormode function
user_data.dive_numbers = dive_numbers;
user_data.Profiles_i = unique([Dive_i; Climb_i]); % initialize this common bit

ctd_depth = results.ctd_depth; % This is the averaged! depth of the bin, not the bin depth
depth = results.depth; % the depth bins
ctd_depth = repmat(depth,1,size(ctd_depth,2)); % expand to size of data arrays

temperature = results.temperature;
salinity = results.salinity;
sigma_t = results.sigma_t;
[depmin,depmax,depdel] = plot_bounds(ctd_depth,0,6000);
[tmin,tmax,tdel] = plot_bounds(temperature,-3,40,[12,2; 6,1; 0,0.5],0.1);
[smin,smax,sdel] = plot_bounds(salinity,19,40,[12,2; 6,1; 0,0.5],0.1);
[sigmin,sigmax,sigdel] = plot_bounds(sigma_t,19,40);


%   T vs. S
fig = figure;
grid on
hold on
plot(salinity(:,Dive_i), temperature(:,Dive_i),  '.b', 'MarkerSize',1)
plot(salinity(:,Climb_i),temperature(:,Climb_i), '.r', 'MarkerSize',1)
axis([smin smax tmin tmax]);
title(titlestring)
xlabel('Salinity [PSU]')
ylabel('Temperature [\circC]')
% Add density contours
sgrid = [smin:(smax-smin)/50:smax];
tgrid = [tmin:(tmax-tmin)/50:tmax];
[Sg, Tg] = meshgrid(sgrid, tgrid);
Pg = zeros( length(tgrid), length(sgrid) );
sigma_grid = sw_dens( Sg, Tg, Pg ) - 1000;
sigma_levels = [5:0.5:28];
[C, h] =contour( Sg, Tg, sigma_grid, sigma_levels, 'k' );
clabel(C, h, sigma_levels); % this takes some time to label...

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);
user_data.X = {salinity,'Salinity',0.01};
user_data.Y = {temperature,'Temp',0.01};
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@nc_bin_edit_dc);
if (ts_only)
  return; % bail early after just TS diagram
end


%   Temperature vs depth
fig = figure;
grid on
hold on
plot(temperature(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
plot(temperature(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
set(gca, 'YDir', 'reverse')
axis([tmin tmax depmin depmax]);
title(titlestring)
ylabel('Depth [m]')
xlabel('Temperature [\circC]')
v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel); 
user_data.X = {temperature,'Temp',0.01};
user_data.Y = {ctd_depth,'Depth',0.1};
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@nc_bin_edit_dc);

%   Salinity vs depth
fig = figure;
grid on
hold on
plot(salinity(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
plot(salinity(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
set(gca, 'YDir', 'reverse')
axis([smin smax depmin depmax]);
title(titlestring)
ylabel('Depth [m]')
xlabel('Salinity [PSU]')
v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel); 
user_data.X = {salinity,'Salinity',0.01};
user_data.Y = {ctd_depth,'Depth',0.1};
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@nc_bin_edit_dc);

%   Potential density vs depth
fig = figure;
grid on
hold on
plot(sigma_t(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
plot(sigma_t(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
set(gca, 'YDir', 'reverse')
axis([sigmin sigmax depmin depmax]);
title(titlestring)
ylabel('Depth [m]')
xlabel('\sigma_t [kg/m^3]')
v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel); 
user_data.X = {sigma_t,'Sigma_t',0.01};
user_data.Y = {ctd_depth,'Depth',0.1};
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@nc_bin_edit_dc);

if (isfield(results,'sbe43_dissolved_oxygen'))
	oxygen = results.sbe43_dissolved_oxygen;
	[oxymin,oxymax,oxydel] = plot_bounds(oxygen,0,400);
	%   Oxygen vs depth
	figure
	grid on
	hold on
	plot(oxygen(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
	plot(oxygen(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
	set(gca, 'YDir', 'reverse')
	axis([oxymin oxymax depmin depmax]);
	title(titlestring)
	ylabel('Depth [m]')
	xlabel('Dissolved Oxygen SBE-43 [\mumol/kg]')
	v = axis;
	ytxt = v(4) + 0.1*(v(4)-v(3));
	xtxt = v(1);
	text(xtxt , ytxt, timelabel); 
	user_data.X = {sigma_t,'Oxygen',0.01};
	user_data.Y = {ctd_depth,'Depth',0.1};
	set(gcf,'UserData',user_data);
	dcm_obj = datacursormode(fig);
	set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@nc_bin_edit_dc);
end

optode_present = 0;
if (isfield(results,'aanderaa3830_dissolved_oxygen'))
	oxygen = results.aanderaa3830_dissolved_oxygen;
	optode_present = 1;
elseif (isfield(results,'aanderaa4330_dissolved_oxygen'))
	oxygen = results.aanderaa4330_dissolved_oxygen;
	optode_present = 1;
end
if (optode_present)
	[oxymin,oxymax,oxydel] = plot_bounds(oxygen,0,400);
	%   optode Oxygen vs depth
	fig = figure;
	grid on
	hold on
	plot(oxygen(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
	plot(oxygen(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
	set(gca, 'YDir', 'reverse')
	axis([oxymin oxymax depmin depmax]);
	title(titlestring)
	ylabel('Depth [m]')
	xlabel('Dissolved Oxygen Aanderaa Optode [\mumol/kg]')
	v = axis;
	ytxt = v(4) + 0.1*(v(4)-v(3));
	xtxt = v(1);
	text(xtxt , ytxt, timelabel); 
	user_data.X = {sigma_t,'Oxygen',0.01};
	user_data.Y = {ctd_depth,'Depth',0.1};
	set(gcf,'UserData',user_data);
	dcm_obj = datacursormode(fig);
	set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@nc_bin_edit_dc);
end

if (isfield(results,'eng_wlbb2f_redCount') && ...
	isfield(results,'eng_wlbb2f_blueCount') && ...
	isfield(results,'eng_wlbb2f_fluorCount'))

	red_scttr = results.eng_wlbb2f_redCount;
	blue_scttr = results.eng_wlbb2f_blueCount;
	fluor = results.eng_wlbb2f_fluorCount;
	
	% prep for wetlabs
	[redmin,redmax,reddel] = plot_bounds(red_scttr,50,550);
	[bluemin,bluemax,bluedel] = plot_bounds(blue_scttr,50,300);
	[fluormin,fluormax,fluordel] = plot_bounds(fluor,0,500);
	
	%   Red Scatter vs depth
  fig = figure;
  grid on
  hold on
  plot(red_scttr(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
  plot(red_scttr(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
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
  grid on
  hold on
  plot(blue_scttr(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
  plot(blue_scttr(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
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
  grid on
  hold on
  plot(fluor(:,Dive_i), ctd_depth(:,Dive_i), '.b', 'MarkerSize',1)
  plot(fluor(:,Climb_i),ctd_depth(:,Climb_i),'.r', 'MarkerSize',1)
  set(gca, 'YDir', 'reverse')
  axis([fluormin fluormax depmin depmax]);
  title(titlestring)
  ylabel('Depth [m]')
  xlabel('Fluorescence [counts]')
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel); 
end
