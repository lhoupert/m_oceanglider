% 
% Copyright (c) 2006-2012, 2014 by University of Washington.  All rights reserved.
%
% This file contains proprietary information and remains the 
% unpublished property of the University of Washington. Use, disclosure,
% or reproduction is prohibited.
% this version has datacursormode so you can find out points, etc. for directives
show_original_data = 1;
show_depth_markers = 0; % show depths on plot
% skip_density_contours = 1;
show_legend = 0;
titlestring = sprintf('%s %s', underscore(mission_title), timestr); % UTC
timelabel = datestr(clock);

%
% get derived/corrected quantities
%
% temp - full measured temperature
% salin_raw -  initial salinity guess
% salin - corrected salinity
compute_dive_climb_indices
if (~exist('tmin','var'))
  [tmin,tmax,tdel] = plot_bounds(temp,-3,40,[12,2; 6,1; 0,0.5],0.1);
  [smin,smax,sdel] = plot_bounds([salin_raw(valid_i);salin],19,40,[12,2; 6,1; 0,0.5],0.1);
end
not_shown_i = find(temp > tmax | temp < tmin | salin > smax | salin < smin);
if (length(not_shown_i))
  fprintf(1,'WARNING: %d points not shown on TS diagram!\n',length(not_shown_i));
end

% Plot T vs. S
fig = figure('Name',sprintf('%s T vs. S',base_file));
set(gca, 'XLim', [smin,smax], 'XTick', [smin:sdel:smax], ...
         'YLim', [tmin,tmax], 'YTick', [tmin:tdel:tmax]);
hold on;

if (~exist('skip_density_contours','var') | ~skip_density_contours)
  v = axis;
  smin = v(1);
  smax = v(2);
  tmin = v(3);
  tmax = v(4);
  s_del = (smax-smin)/50;
  t_del = (tmax -tmin)/50;
  sgrid = [smin:s_del:smax];
  tgrid = [tmin:t_del:tmax];
  [Sg, Tg] = meshgrid(sgrid, tgrid);
  Pg = zeros( length(tgrid), length(sgrid) );
  sigma_grid = sw_dens( Sg, Tg, Pg ) - 1000;
  sigma_levels = [floor(min(min(sigma_grid))):0.5:ceil(max(max(sigma_grid)))];
  [Cg, h] =contour( Sg, Tg, sigma_grid, sigma_levels, 'k' ); % don't use C; clobbers oyxgen calibration constant!!
  clabel(Cg, h, sigma_levels);
else
  grid on;
end


ts_marker_size = 2;
ltxt = {}; % only add items to ltxt if there is a corresponding set of data...
if (show_original_data)
  if (dive_i) plot( salin_raw(dive_i), temp(dive_i), 'vb', 'MarkerSize', ts_marker_size); ltxt(end+1) = {'dive original'}; end
  if (climb_i) plot( salin_raw(climb_i), temp(climb_i), '^r', 'MarkerSize', ts_marker_size); ltxt(end+1) = {'climb original'}; end
end

bad_i = intersect(uncorrectable_i,dive_i);
if (bad_i) plot( salin_raw(bad_i), temp(bad_i), 'vk','MarkerSize', 4); ltxt(end+1) = {'dive uncorrectable'}; end
bad_i = intersect(uncorrectable_i,climb_i);
if (bad_i) plot( salin_raw(bad_i), temp(bad_i), '^k','MarkerSize', 4); ltxt(end+1) = {'climb uncorrectable'}; end

% if (dive_i_TS) plot( salin(dive_i_TS), temp(dive_i_TS), 'vg', 'MarkerSize', ts_marker_size); ltxt(end+1) = {'dive interpolated'}; end
% if (climb_i_TS) plot( salin(climb_i_TS), temp(climb_i_TS), '^g', 'MarkerSize', ts_marker_size); ltxt(end+1) = {'climb interpolated'}; end

% Connect the dots over the corrected data (gaps were uncorrectable)
i_TS_segs = [find(diff(dive_i_corrected) > 1)'; length(dive_i_corrected)];
if (i_TS_segs)
  start_i = 1;
  for last_i = i_TS_segs'
    s_start_i = dive_i_corrected(start_i);
    s_end_i = dive_i_corrected(last_i);
    plot(salin(s_start_i:s_end_i),temp(s_start_i:s_end_i),'-c');
    start_i = last_i+1;
  end
end
i_TS_segs = [find(diff(climb_i_corrected) > 1)'; length(climb_i_corrected)];
if (i_TS_segs)
  start_i = 1;
  for last_i = i_TS_segs'
    s_start_i = climb_i_corrected(start_i);
    s_end_i = climb_i_corrected(last_i);
    plot(salin(s_start_i:s_end_i),temp(s_start_i:s_end_i),'-m');
    start_i = last_i+1;
  end
end

% Now show places where we interpolated
i_TS_segs = [find(diff(dive_i_TS) > 1); length(dive_i_TS)];
if (i_TS_segs)
  start_i = 1;
  for last_i = i_TS_segs'
    s_start_i = dive_i_TS(start_i);
    s_end_i = dive_i_TS(last_i);
    plot(salin(s_start_i:s_end_i),temp(s_start_i:s_end_i),'-k');
    plot(salin(s_start_i),temp(s_start_i),'xk','MarkerSize', 4);
    plot(salin(s_end_i),temp(s_end_i),'xk','MarkerSize', 4);
    start_i = last_i+1;
  end
end
i_TS_segs = [find(diff(climb_i_TS) > 1); length(climb_i_TS)];
if (i_TS_segs)
  start_i = 1;
  for last_i = i_TS_segs'
    s_start_i = climb_i_TS(start_i);
    s_end_i = climb_i_TS(last_i);
    plot(salin(s_start_i:s_end_i),temp(s_start_i:s_end_i),'-k');
    plot(salin(s_start_i),temp(s_start_i),'xk','MarkerSize', 4);
    plot(salin(s_end_i),temp(s_end_i),'xk','MarkerSize', 4);
    start_i = last_i+1;
  end
end

if (show_depth_markers)
  depth_markers = [10:10:100 200:100:1000 2000:1000:6000];
  depth_markers = depth_markers(find(depth_markers <= max(sg_depth_m)));
  depth_m_i = [];
  for depth = depth_markers
    dive_pt_i = dive_i(find(sg_depth_m(dive_i) > depth,1,'first')-1);
    climb_pt_i = climb_i(find(sg_depth_m(climb_i) < depth,1,'first')-1);
    plot([salin(dive_pt_i); salin(climb_pt_i)],[temp(dive_pt_i); temp(climb_pt_i)],'-k');
    max_pt_i = climb_pt_i;
    if (salin(dive_pt_i) > salin(climb_pt_i))
      % max_pt_i = dive_pt_i;
      1;
    end
    text(salin(max_pt_i), temp(max_pt_i),sprintf('%d meters',depth));
  end
end

x_label  =  'Salinity (psu)';
xlabel(x_label, 'FontSize', 10); 
y_label = 'Temperature (\circC)';
ylabel(y_label, 'FontSize', 10);
title(strcat(base_file, ' : ', titlestring), 'FontSize', 12);
if show_legend legend(ltxt,0); end;

user_data = ([]); % initialize
% Users prefer to latch to the raw points
user_data.salin = salin_raw;
user_data.temp = temp_raw;
user_data.depth = sg_depth_m;
user_data.directives = directives;
user_data.tolerance = 0.05;
set(gcf,'UserData',user_data);
dcm_obj = datacursormode(fig);
set(dcm_obj,'Enable','on','SnapToDataVertex','off','Updatefcn',@ts_diagram_dc);

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);
if (exist('calibcomm','var')) % might be GPCTD, which has no calibcomm
  xtxt = v(1) + 0.6*(v(2)-v(1));
  text(xtxt, ytxt, calibcomm);
end