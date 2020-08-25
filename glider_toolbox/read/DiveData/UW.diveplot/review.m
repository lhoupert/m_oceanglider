% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% map over mat files and use waterfall to display dives of interest with some editing support
clear_sg_config_constants % clear the top workspace
sg_calib_constants % supplies mission_title, calibcomm
sg_config_constants
qc_declarations
% 0 - always show
% 1 - # total interp_salinity points
% 2 - # anomolies (no bubbles)
% 3 - # stalled points
% 4 - # suspect points
interest_function = 0; % CONTROL
num_interesting_points = 10; % PARAMETER
use_fixed_plot_bounds = 0; % CONTROL
skip_density_contours = 1; % CONTROL for ts_diagram
show_grid = 0; % CONTROL for waterfall
show_original_data = 1; % CONTROL for waterfall
dump_wf_figures = 0; % CONTROL

switch (interest_function)
 case 0
  ifn = 'Selected dives';
  num_interesting_points = 0;
 case 1
  ifn = 'Interpolated points';
  num_interesting_points = 20; % PARAMETER
 case 2
  ifn = 'Cond anomolies';
 case 3
  ifn = 'Stalled points';
 case 4
  ifn = 'Suspect thermal inertia points';
end


fprintf(1,'%s: Looking for dives with more than %d interesting points\n',ifn, num_interesting_points);
fprintf(1,'SG%s: ', id_str);
runs = ask_which_runs(available_profiles());
nruns = length(runs);
interesting_dives = [];
% if users gives CR waterfall = [] == 0
waterfall = input('Waterfall (1) or TS diagrams (0)?');

% support for waterfall ts plots
wffig = 0; % force a new figure...
timelabel = sprintf('Created : %s',datestr(clock));
titlestring = sprintf('SG%s: %s (> %d points) in %s',id_str,ifn,num_interesting_points, underscore(mission_title));
timestr = ''; % need to read log to get log_dn and convert...save in mat
dpp = 15; % PARAMETER how many dives per page?
color = ['r'; 'g'; 'b'];
% how much to offset each dive from one another
temp_offset = 0; % WAS 0.5; % PARAMETER
salin_offset = 1; % PARAMETER
% how close to make the dive and climb
% 0 overlaps dive and climb; what ts_diagram would show...
separate_dive_climb = 0; % PARAMETER
separate_dive_climb = salin_offset/4; % PARAMETER

% initialize for waterfall
min_t = -3;
max_t = 40;
min_s = 19;
max_s = 40;
tmin_all = max_t; tmax_all = min_t; smin_all = max_s; smax_all = min_s;

% Figures take a lot of memory
% For each figure report remaining percentage
% so user can throttle use
available_mem = avail_mem(1);
% loop on runs
ir = 1;
clear tmin tmax tdel smin smax sdel; % force recomputation of bounds
for ix = 1:nruns;
  irun = runs(ix);

  id_dive = sprintf('%s%04d',id_str,irun);
  base_file = sprintf('p%s',id_dive);
  [loginfo,eng,results] = get_dive_data(id_dive);
  if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
    fprintf(1,'Processing error or skipped profile for %d; skipping\n',irun);
    continue;
  end
  unpack_data;

  if (isempty(dive_i_corrected) || isempty(climb_i_corrected))
      continue; % missing data
  end
  switch (interest_function)
   case 0
    points_i = interpolated_i; % show interpolated points, if any
   case 1
    points_i = interpolated_i;
   case 2
    points_i = [directives.fixable_anomalies; directives.uncorrectable_anomalies];
   case 3
    % remove apogee, top of dive and climb to get for 'flight stalls'
    ignore_i = sort(unique([[1:directives.start_of_dive]'; directives.apogee; [directives.end_of_climb:mp]' ; directives.shallow_bubbles]));
    points_i = setdiff(directives.stalled,ignore_i);
   case 4
    points_i = directives.suspect_thermal_inertia;
  end
    
  if (length(points_i) >= num_interesting_points)
    interesting_dive = 1;
  else
    interesting_dive = 0;
  end

  if (interesting_dive)
    interesting_dives = [interesting_dives ; irun];
    divelabel = sprintf('%d/%d',irun,length(points_i));

    if (waterfall)
      % open-coded waterfall_tsv
      ir_offset = mod(ir-1,dpp);
      jcol = mod(ir, length(color)) + 1;
      tcolor = color(jcol);
      dcolor = tcolor; % original
      ccolor = tcolor; % original
      if (1)
        % match ts_diagram
        tcolor = 'b';
        dcolor = 'c';
        ccolor = 'm';
      end

      [tmin,tmax,tdel] = plot_bounds(temp,min_t,max_t);
      [smin,smax,sdel] = plot_bounds(salin,min_s,max_s);
      % add some room on the bottom for the label
      tmin = tmin - tdel; 
      smin = smin - sdel;
      tmaxo = tmax + dpp*temp_offset;
      smaxo = smax + dpp*salin_offset;
      tmin_all = min(tmin,tmin_all);
      tmax_all = max(tmaxo,tmax_all);
      smin_all = min(smin,smin_all);
      smax_all = max(smaxo,smax_all);
      if (wffig)
        figure(wffig); % ensure in focus...
        set(gca,'XLim', [smin_all,smax_all], 'XTick', [smin_all:sdel:smax_all],...
                'YLim', [tmin_all,tmax_all], 'YTick', [tmin_all:tdel:tmax_all]);
      else
        % create new figure
        % report memory remaining before figure

        wffig = figure;
        hold on;
        %DEAD grid on;
        if (show_grid)
          sgrid = [smin:(smaxo-smin)/50:smaxo];
          tgrid = [tmin:(tmaxo-tmin)/50:tmaxo];
          [Sg, Tg] = meshgrid(sgrid, tgrid);
          Pg = zeros( length(tgrid), length(sgrid) );
          sigma_grid = sw_dens( Sg, Tg, Pg ) - 1000;
          sigma_levels = [20:0.5:40];
          grid_color = 'k';
          grid_color = 'c';
          [C, h] =contour( Sg, Tg, sigma_grid, sigma_levels, grid_color );
          clabel(C, h, sigma_levels, 'FontSize', 8);
        end
      
        set(gca, 'XAxisLocation', 'top');
        set(gca,'XLim', [smin_all,smax_all], 'XTick', [smin_all:sdel:smax_all],...
                'YLim', [tmin_all,tmax_all], 'YTick', [tmin_all:tdel:tmax_all]);
        x_label= strcat( 'Salinity (offsets=', num2str(salin_offset), 'PSU &', num2str(temp_offset),'\circC/dive)'); 
        xlabel(x_label, 'FontSize', 10); 
        ylabel('Temperature (\circC)', 'FontSize', 10);

        v = axis;
        xtxt = v(1); % left margin
        ytxt = v(3) - 0.1*(v(4)-v(3));
        text(xtxt , ytxt, timelabel, 'FontSize', 8);
        
        xtxt = v(1); % left margin
        ytxt = v(3) - 0.05*(v(4)-v(3));
        text(xtxt , ytxt, titlestring,'FontSize',8);
        
        current_mem = avail_mem(0);
        fprintf(1,'%d bytes free; %.1f%% space remaining\n',current_mem, (current_mem/available_mem)*100.0);
      end % initialize new figure
      % for this dive...

      T_offset = ir_offset*temp_offset;
      S_offset = ir_offset*salin_offset;
      TS_color = 'k'; % black
      % TS_color = 'r'; % red
      % TS_line_style = 'none';
      TS_line_style = '-'; % original data are ticks, new data are lines
      TS_size = 2;
      if (show_original_data)
        plot( salin_raw(dive_i) + S_offset, ...
              temp(dive_i) + T_offset, ...
              'vb', 'MarkerSize',1);
        % 'Marker', 'v', 'Color', 'b','MarkerSize', 1);
        plot( salin_raw(climb_i) + S_offset + salin_offset*separate_dive_climb, ...
              temp(climb_i) + T_offset + temp_offset*separate_dive_climb, ...
              '^r', 'MarkerSize',1);
        % 'Marker', '^', 'Color', 'r','MarkerSize', 1);
      end
      plot( salin(dive_i_corrected) + S_offset, ...
            temp(dive_i_corrected) + T_offset, ...
            'Marker', 'v', 'Color', dcolor, ...
            'MarkerSize', 1, 'LineStyle', TS_line_style);
      plot( salin(climb_i_corrected) + S_offset + salin_offset*separate_dive_climb, ...
            temp(climb_i_corrected) + T_offset + temp_offset*separate_dive_climb, ...
            'Marker', '^', 'Color', ccolor, ...
            'MarkerSize', 1, 'LineStyle', TS_line_style);
      % BUG if there are several disconnected segments interpolated
      % the black lines look connected but they aren't...so we show as points
      % instead of lines.  Used to use line with option 'LineStyle', '-'
      plot( salin(dive_i_TS) + S_offset, ...
            temp(dive_i_TS) + T_offset, ...
            strcat('v','',TS_color), ...
            'MarkerSize', TS_size);
      plot( salin(climb_i_TS) + S_offset + salin_offset*separate_dive_climb, ...
            temp(climb_i_TS) + T_offset + temp_offset*separate_dive_climb, ...
            strcat('^','',TS_color), ...
            'MarkerSize', TS_size);
      
      text(salin(dive_i_corrected(end)) + S_offset - 0.25, ...
           temp(dive_i_corrected(end))  + T_offset - 0.0, ...
           divelabel,'FontSize', 10, 'Rotation', 90, 'Color', tcolor );

      ir = ir + 1;
      if (ir_offset == dpp-1)
		  if (dump_wf_figures)
			  print('-dpng',sprintf('review_%d',ir));
			  close;
		  end
        wffig = 0; % force a new figure
        tmin_all = max_t; tmax_all = min_t; smin_all = max_s; smax_all = min_s;
      end
    else
      ts_diagram;
      current_mem = java.lang.Runtime.getRuntime.freeMemory();
      fprintf(1,'%d bytes free; %.1f%% space remaining\n',current_mem, (current_mem/available_mem)*100.0);
      if (~use_fixed_plot_bounds)
          clear tmin tmax tdel smin smax sdel; % force recomputation of bounds
      end
    end
  end
end
if (wffig & dump_wf_figures)
	print('-dpng',sprintf('review_%d',ir));
	close;
	wffig = 0;
end
if (length(interesting_dives))
  interesting_dives = interesting_dives';
  fprintf(1, 'Interesting dives: %s\n', succinct_elts(interesting_dives));
else
  fprintf(1, 'No interesting dives found in runs %s\n',succinct_elts(runs));
end
clear tmin tmax tdel smin smax sdel; % force recomputation of (ts_diagram) bounds from later callers
clear skip_density_contours; % important!!