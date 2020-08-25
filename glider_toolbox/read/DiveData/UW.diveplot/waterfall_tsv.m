% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%
%  waterfall_tsv.m
%  Plot sg sequence of offset profiles
% called from a bin_fit*.m variant to accumulate profiles
% globals assumed:
% - the usual data arrays
% - id_str - the id of the glider
% - ir the index of the dive number, starting at 1
% - irun the dive number

if (~exist('dpp'))
  dpp = 15;
end
% to display variants, set the following outside
% you can display up to dnvar variants, where dwvar (display which var)
% indicates which variant we are currently displaying
if (~exist('dnvar'))
  dnvar = 1;
  dwvar = 0; % zero based counter up to dnvar
end
ir_offset = mod(ir-1,dpp);
% this divides the space between variants into dnvar pieces
ir_offset = dnvar*ir_offset + dwvar*(ir_offset/(dnvar + 1));
% separate dives and climbs or overlap?
if (~exist('separate_dive_climb'))
  separate_dive_climb = 0; % overlap
  separate_dive_climb = .5; % separate
end
if (~exist('variant_tag'))
  variant_tag = ''; % assumes you add your own spacing
end

% which charts to make
temp_graph = 0;
salinity_graph = 0;
density_graph = 0;
vv_graph = 0;
ts_graph = 1;

% limit the graphs to the following depth range (top 200 m)
% Sigh, you'd like these to autoscale to the dives we see
% but we need to set these as constants over the whole set of dives
% so the offsets are against constants...
dmin = 0;
dmax = 1000;

tmin = 25;
tdel = 2;
tmax = 30;
if ((temp_graph | ts_graph) & min(temp) < tmin)
  fprintf(1,'Waterfall: Min temp is %.1f!\n', min(temp));
end
if ((temp_graph | ts_graph) & max(temp) > tmax)
  fprintf(1,'Waterfall: Max temp is %.1f!\n', max(temp));
end

smin = 29;
sdel = 1.0;
smax = 36;
if ((salinity_graph | ts_graph) & min(salin) < smin)
  fprintf(1,'Waterfall: Min salin is %.1f!\n', min(salin));
end
if ((salinity_graph | ts_graph) & max(salin) > smax)
  fprintf(1,'Waterfall: Max salin is %.1f!\n', max(salin));
end

sgmin = 23;
sgdel = 1.0;
sgmax = 27;

wmin = -4;
wdel = 4;
wmax = 24;

if (~exist('titlestring'))
  titlestring = sprintf('Seaglider %s %s', id_str,pwd());
end
timelabel = datestr(clock,'mm/dd/yy'); % forget hour stuff
dive_number= irun;
% iter is the number of iterations required to converge for this dive (TS phase)
divelabel = sprintf('%d/%d%s',irun,iter,variant_tag);
% color = ['r'; 'g'; 'b'; 'm'; 'c'];
color = ['r'; 'g'; 'b'];
jcol = mod(ir, length(color)) + 1;

% figure out dive and climb for this dive
iwn = find( pitch_control < 0 );
iwp = find( pitch_control > 0 );

% find where TS correction happened 
iwn_TS = find( pitch_control < 0 & salin_TS ~= salin);
iwp_TS = find( pitch_control > 0 & salin_TS ~= salin);


% plot temperature profiles

temp_offset = tdel;
tmaxo = tmax + temp_offset*dpp;
tlab = 4;
T_offset = ir_offset*temp_offset;

if (temp_graph)
  figure(1);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( temp(iwn) + T_offset, ct_depth_m(iwn), 'Marker', 'v', ...
              'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( temp(iwp) + T_offset + temp_offset*separate_dive_climb, ...
              ct_depth_m(iwp), 'Marker', '^', 'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(ax1, 'XAxisLocation', 'top', 'XColor', 'k', 'XLim', [tmin,tmaxo], ...
           'XTick', [tmin:tdel:tmaxo],'YLim', [dmin,dmax]);
  x_label  =  strcat('Temperature [\circC] (label @ ', ...
                     num2str(tlab), '\circC) ') ;
  y_label  =  'Depth [m]';
  xlabel(x_label, 'FontSize', 10); 
  ylabel(y_label, 'FontSize', 10);
  
  text(tlab + T_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
  
end

% plot salinity profiles
salin_offset = sdel;
smaxo = smax + salin_offset*dpp;
slab = 34.5;
S_offset = ir_offset*salin_offset;

if (salinity_graph)
  figure(2);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( salin(iwn) + S_offset, ct_depth_m(iwn), ...
              'Marker', 'v', 'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( salin(iwp) + S_offset + salin_offset*separate_dive_climb, ...
              ct_depth_m(iwp), 'Marker', '^', 'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( salin_TS(iwn_TS) + S_offset, ct_depth_m(iwn_TS), ...
              'Marker', 'v', 'Color', 'k', ...
              'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( salin_TS(iwp_TS) + S_offset + salin_offset*separate_dive_climb, ...
              ct_depth_m(iwp_TS), 'Marker', '^', 'Color', 'k', ...
              'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(gca,'YDir','reverse', 'XLim', [smin,smaxo], 'YLim', [dmin,dmax], ...
          'XAxisLocation', 'top', 'XTick', [smin:sdel:smaxo] );
  x_label  =  strcat( 'Salinity (label @ ', ...
                      num2str(slab), ')'); 
  y_label  =  'Depth (m)';
  xlabel(x_label, 'FontSize', 10); 
  ylabel(y_label, 'FontSize', 10);
  
  text(slab + S_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
  
end
% plot density profiles

sg_offset = sgdel;
sgmax = sgmax + sg_offset*dpp;
sglab = 27.5;
Sg_offset = ir_offset*sg_offset;

if (density_graph)
  figure(3);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( sigma(iwn) + Sg_offset, ct_depth_m(iwn), ...
              'Marker', 'v', 'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( sigma(iwp) + Sg_offset + sg_offset*separate_dive_climb, ...
              ct_depth_m(iwp), 'Marker', '^', 'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( sigma_TS(iwn_TS) + Sg_offset, ...
              ct_depth_m(iwn_TS), 'Marker', 'v', 'Color', 'k', ...
              'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( sigma_TS(iwp_TS) + Sg_offset + sg_offset*separate_dive_climb, ...
              ct_depth_m(iwp_TS), 'Marker', '^', 'Color', 'k', ...
              'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(gca,'YDir','reverse', 'XLim', [sgmin,sgmaxo], 'YLim', [dmin,dmax], ...
          'XAxisLocation', 'top', 'XTick', [sgmin:sgdel:sgmaxo] );
  x_label  =  strcat( '\sigma_t [kg/m^3] (label @ ', ...
                      num2str(sglab), ' kg/m^3) '); 
  y_label  =  'Depth (m)';
  xlabel(x_label, 'FontSize', 10); 
  ylabel(y_label, 'FontSize', 10);
  
  text(sglab + Sg_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
end
% Plot vertical velocity estimates

wmin = -4;
wdel = 4;
wmax = 24;
w_offset = wdel;
wmax = wmax + w_offset*dpp;
wlab = 0;
W_offset = ir_offset*w_offset;
if (vv_graph)
  figure(4);
  set(gca,'YDir','reverse');
  ax4 = gca;
  set(ax4, 'XLim', [wmin,wmaxo], 'XTick', [wmin:wdel:wmaxo], ...
           'YLim', [dmin,dmax], 'XAxisLocation', 'top' );
  hold on;
  grid on;
  w_h2O = w - w_stdy;
  plot( w_h2O(iwn) + W_offset, ct_depth_m(iwn), strcat('-', ...
                                                    color(jcol)), 'MarkerSize', 1);
  plot( w_h2O(iwp) + W_offset + w_offset*separate_dive_climb, ct_depth_m(iwp), ...
        strcat('-', color(jcol)), 'MarkerSize', 1);
  x_label  = strcat( 'Vertical Velocity [cm/s] (label @', ...
                     num2str(wlab), ' cm/s)' );
  xlabel(x_label, 'FontSize', 10); 
  y_label = 'Depth [m]';
  ylabel(y_label, 'FontSize', 10);
  
  text(wlab + W_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
  
end
if (0)
  % plot oxygen profile
  
  omin = 0;
  odel = 50;
  oxygen_offset = odel;
  omax = 1250 + oxygen_offset*dpp;
  olab = 0;
  O_offset = ir_offset*oxygen_offset;

  figure(5);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( oxygen(iwn) + O_offset, ct_depth_m(iwn), 'Marker', 'v', ...
              'Color', color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( oxygen(iwp) + O_offset +oxygen_offset*separate_dive_climb, ...
              ct_depth_m(iwp), 'Marker', '^', ...
              'Color', color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(ax1, 'XAxisLocation', 'top', 'XColor', 'k', 'XLim', ...
           [omin,omax], 'XTick', [omin:odel:omax],'YLim', [dmin,dmax]);
  xlabel( strcat( 'Dissolved Oxygen [\mum/kg] (label @ ', ...
		num2str(olab), ' \mum/kg)' ) );
  ylabel( 'Depth [m]');
  
  text(olab + O_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
end

if(0)
  % plot fluorometer profile
  
  fmin = 0;
  fmax = 2000;
  fdel = 100;
  fluor_offset = fdel*0.5;
  flab = 100;
  F_offset = ir_offset*fluor_offset;
  
  figure(6);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( fluor(iwn) + F_offset, ct_depth_m(iwn), 'Marker', 'v', ...
              'Color', color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( fluor(iwp) + F_offset + fluor_offset/2 , ct_depth_m(iwp), ...
              'Marker', '^', 'Color', ...
              color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(ax1, 'XAxisLocation', 'top', 'XColor', 'k', 'XLim', ...
           [fmin, fmax], 'XTick', [fmin:fdel:fmax],'YLim', [dmin,dmax]);
  xlabel( strcat( 'Fluorometer Output [counts] (label @', ...
                  num2str(flab), ' counts)' ) );	
  ylabel( 'Depth [m]');
  
  text(flab + F_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
  
  % plot red backscatter profile
  
  red_scttrmin = 0;
  red_scttrmax = 2000;
  red_scttrdel = 100;
  red_scttr_offset = red_scttrdel*0.5;
  red_scttrlab = 100;
  rs_offset = ir_offset*red_scttr_offset;
  
  figure(7);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( red_scttr(iwn) + rs_offset, ct_depth_m(iwn), 'Marker', 'v', ...
              'Color', color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( red_scttr(iwp) + rs_offset + red_scttr_offset*separate_dive_climb, ...
              ct_depth_m(iwp), 'Marker', '^', 'Color', ...
              color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(ax1, 'XAxisLocation', 'top', 'XColor', 'k', 'XLim', ...
           [red_scttrmin, red_scttrmax], 'XTick', ...
           [red_scttrmin:red_scttrdel:red_scttrmax],'YLim', [dmin,dmax]);
  xlabel( strcat( 'Red b_b Output [counts] (label @ ', ...
                  num2str(red_scttrlab), ' counts)' ) );	
  ylabel( 'Depth [m]');
  
  text(red_scttrlab + rs_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
  
  % plot blue backscatter profile
  
  blue_scttrmin = 0;
  blue_scttrmax = 2000;
  blue_scttrdel = 100;
  blue_scttr_offset = blue_scttrdel*0.5;
  blue_scttrlab = 100;
  bs_offset = ir_offset*blue_scttr_offset;
  
  figure(8);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( blue_scttr(iwn) + bs_offset, ct_depth_m(iwn), 'Marker', 'v', ...
              'Color', color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( blue_scttr(iwp) + bs_offset + blue_scttr_offset*separate_dive_climb, ...
              ct_depth_m(iwp), 'Marker', '^', 'Color', ...
              color(jcol), 'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(ax1, 'XAxisLocation', 'top', 'XColor', 'k', 'XLim', ...
           [blue_scttrmin, blue_scttrmax], 'XTick', ...
           [blue_scttrmin:blue_scttrdel:blue_scttrmax],'YLim', [dmin,dmax]);
  xlabel( strcat( 'Blue b_b Output [counts] (label @ ', ...
                  num2str(blue_scttrlab), ' counts)' ) );	
  ylabel( 'Depth [m]');
  
  text(blue_scttrlab + bs_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
  
  % plot buoyancy profiles
  
  buoymin = 0;
  buoymax = 600;
  buoydel = 25;
  buoy_offset = 2*buoydel;
  buoylab = 0;
  Buoy_offset = ir_offset*buoy_offset;
  
  figure(9);
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  hl1 = line( [Buoy_offset Buoy_offset], [dmin dmax], ...
              'Marker', 'none', 'Color', color(jcol),  'LineStyle', '--');
  %	hl1 = line( spd_stdy(iwp) + Buoy_offset, ct_depth_m(iwp), ...
  %		'Marker', 'o', 'Color', color(jcol), ...
  %		'MarkerSize', 1, 'LineStyle', '-');
  hl1 = line( buoy(iwp) + Buoy_offset, ct_depth_m(iwp), ...
              'Marker', '^', 'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', '-');
  ax1 = gca;
  set(gca,'YDir','reverse', 'XLim', [buoymin,buoymax], 'YLim', [dmin,dmax], ...
          'XAxisLocation', 'top', 'XTick', [buoymin:buoydel:buoymax] );
  x_label  =  strcat( 'Climb Buoyancy [g] & Speed [cm/s] (label @ ', ...
                      num2str(buoylab), ' kg/m^3) '); 
  y_label  =  'Depth (m)';
  xlabel(x_label, 'FontSize', 10); 
  ylabel(y_label, 'FontSize', 10);
  
  text(buoylab + Buoy_offset, dmax, divelabel, 'FontSize', 8, ...
       'Rotation', 90, 'Color', color(jcol) );
  v = axis;
  xtxt = v(1);
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(4) + 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(4) + 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
end
% plot temperature vs. salinity

% these use the temp and salinity limits established above
tempTS_offset = 0.5;
S_offset = ir_offset*salin_offset;
T_offset = ir_offset*tempTS_offset;
tmaxo = tmax + tempTS_offset*dpp;
if (ts_graph)
  TS_color = 'k'; % black
  % TS_color = 'r'; % red
  % TS_line_style = 'none';
  TS_line_style = '-';
  TS_size = 2;
  figure(10);
  hold on;
  grid on;
  hl1 = line( salin(iwn) + S_offset, temp(iwn) + T_offset, ...
              'Marker', 'v', 'Color', color(jcol), ...
              'MarkerSize', 1, 'LineStyle', TS_line_style);
  hl1 = line( salin(iwp) + S_offset + salin_offset*separate_dive_climb, ...
              temp(iwp) + T_offset + tempTS_offset*separate_dive_climb, ...
              'Marker', '^', 'Color', color(jcol), ...'
              'MarkerSize', 1, 'LineStyle', TS_line_style);
  % BUG if there are several disconnected segments interpolated
  % the black lines look connected but they aren't...so we show as points
  % instead of lines.  Used to use line with option 'LineStyle', '-'
  hl1 = plot( salin_TS(iwn_TS) + S_offset, ...
              temp(iwn_TS) + T_offset, ...
              strcat('v','',TS_color), ...
              'MarkerSize', TS_size);
  hl1 = plot( salin_TS(iwp_TS) + S_offset + salin_offset*separate_dive_climb, ...
              temp(iwp_TS) + T_offset + tempTS_offset*separate_dive_climb, ...
              strcat('^','',TS_color), ...
              'MarkerSize', TS_size);
  ax1 = gca;
  set(gca,'XLim', [smin,smaxo], 'YLim', [tmin,tmaxo], 'YTick', ...
          [tmin:tdel:tmaxo], 'XAxisLocation', 'top', 'XTick', [smin:sdel:smaxo] );
  x_label  =  strcat( 'Salinity (label @ ', num2str(slab), ', Offsets=', ...
                      num2str(salin_offset), '&', num2str(tempTS_offset), ...
                      '\circC/dive)'); 
  y_label  =  'Temperature [\circC]';
  xlabel(x_label, 'FontSize', 10); 
  ylabel(y_label, 'FontSize', 10);
  
  if(isempty(iwn))
      ti = iwp(end);
  else
      it = iwn(end);
  end
  text(slab + S_offset, temp(ti)-tdel+T_offset, divelabel, ...
       'FontSize', 10, 'Rotation', 90, 'Color', color(jcol) );
  v = axis;
  xtxt = v(1);
  ytxt = v(3) - 0.1*(v(4)-v(3));
  text(xtxt , ytxt, timelabel, 'FontSize', 8);
  
  xtxt = v(1);
  ytxt = v(3) - 0.05*(v(4)-v(3));
  text(xtxt , ytxt, titlestring,'FontSize',8);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  ytxt = v(3) - 0.1*(v(4)-v(3));
  text(xtxt , ytxt, comment,'FontSize',8);
  
  if(ir_offset == 0)
    sgrid = [smin:(smaxo-smin)/50:smaxo];
    tgrid = [tmin:(tmaxo-tmin)/50:tmaxo];
    [Sg, Tg] = meshgrid(sgrid, tgrid);
    Pg = zeros( length(tgrid), length(sgrid) );
    sigma_grid = sw_dens( Sg, Tg, Pg ) - 1000;
    sigma_levels = [20:0.5:40];
    [C, h] =contour( Sg, Tg, sigma_grid, sigma_levels, 'k' );
    clabel(C, h, sigma_levels, 'FontSize', 8);
  end
end  
if (ir_offset == dpp-1)
  wf_command = input('<cr> - close & continue; p - print: ','s');
  doprint = 0;
  if (length(wf_command) & wf_command == 'p')
    doprint = 1;
  end
  if (temp_graph)
    figure_print_close(1,doprint); 	% temperature
  end
  if(salinity_graph)
    figure_print_close(2,doprint); 	% salinity
  end
  if(density_graph)
    figure_print_close(3,doprint); 	% sigma
  end
  if (vv_graph)
    figure_print_close(4,doprint); 	% w
  end
  if(0)
    figure_print_close(5,doprint); 	% oxygen
    figure_print_close(6,doprint);	% F
    figure_print_close(7,doprint) 	% red_scttr
    figure_print_close(8,doprint); 	% blue_scttr
    figure_print_close(9,doprint); 	% buoyancy
  end
  if (ts_graph)
    figure_print_close(10,doprint); 	% T vs S
  end
end