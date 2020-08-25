% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%
% common display routines for the diveplot_X code

%
% plot flight control variables vs. time
%
titlestring = underscore(sprintf('%s %s', mission_title, startstr));
basetitle = sprintf('%s: %s',base_file,titlestring);
mass_vol_comment = sprintf('volmax=%.1fcc, mass=%.3fkg',volmax,mass);
if (mass_comp)
	mass_vol_comment = sprintf('%s, mass_c_o_m_p=%.3fkg',mass_vol_comment,mass_comp);
end
biaslabel = sprintf('VBD bias = %gcc',vbdbias);

cc_per_m3 = 1e6;

figure('Name','Diveplot');

% Plot heading and right axis

timemin = 0;
timedel = 30;
timemax = (floor(sgtime_m(end)/timedel) + 1)*timedel;
maxd = max(vmdepth);
% use 1100 meters for SG so if we go a little over nomimal depth
% we don't jump to 2km scale
if maxd > 1100
  zscl = fix(maxd/100); % 1000 - 6000m
  % round up to the nearest 10
  zscl = (zscl - rem(zscl,10)) + 10;
elseif maxd > 500
  zscl = 10; % 500 - 1000m
elseif maxd > 300
  zscl = 5; % 300 - 500m
elseif maxd > 200
  zscl = 3; % 200 - 300m
elseif maxd > 100
  zscl = 2; % 100 - 200m
else
  zscl = 1; % 0 - 100m
end

dmin = 0;
dmax = zscl*100;
ddel = dmax/10;

ax1 = axes('YAxisLocation', 'right');
set(gca, 'YLim', [0,360], 'YTick', [0:20:360], 'YColor', 'r', ...
         'Xlim', [timemin timemax], 'Xtick', [timemin:timedel:timemax]);
line(sgtime_m, hdg, 'Color', 'r', 'LineStyle', 'none', ...
     'Marker', '+', 'MarkerSize', 3);
y_label = '(+) Heading (\circM)';
ylabel(y_label, 'FontSize', 10); 
hold on;
plot([sgtime_m(1) sgtime_m(length(sgtime_m))], [mag_heading mag_heading], '--r');

% Plot left axis and curves

ax2 = axes('Position', get(ax1,'Position'),...
           'YAxisLocation', 'left', 'Color', 'none');
set(gca, 'YLim', [-100,80], 'YTick', [-100:10:80], 'YColor', 'k');
hold on;
plot(sgtime_m, zpos/zscl, '-r'); 
x_label  =  'Time (minutes)';
xlabel(x_label, 'FontSize', 10); 
y_label  =  base_file;
ylabel(y_label, 'FontSize', 10); 
title(titlestring, 'FontSize', 12);
hold on;
plot(sgtime_m, w, '-b'); 
plot(sgtime_m, hspdg(ctd2sg_i), '-c'); 
plot(sgtime_m(sg_ctd_valid_i), w_stdy(ctd2sg_valid_i), '--b'); 
plot(sgtime_m(sg_ctd_valid_i), hspd_stdy(ctd2sg_valid_i), '--c'); 
plot(sgtime_m, pitch, '-g'); 
plot(sgtime_m, pitch_control*10, '--g');
plot(sgtime_m, roll, '-y'); 
plot(sgtime_m, roll_control, '--y');
plot(sgtime_m(sg_ctd_valid_i), glideangle_stdy(ctd2sg_valid_i)*radpd, '--m'); % convert to degrees
% Ensure that these two use the same scaling factor so they can be compared visually
plot(sgtime_m(sg_ctd_valid_i), buoy(ctd2sg_valid_i)/10, '-m'); 
plot(sgtime_m, vbd/10, '-k'); 
%
v = axis;
axis( [timemin timemax -100 80] );
set(gca, 'Ytick', [-100:10:80], 'Xtick', [timemin:timedel:timemax])
%
% Mark non-zero motor current 
%
ivb = find( vbddiff > 0.5 );
plot(sgtime_m(ivb), zpos(ivb)/zscl, '.k');
%
% Mark negative (port) and positive (starboard) roll changes
%
irn = find( fix(rolldiff) - abs( fix(rolldiff) ) );
irp = find( fix(rolldiff) + abs( fix(rolldiff) ) );
plot(sgtime_m(irn), zpos(irn)/zscl, '<k');
plot(sgtime_m(irp), zpos(irp)/zscl, '>k');
%
% Mark negative (downward) and positive (upward) pitch changes
%
pitchc = fix(10.*pitchdiff);
ipn = find( fix(pitchc) - abs( fix(pitchc) ) );
ipp = find( fix(pitchc) + abs( fix(pitchc) ) );
plot(sgtime_m(ipn), zpos(ipn)/zscl, 'vk');
plot(sgtime_m(ipp), zpos(ipp)/zscl, '^k');

% Draw legend and time stamp

lg = legend(strcat('Height (',num2str(zscl), 'm)'), 'W = dz/dt (cm/s)', ...
       'Horiz Speed(W,pitch) ', ...
       'Vert Speed(buoy,pitch)', ...
       'Horiz Speed(buoy,pitch)', ...
       'Pitch Up (deg)', ...
       'Pitch Control (mm)', ...
       'Vehicle Roll Stbd (deg)', ...
       'Roll Control (deg)', ...
       'Glide Angle (buoy,pitch)', ...
       'Buoyancy (10g)', ...
       'VBD (10 cc)', ...
       'Location','SouthEast');
set(lg,'color','none'); % make transparent

grid on;
plot([sgtime_m(1) sgtime_m(mp)], [intended_pitch intended_pitch], ':g');
plot([sgtime_m(1) sgtime_m(mp)], [-intended_pitch -intended_pitch], ':g');
plot([sgtime_m(1) sgtime_m(mp)], [w_d w_d], ':b');
plot([sgtime_m(1) sgtime_m(mp)], [-w_d -w_d], ':b');
v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

xtxt = v(1) + 0.6*(v(2) - v(1));
text(xtxt, ytxt, mass_vol_comment);

%
% plot track through water
%
rollmin = 5.;
irrn = find(ctd_roll < -rollmin);
irrp = find(ctd_roll > rollmin);

%	[kalman_brg, kalman_spd] = cart2pol(kalman_control_v, kalman_control_u); 
max_disp = hspdg_d*total_dive_time_s; % based on horizontal speed desired
control_disp = max_disp*speed_factor;
model_glide_slope = w_target/hspdg_d;
heading = mag_heading + magvar; % convert to true heading

[control_disp_y, control_disp_x] = pol2cart(heading*radpd, control_disp);
[control_disp_y_right, control_disp_x_right] = ...
    pol2cart((heading + head_errband)*radpd, control_disp);
[control_disp_y_left, control_disp_x_left] = ...
    pol2cart((heading - head_errband)*radpd, control_disp);
[target_disp_y, target_disp_x] = ...
    pol2cart((heading + course_bias)*radpd, control_disp);

[max_disp_y, max_disp_x] = pol2cart(heading*radpd, max_disp);
[max_disp_y_right, max_disp_x_right] = ...
    pol2cart((heading + head_errband)*radpd, max_disp);
[max_disp_y_left, max_disp_x_left] = ...
    pol2cart((heading - head_errband)*radpd, max_disp);

label_angle_offset = 10;
[label_disp_y_bearing, label_disp_x_bearing] = ...
    pol2cart((heading + label_angle_offset)*radpd, ...
             0.1*max_disp);
[label_disp_y_speed, label_disp_x_speed] = ...
    pol2cart((heading - label_angle_offset)*radpd, ...
             0.1*max_disp);

label_angle_offset = 2;
[label_target_disp_y_bearing, label_target_disp_x_bearing] = ...
    pol2cart((heading + course_bias + label_angle_offset)*radpd, ...
             0.75*max_disp);

figure('Name','Path through water');
plot(dx(irrp), dy(irrp), '.g');
hold on
grid on
plot(dx(irrn), dy(irrn), '.r');
plot(dx, dy, '-b')
plot(dx(ivb), dy(ivb), '.k');
% the X components below were [0., foo] but we got CAT args not
% consistent; changing to [0.; foo] worked
plot([0.; max_disp_x], [0.; max_disp_y], '--k')
plot([0.; max_disp_x_right], [0.; max_disp_y_right], '--k')
plot([0.; max_disp_x_left], [0.; max_disp_y_left], '--k')
plot([0.; control_disp_x], [0.; control_disp_y], '-k')
plot([0.; control_disp_x_right], [0.; control_disp_y_right], '-k')
plot([0.; control_disp_x_left], [0.; control_disp_y_left], '-k')
plot([0.; target_disp_x], [0.; target_disp_y], '--m')
axis equal;

label_angle = 90 - heading;
comment = strcat(num2str(heading), '\circT');
text(label_disp_x_bearing, label_disp_y_bearing, comment, ...
     'Rotation', label_angle);

label_target_angle = 90 - (heading + course_bias);
comment = strcat(num2str(heading + course_bias), '\circT');
text(label_target_disp_x_bearing, label_target_disp_y_bearing, ...
     comment, 'Rotation', label_target_angle);

comment = strcat(num2str(hspdg_d), 'm/s for ', ...
                 num2str(dive_duration), 'hr = ', num2str(max_disp/1000), 'km');
text(label_disp_x_speed, label_disp_y_speed, comment, ...
     'Rotation', label_angle);


xlabel('Eastward Displacement Through Water [m]');
ylabel('Northward Displacement Through Water [m]');
title(basetitle, 'FontSize', 12);

v = axis;
if (ave_brg > 270)
  xtxt = v(1) + 0.05*(v(2)-v(1));
  ytxt = v(3) + 0.15*(v(4)-v(3));
  location = 'NorthEast';
elseif ave_brg > 180
  xtxt = v(1) + 0.4*(v(2)-v(1));
  ytxt = v(3) + 0.15*(v(4)-v(3));
  location = 'NorthWest';
elseif ave_brg > 90
  xtxt = v(1) + 0.05*(v(2)-v(1));
  ytxt = v(3) + 0.15*(v(4)-v(3));
  location = 'NorthWest';
else
  xtxt = v(1) + 0.05*(v(2)-v(1));
  ytxt = v(3) + 0.85*(v(4)-v(3));
  location = 'SouthEast';
end
comment = {sprintf('Avg Spd thru water\t= %5.2f m/s @ %5.1f \\circT', ave_spd, ave_brg); ...
           sprintf('Distance thru water\t= %5.2f km', displacement/1000); ...
           sprintf('Max buoy (cmdfile)\t= %g cc', max_buoy); ...
           sprintf('Target w\t= %g m/s', w_target); ...
           sprintf('Model glide slope\t= %6.4f', model_glide_slope); ...
           sprintf('Net glide slope\t= %6.4f', net_glide_slope); ...
          };
text(xtxt, ytxt, comment);
lg = legend('Roll to Right', 'Roll to Left', 'Location', location);
set(lg,'color','none'); % make transparent

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

xtxt = v(1) + 0.6*(v(2) - v(1));
text(xtxt, ytxt, mass_vol_comment);

%
% Plot vertical velocity estimates
%
w_h2o = zeros(mp,1);
w_h2o(sg_ctd_valid_i) = w(sg_ctd_valid_i) - w_stdy(ctd2sg_valid_i);

figure('Name','Vertical Velocities')
% [wmin,wmax,wdel] = plot_bounds([w;w_stdy;w_h2o],-50,50,[100,10; 60,5; 40,4; 20,2; 0,1]);
% scale using only w
[wmin,wmax,wdel] = plot_bounds([w],-50,50,[100,10; 60,5; 40,4; 20,2; 0,1]);

set(gca,'YDir','reverse');
ax4 = gca;
set(ax4, 'XLim', [wmin,wmax], 'XTick', [wmin:wdel:wmax],'YLim', [dmin,dmax]);
hold on;
grid on;
plot( w, vmdepth, '.b');
plot( w_stdy(ctd2sg_valid_i), vmdepth(sg_ctd_valid_i), '-c');
plot( w_h2o(sg_dive_i), vmdepth(sg_dive_i), '-r');
plot( w_h2o(sg_climb_i), vmdepth(sg_climb_i), '-m');
plot([-w_d -w_d], [dmin dmax], '--k')
plot([w_d w_d], [dmin dmax], '--k')
plot([0 0], [dmin dmax], '-k')
plot( w, vmdepth, 'ob', 'MarkerSize', 2);
plot( w_stdy(ctd2sg_valid_i), vmdepth(sg_ctd_valid_i), 'oc', 'MarkerSize', 2);
plot( w_h2o(sg_dive_i), vmdepth(sg_dive_i), 'vr', 'MarkerSize', 2);
plot( w_h2o(sg_climb_i), vmdepth(sg_climb_i), '^m', 'MarkerSize', 2);

if 0
  % display pitch wrt pitch_intended:
  % display hack: offset pitch to use w_desired as pitch_intended
  d_pitch = (intended_pitch - w_d);
  plot(pitch(sg_dive_i)  - d_pitch,vmdepth(sg_dive_i),'.k')
  plot(pitch(sg_climb_i) + d_pitch,vmdepth(sg_climb_i),'.k')
end

if 0
  % Plot against engineering data, not GC data
  % Mark vertical depth for pumps and bleeds (and pitch changes) ala diveplot
  % Note that we might exclude apogee pump and initial bleed as 'bad samples'
  plot( w(ivb), vmdepth(ivb), '.k', 'MarkerSize', 15);
  plot( w(ipn), vmdepth(ipn), 'vk');
  plot( w(ipp), vmdepth(ipp), '^k');
  plot( w(irn), vmdepth(irn), '<k');
  plot( w(irp), vmdepth(irp), '>k');
else
  % Show GC actions at GC start depth and vehicle observed w,
  % which is computed differently from w
  gc_pumps_i = find(loginfo.gc_vbd_secs > 0);
  plot(loginfo.gc_ob_vertv(gc_pumps_i),loginfo.gc_depth(gc_pumps_i),'.r','MarkerSize',15);
  % compute energy (joules) used by these pumps:
  % volt amp secs == watt secs == joules
  % cumsum(loginfo.gc_vbd_volts(gc_pumps_i).*loginfo.gc_vbd_i(gc_pumps_i).*loginfo.gc_vbd_secs(gc_pumps_i))

  gc_pd_i = find(diff(loginfo.gc_pitch_ctl) < 0)+1; % pitch down
  plot(loginfo.gc_ob_vertv(gc_pd_i),loginfo.gc_depth(gc_pd_i),'vk');
  gc_pu_i = find(diff(loginfo.gc_pitch_ctl) > 0)+1; % pitch up
  plot(loginfo.gc_ob_vertv(gc_pu_i),loginfo.gc_depth(gc_pu_i),'^k');

  % sigh, we don't record roll_ctl in GC so we have to infer direction vs center
  % since we use bang-bang control, we can look at AD diff from center directly
  % to know port from starboard turns, since no incremental turns
  % we also compare against overshoot to deal with the noise in the roll system
  if (isfield(loginfo,'gc_roll_ad')) % datasets pre-2004 don't supply this
	  gc_dive_i = 1:gc_pumps_i(1); % C_ROLL_DIVE 
	  gc_climb_i = gc_pumps_i(2):length(loginfo.gc_vbd_secs); % C_ROLL_CLIMB starting at the climb pump
	  r_port_ovr = loginfo.R_PORT_OVSHOOT;
	  r_stbd_ovr = loginfo.R_STBD_OVSHOOT;
	  gc_pt_i = gc_dive_i(find((loginfo.gc_roll_ad(gc_dive_i) - c_roll_dive) > r_port_ovr)); % roll port dive
	  plot(loginfo.gc_ob_vertv(gc_pt_i),loginfo.gc_depth(gc_pt_i),'<k');
	  gc_pt_i = gc_climb_i(find((loginfo.gc_roll_ad(gc_climb_i) - c_roll_climb) > r_port_ovr)); % roll port climb
	  plot(loginfo.gc_ob_vertv(gc_pt_i),loginfo.gc_depth(gc_pt_i),'<k');
	  
	  gc_st_i = gc_dive_i(find((loginfo.gc_roll_ad(gc_dive_i) - c_roll_dive) < -r_stbd_ovr)); % roll stbd dive
	  plot(loginfo.gc_ob_vertv(gc_st_i),loginfo.gc_depth(gc_st_i),'>k');
	  gc_st_i = gc_climb_i(find((loginfo.gc_roll_ad(gc_climb_i) - c_roll_climb) < -r_stbd_ovr)); % roll stbd climb
	  plot(loginfo.gc_ob_vertv(gc_st_i),loginfo.gc_depth(gc_st_i),'>k');
  end
end

x_label  =  'Vertical Velocity [cm/s]';
xlabel(x_label, 'FontSize', 10); 
y_label = 'Depth [m]';
ylabel(y_label, 'FontSize', 10);
title(basetitle, 'FontSize', 12);

lg = legend('Based on Pressure Change', ...
       'Based on Buoyancy and Pitch', ...
       'Upwelling on Descent', ...
       'Upwelling on Ascent', ...
       'Location','SouthWest');
set(lg,'color','none'); % make transparent

v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

xtxt = v(1) + 0.6*(v(2) - v(1));
text(xtxt, ytxt, mass_vol_comment);

%
% plot sg temperature and salinity profiles
%

figure('Name','Temp & Salinity');
[tmin,tmax,tdel] = plot_bounds(temp,-3,40,[12,2; 6,1; 0,0.5]);
[smin,smax,sdel] = plot_bounds(salin,15,40,[12,2; 6,1; 0,0.5]);
hold on;
grid on; % this sets the grid against the first axis
% get the axes to align the S and T data (as well as D)
% must have the same NUMBER (or modulus) of ticks
% make sdel track tdel:
n_ticks = (tmax - tmin)/tdel;
sdel = (smax - smin)/n_ticks;

set(gca,'YDir','reverse');
line( temp(ctd_dive_i), ctd_depth_m(ctd_dive_i), 'Marker', 'v', 'Color', 'r', ...
      'MarkerSize', 2, 'LineStyle', 'none');
line( temp(ctd_climb_i), ctd_depth_m(ctd_climb_i), 'Marker', '^', 'Color', 'm', ...
      'MarkerSize', 2, 'LineStyle', 'none');
ax1 = gca;
set(ax1, 'XAxisLocation', 'top', 'XColor', 'r', ...
         'XLim', [tmin,tmax], 'XTick', [tmin:tdel:tmax],...
         'YLim', [dmin,dmax], 'YTick', [dmin:ddel:dmax]);
% xlabel('Temperature (\circC)', 'FontSize', 10); 
ylabel('Depth (m)', 'FontSize', 10);
% typically starts hi and goes low 
lg = legend('T (\circC) dive', 'T (\circC) climb', ...
       'Location','NorthWest');
set(lg,'color','none'); % make transparent

ax2 = axes('Position', get(ax1,'Position'),...
           'YAxisLocation', 'right',...
           'Color', 'none', 'XColor', 'b');
set(gca,'YDir','reverse', ...
        'XLim', [smin,smax], 'XTick', [smin:sdel:smax],...
        'YLim', [dmin,dmax], 'YTick', [dmin:ddel:dmax]);
line(salin(ctd_dive_i), ctd_depth_m(ctd_dive_i), 'Color', 'b', 'Marker', 'v', ...
     'MarkerSize', 2, 'LineStyle', 'none', 'Parent', ax2);
line(salin(ctd_climb_i), ctd_depth_m(ctd_climb_i), 'Color', 'c', 'Marker', '^', ...
     'MarkerSize', 2, 'LineStyle', 'none', 'Parent', ax2);
% JSB IOP shows potential density....
xlabel('Salinity (PSU)', 'FontSize', 10); 
% salinity starts lo and goes hi
lg = legend('S dive', 'S climb',...
       'Location','SouthEast');
set(lg,'color','none'); % make transparent
title(basetitle, 'FontSize', 12);

v = axis;
ytxt = v(4) + 0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

xtxt = v(1) + 0.6*(v(2)-v(1));
text(xtxt, ytxt, calibcomm);

%
%	Plot T vs S
%
% NOTE: Use tsv_diagram to see *all* the data; this reduces to sg points
figure('Name','T vs. S');
ax3 = gca;
set(ax3, 'XLim', [smin,smax], 'XTick', [smin:sdel:smax], ...
         'YLim', [tmin,tmax], 'YTick', [tmin:tdel:tmax]);
hold on;
grid on;
plot( salin(ctd_dive_i), temp(ctd_dive_i), 'vb', 'MarkerSize', 2);
plot( salin(ctd_climb_i), temp(ctd_climb_i), '^r', 'MarkerSize', 2);
x_label  =  'Salinity (psu)';
xlabel(x_label, 'FontSize', 10); 
y_label = 'Temperature (\circC)';
ylabel(y_label, 'FontSize', 10);
title(basetitle, 'FontSize', 12);
% let it try to find empty space
lg = legend('Descent', 'Ascent', 'Location','Best');
set(lg,'color','none'); % make transparent

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
sigma_levels = [5:0.5:28];
[C, h] =contour( Sg, Tg, sigma_grid, sigma_levels, 'k' );
clabel(C, h, sigma_levels); % this takes some time to label...

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

xtxt = v(1) + 0.6*(v(2)-v(1));
text(xtxt, ytxt, calibcomm);

%
% plot oxygen profile
%
if (sbe43_present | optode_present)
  figure('Name','O2');
  o2 = [solubility_micromoles_per_kg]; % O2 solubility_micromoles_per_kg is always available
  if (sbe43_present)
    o2 = [o2; oxygen];
  end
  if (optode_present)
    o2 = [o2; optode_instrument_oxygen]; % instrument O2 always available
    if (exist('optode_oxygen','var'))
      optode_oxygen_present = 1;
      o2 = [o2;optode_oxygen];      
    else
      optode_oxygen_present = 0;      
    end
  end
  [omin,omax,odel] = plot_bounds(o2,0,375,[0,25]); % umol/kG
  clear o2;
  o2_xlabel = 'Dissolved Oxygen [\muM/kg]';
  hold on;
  grid on;
  set(gca,'YDir','reverse');
  legend_text = {};
  
  if (sbe43_present)
    line( oxygen(oxy_dive_i), oxy_depth_m(oxy_dive_i), 'Marker', 'v', 'Color', 'b', ...
          'MarkerSize', 1, 'LineStyle', '-');
    legend_text(end+1) = {'SBE43 O_2 dive'};
    line( oxygen(oxy_climb_i), oxy_depth_m(oxy_climb_i), 'Marker', '^', 'Color', 'r', ...
          'MarkerSize', 1, 'LineStyle', '-');
    legend_text(end+1) = {'SBE43 O_2 climb'};
  end
  
  if (optode_present)
    line( optode_instrument_oxygen(opt_dive_i), opt_depth_m(opt_dive_i), 'Marker', 'v', 'Color', 'g', ...
          'MarkerSize', 1, 'LineStyle', '--');
    legend_text(end+1) = {'Optode O_2 dive'};
    line( optode_instrument_oxygen(opt_climb_i), opt_depth_m(opt_climb_i), 'Marker', '^', 'Color', 'y', ...
          'MarkerSize', 1, 'LineStyle', '--');
    legend_text(end+1) = {'Optode O_2 climb'};
    if (optode_oxygen_present)
      line( optode_oxygen(opt_dive_i), opt_depth_m(opt_dive_i), 'Marker', 'v', 'Color', 'b', ...
            'MarkerSize', 1, 'LineStyle', '-');
      legend_text(end+1) = {'Corrected Optode O_2 dive'};
      line( optode_oxygen(opt_climb_i), opt_depth_m(opt_climb_i), 'Marker', '^', 'Color', 'r', ...
            'MarkerSize', 1, 'LineStyle', '-');
      legend_text(end+1) = {'Corrected Optode O_2 climb'};
    end
  end
  
  line( solubility_micromoles_per_kg(ctd_dive_i), ctd_depth_m(ctd_dive_i), 'Marker', 'v', 'Color', 'c', ...
        'MarkerSize', 1, 'LineStyle', '-');
  legend_text(end+1) = {'Saturation O_2 dive'};
  line( solubility_micromoles_per_kg(ctd_climb_i), ctd_depth_m(ctd_climb_i), 'Marker', '^', 'Color', 'm', ...
        'MarkerSize', 1, 'LineStyle', '-');
  legend_text(end+1) = {'Saturation O_2 climb'};
  
  ax1 = gca;
  set(ax1, 'XAxisLocation', 'top', 'XColor', 'r', 'XLim', [omin,omax], ...
           'XTick', [omin:odel:omax],'YLim', [dmin,dmax]);
  xlabel( o2_xlabel );
  ylabel( 'Depth (m)');
  lg = legend(legend_text,'Location','SouthEast');
  set(lg,'color','none'); % make transparent
  
  v = axis;
  ytxt = v(4) + 0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel);
  xtxt = v(1) +0.4*(v(2)-v(1));
  text(xtxt , ytxt, titlestring);
  
  xtxt = v(1);
  ytxt = v(3) - 0.05*(v(4)-v(3));
  text(xtxt , ytxt, base_file);
  
  xtxt = v(1) + 0.6*(v(2)-v(1));
  text(xtxt, ytxt, ox_comment);
end

%
% plot backscatter profiles
%
if (wetlabs_present) % got a WETlabs sensor
	% ignore any reference data (old style sensor)
	% report any red/blue/fluor info
	figure('Name','Backscatter');
	scttr = []; % See what red or blue channel data we have
	if (exist('red_scttr','var'))
		scttr = [scttr; red_scttr];
	end
	if (exist('blue_scttr','var'))
		scttr = [scttr; blue_scttr];
	end
	[bmin,bmax,bdel] = plot_bounds(scttr,0,2400,[0,200]);
	clear scttr;
	hold on;
	grid on;
	set(gca,'YDir','reverse');
	
	legend_text = {};
	if (exist('red_scttr','var'))
		line( red_scttr(wet_dive_i), wet_depth_m(wet_dive_i), 'Marker', 'v', 'Color', 'm', ...
			  'MarkerSize', 1, 'LineStyle', '-');
		legend_text(end+1) = {'Red scatter descent'};
		line( red_scttr(wet_climb_i), wet_depth_m(wet_climb_i), 'Marker', '^', 'Color', 'r', ...
			  'MarkerSize', 1, 'LineStyle', '-');
		legend_text(end+1) = {'Red scatter ascent'};
	end
	if (exist('blue_scttr','var'))
		line( blue_scttr(wet_dive_i), wet_depth_m(wet_dive_i), 'Marker', 'v', 'Color', 'c', ...
			  'MarkerSize', 1, 'LineStyle', '-');
		legend_text(end+1) = {'Blue scatter descent'};
		line( blue_scttr(wet_climb_i), wet_depth_m(wet_climb_i), 'Marker', '^', 'Color', 'b', ...
			  'MarkerSize', 1, 'LineStyle', '-');
		legend_text(end+1) = {'Blue scatter ascent'};
	end	
	ax1 = gca;
	set(ax1, 'XLim', [bmin,bmax], ...
			 'XTick', [bmin:bdel:bmax],'YLim', [dmin,dmax]);
	xlabel( 'Backscatter Output [counts]' );
	ylabel( 'Depth [counts]');
	lg = legend(legend_text,'Location','SouthEast');
	set(lg,'color','none'); % make transparent
	
	title(basetitle, 'FontSize', 12);
	
	v = axis;
	ytxt = v(4) + 0.1*(v(4)-v(3));
	xtxt = v(1);
	text(xtxt , ytxt, timelabel);
	
	%
	% plot fluorometer profile
	%
	% IOP site skips showing fluor_temp on the same graph
	% This code assumes you always get a fluorescence channel and associated temp
	show_fluor_temp = 0;
	if (show_fluor_temp)
		[fmin,fmax,fdel] = plot_bounds([fluor;fluor_temp],0,2000,[0,200]);
	else
		[fmin,fmax,fdel] = plot_bounds(fluor,0,2000,[0,200]);
	end
	figure('Name','Fluorometer');
	hold on;
	grid on;
	set(gca,'YDir','reverse');
	
	legend_text = {};
	line( fluor(wet_dive_i), wet_depth_m(wet_dive_i), 'Marker', 'v', 'Color', 'b', ...
		  'MarkerSize', 1, 'LineStyle', '-');
	legend_text(end+1) = {'Descent'};
	line( fluor(wet_climb_i), wet_depth_m(wet_climb_i), 'Marker', '^', 'Color', 'r', ...
		  'MarkerSize', 1, 'LineStyle', '-');
	legend_text(end+1) = {'Ascent'};
	if (show_fluor_temp)
		line( fluor_temp(wet_dive_i), wet_depth_m(wet_dive_i), 'Marker', 'v', 'Color', 'c', ...
			  'MarkerSize', 1, 'LineStyle', '-');
		legend_text(end+1) = {'Temperature descent'};
		line( fluor_temp(wet_climb_i), wet_depth_m(wet_climb_i), 'Marker', '^', 'Color', 'm', ...
			  'MarkerSize', 1, 'LineStyle', '-');
		legend_text(end+1) = {'Temperature ascent'};
	end
	ax1 = gca;
	set(ax1, 'XLim', [fmin,fmax], ...
			 'XTick', [fmin:fdel:fmax],'YLim', [dmin,dmax]);
	xlabel( 'Fluorometer Output [counts]' );
	ylabel( 'Depth [counts]');
	
	lg = legend(legend_text,'Location','SouthEast');
	set(lg,'color','none'); % make transparent
	
	title(basetitle, 'FontSize', 12);
	v = axis;
	ytxt = v(4) + 0.1*(v(4)-v(3));
	xtxt = v(1);
	text(xtxt , ytxt, timelabel);
end % got a wlbb2f

if(ops_charts)
  %
  % find linear fit between pitch and pitch_control
  %
  
  ymax = 60; % PARAMETER max positive measured pitch (degrees)
  ymin = -60; % PARAMETER min positive measured pitch (degrees)
  ytic = 10;
  xmax = 3.5; % PARAMETER max positive pitch control position (cm)
  xmin = -3.5; % PARAMETER max megative pitch control position (cm)
  xtic = 0.5;
  
  ip = find(pitch_control > xmin & pitch_control < xmax);
  X = [ones(size(pitch_control(ip))) pitch_control(ip)];
  A = X\pitch(ip);
  inv_pitch_cnv = 1/loginfo.PITCH_CNV; % 217.39 % AD counts per cm
  c_pitch_imp = c_pitch - ( A(1)/A(2) )*inv_pitch_cnv;
  
  pitch_control_max = (pitch_max_cnts - c_pitch)/inv_pitch_cnv;
  pitch_control_min = (pitch_min_cnts - c_pitch)/inv_pitch_cnv;
  pitch_control_max_imp = ...
      (pitch_max_cnts - c_pitch_imp)/inv_pitch_cnv;
  pitch_control_min_imp = ...
      (pitch_min_cnts - c_pitch_imp)/inv_pitch_cnv;
  
  %
  % plot scatter of pitch with pitch_control
  %
  
  figure('Name','Pitch');
  plot(pitch_control, pitch, '.b');
  hold on;
  grid on;
  plot(pitch_control_intended, pitch_intended, '.r')
  plot([0], [0], '+r', 'MarkerSize', 10)
  axis([xmin xmax ymin ymax]);
  set(gca, 'XTick', [xmin:xtic:xmax], 'YTick', [ymin:ytic:ymax]);
  xlabel('Pitch Control [cm]');
  ylabel('Pitch [\circ]');
  title(basetitle, 'FontSize', 12);
  
  ix = [xmin:0.01*(xmax-xmin):xmax];
  plot(ix, A(1) + A(2)*ix, '-k');
  
  xtxt = xmin + 0.1*(xmax-xmin);
  ytxt = ymin + 0.95*(ymax-ymin);
  comment = strcat( 'C_{PITCH} =', num2str(c_pitch) );
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.9*(ymax-ymin);
  comment = strcat( 'Pitch control range =', ...
                    num2str(pitch_control_min), ' to ', ...
                    num2str(pitch_control_max), ' cm'); 
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.85*(ymax-ymin);
  comment = strcat( 'Implied Pitch Center =', num2str(c_pitch_imp) );
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.8*(ymax-ymin);
  comment = strcat( 'Implied pitch control range =', ...
                    num2str(pitch_control_min_imp), ' to ', ...
                    num2str(pitch_control_max_imp), ' cm'); 
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.75*(ymax-ymin);
  comment = strcat( 'Assumed pitch control gain =', ...
                    num2str(inv_pitch_cnv), ' [counts/cm aft]' );
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.3*(xmax-xmin);
  ytxt = ymin + 0.14*(ymax-ymin);
  comment = strcat( 'Pitch [\circ] =', num2str(A(1)), '\circ +', ...
                    'Pitch Control [cm] *', num2str(A(2)), '\circ/cm' );
  text(xtxt, ytxt, comment);
  
  v = axis;
  ytxt = v(3) -0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel);

end

if(ops_charts)
  %
  % find linear fit between roll and roll_control
  %
  
  inv_roll_cnv = 1/loginfo.ROLL_CNV; % 35.37 AD counts per degree

  X = [ones(size(roll_control)) roll_control];
  B = X\roll;
  c_roll_imp = c_roll - ( B(1)/B(2) )*inv_roll_cnv;

  roll_control_max = (roll_max_cnts - c_roll)/inv_roll_cnv;
  roll_control_min = (roll_min_cnts - c_roll)/inv_roll_cnv;
  roll_control_max_imp = ...
      (roll_max_cnts - c_roll_imp)/inv_roll_cnv;
  roll_control_min_imp = ...
      (roll_min_cnts - c_roll_imp)/inv_roll_cnv;
  
  Xd = [ones(size(roll_control(sg_dive_i))) roll_control(sg_dive_i)];
  Bd = Xd\roll(sg_dive_i);
  c_roll_dive_imp = c_roll_dive - ( Bd(1)/Bd(2) )*inv_roll_cnv;
  
  Xc = [ones(size(roll_control(sg_climb_i))) roll_control(sg_climb_i)];
  Bc = Xc\roll(sg_climb_i);
  c_roll_climb_imp = c_roll_climb - ( Bc(1)/Bc(2) )*inv_roll_cnv;
  
  
  %
  % plot scatter of roll with roll_control
  %
  
  figure('Name','Roll');
  plot(roll_control(sg_dive_i), roll(sg_dive_i), '.b');
  hold on;
  grid on;
  plot(roll_control(sg_climb_i), roll(sg_climb_i), '.r')
  plot([0], [0], '+k', 'MarkerSize', 10)
  ymax = 40;
  xmax = 55;
  xmin = -55;
  ymin = -40;
  axis([xmin xmax ymin ymax]);
  xlabel('Roll Control [\circ]');
  ylabel('Roll [\circ]');
  title(basetitle, 'FontSize', 12);
  % the plot is cluttered...
  lg = legend('Dive', 'Climb', 'Location', 'SouthWest');
  set(lg,'color','none'); % make transparent
  
  ix = [xmin:0.01*(xmax-xmin):xmax];
  plot(ix, B(1) + B(2)*ix, '-g');
  plot(ix, Bd(1) + Bd(2)*ix, '-b');
  plot(ix, Bc(1) + Bc(2)*ix, '-r');
  
  xtxt = xmin + 0.1*(xmax-xmin);
  ytxt = ymin + 0.95*(ymax-ymin);
  comment = strcat( 'C_{ROLL}=', num2str(c_roll) );
  text(xtxt, ytxt, comment);
  xtxt = xmin + 0.4*(xmax-xmin);
  ytxt = ymin + 0.95*(ymax-ymin);
  comment = strcat( 'C_{ROLL}^{dive} =', num2str(c_roll_dive) );
  text(xtxt, ytxt, comment);
  xtxt = xmin + 0.7*(xmax-xmin);
  ytxt = ymin + 0.95*(ymax-ymin);
  comment = strcat( 'C_{ROLL}^{climb} =', num2str(c_roll_climb) );
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.1*(xmax-xmin);
  ytxt = ymin + 0.9*(ymax-ymin);
  comment = strcat( 'Roll control range =', ...
                    num2str(roll_control_min), ' to ', ...
                    num2str(roll_control_max), ' \circ'); 
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.1*(xmax-xmin);
  ytxt = ymin + 0.85*(ymax-ymin);
  comment = strcat( 'Implied C_{ROLL} =', num2str(c_roll_imp) );
  text(xtxt, ytxt, comment);
  xtxt = xmin + 0.4*(xmax-xmin);
  comment = strcat( 'Implied C_{ROLL}^{dive} =', ...
                    num2str(c_roll_dive_imp) );
  text(xtxt, ytxt, comment);
  xtxt = xmin + 0.7*(xmax-xmin);
  comment = strcat( 'Implied C_{ROLL}^{climb} =', ...
                    num2str(c_roll_climb_imp) );
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.1*(xmax-xmin);  
  ytxt = ymin + 0.8*(ymax-ymin);
  comment = strcat( 'Implied roll control range =', ...
                    num2str(roll_control_min_imp), ' to ', ...
                    num2str(roll_control_max_imp), ' \circ'); 
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.75*(ymax-ymin);
  comment = strcat( 'Assumed roll control gain =', ...
                    num2str(inv_roll_cnv), ' [counts/\circ to stbd]' );
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.25*(xmax-xmin);
  ytxt = ymin + 0.14*(ymax-ymin);
  comment = strcat( 'Roll [\circ] =', num2str(B(1)), ...
                    '\circ + Roll Control [\circ] *', num2str(B(2)) );
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.08*(ymax-ymin);
  comment = strcat( 'Roll_{dive} [\circ] =', num2str(Bd(1)), ...
                    '\circ + Roll Control [\circ] *', num2str(Bd(2)) );
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.02*(ymax-ymin);
  comment = strcat( 'Roll_{climb} [\circ] =', num2str(Bc(1)), ...
                    '\circ + Roll Control [\circ] *', num2str(Bc(2)) );
  text(xtxt, ytxt, comment);
  
  v = axis;
  ytxt = v(3) -0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel);
end

if(ops_charts)
  %
  % find linear fit between turn_rate with roll_control_counts
  %
  roll_control_counts_dive = c_roll_dive + roll_control(sg_dive_i)*inv_roll_cnv;
  Xd = [ones(size(roll_control_counts_dive)) roll_control_counts_dive];
  Cd = Xd\turn_rate(sg_dive_i);
  c_roll_turn_dive = -Cd(1)/Cd(2);

  roll_control_counts_climb = c_roll_climb + roll_control(sg_climb_i)*inv_roll_cnv;
  Xc = [ones(size(roll_control_counts_climb)) roll_control_counts_climb];
  Cc = Xc\turn_rate(sg_climb_i);
  c_roll_turn_climb = -Cc(1)/Cc(2);
  
  %
  % plot scatter of turn rate with roll_control
  %
  figure('Name','Roll Rate');
  plot(roll_control(sg_dive_i), turn_rate(sg_dive_i), '.b');
  hold on;
  grid on;
  plot(roll_control(sg_climb_i), turn_rate(sg_climb_i), '.r')
  plot([0], [0], '+k', 'MarkerSize', 10)
  ymax = 2.5;
  xmax = 55;
  xmin = -55;
  ymin = -2.5;
  axis([xmin xmax ymin ymax]);
  xlabel('Roll Control [\circ]');
  ylabel('Turn rate [\circ/s] (starboard positive)');
  title(basetitle, 'FontSize', 12);
  lg = legend('Dive', 'Climb', 'Location','NorthWest');
  set(lg,'color','none'); % make transparent
  
  ix = [xmin:0.01*(xmax-xmin):xmax];
  plot(ix, Cd(1) + Cd(2)*c_roll_dive + Cd(2)*inv_roll_cnv*ix, '-b');
  plot(ix, Cc(1) + Cc(2)*c_roll_dive + Cc(2)*inv_roll_cnv*ix, '-r');
  
  text(xtxt, ytxt, comment);
  xtxt = xmin + 0.4*(xmax-xmin);
  ytxt = ymin + 0.95*(ymax-ymin);
  comment = strcat( 'C_{ROLL}^{dive} =', num2str(c_roll_dive) );
  text(xtxt, ytxt, comment);
  xtxt = xmin + 0.7*(xmax-xmin);
  ytxt = ymin + 0.95*(ymax-ymin);
  comment = strcat( 'C_{ROLL}^{climb} =', num2str(c_roll_climb) );
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.4*(xmax-xmin);
  ytxt = ymin + 0.85*(ymax-ymin);
  comment = strcat( 'Implied C_{ROLL}^{dive} =', ...
                    num2str(c_roll_turn_dive) );
  text(xtxt, ytxt, comment);
  xtxt = xmin + 0.7*(xmax-xmin);
  comment = strcat( 'Implied C_{ROLL}^{climb} =', ...
                    num2str(c_roll_turn_climb) );
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.1*(xmax-xmin);  
  ytxt = ymin + 0.75*(ymax-ymin);
  comment = strcat( 'Assumed roll control gain =', ...
                    num2str(inv_roll_cnv), ' [counts/\circ to stbd]' );
  text(xtxt, ytxt, comment);
  
  xtxt = xmin + 0.1*(xmax-xmin);
  ytxt = ymin + 0.14*(ymax-ymin);
  comment = strcat( 'Dive Turn Rate [\circ/s] =', num2str(Cd(1)), ...
                    '\circ/s + Roll Control [counts] *', num2str(Cd(2)) );
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.1*(ymax-ymin);
  comment = strcat( 'Dive Turn Rate [\circ/s] =', ...
                    num2str(Cd(1)+Cd(2)*c_roll_dive), '\circ/s + Roll Control [\circ]', ...
                    num2str(Cd(2)*inv_roll_cnv) );
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.06*(ymax-ymin);
  comment = strcat( 'Climb Turn Rate [\circ/s] =', num2str(Cc(1)), ...
                    '\circ/s + Roll Control [counts] *', num2str(Cc(2)) );
  text(xtxt, ytxt, comment);
  ytxt = ymin + 0.02*(ymax-ymin);
  comment = strcat( 'Climb Turn Rate [\circ/s] =', ...
                    num2str(Cc(1)+Cc(2)*c_roll_climb), '\circ/s + Roll Control [\circ]', ...
                    num2str(Cc(2)*inv_roll_cnv) );
  text(xtxt, ytxt, comment);
  
  v = axis;
  ytxt = v(3) -0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel);
end

if(ops_charts)
  %
  % plot VBD configuration
  %
  % find 'surface' density in top/shallowest 1.0 meter
  zmax = max(zpos); % shallowest point to handle yoyo dive (normally ~ -0.5 m)
  i_1m = find(zpos > (zmax - 1.0) & zpos < zmax); % PARAMETER 1meter
  density_1m = mean(density_raw(i_1m)); % what is density 1m lower?
  sigma_1m = density_1m - 1000;
  % find apogee depth within 1 meter of the deepest point
  zmin = min(zpos);
  i_apogee = find(zpos < zmin + 1); % PARAMETER density at 1 meter from max depth 
  density_apogee = mean(density_raw(i_apogee));
  sigma_apogee = density_apogee - 1000;
  vbd_1m = cc_per_m3*mass*(1./density_1m - 1./density_apogee); % What vbd throw do we need to go between deep and surface points?
  % ignore this on yoyo dive?
  cc_surf_min = vbd_1m + 150; % PARAMETER 150g (or cc) to get the antenna out of the water
  c_vbd_true = c_vbd + vbdbias/vbd_cnv;

  [sigmamin,sigmamax,sigmatic] = plot_bounds(sigma_t,18,40);
  sigmagrid = [sigmamin: sigmatic: sigmamax];
  % min/max vbd in cc assuming well centered (wrong for DG: 50% more throw)
  vbdmin = -300;
  vbdmax = 700;
  vbdtic = 100;
  vbdgrid = [vbdmin: vbdtic: vbdmax];
  [SIGMA, VBD] = meshgrid(sigmagrid,vbdgrid);
  figure('Name','VBD');
  title( sprintf('%s - buoyancy contoured', basetitle),'FontSize', 12);
  hold on
  set(gca, 'XLim', [sigmamin, sigmamax], 'XTick', sigmagrid, ...
           'YLim', [vbdmin, vbdmax], 'YTick', vbdgrid);
  BUOY = VBD - cc_per_m3*mass./(1000+SIGMA) + vol0 - vbdbias;
  buoy_levels = [-500:50:500];
  [C, h] = contour(SIGMA, VBD, BUOY, buoy_levels, 'c');
  if ~isempty(C) & ~isempty(h) & ~isempty(buoy_levels)
    clabel(C, h, buoy_levels, 'labelspacing', 288);
  end
  grid on;
  [C, h] = contour( SIGMA, VBD, BUOY, [-0.1 0.1], 'b'); % highlight zero buoyancy contour
  [C, h] = contour( SIGMA, VBD, BUOY, [100 200], 'm'); % highlight +100 and +200 cc buoyancy contours
  xlabel('\sigma_t [kg/m^3]');
  ylabel('VBD [cc]');
  
  plot(sigma_apogee, [0], '+r'); % at what signa are we neutral with current C_VBD?
  % plot(sigma_apogee, vbdbias, '+r');
  plot([sigmamin sigmamax], [vbd_min, vbd_min], '-k'); % show max positive buoyancy
  plot([sigmamin sigmamax], [vbd_max, vbd_max], '-k'); % show max negative buoyancy
  plot(sigma_1m, vbd_1m, '+r'); 
  plot(sigma_1m, cc_surf, '+r');
  plot(sigma_1m, cc_surf_min, 'xr');
  
  plot(sigma_t(ctd2sg_valid_i), vbd(sg_ctd_valid_i), '-k'); % plot a black line though measured sigma/vbd trajectory
  scatter(sigma_t(ctd2sg_valid_i), vbd(sg_ctd_valid_i), 10.0, vmtime(sg_ctd_valid_i), 'filled'); % enhance trajectory
  
  text(sigma_apogee + 0.1, 0, strcat('Dive C_{VBD}=', num2str(c_vbd)));
  % text(sigma_apogee + 0.1, vbdbias, strcat('True C_{VBD}=', num2str(c_vbd_true)) );
  text(sigmamax + 0.1, vbd_min, num2str(vbd_max_cnts));
  text(sigmamax + 0.1, vbd_max, num2str(vbd_min_cnts));
  text(sigma_1m + 0.1, cc_surf, strcat('CC_{SM}=', num2str(cc_surf)));
  text(sigma_1m + 0.1, vbd_1m, sprintf('VBD_{1m}=%.0f', vbd_1m));
  text(sigma_1m + 0.1, cc_surf_min, sprintf('CC_{surf min}=%3.0f', cc_surf_min));
  
  v = axis;
  ytxt = v(3) -0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel);
  
  xtxt = v(1) + 0.6*(v(2) - v(1));
  text(xtxt, ytxt, biaslabel);
  ytxt = v(3) -0.07*(v(4)-v(3));
  text(xtxt, ytxt, mass_vol_comment);
end

if(1 & ops_charts)
  %
  % plot VBD configuration
  %
  iref = find(ctd_sg_press == max(ctd_sg_press),1,'first');
  % compute the volume we need to be at apogee given the density insitu at max pressure
  vol_apogee = cc_per_m3*mass/sw_dens(salin_raw(iref), temp_raw(iref), ctd_sg_press(iref));
  kg2g = 1000;
  vol_comp = kg2g*mass_comp./cml_dens(temp_raw(iref), ctd_sg_press(iref));
  vol_comp_ref = kg2g*mass_comp/cml_dens(temp_ref, 0);
  % what is volume of hull at surface, independent of compressee?
  vol_hull = volmax - vol_comp_ref;
  vol_hull_compression_factor = (exp(-abs_compress*ctd_sg_press(iref) + therm_expan*(temp_raw(iref) - temp_ref)));
  vol_delta_neutral = (( volmax - vol_hull + ...
                         vol_hull*vol_hull_compression_factor + ...
                         (vol_comp - vol_comp_ref)) - ....
                       vol_apogee);
  c_vbd_true = fix(vbd_min_cnts - vbd_cnts_per_cc*vol_delta_neutral);

  ctd_sigma = sw_pden(salin_raw, temp_raw, ctd_sg_press, 0) - 1000;
  [sigmamin,sigmamax,sigmatic] = plot_bounds(ctd_sigma,10,30,[0,0.5]);
  sigmadel = sigmatic/10;
  sigmagrid = [sigmamin: sigmadel: sigmamax];
  grid_valid = zeros(length(sigmagrid),1);
  vbdtic = 100;
  vbdmin = -500;
  vbdmax = 500;
  vbdgrid = [vbdmin: vbdtic: vbdmax];
  [SIGMA, VBD] = meshgrid(sigmagrid, vbdgrid);
  for js=1:length(sigmagrid)
    it = find(ctd_sigma < sigmagrid(js) + 0.5*sigmadel & ...
              ctd_sigma >= sigmagrid(js) - 0.5*sigmadel);
    if(size(it) >= 1)
      grid_valid(js) = 1;
    else
      grid_valid(js) = 0;
    end
    temp_av_grid(js) = mean(temp_raw(it));
    pres_av_grid(js) = mean(ctd_sg_press(it));
    density_insitu_av_grid(js) = mean(density_insitu(it));
    sigma_av_grid(js) = mean(ctd_sigma(it));
  end
  ig = find(grid_valid > 0);
  [TEMP, VBDC] = meshgrid(temp_av_grid(ig), vbdgrid);
  [PRES, VBDC] = meshgrid(pres_av_grid(ig), vbdgrid);
  [SIGMC, VBDC] = meshgrid(sigmagrid(ig), vbdgrid);
  VOL_COMP = kg2g*mass_comp./cml_dens(TEMP, PRES);
  [DENSI, VBDC] = meshgrid(density_insitu_av_grid(ig), vbdgrid);
  
  figure
  
  title(basetitle, 'FontSize', 12)
  hold on
  set(gca, 'XLim', [sigmamin, sigmamax], 'XTick', ...
           [sigmamin:sigmatic:sigmamax], ...
           'YLim', [vbdmin, vbdmax], 'YTick', vbdgrid);
  BUOY_nc = VBD + volmax + (c_vbd - vbd_min_cnts)/vbd_cnts_per_cc ...
            - cc_per_m3*mass./(1000+SIGMA);
  BUOY = VBDC + volmax + (volmax - vol_comp_ref)* ...
         (exp(-abs_compress*PRES + therm_expan*(TEMP - temp_ref)) -1) + ...
         VOL_COMP - vol_comp_ref + (c_vbd - vbd_min_cnts)/vbd_cnts_per_cc ...
         - cc_per_m3*mass./DENSI;
  buoy_levels = [-1000:50:1000];
  [C, h] = contour( SIGMA, VBD, BUOY_nc, buoy_levels, 'c');
  clabel(C, h, buoy_levels, 'Color', 'c');
  grid on;
  [C, h] = contour( SIGMA, VBD, BUOY_nc, [-0.1 0.1], 'm');
  [C, h] = contour( SIGMA, VBD, BUOY_nc, [100 200], 'm');
  [C, h] = contour( SIGMC, VBDC, BUOY, buoy_levels, 'b');
  clabel(C, h, buoy_levels, 'Color', 'b');
  [C, h] = contour( SIGMC, VBDC, BUOY, [-0.1 0.1], 'r');
  [C, h] = contour( SIGMC, VBDC, BUOY, [100 200], 'r');
  xlabel('\sigma_\theta [kg/m^3]');
  
  ylabel('VBD [cc]');
  
  plot([sigmamin sigmamax], [vbd_min, vbd_min], '-k');
  plot([sigmamin sigmamax], [vbd_max, vbd_max], '-k');
  plot([sigmamin sigmamax], [0 0], '--k');
  
  plot(ctd_sigma, ctd_vbd, '-k');
  
  text(sigmamax + 0.1, vbd_min, num2str(vbd_max_cnts));
  text(sigmamax + 0.1, vbd_max, num2str(vbd_min_cnts));

  plot(ctd_sigma(iref), 0, '+r')
  plot(ctd_sigma(iref), (c_vbd_true - c_vbd)/vbd_cnts_per_cc, '+r')
  text(ctd_sigma(iref) + sigmatic/10, (c_vbd_true - c_vbd)/vbd_cnts_per_cc, ...
       strcat('C_V_B_D^t^r^u^e=', num2str(c_vbd_true)));
  text(ctd_sigma(iref) + sigmatic/10, 0, ...
       strcat('C_V_B_D=', num2str(c_vbd)));
  
  v = axis;
  ytxt = v(3) -0.1*(v(4)-v(3));
  xtxt = v(1);
  text(xtxt , ytxt, timelabel);
  
  xtxt = v(1) + 0.6*(v(2) - v(1));
  text(xtxt, ytxt, biaslabel);
  ytxt = v(3) -0.07*(v(4)-v(3));
  text(xtxt, ytxt, mass_vol_comment);
end

