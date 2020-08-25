% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% regress_vbd2.m

% This code handles finding first vbdbias (volmax), then, once found, abc
% It also handles finding thermal expansion coefficients and will handle (per CCE)
% the impact of any compressee on vehicle volume and density

%
%   Read sg data files for a specified set of dive cycles and 
%	minimize the rms difference between vertical velocity from 
%	pressure and vertical velocity from a steady flight model 
%	with respect to vbd bias, lift & drag coeffs a, b, c, and pressure hull
%	thermal expansion and compressibility
%
%  functions called:
%    sg_calib_constants.m - establishes glider calibration constants
%    various seawater routines
%    get_dive_data.m - data and results from nc or nc+mat
%    glide_slope.m - calculate glide speed and angle starting from from w and pitch
%    flightvec2.m   - calculate glide speed and angle starting from from buoyancy and pitch
%    w_misfit_func.m  - minimizes w error w.r.t. vbd bias, vbd bias a, b, & c, or therm_expan and abs_compress
%    w_rms_func.m   - calculate w error between measured and flight model

% This version uses a modified version of flightvec (flightvec2) to handle interpolation over bad points

global HIST; % regression history
good_w_rms_fit = 2.7; % PARAMETER cm/s if less than this, fit was ok...

% flightvec2 limit paramaeters
min_rel_q_inc = 0.001; % desired tolerance of relative q change
max_num_iters = 25; % maximum iterations if tolerance is not reached (non-convergence)
% hd_s = -1/4;
gravity = 9.82; % m/s^2

% PARAMETERS used to filter which data are used in regressions
% avoid acceleration because of pumping and bleeding and big rolls
vbddiffmax = 0.5; % cc/s 

% On excluding large rolls:
% Fritz points out that this is a problem only if you are rolling a lot
% between extreme angles but if you roll hard and stay that way you have
% good flight until you roll back (you aren't de/accelerating for sure).  So
% those hard-rolled points ought to be included.  A quick fix would be to
% expand rollmax to 50 degrees, which takes nearly all points, even on
% compass cal dives where we force her way over both sides
rollmax = 10; % degrees CONSIDER: 50
rollmin = -rollmax; % degrees
vbdbias = 0.0; % starting bias offset if none provided in sg_calib_constants

if (exist('regress_dives','var'))
	dives = regress_dives';
	clear regress_dives;
else
	dives = ask_which_runs(available_profiles());
end

% PARAMETERS that control which quantities are optimized
fit_abc = 0;
% To determine lift, drag and implied drag (abc) you need to investigate
% the space of horizontal and vertical speeds so we can best fit the power
% flight eqn.  To do this you need different pitches and buoyancies.
% Actual test is whether the dives cover a range of bouy and pitch desired
% you might be able to do ok with 3 dives but 4 is better
if (length(dives) > 2)
	fit_abc  = 1;
end

% By design/assertion the SG hull absolute compressibility is same as
% seawater and, for DG, the absolute compressibility of was measured.  The
% *only* time you need to run this part of the regression is when you have a
% 'mongrel' set of fairings with extra syntactic foam, etc. that will change
% the hull/fairing absolute compressibility and/or thermal compressibility
fit_exp_comp = 0; % CONTROL DO NOT USE unless you have a completely different hull and fairing design from SG and DG!!
if ~mass_comp % if you have compressee then the hull/syntactic effects are masked? (DG will never compute the fit_exp_comp)
	syntactic_foam_cc = input('How many ccs of added syntactic foam as rings or blocks are present, if any? ');
	% Running the regressions on sg190 SPURS 2013 with and without 3k cc of foam rings changes the 
	% absolute compressibility and thermal expansion factors a small amount from nominal values
	% to 4.6599e-6 and 9.0047e-5 and improves RMSE by ~0.1cm/s.  Hardly worth it.
	if (syntactic_foam_cc)
		if (syntactic_foam_cc > 3000) % PARAMETER cc of HPC-100 at least
			fit_exp_comp = 1; % probably should run the regression
		else
			fprintf(1,'Ignoring effect of syntactic foam.\n');
		end
	end
end
	


%define constants
radpd = pi/180.;
cm2m = 0.01;
m2cm = 100;
mS_per_cm_2_S_per_m = 0.1;
kg2g = 1000;

sg_calib_constants; % get id_str (get_dive_data deals with all the other constants)
vbdbias_0 = vbdbias; % record initial vbdbias (from sg_calib_constants)
%intitialize storage for list of c_vbd values
first_c_vbd = NaN;
dives_used = [];
W = [];
YD = [];
P = [];
T = [];
S = [];
V = [];
PI = [];
deepglider = 0; % assume SG
% loop on dive cycles and concatenate data from all dives
loop_pass = 1;
for ir = dives'
  
	id_dive = sprintf('%s%04d',id_str,ir);
	[loginfo,eng,results] = get_dive_data(id_dive,1);
	if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
		fprintf(1,'Unable to include data for dive %d\n',ir);
		continue;
	end
	unpack_data;
	if (isempty(dive_i))
		fprintf(1,'No points on dive %d; skipping\n',ir);
		continue;
	end
	if (isempty(climb_i))
		fprintf(1,'No points on climb %d; skipping\n',ir);
		continue;
	end
	
	fprintf(1,'p%s\n',id_dive);
	if (isfield(loginfo,'DEEPGLIDER'))
		deepglider = loginfo.DEEPGLIDER;
	end
	sg_epoch_time = results.time;
	ctd_epoch_time = results.ctd_time;
	
	% interpolate SG ops data to ctd time so we can use QC results to filter
	r_time = ctd_epoch_time;
	% CCE says no need to adjust for pitchbias, if any
	r_pitch = interp1d(sg_epoch_time,eng.pitchAng,ctd_epoch_time,'linear'); 
	r_roll = interp1d(sg_epoch_time, eng.rollAng, ctd_epoch_time,'linear');
	if (isnan(first_c_vbd))
		first_c_vbd = loginfo.C_VBD;
		vbd_cnts_per_cc   = 1./loginfo.VBD_CNV;% same for all dives
	end
	% bring vbdc into alignment with first C_VBD
	r_vbd = eng.vbdCC + (loginfo.C_VBD - first_c_vbd)*loginfo.VBD_CNV;
	r_vbd = interp1d(sg_epoch_time,r_vbd,ctd_epoch_time,'linear');
	r_press = interp1d(sg_epoch_time,results.pressure,ctd_epoch_time,'linear'); % at avg_lat
	r_depth = ctd_depth_m;
	
	w = m2cm.*ctr1stdiffderiv(-r_depth,r_time);
	% determine length of data vectors for this dive
	mp = length(r_time);
	
	%define difference variables
	delta_t =  zeros(mp,1);
	delta_t(1:mp-1) = r_time(2:mp) - r_time(1:mp-1);
	delta_t(mp) = delta_t(mp-1);
	
	vbddiff = zeros(mp,1);
	vbddiff(2:mp)  = diff(r_vbd)./delta_t(2:mp);
	
	% should we test CTD_qc and hdm_qc as well and continue of either are bad?
	
	% Use raw or corrected temp/salin?
	% Currently corrected temperature and salinity
	% Use HDM speeds and corrected temp and salinities
	speed = results.speed;
	temperature = results.temperature;
	salinity = results.salinity;
	
	% select data to be used in minimization
	% where the glider is moving steadily (not flare or surface or apogee or stalled -- QC_BAD)
	% vbd isn't changing at too large a rate,
	% and roll isn't too large
	% keep points that are good and where there is steady flight (no accelerations)
	salinity_qc = str2num(results.salinity_qc);
	speed_qc = str2num(results.speed_qc);
	valid_i = find(speed_qc == QC_GOOD & salinity_qc == QC_GOOD & ...
				   abs(fix(vbddiff)) < vbddiffmax & ...
				   r_roll < rollmax & r_roll > rollmin);
	if (isempty(valid_i))
		fprintf(1,'No valid points in dive %d; skipping\n',ir);
		continue;
	end
	
	loop_pass = loop_pass+1;
	
	% concatenate variables from current dive cycle with those from previous dive  
	dive_climb = zeros(length(valid_i),1);
	dive_climb(find(valid_i >= dive_i(1)  & valid_i <= dive_i(end)))  =  double(ir);
	dive_climb(find(valid_i >= climb_i(1) & valid_i <= climb_i(end))) = -double(ir);
	dives_used = [dives_used; dive_climb];
	W = [W ; w(valid_i)];
	YD = [YD ; r_time(valid_i)];
	P = [P ; r_press(valid_i)];
	T = [T ; temperature(valid_i)];
	S = [S; salinity(valid_i)];
	V = [V ; r_vbd(valid_i)];
	PI = [PI ; r_pitch(valid_i)];
end


% rename arrays
w = W;
yday = YD;
press = P;
temp = T;
salin = S;
vbd = V;
pitch = PI;
clear W YD P T S V PI;

mp = length(w);
if (mp == 0)
	fprintf(1,'No data available!\n');
	return;
end

% ensure density in kg/m^3 if they gave it g/cc
if (rho0 < 2)
	rho0 = rho0*1000;
end

% ensure mass in kg if they gave it in grams
if (mass > 1000)
	mass = mass / 1000;
end

%calculate total volume 
%determine the portion of volume not associated with changes in vbd
%(the volume at neutral)
vol0 = volmax + (first_c_vbd-loginfo.VBD_MIN)*loginfo.VBD_CNV;

%potential density
density = sw_dens0(salin,temp);
density = real(density);
%in situ density (includes seawater compressibility)
density_insitu = sw_dens(salin,temp,press);
density_insitu = real(density_insitu);

%display initial values and misfit
fprintf(1,'Initial parameter values:\n');
fprintf(1,'$C_VBD,%4d\n', round(first_c_vbd));
fprintf(1,'Given vbdbias = %5.4fcc and volmax = %6.1fcc\n', vbdbias_0,volmax);
fprintf(1,'Hydrodynamic parameters: hd_a = %6.5e; hd_b = %6.5e;hd_c = %6.5e\n', hd_a, hd_b, hd_c);
fprintf(1,'Hull parameters: abs_compress = %5.4e; therm_expan  = %5.4e\n', abs_compress, therm_expan);
% fprintf(1,'temp_ref     = %5.4e; %% reference temperature\n', temp_ref);


% Our goal is to run all the regressions possible.
% first we always determine vbdbias
% if fit_abc is possible we do it, else assume abc are unchanged
% finally fit_exp_comp, if desired

% A note on the fields:
% w_rms_func only uses use_exp_comp
% w_misfit_func uses fit_abc and fit_exp_comp before calling w_rms_func
vars_need_fit = struct('fit_abc',0,'fit_exp_comp',0,'use_exp_comp',1,...
                       'vbd',vbd,'density',density,'density_insitu',density_insitu,'pitch',pitch,'yday',yday,...
                       'press',press,'w',w,'vol0',vol0,'mass',mass,'glider_length',glider_length,'rho0',rho0,...
                       'temp',temp,'hd_s',hd_s,'min_rel_q_inc',min_rel_q_inc,'max_num_iters',max_num_iters,...
                       'temp_ref',temp_ref,'mass_comp',mass_comp);

% do we do initial for each phase?
% Where do we start?
[w_rms_initial,ignore,ignore,ignore,ignore,w_sel,w_stdy_sel] = w_rms_func(vbdbias_0,hd_a,hd_b,hd_c,abs_compress,therm_expan,temp_ref,vars_need_fit);

%remember these for use in plotting later
w_init = w_sel;
w_stdy_init = w_stdy_sel;
fprintf(1,'Initial RMS diff. between model and actual w is %5.4f cm/s.\n', w_rms_initial);

display_opt = 'off'; % 'iter'
% intiialize regression results
vbdbias = 0;
a = hd_a; b = hd_b; c = hd_c;
abs_compress_opt = abs_compress; therm_expan_opt = therm_expan; temp_ref_opt = temp_ref;


% First regress vbdbias
fprintf(1,'Determining vbdbias...\n');
opt = optimset('maxfunevals',1000,'display',display_opt,'tolfun',1e-3);
vars_need_fit.fit_abc = 0;
vars_need_fit.fit_exp_comp = 0;
x_0 = vbdbias_0;
flight_consts = [hd_a hd_b hd_c];
vol_consts = [abs_compress therm_expan temp_ref];

% Find the value of vbdbias that minimizes w_rms_final within half the total throw of a DG (1200cc)
% clearly if we have dives thay flew we'll find something much smaller
% See OKMC Mar13 sg177 23:25
% There were few to no pumps on the climb and initial volamx
% was incorrect by ~80cc which lead it to compute very slow speeds
% original fminsearch happened to find a local minimum vbdbias at ~0cc but this was incorrect
HIST = []; % initialize regression history
% The @(x) w_misfit_func(x,v1,v2,v3) form is a function (lambda closure over v1,v2,v3) using var x, which is initialized at x_0
[x,w_rms_final,exitflag] = fminbnd(@(x) w_misfit_func(x,vars_need_fit,flight_consts,vol_consts),-600,600,opt);

% TODO if find_regress_dives sets bad_dives, report new volmax and request reprocess before continuing

% TODO if exitflag is bad, exit
vbdbias = x; % update our result...
fprintf(1,'Found vbdbias: %5.4fcc\n',vbdbias); %DEBUG

% TODO set fit_abc only if it appears we have a good spread of dives
if (fit_abc)
	fprintf(1,'Fitting hydrodynamic parameters...\n');
	%if fitting a.b.&c use more evaluations of misfit
	opt = optimset('maxfunevals',2000,'display',display_opt,'tolfun',1e-3);
	vars_need_fit.fit_abc = 1;
	vars_need_fit.fit_exp_comp = 0;
	% (didn't check if this requires all parameters to meet this tolerance)
	% TODO CCE recommends starting hd_b at 0.02 always since it regresses well
	% use our new vbdbais vbdbais to determine abc
	x_0 = [vbdbias hd_a 0.02 hd_c];
	%DEAD x_0 = [vbdbias hd_a hd_b hd_c]; % DEAD
	vol_consts = [abs_compress therm_expan temp_ref];
	HIST = []; % initialize regression history
	[x,w_rms_final,exitflag] = fminsearch(@(x) w_misfit_func(x,vars_need_fit,[],vol_consts),x_0,opt);
	vbdbias = x(1); a = x(2); b = x(3); c = x(4);
	%DEBUG fprintf(1,'Found vbdbias = %5.4fcc and hydrodynamic parameters: a = %6.5e; b = %6.5e;c = %6.5e\n', vbdbias, a, b, c);
else
	fprintf(1,'Skipping fit of hydrodynamic parameters due to insufficient dives...\n');
end

if (fit_exp_comp) 
	fprintf(1,'Determining compressibility and thermal expansion factors for hull and fairing...\n');
	opt = optimset('maxfunevals',2000,'display',display_opt,'tolfun',1e-3,'tolx',1e-8);
	vars_need_fit.fit_abc = 0;
	vars_need_fit.fit_exp_comp = 1;
	x_0 = [abs_compress therm_expan];
	% Use our best (hopefully newly found) parameters
	vbdbias_const = vbdbias;
	flight_consts = [a b c];
	vol_consts = [vbdbias_const temp_ref];
	HIST = []; % initialize regression history
	[x,w_rms_final,exitflag,fms_output] = fminsearch(@(x) w_misfit_func(x,vars_need_fit,flight_consts,vol_consts),x_0,opt);
	abs_compress_opt = x(1); therm_expan_opt = x(2); temp_ref_opt = temp_ref;
end
%% At this point whatever regression we ran leaves the following variables with 'best' values
% vbdbias (to apply to volmax)
% a,b,c
% abs_compress_opt, therm_expan_opt, temp_ref_opt (never changed)

%display final values and misfit
final_dives = unique(abs(dives_used));
dives = succinct_elts(final_dives);
titlestring = sprintf('SG%03d %s [%s]',eng.glider, underscore(mission_title),dives);
timelabel = datestr(clock);

% need some intermediates for plotting
[w_rms_fit,vbdc,buoy,hspd_stdy,w_stdy,w_sel,w_stdy_sel] = w_rms_func(vbdbias,a,b,c,abs_compress_opt,therm_expan_opt,temp_ref_opt,vars_need_fit);
% ensure these are the same
w_rms_fit = w_rms_final;
implied_volmax = (volmax - vbdbias);
%implied rho neutral with current C_VBD
rho_neutral_implied = mass/(1e-6*(vol0-vbdbias));

%implied c_vbd needed to have correct total vol 
implied_c_vbd_volmax = first_c_vbd+vbd_cnts_per_cc*vbdbias;

fprintf(1,'\nFinal RMS diff. between model and actual w is %5.4f cm/s.\n', w_rms_fit);
rms_dives = [];
for ir = final_dives'
	dive_i = find(abs(dives_used) == ir);
	rms_dive = rms(w_sel(dive_i) - w_stdy_sel(dive_i));
	%DEAD fprintf(1,'%d %5.4fcm/s\n',ir,rms_dive);
	rms_dives = [rms_dives; rms_dive];
end

if (w_rms_fit > good_w_rms_fit)
	fprintf(1,'NOTE: Regression did not yield a good set of parameters.\n');
	better_dives = final_dives(find(rms_dives <= good_w_rms_fit));
	if (length(better_dives) > 4)
		fprintf(1,'Try the following set of dives instead: %s\n',succinct_elts(better_dives));
	end
else
	%DEBUG fprintf(1,'Implied neutral density at current $C_VBD is %4.3f (g/cc)\n', rho_neutral_implied);
	%DEBUG fprintf(1,'Original volmax of %6.1fcc required a vbdbias of %5.4fcc\n',volmax, vbdbias);
	ratty_dives = final_dives(find(rms_dives > good_w_rms_fit+1.0));
	if (~isempty(ratty_dives))
		fprintf(1,'NOTE: While the overall RMS was good, the following dives were not terribly good: %s.\n',succinct_elts(ratty_dives));
	end
	fprintf(1,'Recommended $C_VBD,%4d\n', fix(implied_c_vbd_volmax));
	
	fprintf(1,'Add the following lines to the bottom of sg_calib_constants.m:\n');
	if (vbdbias_0) % they set it in sg_calib_constants.m; should be removed
		fprintf(1,'NOTE: Ensure vbdbias setting is removed or is set to 0!\n');
	end
	fprintf(1,'%% %s RMS=%5.4f cm/s Dives: %s\n',timelabel,w_rms_fit,dives);
end

fprintf(1,'volmax = %6.1f;\n',implied_volmax);
% TODO if vbdbias_0 is non-zero emit warning to remove vbdbias
if (fit_abc == 1)
	fprintf(1,'hd_a = %6.5e;\n', a);
	fprintf(1,'hd_b = %6.5e;\n', b);
	fprintf(1,'hd_c = %6.5e;\n', c);
end
if (fit_exp_comp)
	fprintf(1,'abs_compress = %5.4e;\n', abs_compress_opt);
	fprintf(1,'therm_expan  = %5.4e;\n', therm_expan_opt);
end

% BUG: do we actually ever change this?
% DEBUG fprintf(1,'temp_ref     = %5.4e;\n', temp_ref_opt);

ref_spd = 0.25; % PARAMETER m/s reference speed
av_density = rho0;
ref_q = 0.5*av_density*ref_spd*ref_spd;
c2aa = c/(a*a);
performfactor = 1/(b*c2aa*( ref_q^(hd_s) ) );

initial_abc_text = sprintf('a = %.4g; b = %.4g; c = %.4g',hd_a,hd_b,hd_c);
final_abc_text   = sprintf('a = %.4g; b = %.4g; c = %.4g',a,b,c);
abc_lambda_text  = sprintf('%s; \\Lambda_0_._2_5 = %.4f',final_abc_text,performfactor);
pitch_bias_text = sprintf('Pitch bias = %.4g\\circ; VBD bias = %.4gcc; RMS w-w_m_o_d_e_l = %.4g cm/s',pitchbias,vbdbias,w_rms_fit);

% return;

% TODO first 3 plots only if abc?
%-------------------------------------------------------------
% scatter plot of w with w_stdy for w_stdy computed with both initial 
% and final vbdbias values

figure;
ymax = 25; % cm/s
xmax = ymax;
xmin = -xmax;
ymin = -ymax;
plot([xmin, xmax], [ymin, ymax], '-k');
hold on;
grid on;
axis([xmin xmax ymin ymax]);
set(gca, 'Xtick', [xmin:5:xmax], 'Ytick', [ymin:5:ymax]);
axis square;
plot(w_stdy_init, w_init, '.b');
plot(w_stdy_sel, w_sel, '.r');
xlabel('Steady Model Vertical Velocity [cm/s]');
ylabel('Vertical Velocity from Pressure [cm/s]');
title(titlestring);
% num2str(x) formats using, effectively, %-11.4g
xtxt = -10;
text(xtxt, -15, sprintf('VBD bias = %.4gcc RMS w-w_m_o_d_e_l = %.4gcm/s\n%s', ...
						vbdbias_0, w_rms_initial, initial_abc_text), 'Color', 'b');
text(xtxt, -20, sprintf('VBD bias = %.4gcc RMS w-w_m_o_d_e_l = %.4gcm/s\n%s', ...
						vbdbias,   w_rms_fit,     final_abc_text),    'Color', 'r');
mass_comp_str = '';
if (mass_comp)
	mass_comp_str = sprintf('mass_c_o_m_p = %.4gg',1000*mass_comp);
end
text(xtxt,  20, sprintf('Assumed mass = %6.1fg%s\nAssumed volmax = %.4g liters\nImplied volmax = %.4g liters\nImplied neutral density @ current C_V_B_D = %.4g kg/m^3\nInitial C_V_B_D = %.4g A/D counts',...
						1000*mass,mass_comp_str,0.001*volmax,0.001*(volmax - vbdbias),rho_neutral_implied,first_c_vbd),...
	 'Color', 'k');	 

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

% Charlie's additional plots for debugging and understanding the quality of the data and results

%
% Plot w - w_stdy in buoyancy-pitch space
%
figure;
xmax = 60; % PARAMETER pitch (degrees)
ymax = 400; % PARAMETER buoyancy (cc)
xmin = -xmax;
ymin = -ymax;

[P, B] = meshgrid(0:1:xmax, 2:2:ymax);
Time = B; % doesn't matter as long as same 2D size 
Umag0 = 10*ones(size(P));
[Umag, Thdeg, ignore] = flightvec2(B, P, Time, glider_length, min_rel_q_inc, max_num_iters, [a, b, c], rho0);
Wgrd = Umag.*sin(Thdeg*pi/180);
w_levels = [2:2:20];
[C, h] = contour(P, B, Wgrd, w_levels, '-b');
set(gca, 'Xlim', [-xmax,xmax], 'XTick', [-xmax:10:xmax], ...
         'Ylim', [-ymax,ymax], 'Ytick', [-ymax:100:ymax]);
clabel(C, h, w_levels);
hold on
grid on
[P, B] = meshgrid(-xmax:1:0, -ymax:2:-2);
Umag0 = 10*ones(size(P));
[Umag, Thdeg, ignore] = flightvec2(B, P, Time, glider_length, min_rel_q_inc, max_num_iters, [a, b, c], rho0);
Wgrd = Umag.*sin(Thdeg*pi/180);
w_levels = [-20:2:-2];
[C, h] = contour(P, B, Wgrd, w_levels, '-b');
clabel(C, h, w_levels);

wmax = 3; % PARAMETER cm/s difference
wmin = -wmax;
wdel = 1;
ncolors = (wmax-wmin)/wdel;
colormap(prism(ncolors)); % show color bar for +/-3cm/s w differences
h = c_line(pitch, buoy, w - w_stdy, 'o');
set(h, 'MarkerSize', 1);
caxis([wmin wmax]);
colorbar;
legend('W','Location','Southeast')
title(titlestring);
xlabel('Corrected Pitch [\circ]');
ylabel('Corrected Buoyancy [g]');

text(-55, 300, sprintf('%s\nReference C_V_B_D = %4d A/D counts\nAssumed volmax = %.4g liters',...
								pitch_bias_text, first_c_vbd, 0.001*volmax), 'Color', 'k');

text( 0, -100, sprintf('%s\nImplied volmax = %.4g liters\n', ...
						final_abc_text,0.001*(volmax - vbdbias)), 'Color', 'k');

bar_label = 'w_p_r_e_s_s_u_r_e - w_s_t_e_a_d_y _m_o_d_e_l [cm/s]';
text(1.3*xmax, -0.5*ymax, bar_label, 'Rotation', 90);

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);


% set arrays for model performance plot


if (deepglider)
	umax = 0.8; % PARAMETER m/s speed (DG has bigger VBD throw)
else
	umax = 0.4; % PARAMETER m/s speed
end
uu = ( umax*0.01:umax*0.01:umax )';
ww = ( -umax/2:umax*0.01:umax/2 )';
[UU, WW] = meshgrid(uu, ww);
qmag = 0.5*av_density*( UU.*UU + WW.*WW );
q_s_inv = qmag.^(-hd_s);
t2 = (WW.*WW)./(UU.*UU);
th = atan2(WW, UU);

% calculate model power, buoyancy force, and attack angle for non-zero velocity components

pwr = zeros(length(ww), length(uu));
force = pwr;
attack = pwr;
param = pwr;

non0 = find(WW~=0);
param(non0) = sqrt( 1. - ( 4*b*c/(a*a) )*( 1./(t2(non0).*q_s_inv(non0)) ) );
param(non0) = (1 - real(param(non0)) ).*( 1 - ... 
                                          (2/pi)*atan2( imag(param(non0)), real(param(non0)) ) ); 
force(non0) = glider_length*glider_length*a*a.*qmag(non0).*tan(th(non0)).*param(non0)./ ... 
    (2.*c*cos(th(non0)));
pwr(non0) = force(non0).*WW(non0);
force = 1000*force/gravity;
attack = -a*tan(th).*param*0.5/c;

qmag0 = 0.5*av_density*ww.*ww;
force0 = ( qmag0.^(1 + hd_s) )*glider_length*glider_length*b;
pwr0 = ww.*force0;
force0 = 1000*force0/gravity;
attack0 = zeros(length(ww),1);

% calculate model power, buoyancy force, and attack angle for purely vertical motion

force1 = [force0, force];
pwr1 = [pwr0, pwr];
attack1 = [attack0, attack];
uu = [0; uu];

% add row of zeros corresponding to vanishing vertical motion

clear force pwr attack;
force = force1;
pwr = pwr1;
attack = attack1;

% the calls to contour() are fast while those to clabel() are *very* slow.
figure;
set (gca, 'DefaulttextFontSize', 10)
set (gca, 'DefaulttextColor', 'k')
title(titlestring);
hold on;

colormap(prism(ncolors));
h = c_line(0.01*hspd_stdy, 0.01*w_stdy, w - w_stdy, 'o');
set(h, 'MarkerSize', 1);
caxis([wmin wmax]);
colorbar;

h_con = zeros(3,1); % we want to skip one of the contours below in the legend

buoy_levels = [-700:20:700]; % PARAMETER cc
[C, h] = contour(uu, ww, force, buoy_levels, '--r' );
h_con(1) = h;
set (gca, 'DefaulttextFontSize', 7)
set (gca, 'DefaulttextColor', 'r')
clabel(C, h, buoy_levels);

grid on;
hold on;
axis square;
set (gca, 'Xlim', [0, umax], 'Xtick', [0:0.05:umax], ...
          'Ylim', [-0.5*umax,0.5*umax], 'Ytick', [-0.5*umax:0.05:0.5*umax]);
xlabel( 'Horizontal Speed (m/s)');
ylabel( 'Vertical Speed (m/s)');

pwr_levels = [0.01:0.01:0.1];
[C, h] = contour(uu, ww, pwr, pwr_levels, '-k');
h_con(2) = h;
set (gca, 'DefaulttextFontSize', 7)
set (gca, 'DefaulttextColor', 'k')
clabel(C, h, pwr_levels);

pwr_levels = [0.2:0.1:1];
[C, h] = contour(uu, ww, pwr, pwr_levels, '-k'); % skip this contour in the legend (same as above)
clabel(C, h, pwr_levels);

attack_levels = [-10:0.5:10];
[C, h] = contour(uu, ww, attack, attack_levels, ':b');
h_con(3) = h;
set (gca, 'DefaulttextColor', 'b')
clabel(C, h, attack_levels);

set (gca, 'DefaulttextFontSize', 10);
set (gca, 'DefaulttextColor', 'k');

text(0.35*umax, -0.62*umax, abc_lambda_text);
text(0.1,0, pitch_bias_text, 'Color', 'k');

bar_label = 'w_p_r_e_s_s_u_r_e - w_s_t_e_a_d_y _m_o_d_e_l [cm/s]';
text(1.25*umax, -0.25*umax, bar_label, 'Rotation', 90);

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

legend(h_con, 'Buoyancy (g)', 'Power (W)', 'Attack angle (\circ)', 3);

figure

jj = 1:length(w);
plot(jj, w)
hold on
grid on
plot(jj, w_stdy, '-c')
plot(jj, vbdc/10, '-k')
plot(jj, buoy/10, '-m')
plot(jj, pitch, '-g')
plot(jj, w-w_stdy, '-r')
legend('w', 'w_s_t_d_y', 'vbdc/10', 'buoy/10', ...
       'pitch', 'w -w_s_t_d_y', 0)
title(titlestring)
xlabel('Data record number')

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

xtxt = v(1) + 0.3*(v(2)-v(1));
ytxt = v(3) - 0.1*(v(4)-v(3));
text(xtxt, ytxt, abc_lambda_text);

xtxt = v(1) + 0.05*(v(2)-v(1));
ytxt = v(3) + 0.95*(v(4)-v(3));
text(xtxt, ytxt, pitch_bias_text, 'Color', 'k');

