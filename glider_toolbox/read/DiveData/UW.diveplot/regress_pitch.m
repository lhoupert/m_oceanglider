% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

%  regress_pitch.m

%  Read sg data files for a specified set of dive cycles and 
%	find the pitch control bias, as well as pitch_vbd_shift
%
%  functions called:
%    sg_calib_constants.m - establishes glider calibration constants
%    various seawater routines
%    get_dive_data.m - data and results from nc or nc+mat
%


full_pitch_models = 0; % use experimental pitch models
display_fit = 1;
depthmin = 15; % PARAMTER [m] at least 15 down (below flare)
dwdtmax = 1; % PARAMETER [cm/s^2] vertical acceleration of vehicle is < 1/cm/s^2
vbddiffmin = 0.5; % PARAMATER [cc] vbd isn't moving more than 0.5cc (avoid pumps)
rollmin = 55; % PARAMETER [deg] we are not rolled over more than 55deg

% Display bounds
displ_pitch_max = 60; % PARAMETER max positive measured pitch (degrees)
displ_pitch_min = -60; % PARAMETER min positive measured pitch (degrees)
displ_pitch_tick = 10;
% this is the linear range
linear_pitch_control_max = 2.5; % PARAMETER max positive pitch control position (cm)
linear_pitch_control_min = -2.5; % PARAMETER max megative pitch control position (cm)
pitch_control_tick = 0.5;

timelabel = datestr(clock);

%sg_calib_constants; % get id_str % billr: this blows up the compiled code
if (exist('regress_dives','var'))
	dives = regress_dives;
	clear regress_dives;
else
	dives = ask_which_runs(available_profiles());
end


%define constants
radpd = pi/180.;
cm2m = 0.01;
m2cm = 100;
mS_per_cm_2_S_per_m = 0.1;
kg2g = 1000;

c_vbd_ref = NaN;
av_apogee_density = 0;
dives_used = [];
dive_pitch_desired = [];
V = [];
PT = [];
PI = [];
PC = [];
RC = [];
BU = [];
VC = [];
TR = [];
GPC = [];
GPI = [];
MOD = [];
%
% loop on dive cycles
%
id_str = settings.seagliderID; %billr: this replaces the sg_calib_constants line above
for ir = dives'
  
	id_dive = sprintf('%03d%04d',id_str,ir);
	[loginfo,eng,results] = get_dive_data(id_dive,1); % need results: gsm speed, density, etc.
	if (isfield(results,'processing_error') || ...
		isfield(results,'skipped_profile')  || ...
		~isfield(loginfo, 'MHEAD_RNG_PITCHd_Wd'))
		fprintf(1,'Unable to include data for dive %d\n', ir);

		continue;
	end
	unpack_data;
	fprintf(1, 'p%s\n', id_dive)
	% DEBUG fprintf(1,'Pitch intended: %f\n',abs(loginfo.MHEAD_RNG_PITCHd_Wd(3))); % DEBUG
	
	c_pitch = loginfo.C_PITCH;
	c_roll_dive = loginfo.C_ROLL_DIVE;
	c_roll_climb = loginfo.C_ROLL_CLIMB;
	c_vbd = loginfo.C_VBD;
	if (isnan(c_vbd_ref))
		% we use the last vbdbias
		c_vbd_ref = c_vbd;
		c_pitch_ref = c_pitch;
		c_roll_dive_ref = c_roll_dive;
		c_roll_climb_ref = c_roll_climb;
	end

	if (c_vbd ~= c_vbd_ref)
		fprintf(1,'Centers changed; unable to include data for dive %d\n', ir);
		continue;
	end
	if 0
	if (c_vbd ~= c_vbd_ref | ...
		c_pitch ~= c_pitch_ref | ...
		c_roll_dive ~= c_roll_dive_ref | ...
		c_roll_climb ~= c_roll_climb_ref)
		fprintf(1,'Centers changed; unable to include data for dive %d\n', ir);
		continue;
	end
	end
	
	if (isfield(loginfo,'PITCH_DBAND'))
		pitch_dband = loginfo.PITCH_DBAND;
	else
		pitch_dband = 0.1; % (cm) default value for old glider missions
	end
	pitch_ovshoot = loginfo.P_OVSHOOT;
	pitch_dband = max(pitch_dband,pitch_ovshoot);
	
	pitch_min_cnts = loginfo.PITCH_MIN;
	pitch_max_cnts = loginfo.PITCH_MAX;
	
	roll_min_cnts = loginfo.ROLL_MIN;
	roll_max_cnts = loginfo.ROLL_MAX;
	
	vbd_min_cnts = loginfo.VBD_MIN;
	%DEAD vbd_max_cnts = loginfo.VBD_MAX; % UNUSED
	
	vbd_cnts_per_cc   = 1./loginfo.VBD_CNV;% same for all dives
	roll_cnts_per_deg = 1./loginfo.ROLL_CNV; % 35.37 AD counts per degree
	pitch_cnts_per_cm = 1./loginfo.PITCH_CNV; % 217.39 AD counts per cm
	
	gc_start = loginfo.gc_st_secs - loginfo.log_epoch_start;
	gc_end = loginfo.gc_end_secs - loginfo.log_epoch_start;
	gc_vbd_ctl = loginfo.gc_vbd_ctl;
	gc_vbd_secs = loginfo.gc_vbd_secs;
	gc_pitch_ctl = loginfo.gc_pitch_ctl; % pitch intended values
	apogee_control_index = find(gc_vbd_ctl == 0); % where VBD moved to neutral
	
	% What, nominally, the glider should have done for pitch vs what she
	% actually did in the GCs above
	ngc = length(gc_pitch_ctl); % this is number of GCs
	
	gc_pitch_intended = zeros(1,ngc); % the pitch desired at the start of each gc
	pitch_desired = loginfo.MHEAD_RNG_PITCHd_Wd(3);
	dive_pitch_desired = [dive_pitch_desired; pitch_desired]; % a log of the pitch_desired
	gc_pitch_intended(1:apogee_control_index-1) = pitch_desired;
	gc_pitch_intended(apogee_control_index) = loginfo.APOGEE_PITCH; 
	gc_pitch_intended(apogee_control_index+1:ngc) = -pitch_desired;
	
	pitch_intended = zeros(sg_np,1); % the pitch desired at each data point
	model = zeros(sg_np,1); % which model to fit

	% find ctl_i for each GC, the points from the start of GC(gc) to the start of GC(gc+1)
	% compute properties of that GC stretch
	gc_end_i = 1; % data point where GC(gcp) started
	for gc = 2:ngc
		gc_start_i = find(time >= gc_start(gc),1,'first'); % at what datapoint did GG(gc) start?
		ctl_i = gc_end_i:gc_start_i-1; % these were the points responding to GC(gc-1)
		gcp = gc - 1; % the GC that controls these points
		pitch_intended(ctl_i) = gc_pitch_intended(gcp); % the pitch_intended

		% which model should this data be used to fit?
		if (gcp < apogee_control_index)
			model(ctl_i) = 1; % dive
		elseif (gcp == apogee_control_index)
			model(ctl_i) = 2; % apogee model
		else
			model(ctl_i) = 3; % climb		  
		end
		gc_end_i = gc_start_i;
	end
	ctl_i = gc_end_i:sg_np;
	pitch_intended(ctl_i) = gc_pitch_intended(gc);
	model(ctl_i) = 3; % climb		  
	
	%
	% read sg data
	%
	vmtime = eng.elaps_t;
	mp = length(vmtime);
	% depth in meters, as computed onboard by SG
	% using pressure gauge and linear conversion
	vmdepth  = eng.depth;      
	hdg = eng.head;
	pitch = eng.pitchAng;
	roll = eng.rollAng;
	pitch_control = eng.pitchCtl;
	roll_control = eng.rollCtl;
	vbd = eng.vbdCC;
	
	
	% compute heading difference and turn rate
	hdgdiff = zeros(mp,1);
	hdgdiff(2:mp) = diff(hdg);
	hdgdiff = mod(hdgdiff,360);
	in_hdg = find(hdgdiff > 180);
	hdgdiff(in_hdg) = hdgdiff(in_hdg) - 360;
	hdg_wrapped = hdg(1) + cumsum(hdgdiff);
	turn_rate = ctr1stdiffderiv(hdg_wrapped, vmtime);
	
	%
	% evaluate track speeds in cm/s (positions are in m) 
	%
	delta_t =  zeros(mp,1);
	delta_t(1:mp-1) = vmtime(2:mp) - vmtime(1:mp-1);
	delta_t(mp) = delta_t(mp-1);
	% where did we pump (and bleed)?
	vbddiff = zeros(mp,1);
	vbddiff(2:mp)  = diff(vbd)./delta_t(2:mp);
	
	zpos = -vmdepth;
	w  = 100.*ctr1stdiffderiv(zpos, vmtime);
	
	% We need density at apogee for pitch_vbd_shift calc
	% density_raw from salinity_raw, temperature_raw
	r_salin_raw = interp1d(ctd_epoch_time,results.salinity_raw,    sg_epoch_time,'linear');
	r_temp_raw  = interp1d(ctd_epoch_time,results.temperature_raw, sg_epoch_time,'linear');
	r_press     = results.pressure; % already at SG points
	r_density_raw_insitu = sw_dens(r_salin_raw, r_temp_raw,r_press);
	% BETTER: find deepest point in z and use that
	% where did we pump first (to neutral)? hence assumed deepest location
	ivbdp = find( vbddiff > vbddiffmin );
	% what was potential density there?  based on raw salin/temp
	density0 = r_density_raw_insitu(ivbdp(1));
	
	vbdc = vbd - vbdbias;
	vol0 = volmax + (c_vbd - vbd_min_cnts)/vbd_cnts_per_cc;
	vol = (vol0 + vbdc).*exp( -abs_compress*r_press + therm_expan*(r_temp_raw - temp_ref) );
	vol_comp = zeros(sg_np,1);
	if mass_comp
		kg2g = 1000;
		vol_comp_ref = kg2g*mass_comp/cml_dens(temp_ref, 0);
		vol_comp = kg2g*mass_comp./cml_dens(r_temp_raw, r_press);
		vol_comp = vol_comp - vol_comp_ref;
	end
	vol = vol + vol_comp;
	
	% TODO interpolate valid results.density to sg_epoch_time from
	% ctd_epoch_time being careful of NaNs
	valid_i = find(~isnan(results.density));
	r_density_insitu = interp1d(ctd_epoch_time(valid_i), results.density(valid_i), sg_epoch_time,'linear');
	r_density_insitu = r_density_raw_insitu; % what the glider can compute
	buoy = kg2g*(-mass+r_density_insitu.*vol*(cm2m)^3);
	
	spdg = interp1d(ctd_epoch_time,results.speed_gsm,sg_epoch_time,'linear');
	
	
	% Compute iglide_stdy, where are we actually flying well?  do the regressions from that data
	% where are we negative and not stalled (dive?) uses results.speed_gsm
	ivbd = find( vbdc < 0 & spdg > 0);
	% where are we moving at all
	ispp = find( spdg > 0);
	% start of dive to end of climb where we aren't stalled
	istdy = ( ivbd(1):ispp(length(ispp)) )';
	
	%
	% select data to be used in minimization 
	%
	
	dwdt = ctr1stdiffderiv(w,vmtime); % were are we accelerating vertically
	
	% where are we gliding post flare?  no acceleration, turning or pumping
	iglide0 = find( abs(dwdt) < dwdtmax & ...
					abs(vbddiff) < vbddiffmin & ...
					abs(roll) < rollmin & vmdepth > depthmin);
	% include points after start of dive 
	iglideindex = find( iglide0 >= istdy(1) );
	iglide = iglide0(iglideindex);
	
	%	[istdy(1), istdy(length(istdy)), length(istdy)]
	%	[iglide0(1), iglide0(length(iglide0)), length(iglide0)]
	%	[iglide(1), iglide(length(iglide)), length(iglide)]
	
	%	iglide_stdy = iglide - istdy(1) + 1;
	
	in = zeros(mp,1);
	in_steady = in;
	in_steady(istdy) = 1;
	in_glide = in;
	in_glide(iglide) = 1;
	% where are we both steady and gliding?
	iglide_stdy = find(in_steady > 0 & in_glide > 0);
	
	iglide_stdy = 1:sg_np; % limit below
	
	%
	% adjust pitch and roll control to be relative to common reference values
	%
	pitch_control = pitch_control + (c_pitch - c_pitch_ref)/pitch_cnts_per_cm;
	iwn = find( pitch_control < 0 );
	iwp = find( pitch_control > 0 );
	roll_control(iwn) = roll_control(iwn) + (c_roll_dive  - c_roll_dive_ref)/roll_cnts_per_deg;
	roll_control(iwp) = roll_control(iwp) + (c_roll_climb - c_roll_climb_ref)/roll_cnts_per_deg;
	
	roll_control_counts = zeros(mp,1);
	roll_control_counts(iwn) = c_roll_dive  + roll_control(iwn)*roll_cnts_per_deg;
	roll_control_counts(iwp) = c_roll_climb + roll_control(iwp)*roll_cnts_per_deg;
	
	dives_used = [dives_used; ir];
	
	V = [V ; vbd(iglide_stdy)];
	PT = [PT ; pitch(iglide_stdy)];
	PI = [PI ; pitch_intended(iglide_stdy)];
	PC = [PC ; pitch_control(iglide_stdy)];
	RC = [RC ; roll_control(iglide_stdy)];
	BU = [BU ; buoy(iglide_stdy)];
	VC = [VC ; vol_comp(iglide_stdy)];
	TR = [TR; turn_rate(iglide_stdy)];
	GPC = [GPC ; gc_pitch_ctl];
	GPI = [GPI ; gc_pitch_intended'];
	MOD = [MOD ; model(iglide_stdy)];
	av_apogee_density = av_apogee_density + density0;
end; % end dive cycle loop

%
% rename arrays
%
vbd = V;
pitch = PT;
pitch_intended = PI;
pitch_control = PC;
roll_control = RC;
buoy = BU;
vol_comp = VC;
turn_rate = TR;
gc_pitch_ctl = GPC;
gc_pitch_intended = GPI;
model = MOD;

n_dives_used = length(dives_used);
av_apogee_density = av_apogee_density/n_dives_used;
% are we likely degenerate in pitch_desired?
% this is an approximation.  in fact, if there is a sequence of dives near one another 
% in pitch_desired, then the sorted diff might show small differences but the difference
% from the first to the last is greater than our threshold.  so we need to find
% the stretches between the big breaks and see how much change there is, divided by
% the thresold.
pitches_desired = unique(abs(dive_pitch_desired));
n_pitch_desired = length(find(diff(pitches_desired) > 1)); % use 1 degree
n_pitch_desired = n_pitch_desired + 1; % how many actual different pitch_desired there are, unless dive_pitch_desired is empty
min_pitch_desired = min(pitches_desired); % PARAMETER range of usable pitch_intended
max_pitch_desired = max(pitches_desired);

mp = length(vbd);
if (mp == 0)
  fprintf(1,'No data available!\n');
  return;
end
dives_used = succinct_elts(dives_used);
titlestring = sprintf('Seaglider %s %s [%s]', id_str, underscore(mission_title), dives_used);

iwn = find( pitch_control < 0 ); % dive points
iwp = find( pitch_control > 0 ); % climb points
z_pc = zeros(size(pitch_control));

%
% find linear fit between pitch and pitch_control
%
if (full_pitch_models & display_fit)
	full_legend_text = {};
	figure
	hold on
	grid on
	title(titlestring)
	plot(pitch,'.b'); full_legend_text{end+1} = 'Observed';
end
if (full_pitch_models)
	% regress pitch against pitch control, buoy and vbd oil shift and roll control and pitch_intended if possible
	% NOTE: if the set of dives has nearly the same pitch_intended then the solutions below
	% tend to be rank deficient.
	
	% dive model
	% quite often you will get stray points from subsequent GCs because the timing vagaries of collecting data points
	% so sometimes during 'dives' you'll see the initial points of pumps for apogee
	ip_dive = find(model == 1 & pitch_control > linear_pitch_control_min & vbd < 0);
	AA = pitch_regress(ip_dive,pitch,pitch_control,buoy,vbd,roll_control,pitch_intended,vol_comp);
	AA_dive = AA;
	% models to use if the CT is temporarily plugged...ignore: just fall back on std model of $PITCH_GAIN and $PITCH_VBD_SHIFT
	%DEBUG AA_no_trust_C  = pitch_regress(ip_dive,pitch,pitch_control,z_pc,vbd,roll_control,pitch_intended,vol_comp);
	%DEBUG AA_no_trust_CT = pitch_regress(ip_dive,pitch,pitch_control,z_pc,vbd,roll_control,pitch_intended,z_pc); % really only needed if mass_comp > 0
	ip_dive_x = [ones(size(ip_dive)) pitch_control(ip_dive) buoy(ip_dive) vbd(ip_dive) roll_control(ip_dive) pitch_intended(ip_dive) vol_comp(ip_dive)];
	pitch_d = ip_dive_x*AA_dive;
	rmse_d = rms(pitch(ip_dive) - pitch_d);
	if (display_fit)
		plot(ip_dive, pitch_d,'vg'); full_legend_text{end+1} = 'Dive model';
	end
	
	% apogee model
	ip_apo = find(model == 2);
	AA = pitch_regress(ip_apo,pitch,pitch_control,buoy,vbd,roll_control,pitch_intended,vol_comp);
	AA_apogee = AA;
	ip_apo_x = [ones(size(ip_apo)) pitch_control(ip_apo) buoy(ip_apo) vbd(ip_apo) roll_control(ip_apo) pitch_intended(ip_apo) vol_comp(ip_apo)];
	pitch_a = ip_apo_x*AA_apogee;
	rmse_a = rms(pitch(ip_apo) - pitch_a);
	if (display_fit)
		plot(ip_apo, pitch_a,'dm'); full_legend_text{end+1} = 'Apogee model';
	end
	
	% climb model
	ip_climb = find(model == 3 & pitch_control < linear_pitch_control_max & vbd > 0);
	AA = pitch_regress(ip_climb,pitch,pitch_control,buoy,vbd,roll_control,pitch_intended,vol_comp);
	AA_climb = AA;
	ip_climb_x = [ones(size(ip_climb)) pitch_control(ip_climb) buoy(ip_climb) vbd(ip_climb) roll_control(ip_climb) pitch_intended(ip_climb) vol_comp(ip_climb)];
	pitch_c = ip_climb_x*AA_climb;
	rmse_c = rms(pitch(ip_climb) - pitch_c);
	if (display_fit)
		plot(ip_climb, pitch_c,'^g'); full_legend_text{end+1} = 'Climb model';
	end
	if (0) % display
		composite_pitch_gains = [AA_dive(2) AA_apogee(2) AA_climb(2)];
		composite_pitch_gains = [AA_dive(2)  AA_climb(2)] % just dive and climb
		mean_pitch_gains = mean(composite_pitch_gains)
	end
end

ip = find(pitch_control > linear_pitch_control_min & pitch_control < linear_pitch_control_max);
% The original model
X = [ones(size(pitch_control(ip))) pitch_control(ip)];
% regress pitch against pitch control
A = X\pitch(ip);
rmse_o = rms(pitch(ip)- X*A);

% Jason's model
% AA = pitch_regress(ip,pitch,pitch_control,buoy,vbd,roll_control,pitch_intended,vol_comp);
AA = pitch_regress(ip,pitch,pitch_control,buoy,vbd,z_pc,z_pc,z_pc);
% JG just used pitch_control buoy vbd, not roll_control or pitch_intended, let alone by dive phases
XX_jg = [ones(size(pitch_control(ip))) pitch_control(ip) buoy(ip) vbd(ip)];
AA_jg = XX_jg\pitch(ip);
rmse_jg = rms(pitch(ip) - XX_jg*AA_jg);
if (full_pitch_models & display_fit)
	plot(ip, XX_jg*AA_jg,'xk'); full_legend_text{end+1} = 'JG model';
end

if (full_pitch_models)
	% A single composite model -- dive/apogee/climb in one model as our replacement for JG's
	AA = pitch_regress(ip,pitch,pitch_control,buoy,vbd,roll_control,pitch_intended,vol_comp);
	ip_x = [ones(size(ip)) pitch_control(ip) buoy(ip) vbd(ip) roll_control(ip) pitch_intended(ip) vol_comp(ip)];
	rmse_composite = rms(pitch(ip) - ip_x*AA)
	if (display_fit)
		plot(ip, ip_x*AA,'oc'); full_legend_text{end+1} = 'Composite model';
	end
end

if (full_pitch_models & display_fit)
	% end of full_pitch data figure...
	plot(pitch_intended,'.r'); full_legend_text{end+1} = 'Intended';
	xlabel('Data point')
	ylabel('Pitch Degrees');
	lg = legend(full_legend_text,'Location','SouthWest');
	set(lg,'color','none'); % make transparent
	set(gca,'YLim',[-60 40]) % degrees
	text(length(pitch)/5,-50,sprintf('RMSE_d_i_v_e=%.3f RMSE_c_l_i_m_b=%.3f RMSE_o_r_i_g=%.3f RMSE_j_g=%.3f',rmse_d,rmse_c,rmse_o,rmse_jg));
end

c_pitch_imp = c_pitch_ref - ( AA(1)/AA(2) )*pitch_cnts_per_cm;
pitch_vbd_shift_imp =-(AA(4) + AA(3)*av_apogee_density/kg2g)/AA(2);
% display pitch_gain = AA(2)

% what pitch_control the vehicle is capable of assuming different c_pitches
pitch_control_max = (pitch_max_cnts - c_pitch_ref)/pitch_cnts_per_cm;
pitch_control_min = (pitch_min_cnts - c_pitch_ref)/pitch_cnts_per_cm;
pitch_control_max_imp = (pitch_max_cnts - c_pitch_imp)/pitch_cnts_per_cm;
pitch_control_min_imp = (pitch_min_cnts - c_pitch_imp)/pitch_cnts_per_cm;

%
% plot scatter of pitch against pitch_control
%

figure;
hold on;
grid on;
legend_text = {};

plot(pitch_control, pitch, '.b'); legend_text{end+1} = 'Observed';
plot(pitch_control(ip), X*A, 'oc'); legend_text{end+1} = 'Predicted (pitch-only)'; % original
if (full_pitch_models)
	if (1)
		plot(pitch_control(ip_dive), ip_dive_x*AA_dive,   'vg'); legend_text{end+1} = 'Predicted dive';
		plot(pitch_control(ip_apo),  ip_apo_x*AA_apogee,  'dm'); legend_text{end+1} = 'Predicted apogee';
		plot(pitch_control(ip_climb),ip_climb_x*AA_climb, '^g'); legend_text{end+1} = 'Predicted climb';
	else
		plot(pitch_control(ip), ip_x*AA,'og'); legend_text{end+1} = 'Predicted (composite)';
	end
end
plot(pitch_control(ip), XX_jg*AA_jg, 'xk'); legend_text{end+1} = 'Predicted (JG)';
plot(gc_pitch_ctl, gc_pitch_intended, '.r'); legend_text{end+1} = 'Intended';
plot([0], [0], '+r', 'MarkerSize', 10) % the zero point to help see how off center the prediction is
axis([linear_pitch_control_min linear_pitch_control_max displ_pitch_min displ_pitch_max]);
set(gca, 'XTick', [linear_pitch_control_min:pitch_control_tick:linear_pitch_control_max], 'YTick', [displ_pitch_min:displ_pitch_tick:displ_pitch_max]);
xlabel('Pitch Control [cm]');
ylabel('Pitch [\circ]');
title(titlestring, 'FontSize', 12);
lg = legend(legend_text, 'Location','SouthEast');
set(lg,'color','none'); % make transparent

ix = [linear_pitch_control_min:0.01*(linear_pitch_control_max-linear_pitch_control_min):linear_pitch_control_max]; % show predictions within PC linear response range
plot(ix, A(1) + A(2)*ix, '--k'); % pitch vs. pitch_control only
plot(ix, AA(1) + AA(2)*ix, '-k'); % pitch vs. pitch_control neglecting other terms

xtxt = linear_pitch_control_min + 0.1*(linear_pitch_control_max-linear_pitch_control_min);
ytxt = displ_pitch_min + 0.95*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Reference C_P_I_T_C_H =', num2str(c_pitch_ref) );
text(xtxt, ytxt, comment);
ytxt = displ_pitch_min + 0.9*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Pitch control range =', ...
                  num2str(pitch_control_min), ' to ', ...
                  num2str(pitch_control_max), ' cm'); 
text(xtxt, ytxt, comment);
ytxt = displ_pitch_min + 0.85*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Implied C_P_I_T_C_H =', num2str(c_pitch_imp) );
text(xtxt, ytxt, comment);
ytxt = displ_pitch_min + 0.8*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Implied pitch control range =', ...
                  num2str(pitch_control_min_imp), ' to ', ...
                  num2str(pitch_control_max_imp), ' cm'); 
text(xtxt, ytxt, comment);
ytxt = displ_pitch_min + 0.75*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Assumed pitch control gain =', ...
                  num2str(pitch_cnts_per_cm), ' [counts/cm aft]');
text(xtxt, ytxt, comment);

% $PITCH_GAIN from log file
ytxt = displ_pitch_min + 0.70*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Reference pitch gain =', ...
                  num2str(logInfo.PITCH_GAIN), ' [counts/cm aft]');
text(xtxt, ytxt, comment);

xtxt = linear_pitch_control_min + 0.4*(linear_pitch_control_max-linear_pitch_control_min);
ytxt = displ_pitch_min + 0.28*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Implied PITCH GAIN =', num2str(AA(2)), ...
                  '\circ/cm' );
text(xtxt, ytxt, comment);

xtxt = linear_pitch_control_min + 0.4*(linear_pitch_control_max-linear_pitch_control_min);
ytxt = displ_pitch_min + 0.21*(displ_pitch_max-displ_pitch_min);
comment = sprintf( 'Implied PITCH VBD SHIFT = %fcm/cc @ density of %.2fg/cc',pitch_vbd_shift_imp,av_apogee_density);
text(xtxt, ytxt, comment);

xtxt = linear_pitch_control_min + 0.35*(linear_pitch_control_max-linear_pitch_control_min);
ytxt = displ_pitch_min + 0.14*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Pitch [\circ] =', num2str(A(1)), '\circ +', ...
                  'Pitch Control [cm] *', num2str(A(2)), '\circ/cm' );
text(xtxt, ytxt, comment);

xtxt = linear_pitch_control_min + 0.35*(linear_pitch_control_max-linear_pitch_control_min);
ytxt = displ_pitch_min + 0.07*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Pitch [\circ] =', num2str(AA(1)), '\circ +', ...
                  'Pitch Control [cm] *', num2str(AA(2)), '\circ/cm +' );
text(xtxt, ytxt, comment);
xtxt = linear_pitch_control_min + 0.3*(linear_pitch_control_max-linear_pitch_control_min);
ytxt = displ_pitch_min + 0.02*(displ_pitch_max-displ_pitch_min);
comment = strcat( 'Buoyancy [g] *', num2str(AA(3)), '\circ/g +', ...
                  'VBD [cc] *', num2str(AA(4)), '\circ/cc +', ...
				  'roll [cm] *', num2str(AA(5)), '\circ/cm' ...
				  );
text(xtxt, ytxt, comment);

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

