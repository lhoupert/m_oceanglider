% /Users/jsb/Seaglider/software/trunk/basestation2/matlab/regress_roll.m, Fri Aug  9 08:07:01 2013, Edit by jsb
% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

%  regress_roll.m

%  Read sg data files for a specified set of dive cycles and 
%	find the roll control bias
%
%  functions called:
%    sg_calib_constants.m - establishes glider calibration constants
%    various seawater routines
%    get_dive_data.m - data and results from nc or nc+mat
%
%DEAD full_pitch_models = 0; % use experimental pitch models
display_fit = 0;
depthmin = 15; % PARAMTER [m] at least 15 down (below flare)
dwdtmax = 1; % PARAMETER [cm/s^2] vertical acceleration of vehicle is < 1/cm/s^2
vbddiffmin = 0.5; % PARAMATER [cc] vbd isn't moving more than 0.5cc (avoid pumps)
rollmin = 55; % PARAMETER [deg] we are not rolled over more than 55deg

roll_control_max = 4000; % max roll AD counts 
roll_control_min = 0; % min roll AD counts

timelabel = datestr(clock);

%sg_calib_constants; % get id_str %billr: causes failure in compiled version
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

c_pitch_ref = NaN;
dives_used = [];
dive_pitch_desired = [];
PC = [];
RL = [];
RC = [];
RCC = [];
TR = [];
%
% loop on dive cycles
%
id_str = settings.seagliderID; %billr: this replaces the sg_calib_constants line
for ir = dives'
	
	id_dive = sprintf('%03d%04d',id_str,ir);
	[loginfo,eng,results] = get_dive_data(id_dive,1); % need results: gsm speed, density, etc.
	if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
		fprintf(1,'Unable to include data for dive %d\n', ir);
		continue;
	end
	unpack_data;
	fprintf(1, 'p%s\n', id_dive)
	
	c_pitch = loginfo.C_PITCH;
	c_roll_dive = loginfo.C_ROLL_DIVE;
	c_roll_climb = loginfo.C_ROLL_CLIMB;
	c_vbd = loginfo.C_VBD;
	if (isnan(c_pitch_ref))
		c_pitch_ref = c_pitch;
		c_roll_dive_ref = c_roll_dive;
		c_roll_climb_ref = c_roll_climb;
	end
	if 0
	if (c_pitch ~= c_pitch_ref | ...
		c_roll_dive ~= c_roll_dive_ref | ...
		c_roll_climb ~= c_roll_climb_ref)
		fprintf(1,'Centers changed; unable to include data for dive %d\n', ir);
		continue;
	end
	end
	roll_min_cnts = loginfo.ROLL_MIN;
	roll_max_cnts = loginfo.ROLL_MAX;
	vbd_min_cnts = loginfo.VBD_MIN;
	
	vbd_cnts_per_cc   = 1./loginfo.VBD_CNV;% same for all dives
	roll_cnts_per_deg = 1./loginfo.ROLL_CNV; % 35.37 AD counts per degree
	pitch_cnts_per_cm = 1./loginfo.PITCH_CNV; % 217.39 AD counts per cm
	
	%
	% read sg data
	%
	vmtime = eng.elaps_t;
	mp = length(vmtime);
	% depth in meters, as computed onboard by SG
	% using pressure gauge and linear conversion
	vmdepth  = eng.depth;      
	hdg = eng.head;
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
	
	PC = [PC ; pitch_control(iglide_stdy)];
	RL = [RL ; roll(iglide_stdy)];
	RC = [RC ; roll_control(iglide_stdy)];
	RCC = [RCC ; roll_control_counts(iglide_stdy)];
	TR = [TR; turn_rate(iglide_stdy)];
end; % end dive cycle loop

%
% rename arrays
%
pitch_control = PC;
roll = RL;
roll_control = RC;
roll_control_counts = RCC;
turn_rate = TR;

n_dives_used = length(dives_used);
mp = length(roll);
if (mp == 0)
	fprintf(1,'No data available!\n');
	return;
end
dives_used = succinct_elts(dives_used);
titlestring = sprintf('Seaglider %s %s [%s]', id_str, underscore(mission_title), dives_used);

iwn = find( pitch_control < 0 ); % dive points
iwp = find( pitch_control > 0 ); % climb points

%
% find linear fit between roll and roll_control_counts
%

Xd = [ones(size(roll_control_counts(iwn))) roll_control_counts(iwn)];
Bd = Xd\roll(iwn);
c_roll_dive_imp = -Bd(1)/Bd(2);

Xc = [ones(size(roll_control_counts(iwp))) roll_control_counts(iwp)];
Bc = Xc\roll(iwp);
c_roll_climb_imp = -Bc(1)/Bc(2);

roll_control_max_dive = (roll_max_cnts - c_roll_dive_ref)/roll_cnts_per_deg;
roll_control_min_dive = (roll_min_cnts - c_roll_dive_ref)/roll_cnts_per_deg;
roll_control_max_climb = (roll_max_cnts - c_roll_climb_ref)/roll_cnts_per_deg;
roll_control_min_climb = (roll_min_cnts - c_roll_climb_ref)/roll_cnts_per_deg;
roll_control_max_dive_imp = (roll_max_cnts - c_roll_dive_imp)/roll_cnts_per_deg;
roll_control_min_dive_imp = (roll_min_cnts - c_roll_dive_imp)/roll_cnts_per_deg;
roll_control_max_climb_imp = (roll_max_cnts - c_roll_climb_imp)/roll_cnts_per_deg;
roll_control_min_climb_imp = (roll_min_cnts - c_roll_climb_imp)/roll_cnts_per_deg;


%
% plot scatter of roll with roll_control_counts
%
roll_max = 40; % max pitch degrees
roll_min = -40; % min pitch degrees

figure;
plot(roll_control_counts(iwn), roll(iwn), '.b');
hold on;
grid on;
plot(roll_control_counts(iwp), roll(iwp), '.r');
plot([0], [0], '+r', 'MarkerSize', 10)
axis([roll_control_min roll_control_max roll_min roll_max]);
xlabel('Roll Control [counts]');
ylabel('Roll [\circ]');
title(titlestring, 'FontSize', 12);
lg = legend('Dive', 'Climb', 'Location', 'SouthEast');
set(lg,'color','none'); % make transparent

ix = [roll_control_min:0.01*(roll_control_max-roll_control_min):roll_control_max];
plot(ix, Bd(1) + Bd(2)*ix, '-b');
plot(ix, Bc(1) + Bc(2)*ix, '-r');

xtxt = roll_control_min + 0.4*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.95*(roll_max-roll_min);
comment = strcat( 'Reference C_R_O_L_L^d^i^v^e =', num2str(c_roll_dive_ref) );
text(xtxt, ytxt, comment);
xtxt = roll_control_min + 0.7*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.95*(roll_max-roll_min);
comment = strcat( 'Reference C_R_O_L_L^c^l^i^m^b =', num2str(c_roll_climb_ref) );
text(xtxt, ytxt, comment);

xtxt = roll_control_min + 0.1*(roll_control_max-roll_control_min);    
ytxt = roll_min + 0.85*(roll_max-roll_min);
comment = strcat( 'Dive roll control range =', ...
                  num2str(roll_control_min_dive), ' to ', ...
                  num2str(roll_control_max_dive), ' \circ'); 
text(xtxt, ytxt, comment);

xtxt = roll_control_min + 0.4*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.90*(roll_max-roll_min);
comment = strcat( 'Implied C_R_O_L_L^d^i^v^e =', num2str(c_roll_dive_imp) );
text(xtxt, ytxt, comment);
xtxt = roll_control_min + 0.7*(roll_control_max-roll_control_min);
comment = strcat( 'Implied C_R_O_L_L^c^l^i^m^b =', num2str(c_roll_climb_imp) );
text(xtxt, ytxt, comment);

xtxt = roll_control_min + 0.1*(roll_control_max-roll_control_min);  
ytxt = roll_min + 0.80*(roll_max-roll_min);
comment = strcat( 'Implied dive roll control range =', ...
                  num2str(roll_control_min_dive_imp), ' to ', ...
                  num2str(roll_control_max_dive_imp), ' \circ'); 
text(xtxt, ytxt, comment);
ytxt = roll_min + 0.75*(roll_max-roll_min);
comment = strcat( 'Implied climb roll control range =', ...
                  num2str(roll_control_min_climb_imp), ' to ', ...
                  num2str(roll_control_max_climb_imp), ' \circ'); 
text(xtxt, ytxt, comment);    

ytxt = roll_min + 0.70*(roll_max-roll_min);
comment = strcat( 'Assumed roll control gain =', ...
                  num2str(roll_cnts_per_deg), ' [counts/\circ to stbd]' );
text(xtxt, ytxt, comment);

xtxt = roll_control_min + 0.25*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.08*(roll_max-roll_min);
comment = strcat( 'Roll_d_i_v_e [\circ] =', num2str(Bd(1)), ...
                  '\circ + Roll Control [counts] /', num2str(1/Bd(2)), ...
                  '[counts/\circ]' );
text(xtxt, ytxt, comment);
ytxt = roll_min + 0.02*(roll_max-roll_min);
comment = strcat( 'Roll_c_l_i_m_b [\circ] =', num2str(Bc(1)), ...
                  '\circ + Roll Control [counts] /', num2str(1/Bc(2)), ...
                  '[counts/\circ]' );
text(xtxt, ytxt, comment);

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);

%
% find linear fit between turn_rate with roll_control_counts
%
Xd = [ones(size(roll_control_counts(iwn))) roll_control_counts(iwn)];
Cd = Xd\turn_rate(iwn);
c_roll_turn_dive = -Cd(1)/Cd(2);
Xc = [ones(size(roll_control_counts(iwp))) roll_control_counts(iwp)];
Cc = Xc\turn_rate(iwp);
c_roll_turn_climb = -Cc(1)/Cc(2);

%
% plot scatter of turn rate with roll_control_counts
%
roll_max = 2.5; % max stbd turn rate deg/s
roll_min = -2.5; % max port turn rate deg/s

figure;
plot(roll_control_counts(iwn), turn_rate(iwn), '.b');
hold on;
grid on;
plot(roll_control_counts(iwp), turn_rate(iwp), '.r')
plot([0], [0], '+k', 'MarkerSize', 10)
axis([roll_control_min roll_control_max roll_min roll_max]);
xlabel('Roll Control [counts]');
ylabel('Turn rate [degrees/s] (starboard positive)');
title(titlestring, 'FontSize', 12);
lg = legend('Dive', 'Climb', 'Location', 'SouthEast');
set(lg,'color','none'); % make transparent

ix = [roll_control_min:0.01*(roll_control_max-roll_control_min):roll_control_max];
plot(ix, Cd(1) + Cd(2)*ix, '-b');
plot(ix, Cc(1) + Cc(2)*ix, '-r');

%DEADxtxt = roll_control_min + 0.1*(roll_control_max-roll_control_min);
%DEADytxt = roll_min + 0.85*(roll_max-roll_min);
%DEADcomment = strcat( 'Dive Cycles =', runs_txt);
%DEADtext(xtxt, ytxt, comment, 'Color', 'k');

xtxt = roll_control_min + 0.4*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.95*(roll_max-roll_min);
comment = strcat( 'C_R_O_L_L^d^i^v^e reference=', num2str(c_roll_dive_ref) );
text(xtxt, ytxt, comment);
xtxt = roll_control_min + 0.7*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.95*(roll_max-roll_min);
comment = strcat( 'C_R_O_L_L^c^l^i^m^b reference=', num2str(c_roll_climb_ref) );
text(xtxt, ytxt, comment);

xtxt = roll_control_min + 0.4*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.9*(roll_max-roll_min);
comment = strcat( 'Implied C_R_O_L_L^d^i^v^e =', num2str(c_roll_turn_dive) );
text(xtxt, ytxt, comment);
xtxt = roll_control_min + 0.7*(roll_control_max-roll_control_min);
comment = strcat( 'Implied C_R_O_L_L^c^l^i^m^b =', num2str(c_roll_turn_climb) );
text(xtxt, ytxt, comment);

xtxt = roll_control_min + 0.1*(roll_control_max-roll_control_min);  
ytxt = roll_min + 0.8*(roll_max-roll_min);
comment = strcat( 'Assumed roll control gain =', ...
                  num2str(roll_cnts_per_deg), ' [counts/\circ to stbd]' );
text(xtxt, ytxt, comment);

xtxt = roll_control_min + 0.2*(roll_control_max-roll_control_min);
ytxt = roll_min + 0.14*(roll_max-roll_min);
comment = strcat( 'Dive Turn Rate [\circ/s] =', num2str(Cd(1)), ...
                  '\circ/s + Roll Control [counts] /', num2str(1/Cd(2)), ...
                  '[counts/\circ/s');
text(xtxt, ytxt, comment);

ytxt = roll_min + 0.06*(roll_max-roll_min);
comment = strcat( 'Climb Turn Rate [\circ/s] =', num2str(Cc(1)), ...
                  '\circ/s + Roll Control [counts] /', num2str(1/Cc(2)), ...
                  '[counts/\circ/s');
text(xtxt, ytxt, comment);

v = axis;
ytxt = v(3) -0.1*(v(4)-v(3));
xtxt = v(1);
text(xtxt , ytxt, timelabel);
