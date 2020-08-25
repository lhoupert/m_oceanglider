% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%
% common display routines for the diveplot_X code

% cheesy trim recommender
% assumes you have run dive_process or have an nc file

pitch_ok = 20; % PARAMETER AD counts close enough on C_PITCH? (~10% of inv_pitch_gain)
cc_vbd_ok = 5; % PARAMETER number of CCs below which we don't have to adjust VBD
fraction_of_change = .5; % the 50% rule

sgc = dir('sg_calib_constants.m');
if (size(sgc,1) == 0)
  % no sg_calib_constants to define id_str so ask
  sgid = input('Glider id: '); 
else
  sg_calib_constants;
  sgid = str2num(id_str);
end
irun = input('Dive cycle number: ');

id_dive = sprintf('%03d%04d',sgid,irun);
base_file = sprintf('p%s',id_dive);
[loginfo,eng,results] = get_dive_data(id_dive);
if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
  fprintf(1,'Processing error or skipped profile\n');
  return;
end
% assumes loginfo, eng, and results are available from caller
unpack_data;

fprintf(1,'\n');
% ----- Trim pitch
% First, look to see if measured pitch is near zero at pitch_ctl of zero
% which indicates the vehicle is level.
% Regress pitch against pitch_ctl to find pitch intercept (predicted pitch at 0 pitch_ctl)
% and pitch slope as the apparent pitch gain on this vehicle.
pitch_control = eng.pitchCtl;
pitch = eng.pitchAng; % measued pitch
inv_pitch_gain = 1/loginfo.PITCH_CNV; % AD counts per cm
pitch_gain = loginfo.PITCH_GAIN;

c_pitch = loginfo.C_PITCH;
pitch_ctl_limit = 3.5; % PARAMETER max pitch control position (cm)
ip = find(pitch_control > -pitch_ctl_limit & pitch_control < pitch_ctl_limit); % remove extreme values
X = [ones(size(pitch_control(ip))) pitch_control(ip)];
A = X\pitch(ip); % regress pitch against pitch_control: A(1) is intercept, A(2) is slope
pitch_gain_imp = A(2);
% NOTE this assumes we take the new gain completely...no fractional adjustment
c_pitch_adjust = fix(( A(1)/ pitch_gain_imp )*-inv_pitch_gain); % note sign
c_pitch_imp = c_pitch + c_pitch_adjust;
if (abs(c_pitch_adjust) < pitch_ok)
  fprintf(1,'Pitch trim is close enough.\n');
  % fall through and deal with VBD
else
  fprintf(1,'Measured pitch response implies changing $C_PITCH by %d AD counts\n',c_pitch_adjust);
  % TODO emit as lines to copy to cmdifle?? e.g., $C_PITCH,xxxx
  fprintf(1,'Recommend changing $C_PITCH from %d to %d\n',c_pitch,fix(c_pitch + c_pitch_adjust*fraction_of_change));
  % TODO recommend this change if > 1 degree or a percentage?
  fprintf(1,'Recommend changing PITCH_GAIN from %.1f to %.1f degrees\n',pitch_gain,A(2))
  % if we are recommending an AD change that would result in more than a 1.5 degree pitch change
  % don't bother with the rest
  % NOTE we use original gain and the entire adjustment
  if (abs(c_pitch_adjust) > 2*inv_pitch_gain/pitch_gain) % PARAMETER [degrees] expected impact on average pitch in degrees
    fprintf(1,'Recommend doing that first before attending to VBD and roll changes!\n')
    % return;
  end
end

fprintf(1,'\n');
% ----- Trim VBD
% The goal, to reduce cost of pumping at apogee and to maintain control at depth,
% is to have the vehicle neutrally buoyant at D_TGT.  Thus we need to match vehicle
% density to the water density insitu at that depth.  This is done by changing C_VBD.
% As the target depth and/or water densities change during the deployment evolution, the 
% VBD trim will likely need continual adjustment.

% If she's moving up or down in the water column, she's too light or heavy respectively.
c_vbd = loginfo.C_VBD;
vbd_cnts_per_cc = 1.0/loginfo.VBD_CNV;
vbd = eng.vbdCC;

if (0)
  % CCE's alternate C_VBD trim calc, assuming mass, mass_comp, and volmax are accurate
  % volmax is the maximum volume at surface, including compressee
  % NOTE this might be below D_TGT (if very heavy) so we'd be ballasting for overshoot
  % of course, if we were light and never made our depth we'd be undershooting
  ctd_sg_press = interp1d(sg_epoch_time,sg_press,ctd_epoch_time,'linear');
  iref = find(ctd_sg_press == max(ctd_sg_press));

  % report measured insitu density at maximum pressure
  % need to compute raw density since rsults.density is not computed over apogee
  density_raw = sw_dens(salin_raw,temp_raw,ctd_sg_press);
  fprintf(1,'Measured water density at %.1fm is %.1fg/m^3.\n', ctd_depth_m(iref),density_raw(iref));
  % compute the volume we need to be at apogee given the density insitu at max pressure
  cc_per_m3 = 1e6;
  vol_apogee = cc_per_m3*mass/sw_dens(salin_raw(iref), temp_raw(iref), ctd_sg_press(iref));
  kg2g = 1000;
  vol_comp = kg2g*mass_comp./cml_dens(temp_raw(iref), ctd_sg_press(iref));
  vol_comp_ref = kg2g*mass_comp/cml_dens(temp_ref, 0);
  % what is volume of hull at surface, independent of compressee?
  vol_hull = volmax - vol_comp_ref;
  vol_hull_compression_factor = (exp(-abs_compress*ctd_sg_press(iref) + therm_expan*(temp_raw(iref) - temp_ref)));
  vol_delta_neutral = (( volmax - vol_hull + ...
                         vol_hull*vol_hull_compression_factor + ...
                         (vol_comp - vol_comp_ref)) - ...
                       vol_apogee);
  c_vbd_true = fix(vbd_min_cnts - vbd_cnts_per_cc*vol_delta_neutral);
  c_vbd_adjust = c_vbd_true - c_vbd; % AD counts
  fprintf(1,'NOTE: Assuming accurate glider mass of %.1fkg, compressee mass of %.2fkg, and volmax of %.1fcc:\n',mass,mass_comp,volmax);
  fprintf(1,'Volume/density calculation implies adjusting $C_VBD by %.1f cc or %d AD counts\n',c_vbd_adjust/vbd_cnts_per_cc,c_vbd_adjust)
  if (c_vbd_true < 0 | c_vbd_true > 4095)
    fprintf(1,'NOTE: Neutral buoyancy impossible to achieve!? Volmax or mass are likely incorrect.\n')
  else
    fprintf(1,'Recommend changing $C_VBD from %d to %d\n',c_vbd,fix(c_vbd + c_vbd_adjust*fraction_of_change));
  end
end

% The usual pilot heuristic to determine how to trim the vehcile is to look
% around apogee, near where she pumps to putative neutral, and see where
% measured w is actually zero (or crossing zero).  Then look at the VBD at
% that point and adjust C_VBD to eliminate that difference to zero.

% There are some caveats: During the first pump to neutral the vehcile is
% pitched down at PITCH_APOGEE so it is actually still traveling down to
% keep up her momentum.  And depending on your pitch center the glider can
% appear to transistion through w zero crossings several times.  Look for
% the last during the apogee pumps only (if none, you have other problems).
% The other thing that can make it confusing is strong up or downwelling at
% depth due to internal waves.

w  = 100.*ctr1stdiffderiv(-eng.depth, sg_time); % observed w (cm/s)
% To try avoid internal waves, etc. look a few meters vertically around the deepest part of dive
max_depth = max(eng.depth);
apo_i = find(eng.depth >= max_depth - 10); % PARAMETER were we turned around at the bottom 10m
% Where was the last place during apogee where we were still descending?
aneg_i = find(sign(w(apo_i)) < 0,1,'last');
if isempty(aneg_i)
  fprintf(1,'Unable to find neutral apogee point?  Missing data?\n');
  return;
else
  % CONSIDER: for better fit, keep going backward for all contiguous negative points
  apo_i = apo_i(aneg_i:end); % eliminate early bits where she was negative or bouncing
  % Where was the first place we went positive after being negative?
  apos_i = find(sign(w(apo_i)) >= 0,1,'first');
  if isempty(apos_i)
    fprintf(1,'Unable to find neutral apogee point?  Missing data?\n');
    return;
  end
  % CONSIDER: for better fit, keep going formward for all contiguous positive points (all remaining)
  apo_i = apo_i(1:apos_i);
end
fit = polyfit(w(apo_i),vbd(apo_i),1); % how does vbd change with w in this range?
% what was vbd when she was stopped vertically?
% note that this is the intercept of the linear fit, e.g., fit(2)
cc_vbd_zw = polyval(fit,[0]); 
neutral_depth_m = mean(eng.depth(apo_i));
if (abs(cc_vbd_zw) < cc_vbd_ok) 
  fprintf(1,'Vertical velocity at %.1fm implies VBD trim is close enough...\n',neutral_depth_m);
else
  % cc_vbd_zw suggests how C_VBD needs to change:  
  % Negative cc_vbd_zw implies she is actually light so we need to increase C_VBD
  c_vbd_adjust = fix(cc_vbd_zw*vbd_cnts_per_cc); % determine adjust in AD counts
  c_vbd_implied = c_vbd + c_vbd_adjust;
  fprintf(1,'Measured vertical velocity at %.1fm implies adjusting $C_VBD by %.1f cc or %d AD counts\n', neutral_depth_m, cc_vbd_zw, c_vbd_adjust);
  fprintf(1,'Recommend changing $C_VBD from %d to %d\n',c_vbd,fix(c_vbd + c_vbd_adjust*fraction_of_change));
  if (0)
    vbd_scale = 10;
    figure
    hold on
    grid on
    plot(w,'b')
    plot(vbd/vbd_scale,'k')
    plot(apo_i,w(apo_i),'.r')
    plot(apo_i,vbd(apo_i)/vbd_scale,'.g')
    legend('w',sprintf('vbd/%d',vbd_scale))
    xlabel('Data point')
    title(sprintf('SG%0d Dive %d (off by %.1f cc)',sgid,irun, cc_vbd_zw))
  end
end

fprintf(1,'\n');
% ----- Trim Roll
roll_control = eng.rollCtl;
roll = eng.rollAng;
inv_roll_cnv = 1/loginfo.ROLL_CNV; % 35.37 AD counts per degree

c_roll_dive = loginfo.C_ROLL_DIVE;
Xd = [ones(size(roll_control(sg_dive_i))) roll_control(sg_dive_i)];
Bd = Xd\roll(sg_dive_i);
c_roll_dive_imp = c_roll_dive - ( Bd(1)/Bd(2) )*inv_roll_cnv;
fprintf(1,'Measured roll on dive implies adjusting $C_ROLL_DIVE from %d to %d\n',c_roll_dive,fix(c_roll_dive_imp));
c_roll_dive_adjust = c_roll_dive_imp - c_roll_dive;
fprintf(1,'Recommend changing $C_ROLL_DIVE from %d to %d\n',c_roll_dive,fix(c_roll_dive + c_roll_dive_adjust*fraction_of_change));

c_roll_climb = loginfo.C_ROLL_CLIMB;
Xc = [ones(size(roll_control(sg_climb_i))) roll_control(sg_climb_i)];
Bc = Xc\roll(sg_climb_i);
c_roll_climb_imp = c_roll_climb - ( Bc(1)/Bc(2) )*inv_roll_cnv;
fprintf(1,'Measured roll on climb implies adjusting $C_ROLL_CLIMB from %d to %d\n',c_roll_climb,fix(c_roll_climb_imp));
c_roll_climb_adjust = c_roll_climb_imp - c_roll_climb;
fprintf(1,'Recommend changing $C_ROLL_CLIMB from %d to %d\n',c_roll_climb,fix(c_roll_climb + c_roll_climb_adjust*fraction_of_change));
% what about roll gain?
