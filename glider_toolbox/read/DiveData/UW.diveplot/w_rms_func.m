% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% calculates mean square difference between glider 
% vertical velocity calculated from pressure and that 
% calculated from the hydrodynamic model for steady flight
function [w_rms,vbdc,buoy,hspd_stdy,w_stdy,w_sel,w_stdy_sel] = w_rms_func(vbdbias,a,b,c,abs_compress,therm_expan,temp_ref,vars_need_fit);
  radpd = pi/180.;
  kg2g = 1000;
  cm2m = 0.01;
  
  use_exp_comp = vars_need_fit.use_exp_comp;
  vbd = vars_need_fit.vbd;
  density = vars_need_fit.density;
  density_insitu = vars_need_fit.density_insitu;
  pitch = vars_need_fit.pitch;
  yday = vars_need_fit.yday;
  press = vars_need_fit.press;
  w = vars_need_fit.w;
  vol0 = vars_need_fit.vol0;
  mass = vars_need_fit.mass;
  glider_length = vars_need_fit.glider_length;
  rho0 = vars_need_fit.rho0;
  temp = vars_need_fit.temp;
  hd_s = vars_need_fit.hd_s;
  min_rel_q_inc = vars_need_fit.min_rel_q_inc;
  max_num_iters = vars_need_fit.max_num_iters;
  temp_ref = vars_need_fit.temp_ref;
  mass_comp = vars_need_fit.mass_comp;
  
  %calculate sg buoyancy
  vbdc = vbd - vbdbias;
  if use_exp_comp
    %(w/ hull thermal expansion and compressibility)
    vol = (vol0 + vbdc).*exp(-abs_compress*press+therm_expan*(temp-temp_ref));
  else
    %(w/o hull thermal expansion or compressibility)
    vol = vol0 + vbdc;
  end

  if mass_comp
    vol_comp_ref = kg2g*mass_comp/cml_dens(temp_ref, 0);
    vol_comp = kg2g*mass_comp./cml_dens(temp, press);
    vol = vol + vol_comp - vol_comp_ref;
  end
  
  if use_exp_comp
    %(w corrected volume and in situ density)
    buoy = kg2g*(-mass+density_insitu.*vol*(cm2m)^3);
  else
    %(assumes change in seawater density with pressure and hull compressibility offset)
    buoy = kg2g*(-mass+density.*vol*(cm2m)^3);
  end
  
  %calculate model horizontal and vertical speed based on steady flight model 
  %using buoyancy and pitch
  %glidespeed_stdy and glideangle_stdy describe the behavior of the sg in a zero current environment
  %glidespeed_stdy is in cm/s and glideangle_stdy in degrees
  [glidespeed_stdy,glideangle_stdy,flightvec_iters] = flightvec2(buoy,pitch,yday,glider_length,min_rel_q_inc,max_num_iters,[a b c],rho0);
  %compute the horizontal speed and vertical velocity from the steady flight model
  hspd_stdy = glidespeed_stdy.*cos(glideangle_stdy*radpd);
  w_stdy = glidespeed_stdy.*sin(glideangle_stdy*radpd);
  
  % calculate misfit between w and w_stdy
  inds_w_finite = find(isfinite(w) & isfinite(w_stdy));
  w_sel = w(inds_w_finite);
  w_stdy_sel = w_stdy(inds_w_finite);
  
  % calculate misfit between w and w_stdy for selected data 
  w_s = sum((w_sel - w_stdy_sel).^2);
  if length(w_sel) > 0
    w_rms = sqrt(w_s/length(w_sel));
  else
    w_rms = NaN;
  end
  