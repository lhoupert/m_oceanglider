% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% function w_misfit_func(x);
function w_rms = w_misfit_func(x,vars_need_fit,flight_consts,vol_consts);
  global HIST;
  % Unpack arguments from x according to task at hand
  fit_abc = vars_need_fit.fit_abc;
  fit_exp_comp = vars_need_fit.fit_exp_comp;
  
  if fit_abc == 0 & fit_exp_comp == 0 %optimize only vbdbias
    vbdbias = x(1);
    a = flight_consts(1); b = flight_consts(2); c = flight_consts(3);
    abs_compress = vol_consts(1); therm_expan = vol_consts(2); temp_ref = vol_consts(3);
  elseif fit_abc == 1 & fit_exp_comp == 0 %optimize vbdbias, a, b, and c
    vbdbias = x(1);
    a = x(2); b = x(3); c = x(4);
    abs_compress = vol_consts(1); therm_expan = vol_consts(2); temp_ref = vol_consts(3);
  elseif fit_exp_comp == 1 & fit_abc == 0 %optimize therm_expan and abs_compress
    vbdbias = vol_consts(1); % vbdbias_const
    a = flight_consts(1); b = flight_consts(2); c = flight_consts(3);
    abs_compress = x(1); therm_expan = x(2); temp_ref = vol_consts(2); % vol_const
  end
  
  [w_rms,ignore,ignore,ignore,ignore,ignore,ignore] = w_rms_func(vbdbias,a,b,c,abs_compress,therm_expan,temp_ref,vars_need_fit);
  
  HIST = [HIST; [w_rms vbdbias x]];
