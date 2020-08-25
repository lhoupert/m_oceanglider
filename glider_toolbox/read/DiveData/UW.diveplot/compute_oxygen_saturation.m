% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%
function [oxygen_sat_seawater,oxygen_sat_fresh_water,oxygen_sat_salinity_adjustment] = compute_oxygen_saturation(temp,salin)
% Compute (max) oxygen saturation (solubility) (uM/kg)
% at standard air and pressure (1013hPa) for fresh and seawater at given temperature and salinity
%    Garcia and Gordon, 'Oxygen solubility in seawater: Better fitting equations'
%    Limnol. Oceanogr. 37(6), 1992 1307-1312
%    Equation (8) page 1310
  
% Constants for calculation of Oxygen saturation,
% which depends on temperature and salinity only
  
% Note: Argo uses the Benson and Krause constants (first column, Table 1)
% while Aanderaa uses the 'combined fit' constants
  if 0
    % 'combined fit'
    oA_micromolar = [2.00856; 3.22400; 3.99063; 4.80299; 0.978188; 1.71096];
    oB_micromolar = [6.24097E-3; 6.93498E-3; 6.90358E-3; 4.29155E-3];
    oC0_micromolar = -3.11680E-7;
  else
    % Benson & Krause constants, which Argo prefers, and which Garcia and Gordon recommend (last para)
    oA_micromolar = [2.009070; 3.220140; 4.050100; 4.944570; -0.256847; 3.887670];
    oB_micromolar = [-6.24523E-3; -7.37614E-3; -1.0341E-2; -8.17083E-3];
    oC0_micromolar = -4.88682E-7;
  end
  
  temp_s = log( (298.15 - temp)./(273.15 + temp) ); % solubility temperature
  
  % oxygen solubility for fresh water in ml/l = cm^3/dm^3
  oxygen_sat_fresh_water = exp( polyval( flipud(oA_micromolar), temp_s) );  % [ml/L]
  oxygen_sat_salinity_adjustment = exp(salin .* polyval(flipud(oB_micromolar),temp_s) + oC0_micromolar*salin.*salin); % percentage
  oxygen_sat_seawater = oxygen_sat_salinity_adjustment .* oxygen_sat_fresh_water; % [ml/L]
  o2_molar_mass = 44.614; %  uM/kg * L/ml of oxygen at standard pressure and temperature
  oxygen_sat_fresh_water = oxygen_sat_fresh_water*o2_molar_mass; % [uM/kg]
  oxygen_sat_seawater = oxygen_sat_seawater*o2_molar_mass; % [uM/kg]
  