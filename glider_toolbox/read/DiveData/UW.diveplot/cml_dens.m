% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

%
%   Calculates compressee/expansee density [g/cc] 
%   Uses a polynomial fit based on environmental chamber and sound speed measurements
%
%   temperature [deg C], pressure [dbar, 0=sea surface]
%   A are polynomial coefficients for density as a function of temperature
%   at atmospheric pressure (p=0), nominally size(A) = [1 2]
%   B are coefficients for the derivative of density with respect to
%   pressure, nominal size(B) = [12 1] for a fit with up to quadratic terms
%   in temperature and cubic terms in pressure
%

function density = cml_dens(temperature, pressure, A, B, fluidtype)

P1 = pressure;
P2 = P1.*P1;
P3 = P1.*P2;
P4 = P2.*P2;
% Do this once...
% P1 = P1/1;
P2 = P2/2;
P3 = P3/3;
P4 = P4/4;

T1 = temperature;
T2 = T1.*T1;

% from environmental chamber density measurements
    
A_red = [-1.11368512753634 796.461657048578];
A_blue = [-1.07743970211476 790.951260319579];
    
% from PTV sound velocity measurements 10/23/12 and 02/06/13
    
B_red_T2P3 = [0.0102261734465205; 8.47816063186125e-05; 4.47694544079322e-07; ...
			  -1.51575213805701e-06; -2.75610973002708e-08; -5.50259771402447e-11; ...
			  2.03900944994636e-10; 5.90718779795042e-12; -5.66873780087388e-14; ...
			  -1.27260138364051e-14; -5.57122868190274e-16; 1.17249605451857e-17];

B_red_T2P2 = [0.0102052829145449; 8.52182108882249e-05; 4.34927182961885e-07; ...
			  -1.30186206661706e-06; -3.03705760249538e-08; 2.88293344499584e-10; ...
			  9.52846703487369e-11; 4.45151822732093e-12; -1.00703879876029e-13];

B_blue_T2P3 = [0.0102728079859404; 7.87130783574905e-05; 6.0665972076054e-07; ...
			   -1.47915664017977e-06; -3.49451391145939e-08 ;2.20151163102027e-10; ...
			   1.91632424033919e-10; 1.04377013186988e-11; -2.33887126509371e-13; ...
			   -1.17991005158086e-14; -1.12934912904199e-15; 3.47861909690662e-17];

B_blue_T2P2 = [0.0102684042928214; 7.69524163981727e-05; 6.48651329885882e-07; ...
			   -1.30702391855726e-06; -2.79767681939615e-08; 2.16791031520607e-10; ...
			   9.7121608458588e-11; 4.28505458495888e-12; -9.9768459288306e-14];

if (nargin < 3)   % default to red (low purity) fluid, cubic model
    A = A_red;
    B = B_red_T2P3;
end
if (nargin == 5)
    if(fluidtype(1) == 'r')
        A = A_red;
        B = B_red_T2P3;
    elseif(fluidtype(1) == 'b') % blue (reagent grade) fluid, cubic model
        A = A_blue;
        B = B_blue_T2P3;
    end
end

rho0 = polyval(A, T1);

% analytic integration of d(rho)/dp from its polynomial coefficients
drho = B(1)*P1  + B(2)*T1.*P1  + B(3)*T2.*P1  + ...
       B(4)*P2  + B(5)*T1.*P2  + B(6)*T2.*P2  + ...
       B(7)*P3  + B(8)*T1.*P3  + B(9)*T2.*P3  ;
   if(length(B) == 12)
       drho = drho + B(10)*P4 + B(11)*T1.*P4 + B(12)*T2.*P4  ;
   end
   
rho = rho0 + drho; % compressee density [kg/m^3]
rho = rho/1000; % [g/cc]  (1E3 g/kg)/(1E6 cc/m^3)
density = reshape(rho, size(temperature)); 
