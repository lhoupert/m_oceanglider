% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%
function [theta, spdg] = glide_slope(w, pitch, a, b, c, s, rho0)
%
% solve for glide slope angle theta[radians] in terms of w & pitch based on model
% spdg is in cm/s; w is expected in cm/s
%
% NOTE: this version differs from trunk/matlab/glide_slope in taking s (hd_s)

radpd = pi/180;
cm2m = 0.01;

%  make first guess ---------------------------------
%	tic
Phi = radpd*pitch;
th  = Phi;
phi = Phi;

% find non-stall values ---------------------------------

u_prelim = zeros(size(pitch));
non0pitch = find(pitch);
u_prelim(non0pitch) = w(non0pitch)./sin(th(non0pitch));
cy = a*a*((rho0/2)^(-s));
pfac = ( cy*tan(th).*tan(th).*...
         sqrt( cm2m*abs( u_prelim ) ) )./(4*b*c);

mask =  find(pfac <= 1);
valid = find(pfac > 1);

% iterate and hope it converges -----------------------

f=1; % counter
test=1;

while ( ~isempty( find( abs(test) > 0.0001) )  )
  
  fac = zeros(size(th));
  fac(valid) = 4*b*c./(cy*tan(th(valid)).*tan(th(valid)).* ...
                       sqrt(cm2m*w(valid)./sin(th(valid))) );
  th_old=th;
  alp =  -0.5*a*tan(th).*(1 - sqrt(1 - fac) )/c;
  th =  phi -alp*radpd;
  
  test_old=test;
  test = th_old - phi + alp*radpd;
  test(mask)=0;
  mask2 = find( imag(test) ~= 0.0 );
  test(mask2)=0;
  f=f+1; 
  if( f > 20), test=0; end
end

%toc

th(mask)=Phi(mask);
th(mask2)=Phi(mask2);

theta=th;

spdg = zeros(size(th));
spdg(non0pitch) = abs(w(non0pitch)).*sqrt(1 + ...
  1./( tan(theta(non0pitch)).*tan(theta(non0pitch)) ) );
