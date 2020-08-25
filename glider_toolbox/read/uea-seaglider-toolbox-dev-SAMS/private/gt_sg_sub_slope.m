function data = gt_sg_sub_slope(data)
% 
% data = gt_sg_sub_slope(data)
%
% Wrapper for the UW glide_slope.m routines which calculates a first
% approximation of glide slope and speed from pressure-derived vertical
% velocities, pitch and hydrodynamic parameters.
%
% All credits to original authors at the University of Washington
%
% B.Y.QUESTE Feb 2015
%

gt_sg_sub_echo({'Calculating naive estimate of glider flight based on pressure change.',...
    'Output in flight substructure as "glide_spd" and "glide_slope".'});

for istep = [data.eng.dive]
    data.flight(istep).glide_vert_spd = -gt_sg_sub_diff({data.hydrography(istep).depth,data.eng(istep).elaps_t}); % m.s-1
    [data.flight(istep).glide_slope, data.flight(istep).glide_spd, data.flight(istep).glide_horz_spd]... % degrees, m.s-1 ?, m.s-1 ?
        = glide_slope(...
        data.flight(istep).glide_vert_spd,... % m.s-1
        data.eng(istep).pitchAng,... % degrees
        data.gt_sg_settings.hd_a,...
        data.gt_sg_settings.hd_b,...
        data.gt_sg_settings.hd_c,...
        data.gt_sg_settings.rho0);
end
end

% ********************************************************************************
% Calculate the theoretical glide angle and speed baseed upon the input
% parameters.
% ********************************************************************************
function [thetadeg, spdg, hspdg] = glide_slope(w, pitch, a, b, c, rho0)

%
% solve for glide slope angle theta[radians] in terms of w & pitch based on model
%

radpd = pi/180;
%  make first guess ---------------------------------

Phi=radpd*pitch;
th = Phi ;
phi=Phi ;

% find non-stall values ---------------------------------

u_prelim = zeros(size(pitch));
non0pitch = find(pitch);
u_prelim(non0pitch) = w(non0pitch)./sin(th(non0pitch));

pfac1 = a*a*sqrt(sqrt(rho0/2));
pfac2 = pfac1*tan(th).*tan(th).*sqrt(0.01*abs(u_prelim));

pfac = pfac2./(4*b*c);

mask =  find(pfac <= 1);
valid = find(pfac > 1);

% iterate and hope it converges -----------------------

f=1;
test=1;
th_rec=0;

while ( ~isempty( find( abs(test) > 0.0001) )  )
  
  fac = zeros(size(th));
  fac(valid) = 4*b*c./(a*a*sqrt(sqrt(rho0/2)) ...
    *tan(th(valid)).*tan(th(valid)).* ...
    sqrt(0.01*w(valid)./sin(th(valid))) );
  th_old=th;
  alp =  -0.5*a*tan(th).*(1 - sqrt(1 - fac) )/c;
  th =  phi -alp*radpd;
  
  test_old=test;
  test = th_old - phi + alp*radpd;
  test(mask)=0;
  mask2 = find( imag(test) ~= 0.0 );
  test(mask2)=0;
  f=f+1; 
  if( f > 20), test=0;,end
end

%toc

th(mask)=Phi(mask);
th(mask2)=Phi(mask2);

theta=th;

spdg = zeros(size(th));
spdg(non0pitch) = abs(w(non0pitch)).*sqrt(1 + ...
  1./( tan(theta(non0pitch)).*tan(theta(non0pitch)) ) );
hspdg = spdg.*cos(theta);

thetadeg = rad2deg(theta);
end