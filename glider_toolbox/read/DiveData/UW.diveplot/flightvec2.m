% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

% [umag,thdeg] = flightvec2(bu,ph,vmtime,xl,min_rel_q_inc,max_num_iters,flight_consts,rho0)
%
% solves unaccelerated flight equations iteratively 
% for along glideslope speed umag and glideangle thdeg
% umag and thdeg describe flight in the presence of zero current
%
% inputs are: the net buoyancy (bu), pitch in degrees (ph), 
%             time since start of dive (vmtime), sg fairing length (xl),
%             min_rel_q_inc: the desired tolerance for relative q change,
%                            between successive iterations,
%             max_num_iters: the maximum number of iterations if 
%                            min_rel_q_inc is not reached,
%             flight_consts = [hd_a hd_b hd_c], and rho0 a ref. density

function [umag,thdeg,j] = flightvec2(bu,ph,vmtime,xl,min_rel_q_inc,max_num_iters,flight_consts,rho0);

gravity = 9.82;
g2kg = 0.001;
m2cm = 100;

%hd_a is in 1/deg. units, hd_b has dimensions q^(1/4), hd_c is in 1/deg.^2 units
hd_a = flight_consts(1); hd_b = flight_consts(2); hd_c = flight_consts(3);
hd_s = -1/4;
%convert net buoyancy to force dimensions (N units)
buoyforce = g2kg*gravity*bu;

%initialize flight variables
%initial glideslope set to 1 or -1 depending on whether sg is pos. or neg. buoyant
%if sign(buoyforce)=0, that value will not be a member of valid below
th = (pi/4)*sign(buoyforce);
%initial dynamic pressure q for vertical flight (always a positive quantity)
%q is determined from the drag eqn. (Eqn. 2 in Eriksen et al.) by assuming alpha = 0
q = (abs(buoyforce.*sin(th))/(xl*xl*hd_b)).^(1/(1+hd_s)); 	

%initialize variables used in while loop
param = NaN*ones(size(bu));
alpha = NaN*ones(size(bu));
thdeg = NaN*ones(size(bu));

%this loop iterates the dynamic pressure (q) to get glidespeed 
%and the attack angle (alpha) to get glideangle 
%buoyancy is taken to be an accurately known quantity and is not iterated
j = 1;
test = (1+0.1)*min_rel_q_inc*ones(size(q)); %force loop to run at least once
while(~isempty(find(abs(test) > min_rel_q_inc)) & j <= max_num_iters)   
   %param_inv is the reciprocal of the term subtracted 
   %from 1 under the radical in the solution equations (Eqs. 7-8)
   %here param_inv is estimated using the current th
   param_inv = hd_a^2*tan(th).^2.*q.^(-hd_s)/(4*hd_b*hd_c);
   
   %valid solutions exist for param_inv >= 1
   %also only consider points where buoyancy and pitch are both up or both down
   valid = find(param_inv >= 1 & sign(bu).*sign(ph) > 0);

   %th and q are not zero for valid points, since otherwise param_inv would = 0
   param(valid) = 1./param_inv(valid);
   %hold current q in q_old
   q_old = q;
   %Eq. 7 in the reference, obtained using the quadratic formula ...
   %q^(1/4) considered to vary slowly compared to q
   q(valid) = buoyforce(valid).*sin(th(valid)).*q(valid).^(-hd_s)./(2*xl^2*hd_b) ...
       .*(1 + sqrt(1-param(valid)));

   % NOTE that once a point stalls, it can never recover
   % even though our guess at th and q might be off...
   % the problem is we don't know how to relax either
   % locally to possibly ensure or encourage a lasting
   % non-stall solution.
   stall = setdiff([1:length(q)],valid);
   q(stall) = NaN;

   %calculate relative difference to see if next iteration is necessary
   test = (q-q_old)./q;
   %Eq. 8 in the reference       
   alpha(valid) = -hd_a*tan(th(valid))./(2*hd_c).*(1 - sqrt(1-param(valid)));
   alpha(stall) = NaN;
   %glideangle is the difference between pitch and angle of attack
   thdeg(valid) = ph(valid) - alpha(valid);
   thdeg(stall) = NaN;   
   
   %switch theta to radians
   th = thdeg*pi/180;
   
   %increase iteration number
   j = j+1;
end


[m,n] = size(vmtime);
if m == 1 | n == 1 % vector?

  %for values of th and q near apogee for which model couldn't be used,
  %linearly interpolate to determine these values
  good_inds = find(isfinite(vmtime) & isfinite(thdeg));
  if length(good_inds) >= 2
    thdeg = interp1(vmtime(good_inds),thdeg(good_inds),vmtime,'linear');
    %force leading and trailing values where NaNs occur since interpolation wasn't possible
    %to be 0 ... otherwise extreme values may be produced here
    good_thdeg_inds = find(isfinite(thdeg));
    if good_thdeg_inds(1) ~= 1
      thdeg(1:good_thdeg_inds(1)-1) = 0;
    end
    if good_thdeg_inds(end) ~= length(thdeg)
      thdeg(good_thdeg_inds(end)+1:end) = 0;
    end
  end
  
  good_inds_2 = find(isfinite(vmtime) & isfinite(q));
  if length(good_inds_2) >= 2
    q = interp1(vmtime(good_inds_2),q(good_inds_2),vmtime,'linear');
    good_q_inds = find(isfinite(q));
    if good_q_inds(1) ~= 1
      q(1:good_q_inds(1)-1) = 0;
    end
    if good_q_inds(end) ~= length(q)
      q(good_q_inds(end)+1:end) = 0;
    end
  end
end
%glidespeed in cm/s
umag = m2cm*sqrt(2*q/rho0);
