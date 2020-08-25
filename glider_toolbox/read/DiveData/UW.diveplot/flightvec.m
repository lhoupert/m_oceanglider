% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%

%	function flightvec(bu, ph, xl, a, b, c, rho0, s)
%		solves unaccelerated flight equations iteratively 
%		for speed magnitude umag (cm/s) and glideangle thdeg (radians)

function [ umag, thdeg ] = flightvec( bu, ph, xl, a, b, c, rho0, s );
  
  gravity = 9.82;
  g2kg = 0.001;
  m2cm = 100;
  %convert net buoyancy to force dimensions (N units)
  buoyforce = g2kg*gravity*bu;
  
  %initialize flight variables
  %initial glideslope set to 1 or -1 depending on whether sg is pos. or neg. buoyant
  %if sign(buoyforce)=0, that value will not be a member of valid below
  th = (pi/4)*sign(buoyforce);

  %initial dynamic pressure q for vertical flight (always a positive quantity)
  %q is determined from the drag eqn. (Eqn. 2 in Eriksen et al.) by assuming alpha = 0
  q = ( sign(bu).*buoyforce/(xl*xl*b) ).^(1/(1+s)); 	% dynamic pressure for vertical flight
  
  param = NaN*ones(size(bu));
  alpha = NaN*ones(size(bu));
  %DEAD valid = find( bu ~= 0 & sign(bu).*sign(ph) > 0 );
  %DEAD umag = zeros(size(bu));
  thdeg = NaN*ones(size(bu));
  
  tol = 0.001;
  j = 0;
  q_old = zeros(size(bu));
  while(~isempty( find( abs( (q-q_old)./q ) > tol )) & j <= 15);
    % trace_array(sprintf('q_%d', j), q)
    %param_inv is the reciprocal of the term subtracted 
    %from 1 under the radical in the solution equations (Eqs. 7-8)
    %here param_inv is estimated using the current th
    param_inv = a*a*tan(th).*tan(th).*q.^(-s)/(4*b*c);
    
    %valid solutions exist for param_inv >= 1
    %also only consider points where buoyancy and pitch are both up or both down
    valid = find( param_inv > 1 & sign(bu).*sign(ph) > 0);	% valid solutions for param < 1
    
    %DEAD param(valid) = 4*b*c./(a*a*tan(th(valid)).*tan(th(valid)).*q(valid).^(-s));
    %th and q are not zero for valid points, since otherwise param_inv would = 0
    param(valid) = 1./param_inv(valid);

    % NOTE that once a point stalls, it can never recover
    % even though our guess at th and q might be off...
    % the problem is we don't know how to relax either
    % locally to possibly ensure or encourage a lasting
    % non-stall solution.
    stall = setdiff([1:length(q)],valid);
    q(stall) = NaN;
    
    %hold current q in q_old
    q_old = q;

    %Eq. 7 in the reference, obtained using the quadratic formula ...
    %q^(1/4) considered to vary slowly compared to q
    q(valid) = ( buoyforce(valid).*sin(th(valid))./(2*xl*xl*b*q(valid).^(s)) ).*  ...
        (1 + sqrt(1-param(valid)));
    %Eq. 8 in the reference       
    alpha(valid) = ( -a*tan(th(valid))./(2*c) ).*(1 - sqrt(1-param(valid)));
    alpha(stall) = NaN;

    thdeg(valid) = ph(valid) - alpha(valid);
    thdeg(stall) = NaN;
    
    th = thdeg*pi/180;
    % trace_array(sprintf('th_%d', j), th)
    % trace_array(sprintf('param_inv_%d', j), param_inv)
    j = j+1;
    
  end;
  %glidespeed in cm/s
  umag = m2cm*sqrt( 2*q/rho0 );


