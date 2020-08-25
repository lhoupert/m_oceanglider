% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%
function dxdt = ctr1stdiffderiv(x, t)
%
%	Evaluates the centered first difference derivative of x on grid t
%	Endpoints are evaluated using a one-sided difference
%
	n = length(x);
	dxdt = zeros(n,1);
	dxdt(2:n-1) = (x(3:n) - x(1:n-2))./(t(3:n) - t(1:n-2));
	dxdt(1) = (x(2) - x(1))/(t(2) - t(1));
	dxdt(n) = (x(n) - x(n-1))/(t(n) - t(n-1));
