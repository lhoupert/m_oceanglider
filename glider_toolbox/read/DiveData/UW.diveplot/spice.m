% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%
function spiciness = spice(p,t,s)
% function spiciness = spice(p,t,s)
% adapted from algorithm developed by
% P. Flament.
%
% SCKennan(Dec92)

B(1,1) = 0;
B(1,2) = 7.7442e-001;
B(1,3) = -5.85e-003;
B(1,4) = -9.84e-004;
B(1,5) = -2.06e-004;

B(2,1) = 5.1655e-002;
B(2,2) = 2.034e-003;
B(2,3) = -2.742e-004;
B(2,4) = -8.5e-006;
B(2,5) = 1.36e-005;

B(3,1) = 6.64783e-003;
B(3,2) = -2.4681e-004;
B(3,3) = -1.428e-005;
B(3,4) = 3.337e-005;
B(3,5) = 7.894e-006;

B(4,1) = -5.4023e-005;
B(4,2) = 7.326e-006;
B(4,3) = 7.0036e-006;
B(4,4) = -3.0412e-006;
B(4,5) = -1.0853e-006;

B(5,1) = 3.949e-007;
B(5,2) = -3.029e-008;
B(5,3) = -3.8209e-007;
B(5,4) = 1.0012e-007;
B(5,5) = 4.7133e-008;

B(6,1) = -6.36e-010;
B(6,2) = -1.309e-009;
B(6,3) = 6.048e-009;
B(6,4) = -1.1409e-009;
B(6,5) = -6.676e-010;


[r,c] = size(t);
sp = zeros(r,c);
s = s - 35.*ones(r,c);
T = 1.*ones(r,c);
for i = 1:6
    S = ones(r,c);
    for j = 1:5
        sp = sp + B(i,j).*T.*S;
        S = S.*s;
      end
    T = T.*t;
  end
spiciness = sp;
