% 
% Copyright (c) 2012 by University of Washington.  All rights reserved. Confidential
%
function [date_num] = unix_to_datenum(unix_t)
% matlab does not support ctime, asctime, etc. so...
% convert from UNIX epoch secs to serial date time
% dn 719529 is the number of days from 1 Jan 0000 to 1 Jan 1970
% => number of days/year = 365.2431472081219
  unix_start_dn = datenum([1970, 1, 1, 0, 0, 0]); % for computing secs since unix epoch
  secs_per_day = 86400.0;
  date_num = floor(unix_t/secs_per_day) + ... 
      rem(unix_t,secs_per_day)/ secs_per_day + ...
      unix_start_dn; 
  