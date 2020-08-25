% 
% Copyright (c) 2012 by University of Washington.  All rights reserved. Confidential
%
function [unix_t] = datenum_to_unix(date_num)
  unix_start_dn = datenum([1970, 1, 1, 0, 0, 0]); % for computing secs since unix epoch
  secs_per_day = 86400.0;
  unix_t = (date_num - unix_start_dn)*secs_per_day;