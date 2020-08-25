% 
% Copyright (c) 2012 by University of Washington.  All rights reserved. Confidential
%
function [interp_data_v] = interp1d(first_time_v,data_v,second_time_v,kind)
  % like interp1 but handles non-overlapping time vectors
  if (nargin == 3)
    kind = 'linear';
  end
  % add 'nearest' data item to the ends of data and first_time_s_v
  if (second_time_v(1) < first_time_v(1))
    data_v = [data_v(1); data_v];
    first_time_v = [second_time_v(1); first_time_v];
  end
  if (second_time_v(end) > first_time_v(end))
    data_v = [data_v; data_v(end)];
    first_time_v = [first_time_v; second_time_v(end)];
  end
  % If data_v contains NaN matlab warns about undefined values but does the right thing
  interp_data_v = interp1(first_time_v,data_v,second_time_v,kind);