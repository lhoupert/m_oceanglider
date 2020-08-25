% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%
% determine which valid points in the sensor data stream are parts of dive and climb
% should match to the sensor depth!!
function [sen_depth,sen2sg_i,sg_sen_valid_i,sen2sg_valid_i,sen_dive_i,sen_climb_i] = compute_sensor_sg_indices(sen_time,sen_qc,sg_time,sg_depth,sg_dive_i,sg_climb_i)
  sen_depth = interp1d(sg_time,sg_depth,sen_time,'linear'); % depth of readings (for plots)
  sg_np = length(sg_time);
  sen_np = length(sen_time);
  sen2sg_i = interp1d(sen_time,[1:sen_np]',sg_time,'nearest'); % nearest sensor point to a SG point

  sg_uncorrectable_i = bad_qc(sen_qc(sen2sg_i)); % nearest uncorrectable SG points
  sg_sen_valid_i = setdiff([1:sg_np],sg_uncorrectable_i); % which sg points have valid sensor points?
  sen2sg_valid_i = sen2sg_i(sg_sen_valid_i); % nearest VALID sensor point to a SG point

  sen_uncorrectable_i = bad_qc(sen_qc);
  sen_valid_i = setdiff([1:sen_np],sen_uncorrectable_i); % which are the valid sensor points?
  sen_dive_i  = sen_valid_i(find(sen_time(sen_valid_i) <= sg_time(sg_dive_i(end))));
  sen_climb_i = sen_valid_i(find(sen_time(sen_valid_i) >= sg_time(sg_climb_i(1))));
