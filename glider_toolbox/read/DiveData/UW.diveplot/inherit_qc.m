% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%

% An important note about the programmatic use of qc values
% Essesntially there are two phases in the life of a QC tag in a vector
% The first is 'imperative' wherein the code decides something is, say, QC_BAD
% or needs to be QC_INTERPOATED.  The value is asserted and then a later bit of code
% looks at the tags and does something to the associated value (sets it to NaN or
% actually does the interpolation).
% The second is 'declarative' where a change in one variable implies a change in another.
% The classic example is interpolating temp spikes will cause a change in some temp values
% which implies a change in the corresponding derived salinity value. 
% Of course, at the end, you want to declare that salinity was interpolated
% but you don't want to mark it (for) interpolation during the iperative phase.
% For this reason, use inherit_qc after all the calculations are settled down
% based on the imperative phases.
% If you don't, then you will be double interpolating in salinity!!
function [to_qc_v] = inherit_qc(from_qc_v,to_qc_v,from_data_type,to_data_type)
  qc_declarations;
  reason = sprintf('changed %s implies changed %s', from_data_type,to_data_type);
  non_good_i = find(from_qc_v ~= QC_GOOD);
  for qc_tag = sort(unique(from_qc_v(non_good_i)))'
    qc_i_v = find(from_qc_v == qc_tag);
    to_qc_v = assert_qc(qc_tag,to_qc_v,qc_i_v,reason);
  end
