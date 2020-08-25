%
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%
function [indices] = bad_qc(qc_v)
  qc_declarations;
  % include both truely and probably bad points plus the unsampled points
  indices = find(qc_v == QC_BAD | qc_v == QC_PROBABLY_BAD | qc_v == QC_UNSAMPLED);

