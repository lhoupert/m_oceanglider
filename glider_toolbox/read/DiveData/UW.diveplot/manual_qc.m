% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%
function [qc_v,indices] = manual_qc(directives,fn,assertion_tag,qc_tag,qc_v,data_tag)
  o_indices = find(qc_v == qc_tag);
  eval(sprintf('directives.%s = o_indices;', assertion_tag)); % python: setattr(directives,assertion_tag,o_indices)
  indices = drv_eval(directives,fn);
  changed_i = setdiff(o_indices,indices);
  if (~isempty(changed_i))
    qc_v = assert_qc(1,qc_v,changed_i,sprintf('%s QC reset manually',data_tag)); % hardwired QC_GOOD
  end
  changed_i = setdiff(indices,o_indices);
  if (~isempty(changed_i))
    qc_v = assert_qc(qc_tag,qc_v,changed_i,sprintf('%s QC set manually',data_tag));
  end
  % since we don't return directives, this is useless
  %DEAD eval(sprintf('directives.%s = indices;', assertion_tag)); % update
end
