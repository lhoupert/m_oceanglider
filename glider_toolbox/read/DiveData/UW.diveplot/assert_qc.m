% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%
function [qc_v] = assert_qc(qc_tag,qc_v,indices_i,reason)
  qc_declarations;
  if (qc_tag == QC_NO_CHANGE) % QC_NO_CHANGE
    return % nothing to do
  end
  trump_qc_v = [];
  switch qc_tag
   case QC_INTERPOLATED
    trump_qc_v = [QC_PROBABLY_BAD; QC_BAD];
   case QC_PROBABLY_BAD
    trump_qc_v = [QC_BAD];
   case {QC_BAD,QC_UNSAMPLED}
    % no trump
   otherwise
    fprintf(1,'No QC preference order for %s!\n', qc_tag_names{qc_tag+1});
    % no trump
  end
  trump_qc_v = [trump_qc_v ; qc_tag]; % skip those already set
  old_qc = qc_v(indices_i);
  already_set_i = [];
  for trump_qc = trump_qc_v'
    set_i = indices_i(find(old_qc == trump_qc));
    if (~isempty(set_i))
      already_set_i = [already_set_i; set_i];
    end
  end
  changed_i = setdiff(indices_i,already_set_i);
  if (~isempty(changed_i))
    qc_v(changed_i) = qc_tag;
    trace_array(sprintf('QC %s -> %d',reason,qc_tag), changed_i);
    fprintf(1,'Changed (%d/%d) %s to %s because %s\n',length(changed_i),length(qc_v),...
            succinct_elts(changed_i),qc_tag_names{qc_tag+1},reason);
  end
end