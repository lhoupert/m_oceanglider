% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

% dive numbers of paired log and eng files
% succinct_elts(available_runs())
function [runs] = available_runs(id_str,report_mismatch)
  if (nargin == 0 | id_str == 0)
    sg_calib_constants;
  end
  if (nargin < 2)
    report_mismatch = 0;
  end
  files = dir(sprintf('p%s*.*',id_str));
  log_files = [];
  eng_files = [];
  for file = files'
    file_name = file.name;
    if (file.bytes > 0) % sometimes log and eng files can be zero bytes...
      fields = sscanf(file_name,'p%3d%4d.%s');
      % fields 3:end are the ascii encodings of each character
      if (strfind(file_name,'log'))
        log_files = [log_files ; fields(2)];
      elseif (strfind(file_name,'eng'))
        eng_files = [eng_files ; fields(2)];
      end
    end
  end
  if (report_mismatch)
    % report missing files
    eng_missing_log = setdiff(eng_files,log_files);
    if (~isempty(eng_missing_log))
      fprintf(1,'Eng files missing log files: %s\n',mat2str(eng_missing_log));
    end
    log_missing_eng = setdiff(log_files,eng_files);
    if (~isempty(log_missing_eng))
      fprintf(1,'Log files missing eng files: %s\n',mat2str(log_missing_eng));
    end
  end
  
  runs = intersect(log_files, eng_files)';
  runs = setdiff(runs,[0]); % skip dive 0
