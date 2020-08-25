% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

% find any profile file (log, eng, nc or mat), rather than original just log/eng files (see available_runs())
% succinct_elts(available_profiles())
function [runs] = available_profiles(data_ext, id_str)
%
% data_ext - Specifric data extensions to look for ('nc','log','eng' or
% 'mat'). If none is specified, then all valid extensions will be used in
% search. Specify as cell array object: e.g. {'nc','log'}.
%
% id_str - Seaglider id #
%

%  sg_calib_constants; % billr: doesn't work in compiled code
  
  if exist('data_ext', 'var') == 0
      data_ext = {'nc','log','eng','mat'};
  end
  
  files = dir(sprintf('p%03d*.*',id_str));
  runs = [];
  for file = files'
    file_name = file.name;
    if (file.bytes > 0)
      fields = textscan(file_name,'p%3d%4d.%s');
      for i = 1:size(data_ext,2)
        if (strcmp(fields{3},data_ext{i}))
          runs = [runs ; fields{2}];
          break;
        end
      end
    end
  end
  runs = sort(unique(runs));
