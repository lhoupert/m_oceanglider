% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

% e.g., bin_edit_dir('/Users/jsb/Seaglider/TestData/v66/fin/gliders');
function [dirs] = bin_edit_dir(starting_dir)
  cd(starting_dir) % ensure we are in the starting directory and report it...
  dirs = [];
  files = dir('sg_calib_constants.m');
  bin_file = 'Bin.mat';
  if (length(files)) % looks like a dive directory (and we need sg_calib_constants)
    [fid,message] = fopen(bin_file,'r'); % any Bin to show?
    if (fid ~= -1)
      fclose(fid); % don't need this...
      load(bin_file);
      if (size(Bin))
        runs = sort(unique(Bin(:,1))); % all available in Bin 1 = dive_num_index
        bin_fit_runs = runs;
        fprintf(1,'%s %s\n',cd(),succinct_elts(bin_fit_runs));
        bin_edit;
        if (0)
          % assuming you are only displaying the ts diagram (see bin_edit and premature return...)
          % retitle the single figure
          title(cd()) 
          print;close
        end
      else
        fprintf(1,'%s NO DIVES!\n',cd());        
      end
    end
  end

  files = dir('*');
  for file = files'
    if (strcmp(file.name,'.') || strcmp(file.name,'..'))  
      continue;
    end
    if (file.isdir)
      subdirs = bin_edit_dir(sprintf('%s/%s',starting_dir,file.name));
      % dirs = [dirs; subdirs];
    end
  end
