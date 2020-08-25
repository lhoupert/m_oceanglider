% 
% Copyright (c) 2006-2011, 2013 by University of Washington.  All rights reserved. Confidential
%

% attempts to 'succinctly' describe an array of runs or indices in matlab
function [xruns] = succinct_elts(runs);
  runs = sort(unique(runs)); % make canonical
  xruns = '';
  prefix = '';
  [m,n] = size(runs);
  if (n == 1)
    n = m;
    runs = runs';
  end
  if n
    last_i = 1;
    for break_i = [find(diff(runs) > 1)';length(runs)]'
      nruns =  runs(break_i) - runs(last_i);
      if (nruns == 0)
        xruns = sprintf('%s%s%d',xruns,prefix,runs(break_i));
	  elseif (nruns == 1)
        xruns = sprintf('%s%s%d %d',xruns,prefix,runs(last_i),runs(break_i));
      else
        xruns = sprintf('%s%s%d:%d',xruns,prefix,runs(last_i),runs(break_i));
      end
      last_i = break_i+1;
      prefix = ' ';
    end
  end
