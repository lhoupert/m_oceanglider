% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

function [memory] = avail_mem(gc)
  if (nargin == 0)
    gc = 0;
  end
  if (gc)
    % call gc() (this is java so this is only a hint)
    java.lang.Runtime.getRuntime.gc();
  end
  memory = java.lang.Runtime.getRuntime.freeMemory(); % bytes