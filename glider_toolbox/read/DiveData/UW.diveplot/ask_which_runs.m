% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%

% runs passed are assumed to be 'all available runs'
% script returns a potentially modified list of runs
function [runs] = ask_which_runs(runs)
  xruns = input('Which dives? (CR for all) ','s');
  if (length(xruns) ~= 0)
    if (xruns(1) == '!')
      xruns = xruns(2:end); % strip it...
      sruns = evalin('caller',xruns); % eval and use the results, if possible
    else
      sruns = str2num(xruns);
    end
    sruns = sruns';
  else
    sruns = runs; % use them all
  end
  runs = sort(unique(intersect(sruns,runs)));
