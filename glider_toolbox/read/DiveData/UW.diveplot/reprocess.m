% 
% Copyright (c) 2013, 2014 by University of Washington.  All rights reserved. Confidential
%

% Run basestation reprocess script to reprocess all dive files
% This blocks and can take a long time...
dives = ask_which_runs(available_profiles());
% do not use succinct_elts since basestation code can't process the result
% NOTE: this relies on matlab's conversion of a list [1 2 3] into a string *without brackets* in sprintf!!
system(sprintf('/usr/local/basestation/reprocess.sh %s',dives));
