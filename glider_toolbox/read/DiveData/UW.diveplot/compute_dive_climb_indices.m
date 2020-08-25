% 
% Copyright (c) 2006-2012, 2014 by University of Washington.  All rights reserved.
%
% This file contains proprietary information and remains the 
% unpublished property of the University of Washington. Use, disclosure,
% or reproduction is prohibited.
% assumes you have run unpack_mat/unpack_nc

% some indices in sg space
sg_start_of_climb_i = find(sg_time >= start_of_climb_time,1,'first');
sg_dive_i = [1:sg_start_of_climb_i - 1]';
sg_climb_i = [sg_start_of_climb_i:sg_np]';

% all these values are in ctd space
uncorrectable_i = bad_qc(salin_qc);
valid_i = setdiff([1:ctd_np]',uncorrectable_i);
interpolated_i = find(salin_qc == QC_INTERPOLATED);

dive_i = [1:start_of_climb_i - 1]';
dive_i_corrected = setdiff(dive_i,uncorrectable_i);
dive_i_TS = intersect(dive_i_corrected,interpolated_i); % where we interpolated in TS

climb_i = [start_of_climb_i:ctd_np]';
climb_i_corrected = setdiff(climb_i,uncorrectable_i);
climb_i_TS = intersect(climb_i_corrected,interpolated_i); % where we interpolated in TS

