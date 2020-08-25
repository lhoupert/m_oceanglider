% 
% Copyright (c) 2013 by University of Washington.  All rights reserved. Confidential
%

% like get_dive_data except it only reads nc files and puts all values on results
% and no other data sources are read.  Useful for timeseries and profile files
% from basestation but can be used for per-dive files if desired.
function [results] = get_nc_data(nc_bin_filename)
	
    results = [];
    old_format = 0; % assume the best
	try
        ncid = netcdf.open(nc_bin_filename,'NC_NOWRITE'); % this works or gives errors
        
        % BUG: if any dimname, varname, or attname starts with an '_' matlab
        % hates it and complains about bad input character (following a '.')
        % E.g., ARGO has attibutes like _FillValue
        [ndims,nvars,ngatts,unlimdimID]= netcdf.inq(ncid);
        % intern global attributes into the structure
        global_id = netcdf.getConstant('NC_GLOBAL');
        for gvarid = 0:ngatts-1 % the ids are 0-based
          attname = netcdf.inqAttName(ncid,global_id,gvarid);
          if isempty(regexp(attname(1), '[a-zA-Z]'))
            attname = ['a' attname];
          end
          attvalue = netcdf.getAtt(ncid,global_id,attname);
          % do not need to transpose strings from global attributes
          %DEAD eval(['globals.' attname ' = attvalue;']); % intern
          globals.(attname) = attvalue; % intern
        end

        if (isfield(globals,'base_station_version'))
          if (globals.base_station_version < 2.7)
            old_format = 1; % too early
          end
        else
          old_format = 1; % no base_station_version?
        end
        if (isfield(globals,'file_version'))
          if (globals.file_version < 2.7)
            old_format = 1; % too early
          end
        else
          old_format = 1; % no file_version?
        end
        if old_format
          error(sprintf('%s too old; consider --force to rebuild!\n',nc_bin_filename));
        end
        results.GLOBALS = globals; % cache for later access if desired
        
        % intern each variable into the structure
        columns = {};
        for varid = 0:nvars-1 % the ids are 0-based
          [varname vartype dimids natts] = netcdf.inqVar(ncid,varid);
          value = netcdf.getVar(ncid,varid); % into local scope

		  % .mat files write the centers, max and min log values as doubles (via read_log2)
		  % .nc files write these values as int32, per spec.  A consequence is that computation
		  % of, for examples, c_pitch_imp in the .mat version yields non-integer centers.
		  % We appear to have a choice: coerce .mat to int or .nc to double.
		  % However we are forced to choose the later because matlab has trouble 
		  % combining int32 with vector doubles (e.g., vol0 + vbd). Idiots.
          if (strcmp(class(value),'int32'))
            value = double(value);
          end
		  results.(varname) = value; % let the caller sort it out...
          %NOTYET [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
        end
		netcdf.close(ncid);
		% start_time_string = results.GLOBALS.time_coverage_start
		% start_time = datevec(unix_to_datenum(results.time(1)))
	catch ME
        fprintf(1,'%s\n',ME.message);
	end
