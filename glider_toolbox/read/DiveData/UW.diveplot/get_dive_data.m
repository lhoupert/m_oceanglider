% 
% Copyright (c) 2006-2014 by University of Washington.  All rights reserved. Confidential
%

% the protocol:

% get_dive_data() returns loginfo, eng, and results structures with 'raw'
% data and results from nc files only. Consumers should NOT modify these
% structures but copy into variables to manipulate.  Much of this copying
% into variables happens in unpack_data, which see.

% NOTE: loginfo, NOT log, since that kills the 'log' function!
function [loginfo,eng,results] = get_dive_data(id_dive,require_results)
% works somewhat like python's load_dive_profile_data() 
% except dumps data in matlab structures and in caller environment
% assumes basename looks like pGGGYYYY and we are cd'd to the proper dir
	% id_dive can be a string like 'GGGYYYY' or a dive integer YYYY
		
	% require_results can take 2 values:
	% load raw data and results from nc file, then
	% 1 - make no additional updates from sgc
	% 2 - apply updates from sgc if more-recent or changed, documenting the changes on the history

	if (nargin < 2)
      require_results = 1; % require loading from nc file
    end
	if (require_results == 0)
		require_results = 1; % in case of old calls
	end
    
	if (~ischar(id_dive)) % just an integer dive number?
		% determine the basename from the first to match?
		% no, avoid pt files
		basename = sprintf('p*%04d.*',id_dive);
		files = dir(basename);
		if ~isempty(files)
			for file = files
				basename = files(1).name;
				if (basename(2) == 't')
					continue; % skip self-test files
				end
				ext_i = strfind(basename,'.');
				if ext_i
					basename = basename(1:ext_i-1);
					break;
				end
			end
		end
	else
		basename = sprintf('p%s',id_dive);
    end

    sg_calib_constants_file = 'sg_calib_constants.m';
    sgc_fileinfo = dir(sg_calib_constants_file);
    sgc_exists = (size(sgc_fileinfo,1) == 1);

    gz_filename = sprintf('%s.nc.gz',basename);
    gz_fileinfo = dir(gz_filename);
    if (size(gz_fileinfo,1) == 1)
      system(sprintf('gunzip %s',gz_filename)); % should yield nc_filename
    end

    nc_filename = sprintf('%s.nc',basename);
    nc_fileinfo = dir(nc_filename);
    nc_exists = (size(nc_fileinfo,1) == 1);
    matlab_versions = ver();
    min_required_version = 7.8; % aka r2009a
    if (nc_exists & matlab_versions(1).Version < min_required_version)
      fprintf(1,'Matlab version %.1f required to read %s -- ignoring file.\n',min_required_version,nc_filename);
      nc_exists = 0; % Like it wasn't even there...
    end
    
    loginfo = [];
    results = [];
    eng = [];
    old_format = 0; % assume the best
    % no matter what, and before we evaluate unpack_mat or unpack_nc, evaluate this file
    evalin('caller','qc_declarations');
	results.basename = basename;
    if (nc_exists)
      try
        % a version of ncload that tries to reconstruct the equivalent log and eng files
        % BUT note for gpctd and scicon data the timevars and data fields are NOT the same length
        % let caller decimate to same length if desired
        
        ncid = netcdf.open(nc_filename,'NC_NOWRITE'); % this works or gives errors
        
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
          error(sprintf('%s too old; consider the --force basestation option to rebuild from log and eng files!\n',nc_filename));
        end

        dive_start_time = globals.start_time;
        loginfo.log_epoch_start = dive_start_time;
        start_sdn = unix_to_datenum(dive_start_time);
        loginfo.log_dn = start_sdn;
        start_vec = datevec(start_sdn);
        loginfo.year   = start_vec(1);
        loginfo.month  = start_vec(2);
        loginfo.date   = start_vec(3);
        loginfo.hour   = start_vec(4);
        loginfo.minute = start_vec(5);
        loginfo.second = start_vec(6);
        loginfo.log_file = sprintf('%s.log',basename); % this may or may not exist
        logfinfo.version = globals.seaglider_software_version;
        
        eng.version = globals.seaglider_software_version;
        eng.basestation_version = globals.base_station_version;
        eng.glider = globals.glider;
        eng.mission = globals.mission;
        eng.dive = globals.dive_number;
        eng.start = [loginfo.month, loginfo.date, loginfo.year, loginfo.hour, loginfo.minute, loginfo.second];
        
        results.GLOBALS = globals; % cache for later access if desired
        sgc = []; % collect the sg_calib_constants on a structure in addition to spreading into the caller space
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
          if (isa(value,'int32'))
            value = double(value);
          end
          if (~isempty(regexp(varname,'^eng_')))
            column = varname(5:end);
            columns{end+1} = column;
            %DEAD eval(['eng.' column ' = value;']);
            eng.(column) = value;
          elseif (~isempty(regexp(varname,'^gc_')))
            % put the various GC columns on loginfo; see read_log2
            %DEAD eval(['loginfo.' varname ' = value;']);
            loginfo.(varname) = value;
          elseif (~isempty(regexp(varname,'^log_')))
            tag = varname(5:end);
            if isempty(regexp(tag(1), '[a-zA-Z]'))
              tag = ['a' tag];
            end
            value = value'; % in case of string; NOP if constant
            if (isa(value,'char'))
              % follow read_log2 protocol and break strings with , into an array of values
              % but note that some strings ($GPS1) still have their tags attached, ugh...
              % strip those first
              if (value(1) == '$')
                [tag, remainder] = strtok(strtok(value, '$'), ',');
              else
                remainder = value;
              end
              % Read the values associated with the tag. 
              % NOTE: none of these tags are prefixed with _ so we never have an 'a_tag' tag
              exception_tags = {'TGT_NAME', 'GCHEAD', 'DEVICES', 'SENSORS'};
              exception_fmts = {'%s', '%s', '%s', '%s'};
              
              ii = 0;
              values = [];
              while (~isempty(remainder))
                [tmp, remainder] = strtok(remainder, ',');
                ii = ii + 1;
                jj = strmatch(tag, exception_tags, 'exact');
                if isempty(jj)
                  scan = sscanf(tmp, '%f');
                  if ~isempty(scan)
                    values(ii) = scan;
                  end
                else
                  values{ii} = sscanf(tmp, char(exception_fmts{jj}));
                end
              end
              value = values;
            end
            %DEAD eval(['loginfo.' tag ' = value;']);
            loginfo.(tag) = value;
          elseif (~isempty(regexp(varname,'^sg_cal_')))
            value = value'; % in case of string; NOP if constant
			variable = varname(8:end);
			sgc.(variable) = value;
            evalin('caller', [variable ' = ' mat2str(value) ';']); % use value from local scope
          else
            % gps_ structures go here
            %DEAD eval(['results.' varname ' = value;']); % let the caller sort it out...
            results.(varname) = value; % let the caller sort it out...
          end
          eng.columns = columns;
          %NOTYET [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
        end
      catch
        % be silent
      end
	  if (isfield(results.GLOBALS,'date_modified'))
		  % when results were last computed on the basestation, not when it was copied
		  % also not when the file was first created (date_created) -- we want to know when sg_calib_constants were incorporated last
		  results_datenum = datenum(results.GLOBALS.date_modified,'yyyy-mm-ddTHH:MM:SS');
	  else
		  % This is the local time of the nc file, which could be its copy time...
		  results_datenum = nc_fileinfo.datenum;
	  end
	  sgc.fileinfo = nc_fileinfo;
	  sgc.fileinfo.datenum = results_datenum;
	  results.sg_calib_constants = sgc; % return the constants
      netcdf.close(ncid);
      
      if (~old_format)
        % from read_eng()
        if isfield(eng, 'depth') % ensure we got this field from NC land
          eng.depth = eng.depth / 100.0;
        end
        if (isfield(results,'processing_error') || isfield(results,'skipped_profile'))
          %DEBUG fprintf(1,'Processing error or skipped profile reported from NC file\n');
        else
          % unpack_nc
        end
      end
      
    else
      fprintf(1,'No nc file available for %s!\n',basename);
      results.processing_error = -1; % alert caller
      return
    end

    % update (again) sg_calib_constants, if needed
	if (require_results == 2 & ~isempty(results))
		% Loading a local sg_calib_constants file is a complicated matter:
		% First, some callers don't want it loaded (just the basestation facts ma'am);
		% they should use require_results = 1
		% But for those that do (e.g, diveplot, since pilots use it to check if
		% vbd regressions will have good effects) they should use require_results = 2
		% Even then we only want to load it if there were actual changes.
		% Checking file dates is insufficient....
		if (sgc_exists)
			changed_hx = '';
			sgcf = read_sg_calib_constants();
			sgcr = results.sg_calib_constants;
			changed = 0; % assume no changes
			if (sgcf.fileinfo.datenum > sgcr.fileinfo.datenum)
				% TODO report the dates?
				changed = 'more-recent';
			else

				% The local sg_calib_constants.m file is older than the one
				% used to create/modify the nc file last.  This can happen
                % if the basestation version of sgc is even older but gets
                % incorporated into nc files at a later date. 

				% The local file probably has changes that have not been
                % broadcast yet or it could really be the older file.  Since
                % they requested updated results (2) let them know about any
                % changes

				for field = fieldnames(sgcf)'
                    field = field{1};
                    if (strcmp(field,'fileinfo'))
                       continue; % skip the date info
                    end
					if (isfield(sgcr,field))
						vr = sgcr.(field);
						vf = sgcf.(field);
						same = 0;
						if ischar(vr)
							same = strcmp(vr,vf);
                        elseif isstruct(vr) % what are these doing here?
                            same = 1; % can't compare these; assume no change
						else
							same = (vr == vf);
						end
						if (~same)
							changed_hx = sprintf('%sVariable %s changed from %g to %g in %s\n',...
												 changed_hx,field,vr,vf,sg_calib_constants_file);
							changed = 'changed';
						end
					else
						changed_hx = sprintf('%sVariable %s added to %s\n',...
											 changed_hx,field,sg_calib_constants_file);
						changed = 'changed';
					end
				end
			end
			if changed
				changed_hx = sprintf('%sApplying %s %s dated %s...\n',...
									 changed_hx,changed,sg_calib_constants_file,datestr(sgcf.fileinfo.datenum));
				evalin('caller','sg_calib_constants');
				results.GLOBALS.history = sprintf('%s%s',results.GLOBALS.history,changed_hx);
			end
		end
	end
    
    evalin('caller','sg_config_constants'); % update these constants based on prevailing sg_calib_constant values

    % At this point we know we have a good loginfo and eng structure
    % Finally, fixup various log GPS entries
    
    % modified version of the tail end of read_log2
    % but without the recursive call to read_log2
    % this duplicates computation from read_log2 but otherwise handles 
    % putting GPS info into matlab usable format if coming from nc file
    report_midnight_crossing = 0; % CONTROL DEBUG
    % 'modern' GPS structures look like:
    % DDMMYY,HHMMSS,lat,lon,first_fix_time,hdop,final_fix_time,magvar
    serial_date_time = datenum(loginfo.year, loginfo.month, loginfo.date, ...
                               loginfo.hour, loginfo.minute, loginfo.second);
    loginfo.log_dn = serial_date_time; % save this...
    loginfo.log_epoch_start = datenum_to_unix(loginfo.log_dn);
    dive_mo_str = datestr(serial_date_time, 5); % month
    dive_da_str = datestr(serial_date_time, 7); % day
    dive_yr_str = datestr(serial_date_time, 11);% year(YY)
    gps_yr = str2num(dive_yr_str);
    gps_mo = str2num(dive_mo_str);
    gps_da = str2num(dive_da_str);
    loginfo.gps_date = gps_da*1e4 + gps_mo*100 + gps_yr; % in GPS date format (DDMMYY)
    
    % hdop = horizontal degree of precision
    % error in meters ~= 3.04*hdop + 3.57
    % 4 =>  ~16m error
    % 10 => ~35m error
    % 99 indicates no valid fix acquired (or no response to $GPGGA query)
    valid_hdop = 10; % less than this is ok
    loginfo.old_style_GPS = 0; % assume GPS supplies dates as well as times
    if isfield(loginfo, 'GPS1')
      if length(loginfo.GPS1) < 5
        % ancient style: HHMMSS,lat,lon,secs
        loginfo.old_style_GPS = 1;
        loginfo.GPS1 = [loginfo.gps_date;loginfo.GPS1';2.0;loginfo.GPS1(end);0]'; % 6:hdop;7:final_fix_secs;8:magvar
      end
      if length(loginfo.GPS1) < 8
        loginfo.old_style_GPS = 1;
        loginfo.GPS1 = [loginfo.gps_date;loginfo.GPS1']'; % convert to new style
      end
      t = sscanf(sprintf('%06d', loginfo.GPS1(2)), '%02d%02d%02d');
      d = sscanf(sprintf('%06d', loginfo.GPS1(1)), '%02d%02d%02d'); % DDMMYY
      loginfo.GPS1_dn = datenum([2000 + d(3), d(2), d(1), t(1), t(2), t(3)]);
      loginfo.GPS1_t = datenum_to_unix(loginfo.GPS1_dn);
      degrees = fix(loginfo.GPS1(3)/100);
      minutes = rem(loginfo.GPS1(3),100);
      loginfo.GPS1_lat = degrees + minutes / 60.0;
      loginfo.GPS1_lat_deg = degrees;
      loginfo.GPS1_lat_min = minutes;
      
      degrees = fix(loginfo.GPS1(4)/100);
      minutes = rem(loginfo.GPS1(4), 100);
      loginfo.GPS1_lon = degrees + minutes / 60.0;
      loginfo.GPS1_lon_deg = degrees;
      loginfo.GPS1_lon_min = minutes;
      
      loginfo.GPS1_sec = loginfo.GPS1(5);
      loginfo.GPS1_hdop = loginfo.GPS1(6);
      % old data sets didn't record magvar (or even final_fix_secs)
      if length(loginfo.GPS1) >= 7
        loginfo.GPS1_magvar = loginfo.GPS1(end);
      end
      loginfo.GPS1_valid = (loginfo.GPS1_hdop < valid_hdop & loginfo.GPS1_sec > 0);
    end
    
    if isfield(loginfo, 'GPS2')
      if length(loginfo.GPS2) < 5
        % ancient style
        loginfo.old_style_GPS = 1;
        loginfo.GPS2 = [loginfo.gps_date;loginfo.GPS2';2.0;loginfo.GPS2(end);0]'; % 6:hdop;7:final_fix_secs;8:magvar
      end
      if length(loginfo.GPS2) < 8
        loginfo.old_style_GPS = 1;
        loginfo.GPS2 = [loginfo.gps_date;loginfo.GPS2']'; % convert to new style
      end
      t = sscanf(sprintf('%06d', loginfo.GPS2(2)), '%02d%02d%02d');
      d = sscanf(sprintf('%06d', loginfo.GPS2(1)), '%02d%02d%02d'); % DDMMYY
      loginfo.GPS2_dn = datenum([2000 + d(3), d(2), d(1), t(1), t(2), t(3)]);
      loginfo.GPS2_t = datenum_to_unix(loginfo.GPS2_dn);
      if (loginfo.old_style_GPS)
        % In old log files we didn't record both date and time for the GPS reading
        % instead we recorded the HHMMSS wrt to the start time of the file,
        % which was sync'd w/ the GPS.  But it could happen that we surfaced
        % just before midnight, took GPS1, did surface stuff, then took GPS2
        % (which could still be before midnight! -- see ak/nov03 dive 137),
        % then opened and tagged the new log file start time that was after
        % midnight.  In this case we would see GPS1 as apparently later than
        % the log time (and perhaps also see GPS2 time ahead of log time) by
        % nearly 24hrs.  Back correct GPS1 and perhaps GPS2 time back by 1
        % day's worth of seconds,
        
        % datenum gives fraction of days new format GPS supplies the actual date
        % avoiding this bug
        
        if (loginfo.GPS1_dn > loginfo.log_dn)
          if (report_midnight_crossing && initial_call)
            fprintf(1,'%s Correcting GPS1 for midnight crossing.\n',log_file); % DEBUG
          end
          t = sscanf(sprintf('%06d', loginfo.GPS1(2)), '%02d%02d%02d'); % get the GPS1 time back
          loginfo.GPS1_dn = datenum([2000 + d(3), d(2), d(1), t(1), t(2), t(3)])-1.0;
          loginfo.GPS1(1) = str2num(datestr(loginfo.GPS1_dn,'mmddyy')); % y2k?
          loginfo.GPS1_t = datenum_to_unix(loginfo.GPS1_dn);
        end
        if (loginfo.GPS2_dn > loginfo.log_dn)
          if (report_midnight_crossing && initial_call)
            fprintf(1,'%s Correcting GPS2 for midnight crossing.\n',log_file); % DBEUG
          end
          t = sscanf(sprintf('%06d', loginfo.GPS2(2)), '%02d%02d%02d'); % get the GPS2 time back
          loginfo.GPS2_dn = datenum([2000 + d(3), d(2), d(1), t(1), t(2), t(3)])-1.0;
          loginfo.GPS2(1) = str2num(datestr(loginfo.GPS2_dn,'mmddyy')); % y2k?
          loginfo.GPS2_t = datenum_to_unix(loginfo.GPS1_dn);
        end
      end
      
      degrees = fix(loginfo.GPS2(3)/100);
      minutes = rem(loginfo.GPS2(3),100);
      loginfo.GPS2_lat = degrees + minutes / 60.0;
      loginfo.GPS2_lat_deg = degrees;
      loginfo.GPS2_lat_min = minutes;
      
      degrees = fix(loginfo.GPS2(4)/100);
      minutes = rem(loginfo.GPS2(4), 100);
      loginfo.GPS2_lon = degrees + minutes / 60.0;
      loginfo.GPS2_lon_deg = degrees;
      loginfo.GPS2_lon_min = minutes;
      
      loginfo.GPS2_sec = loginfo.GPS2(5);
      loginfo.GPS2_hdop = loginfo.GPS2(6);
      % 7 is final_fix_sec
      if length(loginfo.GPS2) >= 7
        loginfo.GPS2_magvar = loginfo.GPS2(end);
      end
      loginfo.GPS2_valid = (loginfo.GPS2_hdop < valid_hdop & loginfo.GPS2_sec > 0);
      loginfo.GPS2_valid = (loginfo.GPS2_valid & (loginfo.GPS2_dn > loginfo.GPS1_dn)); % GPS2 changed from GPS1?
    end
    
    if isfield(loginfo, 'GPS')
      % this might be a new GPS structure OR it might have come from reading ahead in read_log2
      % can't assume it is always a correct, new GPS structure!
      if (~isfield(loginfo,'GPS_dn'))
        t = sscanf(sprintf('%06d', loginfo.GPS(2)), '%02d%02d%02d');
        d = sscanf(sprintf('%06d', loginfo.GPS(1)), '%02d%02d%02d'); % DDMMYY
        loginfo.GPS_dn = datenum([2000 + d(3), d(2), d(1), t(1), t(2), t(3)]);
        loginfo.GPS_t = datenum_to_unix(loginfo.GPS_dn);
      end
      degrees = fix(loginfo.GPS(3)/100);
      minutes = rem(loginfo.GPS(3),100);
      loginfo.GPS_lat = degrees + minutes / 60.0;
      loginfo.GPS_lat_deg = degrees;
      loginfo.GPS_lat_min = minutes;
      
      degrees = fix(loginfo.GPS(4)/100);
      minutes = rem(loginfo.GPS(4),100);
      loginfo.GPS_lon = degrees + minutes / 60.0;
      loginfo.GPS_lon_deg = degrees;
      loginfo.GPS_lon_min = minutes;
      
      loginfo.GPS_hdop = loginfo.GPS(6);
      loginfo.GPS_sec = loginfo.GPS1(5);
      if length(loginfo.GPS) >= 7
        loginfo.GPS_magvar = loginfo.GPS(end);
      end
      loginfo.GPS_valid = (loginfo.GPS_hdop < valid_hdop & loginfo.GPS_sec > 0);
      loginfo.GPS_valid = (loginfo.GPS_valid & (loginfo.GPS_dn > loginfo.GPS2_dn)); % GPS changed from GPS2?
    end

    % from read_eng()
    if isfield(eng, 'start')
      eng.start_ts_dn = datenum(eng.start(3),eng.start(1),eng.start(2),eng.start(4),eng.start(5),eng.start(6));
      eng.start_ts = datenum_to_unix(eng.start_ts_dn);
    end
    if ~isfield(eng,'start_ts_dn')
      % this is close but off by several seconds
      eng.start_ts_dn = loginfo.GPS2_dn;
      eng.start_ts = loginfo.GPS2_t;
    end
    
    % Now, wouldn't it be nice to evalin('caller','unpack_data')
    % to assert the variables CCE's scripts typically use but
    % we can't since we haven't returned log, eng, and results yet.
    % our other alternative would be to NOT return them
    % but eval each of the expressions in caller
    % equally ugly...

    if (require_results && isempty(results))
      results.processing_error = 1; % let caller know
    end
    
    