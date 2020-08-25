% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved.
%
% This file contains proprietary information and remains the 
% unpublished property of the University of Washington. Use, disclosure,
% or reproduction is prohibited.

%********************************************************************
% A structure-based Seaglider log reader.
% Returns the structure 'loginfo' with elements names identical to
% the tags in the log file with the following modifications:
% - The leading '$' has been removed.
% - Tags beginning with characters not in [a-zA-Z] 
%   (e.g. '_' and numbers) have a leading character appended.
% Thus, if pad_char = 'a': $DIVE -> loginfo.DIVE, 
% $24V_AH -> loginfo.a24_AH and $_SM_ANGLEo -> loginfo.a_SM_ANGLEo
% Most structure elements will be single values or vectors.
% The exceptions are those listed in 'multi_tags'. These are tags that
% have multple instances in the log file ($GC is the orginal example).
% The elements associated with these tags will be cell arrays. Thus,
% loginfo.GC{i} = <vector of GC values for the ith occurance of the 
% $GC tag in the logfile.
% Some sample output...
% loginfo.KALMAN_Y: [1.0386e+004 796.6000 -625.3000 -1.7215e+004 999.6000]
% loginfo.MHEAD_RNG_PITCHd_Wd: [179 186617 -13.6000 -9.9800]
% loginfo.D_GRID: 2031
% loginfo.T_DIVE2: 167
% loginfo.GC: {1x40 cell}
% loginfo.GC{2}: [54.0000 -1.0600 -146.1000 ... 0.1840 0.1050]
% loginfo.HUMID: 1834
% loginfo.a24V_AH: [23.7000 5.4850]
%
% cml 14 April 2004

% extended to deal with old-style GPS formats and missing GPS fields
function loginfo = read_log2 (log_file, varargin)
  initial_call = (size(varargin,2) == 0); % is this a lookahead call?
  report_midnight_crossing = 0; % CONTROL DEBUG
  % Character for padding troublesome tag names.
  pad_char = 'a';
  
  % A cell array of tag names that can occur more than once in a log file.
  
  multi_tags = {'GC', 'RAFOS'};
  exception_tags = {'TGT_NAME', 'GCHEAD', 'DEVICES', 'SENSORS'};
  exception_fmts = {'%s', '%s', '%s', '%s'};
  
  %*** Need to work in parsing for lat/lon locations.
  
  %---------------------------------------------------------------
  % Bookeeping for multi-lne tags. 
  % 'mc_count' will be used to index instances.
  
  n_mc = length(multi_tags);
  mc_count = zeros(n_mc,1);

  % Open the log file.
  
  fp = fopen(log_file, 'r');
  if fp == -1
    % non-existent log file
    loginfo = struct([]);
    return
  end
  
  line = fgetl(fp);
  if (line == -1)
    % empty log file
    fclose(fp);
    loginfo = struct([]);
    return
  end
  loginfo.log_file = log_file; % save original log filename...

  if strncmp(line, 'version:', 8) == 1
    fclose(fp);
    [hdr, fp] = read_header(log_file);	% read_header opens the file and returns the gently used file pointer
    if isfield(hdr, 'start')
      loginfo.month = hdr.start(1);
      loginfo.date = hdr.start(2);
      loginfo.year = hdr.start(3);
      loginfo.hour = hdr.start(4);
      loginfo.minute = hdr.start(5);
      loginfo.second = hdr.start(6);
    end
    if isfield(hdr, 'version')
      loginfo.version = hdr.version;
    end
  else
    if (length(strfind(line,'$SM_CCo')))
      % for one deployment ps/aug04/sg016 we wrote the old SM_CC data afer we opened the new logfile
      % drop it for the moment, parse it and put on another tag and copy it below
      line = fgetl(fp);
    end
    tmp = sscanf(line, '%2d %2d %3d %2d %2d %2d');
    
    if length(tmp) == 6 & tmp(1) > 0 & tmp(1) <= 12 & tmp(2) > 0 & tmp(2) <= 31
      loginfo.month = tmp(1);
      loginfo.date = tmp(2);
      loginfo.year = 1900 + tmp(3);
      loginfo.hour = tmp(4);
      loginfo.minute = tmp(5);
      loginfo.second = tmp(6);
    end
  end
  if ~isfield(loginfo, 'version')
    loginfo.version = 62; % ancient
  end
  
  % Assumes that log files begin with a 6-integer time stamp of the
  % form: mm dd yyy (+1900) hh mm ss (UTC). 
  % Read and assign the time stamp.
  
  tag = 0;
  while (tag ~= -1)
    
    
    line = fgetl(fp);
    
    if line == -1
      break;
    end
    
    % if no leading $ then this is probably the date line
    
    if line(1) ~= '$' | strncmp(line, '$ESCAPE', 7) == 1
      continue;
    end
    
    % Strip the leading '$'.
    
    [tag, remainder] = strtok(strtok(line, '$'), ',');
    
    % Read the values associated with the tag. 
    
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
    if (tag ~= -1)
      
      % If the tag begins with something other than a letter, pad
      % the name with 'pad_char' before defining the structure element.
      % MATLAB variables must begin with [a-zA-Z].
      
      if isempty(regexp(tag(1), '[a-zA-Z]'))
        tag = [pad_char tag];
      end
      
      % Is the tag on the multi-line tag list?
      
      jj = strmatch(tag, multi_tags, 'exact');
      if isempty(jj)
        eval (['loginfo.' tag ' = values;']);
      else
        % Special case where there may be multiple lines of the same
        % tag (e.g. $GC)...     
        mc_count(jj) = mc_count(jj) + 1;
        %      eval (['loginfo.' tag '{' num2str(mc_count(jj)) '} = values;']);
        if mc_count(jj) == 1
          eval (['loginfo.' tag '(' num2str(mc_count(jj)) ',:) = values;']);
        else
          base_cols = length( eval (['loginfo.' tag '(1,:);']) );
          cols = length(values);
          if cols == base_cols
            eval (['loginfo.' tag '(' num2str(mc_count(jj)) ',:) = values;']);
          end
        end
      end
      
    end
  end
  
  fclose(fp);
  
  %**********************************************************************
  % some useful derived quantities here -- make sure to protect with isfield
  

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

  % drop GC and replace with gc_XXX entries from GCHEAD
  % if not GCHEAD, use default
  if (~isfield(loginfo,'GCHEAD'))
    loginfo.GCHEAD = {'st_secs','pitch_ctl','vbd_ctl','depth','ob_vertv','data_pts','end_secs','pitch_secs','roll_secs','vbd_secs','vbd_i','gcphase'};
  end
  for f = 1:size(loginfo.GCHEAD,2)
    fn = loginfo.GCHEAD{f};
    eval(sprintf('loginfo.gc_%s = loginfo.GC(:,%d);',fn,f));
  end
  loginfo = rmfield(loginfo,{'GC'});
  % nc saves these in epoch times (the other secs are elapsed times for motors)
  loginfo.gc_st_secs = loginfo.gc_st_secs   + loginfo.log_epoch_start;
  loginfo.gcend_secs = loginfo.gc_end_secs + loginfo.log_epoch_start;
  
  % hdop = horizontal degree of precision
  % error in meters ~= 3.04*hdop + 3.57
  % 4 =>  ~16m error
  % 10 => ~35m error
  % 99 indicates no valid fix acquired (or no response to $GPGGA query)
  valid_hdop = 10; % less than this is ok
  loginfo.missing_GPS = ~isfield(loginfo, 'GPS');
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

  if (loginfo.missing_GPS && initial_call)
    % This is an old-style logfile which did not have a $GPS field entry
    % try to get it from the next logfile
    % assume the log_file name is 'pNNNDDDD.log'
    run = str2num(log_file(5:8));
    next_logfile = sprintf('%s%04d.log',log_file(1:4),run+1);
    next_loginfo = read_log2(next_logfile,1); % lookahead call, but only once
    % next_loginfo could be an empty structure if file does not exist
    if (~isfield(next_loginfo,'GPS1'))
      % DEBUG fprintf(1,'Unable to read GPS1 field from %s!\n',next_logfile);
      % Compose an untrusted location based on last fix for this dive...
      % this is what the glider does anyway...
      GPS = loginfo.GPS2; % punt w/ our last location
      GPS(5) = loginfo.T_GPS*60 + 1; % 'timed out' or -99
      GPS_dn = loginfo.GPS2_dn;
      GPS_t = loginfo.GPS2_t;
    else
      GPS = next_loginfo.GPS1;
      GPS_dn = next_loginfo.GPS1_dn;
      GPS_t = next_loginfo.GPS1_t;
    end
    % add a date value, time value already exists....
    loginfo.GPS = GPS;
    loginfo.GPS_dn = GPS_dn;
    loginfo.GPS_t = GPS_t;
    % copy SM_CCo here if you preserve it above
    % fallthrough
  end
  
  if isfield(loginfo, 'GPS')
    % this might be a new GPS structure OR it might have come from reading ahead
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
  