% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved.
%
% This file contains proprietary information and remains the 
% unpublished property of the University of Washington. Use, disclosure,
% or reproduction is prohibited.

% read the pair of gpctd files
% returns a struct w/ serveral fields
function [gpctd] = read_gpctd(id_dive)
% read commented headers
  file_tags = {'a','b'};
  first_epoch_time = 0;
  sampling_interval = 1;
  pressure_i = 0; % index of Pressure column
  cond_i = 0; % index of Conductivity column
  gpctd.Time = [];
  unix_start_dn = datenum([1970, 1, 1, 0, 0, 0]); % for computing secs since unix epoch
  for i = 1:size(file_tags,2)
    file_name = sprintf('ppc%s%s.eng',id_dive,char(file_tags{i}));
    fp = fopen(file_name, 'r');
    if fp == -1
      fprintf(1,'Unable to find %s!\n',file_name);
      continue;
    end

    % parse the header and find the i'th start time and sample interval
    this_epoch_time = 0;
    while 1,
      line = fgetl(fp);
      if line == -1
        break;
      end
      line = strtok(line, '%'); % strip leading % sign 
      [tag, value] = strtok(line, ' ');
      if isempty(tag)
        break; % a data line
      end
      tag = char(tag);
      switch char(tag)
       case 'cast'
        values = textscan(char(value),'%d %d %s %d %d:%d:%d samples %d to %d, int = %d, %s');
        cast_num = values{1};
        epoch_time = datenum(sprintf('%d %s %d %d:%d:%d',values{2},char(values{3}),values{4},values{5},values{6},values{7}));
        % epoch_time = epoch_time - unix_start_dn; % bring into alignment with read_log2 times
        if (cast_num == 1)
          first_epoch_time = epoch_time;
        end
        if (cast_num == i)
          this_epoch_time = epoch_time;
          sampling_interval = eval(num2str(values{10})); % gotta be a better way
        end
       
       case 'columns:'
        [tmp, remainder] = strtok(value, ',');
        ii = 0;
        columns = [];
        while ~isempty(remainder) & ~isempty(tmp)
          ii = ii + 1;
          columns{ii} = char(tmp);
          [tmp, remainder] = strtok(remainder, ',');
        end
        if length(tmp) > 0
          ii = ii + 1;
          columns{ii} = char(tmp);
        end
        % prepare tag names and find 'Pressure', if any
        for jj = 1:ii
          tmp = columns{jj};
          [s,v] = strtok(tmp, '.');
          if ~isempty(v)
            s = v(2:end); % use the last tag, stripping the .
          end
          columns{jj} = s;
          switch s
           case 'Pressure'
            pressure_i = jj;
           case 'Cond'
            cond_i = jj;
           otherwise
          end
        end
       case 'data:'
        break; % done parsing...data follows so use load below
       
       otherwise
        % who knows that this header line is....
      end
    end
    fclose(fp); % done w/ header
    
    data = load(file_name);
    [n, m] = size(data);
    if (m ~= length(columns))
      error ('read_gpctd.m: number of data columns and tags does not match.');
    end
    valid_i = [1:n];
    gpctd.Time = [gpctd.Time; [this_epoch_time + (valid_i - 1)*sampling_interval/86400]']; % fractions of a day
    for ii = 1:m
      s = char(columns{ii});
      if (~isfield(gpctd,s))
        eval([sprintf('gpctd.%s = [];',s)]);
      end
      eval([sprintf('gpctd.%s = [gpctd.%s; data(valid_i,ii)];',s,s)]); % append
    end
  end
  %DEAD gpctd.elapsedTime = (gpctd.Time - first_epoch_time)*86400; % elapsed secs
  %DEAD gpctd.elapsedTime = fix(gpctd.elapsedTime); % integer elapsed secs?