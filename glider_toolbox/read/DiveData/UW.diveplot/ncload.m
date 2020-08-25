% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

% Read an arbitrary netcdf file and return its variables as a matlab structure
% Requires MATLAB 2008b at least for the netcdf primitives
function ncinfo = ncload(ncdfile)

  ncid = netcdf.open(ncdfile,'NC_NOWRITE'); % this works or gives errors
  % BUG: if any dimname, varname, or attname starts with an '_' matlab
  % hates it and complains about bad input character (following a '.')
  % E.g., ARGO has attibutes like _FillValue
  [ndims,nvars,ngatts,unlimdimID]= netcdf.inq(ncid);

  % intern global attributes into the structure
  global_id = netcdf.getConstant('NC_GLOBAL');
  for gvarid = 0:ngatts-1 % the ids are 0-based
    attname = netcdf.inqAttName(ncid,global_id,gvarid);
    % For strings, see note below about possible transposition required
    attvalue = netcdf.getAtt(ncid,global_id,attname);
    eval(['ncinfo.GLOBALS.' attname ' = attvalue;']); % intern
  end
  
  % intern each dimension into the structure
  for dimid = 0:ndims-1 % the ids are 0-based
    [dimname, dimlen] = netcdf.inqDim(ncid,dimid);
    eval(['ncinfo.DIMENSIONS.' dimname ' = dimlen;']); % intern
  end

  % intern each variable into the structure
  for varid = 0:nvars-1 % the ids are 0-based
    [varname vartype dimids natts] = netcdf.inqVar(ncid,varid);
    value = netcdf.getVar(ncid,varid); % into local scope
    % Quite often char type (2) arrays need transposition to be read
    % but it is complicated depending on the number of dimensions.  
    % Leave this to caller? Or handle the [n x 1] case for them?
    % value = value'; % Leave to caller...
    % in addition, to convert _qc character arrays use str2num() on the array
    eval(['ncinfo.' varname ' = value;']); % intern
    if (0)
      for attid = 0:natts-1 % the ids are 0-based
        attname = netcdf.inqAttName(ncid,varid,attid);
        attvalue = netcdf.getAtt(ncid,varid,attname);
        eval(['ncinfo.' varname '_ATTRIBUTES.' attname ' = attvalue;']); % intern
      end
    end
  end
  netcdf.close(ncid);
end
