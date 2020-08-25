% *************************************************************************
%
% *************************************************************************
function [calibData] = readCalibConstants(calibFilename)
  calibData = [];
  strCalibData = [];
  nCommentCount = 0;
  
  fid = fopen(calibFilename);
  if (fid > 0)
    tline = fgets(fid);
    i = 0;
    while ischar(tline)
      if (tline(1) ~= '%')
        i = i + 1;
        strCalibData{i} = tline;
      else
        nCommentCount = nCommentCount + 1;
      end
      
      tline = fgets(fid);
    end
    fclose(fid);
    
    if (isempty(strCalibData) == 0)
      c = char(cellstr(cell2mat(strCalibData)));
      pos = findstr(c, 'id_str');
      t = c(pos:end);
      sq = findstr(t, '''');
      id_str = t(sq+1:sq(2)-1);
    end
  end
  
  for i=1:length(strCalibData)
    try
      strEval = ['calibData.' char(strCalibData(i))];

      % SOMETIMES THE END-USERS OR PILOTS FORGET TO PUT A SEMI-COLON AT THE END OF THE LINE IN THE CALIBRATION FILE. IF THE THE
      % SEMI-COLON IS MISSING, MATLAB WILL PRINT OUT TO THE COMMAND WINDOW THE ENTIRE DATA STRUCTURE FOR EACH AND EVERY INSTANCE
      % WHERE THE SEMI-COLON WAS OMITTED. THIS CAN RESULT IN MISSING TRUE ANOMALIES.  THE FOLLOWING LINES OF CODE ADD A SEMI-COLON
      % TO THE END OF THE LINE IF IT IS MISSING BEFORE EXECUTING (EVAL'ing) THE LINE.
      if strEval(end) ~= ';'
          strEval = strcat(strEval, ';');
      end
      eval(strEval);
    catch
    end
  end
end
