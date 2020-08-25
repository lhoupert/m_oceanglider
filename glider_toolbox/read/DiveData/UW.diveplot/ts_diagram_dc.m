% 
% Copyright (c) 2006-2011, 2013 by University of Washington.  All rights reserved. Confidential
%
function txt = ts_diagram_dc(h_obj,event_obj)
  eventdata = get(event_obj);
  user_data = get(gcbf,'UserData'); % a struct
  pos = eventdata.Position; % ensured to exist
  s= pos(1); % salin
  t= pos(2); % temp
  salin = user_data.salin;
  temp = user_data.temp;
  depth = user_data.depth;
  directives = user_data.directives;
  ttol = user_data.tolerance;
  txt = {};
  txt{end+1} = sprintf('Point: %.2fpsu %.2fC',s,t); % debug...
  try
    % widen tolerances until we get points
    points_i = [];
    while (isempty(points_i))
      stol = 10*ttol;
      points_i = find(salin >= s-stol & salin <= s+stol & ...
                      temp  >= t-ttol & temp  <= t+ttol);
      ttol = ttol * 2;
    end
    
    dive_i  = intersect(points_i,[directives.start_of_dive:directives.start_of_climb]);
    climb_i = intersect(points_i,[directives.start_of_climb:directives.end_of_climb]);
    % find points on dive and climb
    dive_point_i = 0;
    if (~isempty(dive_i))
      dive_point_i = dive_i(1);
    end
    climb_point_i = 0;
    if (~isempty(climb_i))
      climb_point_i = climb_i(1);
    end

    interp_salinity = directives.interpolated_salinity;
    stable_ts = setdiff([1:length(salin)],interp_salinity);
    if (dive_point_i)
      interp = '';
      is = intersect(interp_salinity,dive_point_i);
      if (is)
        stable_avant_i = stable_ts(find(stable_ts <= is));
        stable_apres_i = stable_ts(find(stable_ts >= is));
        interp = sprintf(' interp [%d:%d]',stable_avant_i(end),stable_apres_i(1));
      end
      txt{end+1} = sprintf('D: %3d %.1fm%s',dive_point_i,depth(dive_point_i),interp);
    end
    if (climb_point_i)
      interp = '';
      is = intersect(interp_salinity,climb_point_i);
      if (is)
        stable_avant_i = stable_ts(find(stable_ts <= is));
        stable_apres_i = stable_ts(find(stable_ts >= is));
        interp = sprintf(' interp [%d:%d]',stable_avant_i(end),stable_apres_i(1));
      end
      txt{end+1} = sprintf('C: %3d %.1fm%s',climb_point_i,depth(climb_point_i),interp);
    end
  catch
    1; % pass
  end
  