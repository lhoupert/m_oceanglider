% 
% Copyright (c) 2006-2011, 2013 by University of Washington.  All rights reserved. Confidential
%
function txt = nc_bin_edit_dc(h_obj,event_obj)
  eventdata = get(event_obj);
  user_data = get(gcbf, 'UserData');
  pos = eventdata.Position; % ensured to exist
  x= pos(1); % X
  x_info = user_data.X;
  x_values = x_info{1};
  x_label = x_info{2};
  x_tol = x_info{3};
  
  y = pos(2); % Y
  y_info = user_data.Y;
  y_values = y_info{1};
  y_label = y_info{2};
  y_tol = y_info{3};

  Profiles_i = user_data.Profiles_i;
  dive_numbers = user_data.dive_numbers;
  % assumes Bin_i are the indices in Bin of the set of elements displayed
  try
	  [values_i,dives_i] = find(x_values >= x-x_tol & x_values <= x+x_tol &...
								y_values >= y-y_tol & y_values <= y+y_tol);
	  % dives_i has the indices of *all* the dives with those values
	  % restrict to the ones the user specified (Profiles_i)
	  dives = unique(dive_numbers(intersect(Profiles_i,unique(dives_i))));
	  dives = succinct_elts(dives'); % make a nice string
	  evalin('base',sprintf('selected_dives = [%s]'';',dives)); % make available
  catch
	  dives = 'Cannot determine dives! UserData lost?';
  end
  txt = {sprintf('%s: %.2f %s: %.2f',x_label,x,y_label,y); dives};
  