% 
% Copyright (c) 2006-2011, 2013 by University of Washington.  All rights reserved. Confidential
%
function txt = bin_edit_dc(h_obj,event_obj)
  eventdata = get(event_obj);
  user_data = get(gcbf, 'UserData');
  pos = eventdata.Position; % ensured to exist
  x= pos(1); % X
  x_info = user_data.X;
  x_index = x_info{1};
  x_label = x_info{2};
  x_tol = x_info{3};
  
  y = pos(2); % Y
  y_info = user_data.Y;
  y_index = y_info{1};
  y_label = y_info{2};
  y_tol = y_info{3};

  Bin = user_data.Bin_data.Bin;
  Bin_i = user_data.Bin_data.Bin_i;
  % assumes Bin_i are the indices in Bin of the set of elements displayed
  try
    dives = unique(Bin(Bin_i(find(Bin(Bin_i,x_index) >= x-x_tol & Bin(Bin_i,x_index) <= x+x_tol &...
                                  Bin(Bin_i,y_index) >= y-y_tol & Bin(Bin_i,y_index) <= y+y_tol)),1));
    dives = succinct_elts(dives'); % make a nice string
    evalin('base',sprintf('selected_dives = [%s]'';',dives)); % make available
  catch
    dives = 'Cannot determine dives! UserData lost?';
  end
  txt = {sprintf('%s: %.2f %s: %.2f',x_label,x,y_label,y); dives};
