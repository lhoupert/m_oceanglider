% 
% Copyright (c) 2006-2013 by University of Washington.  All rights reserved. Confidential
%

% [tmin,tmax,tdel] = plot_bounds(temp_m,-3,40,[12,2; 6,1; 0,0.5]);
% [smin,smax,sdel] = plot_bounds([salin_raw;salin],19,40,[12,2; 6,1; 0,0.5]);
% [wmin,wmax,wdel] = plot_bounds([w;w_stdy;w_h2o],-24,16,[100,10; 60,5; 40,4; 20,2; 0,1]);

function [pmin,pmax,pdel] = plot_bounds(data,pmin_valid,pmax_valid,varargin)
	data = reshape(data,prod(size(data)),1);
  real_i = find(imag(data) == 0 & ~isnan(data));
  if (isempty(data(real_i))) 
      % e.g., QPE May09 sg167 dive 328: wetlabs present but no data recorded
      pmin = pmin_valid;
      pmax = pmax_valid;
      pdel = 1;
      return;
  end
  data_min = min(data(real_i));
  data_max = max(data(real_i));
  if (size(varargin,2) == 2)
    m = varargin{2}; % round to given bound
  else
    m = 1; % round to nearest integer
  end
  y = data_min/m; q = floor(y); data_min = q*m;
  y = data_max/m; q = ceil (y); data_max = q*m;
  pdel = 1; % default
  if (data_min == data_max)
    data_max = data_max + pdel;
  end
  if (data_min > pmax_valid || data_max < pmin_valid)
      % all the data is out of range...trust the valid bounds
      pmin = pmin_valid;
      pmax = pmax_valid;
  else
    pmin = max(pmin_valid, data_min);
    pmax = min(pmax_valid, data_max);
  end
  % caller can supply an array of [range delta] pairs, e.g.,
  %      [12,2; 6,1; 0,0.5]
  % with ranges in descending order which indicate if the actual
  % range is greater than the specified range, use the associated delta
  if (size(varargin,2) >= 1)
    delta_ranges = varargin{1};
    [m,n] = size(delta_ranges);
    pdiff = pmax - pmin;
    for i = 1:m
      if (pdiff > delta_ranges(i,1))
        pdel = delta_ranges(i,2);
        break;
      end
    end
  end
