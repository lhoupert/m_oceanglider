% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%

% linear interpolation of y wrt x using the pre and post points of regions in interp_points_i
% optional anchors are legit points to interp from that are part of the interp range
% if no anchors are given then the points just adjacent to a range are used, if they are good
% qc is the associated qc variable, qc_tag to use if there are issues with the interpolation
function [y_interp,y_qc,x] = interpolate_data(y,x,interp_points_i,type,directives,y_qc,qc_tag)
  y_interp = y;
  if (isempty(interp_points_i))
    return; % nothing to do
  end
  np = length(y);
  qc_declarations;

  % anchors are adjacent points, unless QC_BAD
  % CONSIDER: find all QC_GOOD points and then, for each interp run, finding the nearest good point as anchors
  % could interp over bad points between good points
  % This interpolation routine assumes that the x(interp_i) points are monotonically increasing and bounded by the anchors
  % this is the case w/ time but not temp (in TS interp).
  last_i = 1;
  for break_i = [find(diff(interp_points_i) > 1);length(interp_points_i)]'
    pre_index = interp_points_i(last_i);
    post_index = interp_points_i(break_i);
    ip_i = [pre_index:post_index]; % before extension to anchors for messages
    pre_index = max(pre_index - 1,1);
    post_index = min(post_index + 1,np);
    interp_i = [pre_index:post_index]; % which points to recompute
    anchors_i = [pre_index;post_index];
    bad_anchors_i = bad_qc(y_qc(anchors_i));
    if (~isempty(bad_anchors_i))
      if (~isempty(find(y_qc(ip_i) == QC_INTERPOLATED)))
        % only suggest this if there were points that needed work (they could all be BAD)
        reason = 'bad interpolation anchors';
        y_qc = assert_qc(qc_tag,y_qc,ip_i,reason);
        drv_suggest(directives,sprintf('bad_%s data_points in_between %d %d %% %s',type,pre_index,post_index,reason));
      end
    else
      coefficents = polyfit(x(anchors_i),y(anchors_i),1); % linear fit amongst the anchors
      if (find(isnan(coefficents) | isinf(coefficents)))
        reason = 'unable to interpolate (NaN/Inf)';
        y_qc = assert_qc(qc_tag,y_qc,ip_i,reason);
        drv_suggest(directives,sprintf('bad_%s data_points in_between %d %d %% %s',type,pre_index,post_index,reason));
      else
        y_interp(interp_i) = polyval(coefficents,x(interp_i));
      end
    end
    last_i = break_i+1;
  end
