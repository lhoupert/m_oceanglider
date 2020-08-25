% 
% Copyright (c) 2006-2011 by University of Washington.  All rights reserved. Confidential
%

% color is base, black dots are overlay pts that are different by threshold
bin_base = input('Base bin file: ','s');
load(bin_base);
Bin_base = Bin;

overlay_bin = input('Overlay bin file: ','s');
load(overlay_bin);
Bin_overlay = Bin;
clear Bin;

salin_diff_threshold = 0.1;
threshold = input(sprintf('Salinity difference threshold (CR = %.2f psu): ',salin_diff_threshold),'s');
if (~strcmp(threshold,''))
  salin_diff_threshold = str2num(threshold);
end
sg_calib_constants
titlestring = sprintf('SG%s %s %s vs. %s (%.2fpsu)', ...
                      id_str, pwd(), ...
                      regexprep(overlay_bin,'\_','-'),...
                      regexprep(bin_base,'\_','-'),...
                      salin_diff_threshold);
% dpp = 15; % dives per plot

[bin_total_points,ncolumns] = size(Bin_base);
bin_runs = unique(Bin_base(:,1));
[binc_total_points,ncolumns] = size(Bin_overlay);
binc_runs = unique(Bin_overlay(:,1));

shared_runs = intersect(bin_runs,binc_runs); % shared
runs = ask_which_runs(shared_runs);
runs = intersect(runs,shared_runs); % what we can actually show...
xruns = succinct_elts(runs);
fprintf(1,'Dives: %s\n',xruns);
nruns = length(runs);
iter = 1; % for waterfall_tsv
comment = '';
ir = 0;
max_diff = 0;
max_diff_run = 0;
for xr = 1:nruns;
  xrun = runs(xr);

  run_i = find(Bin_base(:,1) == xrun); % valid data points for this dive
  ct_depth_m = Bin_base(run_i,6);
  runc_i = find(Bin_overlay(:,1) == xrun); % valid data points for this dive in compare
  ct_depth_m_c = Bin_overlay(runc_i,6); % older
  apogee_i = find(ct_depth_m(1:end-1) > ct_depth_m(2:end));
  ct_depth_m(1:apogee_i) = -ct_depth_m(1:apogee_i); % distinguish dive/climb
  apogee_i = find(ct_depth_m_c(1:end-1) > ct_depth_m_c(2:end));
  ct_depth_m_c(1:apogee_i) = -ct_depth_m_c(1:apogee_i);  % distinguish dive/climb
  
  shared_depths = intersect(ct_depth_m, ct_depth_m_c);
  s_v = ismember(ct_depth_m,shared_depths);
  run_s_i = find(s_v);
  s_v = ismember(ct_depth_m_c,shared_depths);
  runc_s_i = find(s_v);
  run_i = run_i(run_s_i);
  runc_i = runc_i(runc_s_i);
  if (length(run_i) ~= length(runc_i))
      % sometimes one dive will have two entries for bin at depth 0, while
      % the other has only one. consider dropping the first...
      fprintf(1,'Different bin sizes for dive %d -- skipping\n',xrun);
      shared_depths(find(histc(ct_depth_m(run_s_i),shared_depths) > 1));
      shared_depths(find(histc(ct_depth_m_c(runc_s_i),shared_depths) > 1));
      continue;
  end

  ct_depth_m = ct_depth_m(run_s_i); % preserve dive/climb
  salin = Bin_base(run_i,8);
  ct_depth_m_c = ct_depth_m_c(runc_s_i); % preserve dive/climb
  salin_c = Bin_overlay(runc_i,8);
  match_i = find(ct_depth_m_c == ct_depth_m);
  diff_i = match_i(find(abs(salin(match_i) - salin_c(match_i)) > salin_diff_threshold));
  if (diff_i)
    % show this one...
    % used by waterfall to distinguish dives
    ir = ir + 1;
    irun = xrun;
    ct_depth_m = abs(ct_depth_m); % lose dive/climb distinction for display
    salin_TS = salin; % assume base is the same as overlay
    salin_TS(diff_i) = salin_c(diff_i); % except these overlay points, displayed as "interpolations"
    pitch_control = Bin_base(run_i,2); % sic but used to distinguish dive/climb
    temp = Bin_base(run_i,7);
    % bind these just in case those diagrams are enabled...always show base data
    % but really these are ignored...
    sigma = Bin_base(run_i,9);
    oxygen = Bin_base(run_i,10);
    red_scttr = Bin_base(run_i,11);
    blue_scttr = Bin_base(run_i,12);
    fluor = Bin_base(run_i,13);
    if (ncolumns >= 14)
      % optode present
      optode_oxygen = Bin_base(run_i,14);
      oxygen_dphase = Bin_base(run_i,15);
    end
    % can't do VV plot (no w)
    % can't do buoyancy plot (no buoy)
    waterfall_tsv;
  else
    % for those runs not showing differences, record the largest diff
    maxd = abs(salin(match_i) - salin_c(match_i));
    if (maxd > max_diff)
      max_diff = maxd;
      max_diff_run = xrun;
    end
  end
end
if (max_diff > 0)
  fprintf(1,'Max diff of non-displayed dives %.4fPSU on dive %d\n',max_diff,max_diff_run);
end
