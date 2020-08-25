% 
% Copyright (c) 2006-2012 by University of Washington.  All rights reserved. Confidential
%

% For QC:
% flags used by ARGO
QC_NO_CHANGE = 0;
QC_GOOD = 1;
QC_PROBABLY_GOOD = 2;
QC_PROBABLY_BAD = 3; % potentially correctable
QC_BAD = 4;
QC_CHANGED = 5; % explicit manual change
QC_UNSAMPLED = 6; % explicitly not sampled (vs. expected but missing)
QC_INTERPOLATED = 8; % interpolated value
QC_MISSING = 9; % value missing -- instrument timed out

qc_tag_names =  {'QC_NO_CHANGE','QC_GOOD','QC_PROBABLY_GOOD','QC_PROBABLY_BAD','QC_BAD','QC_CHANGED','QC_UNSAMPLED', 'QC_UNKNOWN7','QC_INTERPOLATED','QC_MISSING'};

ARGO_sample_depth_m = 2*25; % 1000m/40 samples, or 2000m/70 samples (spike and gradient tests use 3 samples, hence 2x the typical sampling distance)
% these can be overriden in sg_calib_constants
QC_bound_action = QC_BAD; % PARAMETER what QC flag to use on bound violations
QC_spike_action = QC_INTERPOLATED; % PARAMETER what QC flag to use on spike detection
QC_temp_min = -2.5; % PARAMETER [degC] Carnes; compare global Schmid -2.5 (labsea?) MDP -4.0
QC_temp_max = 43.0; % PARAMETER [degC] Carnes; compare global Schmid 40.0

QC_temp_spike_depth = 500.0; % PARAMETER [m] Carnes; ditto Schmid (db)
QC_temp_spike_shallow = 6.0/ARGO_sample_depth_m; % PARAMETER [degC/m] Carnes 2.0; Schmid 6.0
QC_temp_spike_deep = 2.0/ARGO_sample_depth_m; % PARAMETER [degC/m] Carnes 1.0; Schmid 2.0
% by inspection wa/aug04 see dives 56/57
QC_temp_spike_shallow = 0.05; % PARAMETER [degC/m] 
QC_temp_spike_deep = 0.01; % PARAMETER [degC/m]

QC_cond_min = 0.0; % PARAMETER [S/ml]
QC_cond_max = 10.0; % PARAMETER [S/ml] -- was 40

QC_cond_spike_depth = 500.0; % PARAMETER [m] Carnes
QC_cond_spike_shallow = 0.006; % QC_temp_spike_shallow/9.0; % PARAMETER [S/ml/m] scale by expected warm temp
QC_cond_spike_deep = 0.001; % QC_temp_spike_deep/9.0; % PARAMETER [S/ml/m] scale by expected cold temp

QC_salin_min = 19.0; % PARAMETER [PSU] was 2.0 per Carnes; ditto Schmid but we can't fly in waters that fresh
QC_salin_max = 45.0; % PARAMETER [PSU] Carnes; compare global Schmid 41.0

QC_overall_ctd_percentage = .3; % PARAMETER Carnes 30%
