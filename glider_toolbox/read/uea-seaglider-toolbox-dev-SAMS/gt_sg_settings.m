% these parameters are only used if they are not defined in the glider
% specific sg_calb_const.m file!!!
% otherwise the sg_calib_const.m file is the higher weighted one!

% length = 1.8;
% flight_regression_passes = 0;%3;

% min_flush_speed = 3; % cm/s ?% 0.1
% max_flush_speed = 80;          %  0.5
% 
% 
% tautemp = 0.6; % 0.2 - 1.2
% tol_glide = 0.0001;
% 
% %conv_filter = [1 5 15 35 70 35 15 5 1]; 
% conv_filter = [ 1 6 15 20 15 6 1]; conv_filter = conv_filter ./ sum(conv_filter);
%                  
%                 
% conv_filter_slowsampling = [0.5 1 1 1 1 0.5];
% conv_filter_fast_sampling_runing_mean = [0.5 1 1 1 1 1 1 1 1 1 0.5];
% 
% 
% tconinv_auto = true; % not working


%NEWONES:

% sbe_temp_tau
% sbe_cond_coefficients
%sbe_flush_speed_min
% WHat to do with flagged values: NaN; leave; or interp
% flight_regression_passes
