function data = gt_sg_sub_currents(data)

for istep = intersect([data.eng.dive],[data.log.dive]);
    % Duration, as defined by time between the two GPS coordinates
    data.flight(istep).duration = (data.gps_postdive(istep,3) - data.gps_predive(istep,3)) * 24 * 60 * 60; % seconds
    delta_t =  gt_sg_sub_diff(data.eng(istep).elaps_t);

    %% Estimate of travel from flight model:
    % Calculate mean vec(u) and vec(v) for each sampling interval.
    [vv_s, uu_s] = pol2cart( data.eng(istep).head * pi/180 , data.flight(istep).model_horz_spd); % m.s-1
    % Multiply by duration of sampling interval then add cumulatively to
    % obtain x and y transport from start at each time step.
    x_s = cumtrapz(uu_s.*delta_t); % m
    y_s = cumtrapz(vv_s.*delta_t); % m
    % Export data
    data.flight(istep).xy_flightdistance = [x_s(end), y_s(end)];
    
    %% Travel distance from GPS coordinates
    % Next bit of code is legacy from S. Schmidtko. Unverified.
    data.flight(istep).xy_surfaceings = [(data.gps_postdive(istep,2)-...
        data.gps_predive(istep,2)).*111319.5.*cosd(data.gps_postdive(istep,1)),...
        (data.gps_postdive(istep,1)-data.gps_predive(istep,1)).*111319.5];
        % TODO: verify units
    
    %% Estimate currents
    data.hydrography(istep).w_H2O = data.flight(istep).glide_vert_spd - data.flight(istep).model_vert_spd;
    data.hydrography(istep).DAC_u = (data.flight(istep).xy_surfaceings(1)-data.flight(istep).xy_flightdistance(1))/data.flight(istep).duration;
    data.hydrography(istep).DAC_v = (data.flight(istep).xy_surfaceings(2)-data.flight(istep).xy_flightdistance(2))/data.flight(istep).duration;

    %% Estimate glider's trajectory
    % Glider's travel through the water, not accounting for currents.
    data.flight(istep).trajectory_xy_nocurrents = [x_s', y_s'];
    % Estimated GPS coordinates using GPS surfaceings and DAC
    % DOES NOT ACCOUNT FOR CURVATURE OF THE EARTH
    %data.flight(istep).trajectory_latlon_estimated = [data.gps_predive(istep,1) + ( data.gps_postdive(istep,1) - data.gps_predive(istep,1)) / y_s(end) * y_s',...
    %    data.gps_predive(istep,2) + ( data.gps_postdive(istep,2) - data.gps_predive(istep,2)) / x_s(end) * x_s'];
    
    data.flight(istep).trajectory_latlon_estimated = [...
        interp1(real([y_s(1) y_s(end)]),real([data.gps_predive(istep,1) data.gps_postdive(istep,1)]),real(y_s),'linear','extrap')',...
        interp1(real([x_s(1) x_s(end)]),real([data.gps_predive(istep,2) data.gps_postdive(istep,2)]),real(x_s),'linear','extrap')'];
end
end
