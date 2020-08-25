function data = gt_sg_sub_flightmodelregression(data)
%
%
% 20 dives enough for regression. Check moving 20 dive window for time
% varying parameters.
%
% w2 vs N
%
%
%
% data = gt_sg_sub_flightvec(data)
%
% Wrapper for the UW flightvec.m routines which calculates, using the
% hydrodynamic flight model and buoyancy, the speed and glide angle of the
% glider. This is a better (and iterative) approximation of flight than the
% glider_slope version.
%
% All credits to original authors at the University of Washington
%
% B.Y.QUESTE Apr 2015
%

dives = intersect([data.log.dive],[data.eng.dive]);

% Determine hull length if at all possible
if isfield(data.gt_sg_settings,'length')
    hullLength = data.gt_sg_settings.length;
    hullLengthMessage = {['Using a hull length of ',num2str(hullLength),' m.']};
else
    try
        hullLength = mode([data.log.LENGTH]);
        hullLengthMessage = {['Using a hull length of ',num2str(hullLength),' m as provided by $LENGTH parameter in the log files.']};
    catch
        hullLength = [];
        while any([~isnumeric(hullLength),isnan(hullLength),isempty(hullLength)])
            hullLength = input('Hull length in meters (1.8 for standard, 2.0 for ogive):');
        end
        hullLengthMessage = {['Using a hull length of ',num2str(hullLength),' m.']};
    end
end

gt_sg_sub_echo({'Calculating glider flight parameters based on buoyancy and pitch.',...
    'Output in flight substructure as "model_spd" and "model_slope".',...
    hullLengthMessage{:}});

%% User input
% Regress volmax yes or no?
flight_regression_passes = 'string';

if ~isfield(data.gt_sg_settings,'flight_regression_passes')
    while flight_regression_passes > 10 | flight_regression_passes < 0;
        flight_regression_passes = input('Number of flight model regression passes to do? (min 0, max 10, default 3) :  ');
    end
else
    flight_regression_passes = data.gt_sg_settings.flight_regression_passes;
end

%% If we want to regress, do the following...
if flight_regression_passes ~= 0
    gt_sg_sub_echo({'Beginning regressions. This WILL take a while.',...
        'Maybe go write a paper or something...'});
    
    % Determine cast for gridding in hd_ regressions
    cast_dir = arrayfun(@(x) ((x.elaps_t*0)+1).*x.dive,data.eng,'Uni',0); % Dive no
    for dstep=1:numel(cast_dir)
        cast_dir{dstep}(1:data.hydrography(dstep).max_pressure_index) = -cast_dir{dstep}(1:data.hydrography(dstep).max_pressure_index); % Sign for direction
    end
    
    % Create get_exp handle to make all regressed coefficients of similar
    % units (works better).
    get_exp = @(x) floor(log10(x));
    get_coeff = @(x) x.*10.^-get_exp(x);
        
    run_num = 0;
    plot_vertical_velocities(run_num);
    for regress_step = 1:flight_regression_passes
        
        run_num = run_num+1;        
        gt_sg_sub_echo({['Beginning pass ',num2str(run_num),' of ',num2str(2*flight_regression_passes),' (hydrodynamic parameter regressions).']});
        % REGRESS hd_a, b and c
        % Target is reducing difference between up and down estimates of w_H20
        % 1) make all inputs same order of magnitude
        hd_coeffs = get_coeff([data.gt_sg_settings.hd_a,data.gt_sg_settings.hd_b,data.gt_sg_settings.hd_c]);
        hd_exp = get_exp([data.gt_sg_settings.hd_a,data.gt_sg_settings.hd_b,data.gt_sg_settings.hd_c]);
        % 2) set options
        options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',200);
        % 3) do regression AND reapply exponents
        abc_regressed = fminsearch(@run_abcregression,hd_coeffs,options) .* (10.^hd_exp);
        % 4) joyfully announce it to the world
        gt_sg_sub_echo({'HD_A, B & C regression finished. Proceeding with the new values.',...
            ['        Original values  |  Regressed values'],...
            ['hd_a : ',num2str(data.gt_sg_settings.hd_a),' | ',num2str(abc_regressed(1))],...
            ['hd_b : ',num2str(data.gt_sg_settings.hd_b),' | ',num2str(abc_regressed(2))],...
            ['hd_c : ',num2str(data.gt_sg_settings.hd_c),' | ',num2str(abc_regressed(3))]});
        data.gt_sg_settings.hd_a = abc_regressed(1);
        data.gt_sg_settings.hd_b = abc_regressed(2);
        data.gt_sg_settings.hd_c = abc_regressed(3);
        plot_vertical_velocities(run_num);

        run_num = run_num+1;        
        gt_sg_sub_echo({['Beginning pass ',num2str(run_num),' of ',num2str(2*flight_regression_passes),' (volume parameter regressions).']});
        % REGRESS volmax (= vbdbias), abs_compress, therm_expan
        % Target is reducing w_H20; which is making glide_ and model_ values as
        % close as possible to each other.
        % 1) make all inputs same order of magnitude
        vol_coeffs = get_coeff([data.gt_sg_settings.volmax,data.gt_sg_settings.abs_compress,data.gt_sg_settings.therm_expan]);
        vol_exp = get_exp([data.gt_sg_settings.volmax,data.gt_sg_settings.abs_compress,data.gt_sg_settings.therm_expan]);
        % 2) set options
        options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',200);
        % 3) do regression AND reapply exponents
        volparam_regressed = fminsearch(@run_volparam_regression,vol_coeffs,options) .* (10.^vol_exp);
        % 4) joyfully announce it to the world
        gt_sg_sub_echo({'Volume parameter estimation finished. Proceeding with the new values.',...|
            ['               Original values  |  Regressed values'],...
            ['volmax :       ',num2str(data.gt_sg_settings.volmax),' -> ',num2str(volparam_regressed(1))],...
            ['abs_compress : ',num2str(data.gt_sg_settings.abs_compress),' -> ',num2str(volparam_regressed(2))],...
            ['therm_expan :  ',num2str(data.gt_sg_settings.therm_expan),' -> ',num2str(volparam_regressed(3))]});
        data.gt_sg_settings.volmax = volparam_regressed(1);
        data.gt_sg_settings.abs_compress = volparam_regressed(2);
        data.gt_sg_settings.therm_expan = volparam_regressed(3);
        plot_vertical_velocities(run_num);
        
    end
    
    % Plot difference of up and down casts
    % 1. Grid all up/down estimates of w_H20
    w_H2O = bindata2(...
        [data.flight(dives).glide_vert_spd] - [data.flight(dives).model_vert_spd],... ( = w_H20 )
        abs([cast_dir{dives}]),[data.hydrography(dives).depth],... % (x and y positions)
        [0:max(dives)],[10:4:990]); % (x and y bins)
    w_H2O_diff = bindata2(...
        [data.flight(dives).glide_vert_spd] - [data.flight(dives).model_vert_spd],... ( = w_H20 )
        [cast_dir{dives}],[data.hydrography(dives).depth],... % (x and y positions)
        [-max(dives):max(dives)],[10:4:990]); % (x and y bins)
    % w_H2O = downcasts on the left, up casts on the right, with dives
    % mirrored, so "fold" the matrix in the middle and substract the
    % two sides.
    w_H2O_diff = w_H2O_diff(:,max(dives):-1:1) - w_H2O_diff(:,max(dives)+1:end);
    
    figure
    subplot(2,1,1)
    imagesc([1:max(dives)],[10:4:990],real(w_H2O))
    xlabel('Dive Number'); ylabel('Depth (m)');
    title({'Estimated W_H_2_O'});
    colorbar; caxis([-1 1]*0.05);  
    subplot(2,1,2)
    imagesc([1:max(dives)],[10:4:990],real(w_H2O_diff))
    xlabel('Dive Number'); ylabel('Depth (m)');
    title({'Absolute difference between up- and downcast estimates of W_H_2_O','Ideally, it should appear as random noise dissimilar T&S structures or dive profile shape.'});
    colorbar; caxis([-0.05 0.05]);
end

%% Run the flight model

settings.volmax = data.gt_sg_settings.volmax;
settings.hd_a = data.gt_sg_settings.hd_a;
settings.hd_b = data.gt_sg_settings.hd_b;
settings.hd_c = data.gt_sg_settings.hd_c;
settings.abs_compress = data.gt_sg_settings.abs_compress;
settings.therm_expan = data.gt_sg_settings.therm_expan;
settings.hull_length = hullLength;

run_flightmodel(settings)

%%%%%%%%%%%%%%%%% END OF MAIN FUNCTION %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% START OF MINIMISATION CALLS %%%%%%%%%%%%%%

    function w_rmsd = run_volparam_regression(vol_param)
        % Define settings
        params = vol_param .* (10.^vol_exp);
        settings.volmax = params(1);
        settings.hd_a = data.gt_sg_settings.hd_a;
        settings.hd_b = data.gt_sg_settings.hd_b;
        settings.hd_c = data.gt_sg_settings.hd_c;
        settings.abs_compress = params(2);
        settings.therm_expan = params(3);
        settings.hull_length = hullLength;
        % Run the model
        run_flightmodel(settings)
        % Score function: root mean square difference of the glide and
        % model vertical speeds. This should effectively minimise w_H20 and
        % hopefully also correct slanted w_H20 linked to thermal expansion
        % in high thermal gradient regions.
        w_rmsd = nanmean(abs([data.flight(dives).glide_vert_spd]-[data.flight(dives).model_vert_spd])); % TODO: grid in time to avoid aliasaing of vertical sampling of vertically moving features
    end

    function w_rmsd = run_abcregression(hd_param)
        % Define settings
        settings.volmax = data.gt_sg_settings.volmax;
        params = hd_param .* (10.^hd_exp);
        settings.hd_a = params(1);
        settings.hd_b = params(2);
        settings.hd_c = params(3);
        settings.abs_compress = data.gt_sg_settings.abs_compress;
        settings.therm_expan = data.gt_sg_settings.therm_expan;
        settings.hull_length = hullLength;
        % Run the model
        run_flightmodel(settings)
        % Score function: root mean square difference of the up and down
        % estimates of vertical velocity across all dives. This requires
        % gridding of all the estimates and then an RMSD.
        
        % 1. Grid all up/down estimates of w_H20
        w_H2O = bindata2(...
            [data.flight(dives).glide_vert_spd] - [data.flight(dives).model_vert_spd],... ( = w_H20 )
            [cast_dir{dives}],[data.hydrography(dives).depth],... % (x and y positions)
            [-max(dives):max(dives)],[10:4:990]); % (x and y bins)
        % w_H2O = downcasts on the left, up casts on the right, with dives
        % mirrored, so "fold" the matrix in the middle and substract the
        % two sides.
        w_H2O = w_H2O(:,max(dives):-1:1) - w_H2O(:,max(dives)+1:end);
        % TODO: Do we want to weight bins by number of samples in each bin
        % (could be useful for severely skewed profile depth distributions)
        w_rmsd = nanmean(abs(w_H2O(:)));
    end

%%%%%%%%%%%%%%%% END OF MINIMISATION CALLS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% START OF FLIGHT MODEL ROUTINES %%%%%%%%%%%%%%

    function run_flightmodel(settings)
        for istep = intersect([data.eng.dive],[data.log.dive])
            
            % settings. should have volmax,hd_a,hd_b,hd_c,abs_compress,therm_expan,hull_length
            
            %% Determine glider buoyancy based on volmax and in situ density
            % Calculate volume
            vbd = data.eng(istep).vbdCC; % cm^3
            vol0 = settings.volmax + (data.log(istep).C_VBD - data.gt_sg_settings.vbd_min_cnts)/data.gt_sg_settings.vbd_cnts_per_cc; % cm^3
            vol = (vol0 + vbd) .* ...
                exp(-settings.abs_compress * data.hydrography(istep).pressure + settings.therm_expan * (data.hydrography(istep).temp - data.gt_sg_settings.temp_ref)); % cm^3
            
            % Calculate buoyancy of the glider
            % Buoyancy units
            %   kg + (kg.m-3 * cm^3 * 10^-6)
            % = kg + (kg.m-3 * m^3)
            % = kg
            data.flight(istep).buoyancy =  -data.gt_sg_settings.mass + ( data.hydrography(istep).rho .* vol * 1.e-6 ); % kg
            
            % Calculate glide slope and speed
            [data.flight(istep).model_spd, data.flight(istep).model_slope]... % m.s-1, degrees
                = flightvec(...
                data.flight(istep).buoyancy,... % kg
                data.eng(istep).pitchAng,... % degrees
                settings.hull_length,... % Hull length in meters, excluding antenna
                settings.hd_a,... % rad^-1
                settings.hd_b,... % m^1/4 . kg^1/4 . s^-1/2
                settings.hd_c,... % rad^-2
                data.gt_sg_settings.rho0,...
                istep); % kg.m-3

            % Determine vertical and horizontal speed
            data.flight(istep).model_horz_spd = data.flight(istep).model_spd.*cos(data.flight(istep).model_slope*pi/180);
            data.flight(istep).model_vert_spd = data.flight(istep).model_spd.*sin(data.flight(istep).model_slope*pi/180);
        end
    end

    function [ umag, thdeg ] = flightvec( bu, ph, xl, a, b, c, rho0, istep)
        % ********************************************************************************
        %	function flightvec0(bu, ph, xl, a, b, c, rho0 )
        %		Solves unaccelerated flight equations iteratively
        %		for speed magnitude umag and glideangle thdeg
        % ********************************************************************************
        
        gravity = gsw_grav(data.gps_postdive(istep,1),...
            data.hydrography(istep).pressure); % m.s-2
        
        th = (pi/4)*sign(bu);	% initial slope set to 1 or -1
        
        buoyforce = gravity .* bu; % kg.m.s-2
        
        q = ( sign(bu).*buoyforce./(xl*xl*b) ).^(4/3); 	% dynamic pressure for vertical flight
        alpha = 0.;
        
        tol = 0.001;
        j = 0;
        q_old = zeros(size(bu));
        param = ones(size(bu));
        valid = find( bu ~= 0 & sign(bu).*sign(ph) > 0 );
        umag = zeros(size(bu));
        thdeg = zeros(size(bu));
        
        while(~isempty( find( abs( (q(valid)-q_old(valid))./q(valid) ) > tol )) & j <= 15);
            
            q_old = q;
            param_inv = a*a*tan(th).*tan(th).*q.^(0.25)/(4*b*c);
            
            valid = find( param_inv > 1 & sign(bu).*sign(ph) > 0);	% valid solutions for param < 1
            
            param(valid) = 4*b*c./(a*a*tan(th(valid)).*tan(th(valid)).*q(valid).^(0.25));
            q(valid) = ( buoyforce(valid).*sin(th(valid))./(2*xl*xl*b*q(valid).^(-0.25)) ).*  ...
                (1 + sqrt(1-param(valid)));
            alpha(valid) = ( -a*tan(th(valid))./(2*c) ).*(1 - sqrt(1-param(valid)));
            thdeg(valid) = ph(valid) - alpha(valid);
            
            stall = find( param_inv <= 1 | sign(bu).*sign(ph) < 0 );	% stall if param >= 1
            
            % TODO: REMOVE ME
            stall = [];
            %
            
            q(stall) = 0.;
            thdeg(stall) = 0.;
            
            th = thdeg*pi/180;
            
            j = j+1;
            % TODO: Should this be replaced with actual density.....??
            umag = sqrt( 2*q/rho0 );
            
        end;
        
    end
%%%%%%%%%%%%%% END OF FLIGHT MODEL ROUTINES %%%%%%%%%%%%%%
%%%%%%%%%%%%%% START OF PLOTTING ROUTINES %%%%%%%%%%%%%%
    function plot_vertical_velocities(run_num)
        
        settings.volmax = data.gt_sg_settings.volmax;
        settings.hd_a = data.gt_sg_settings.hd_a;
        settings.hd_b = data.gt_sg_settings.hd_b;
        settings.hd_c = data.gt_sg_settings.hd_c;
        settings.abs_compress = data.gt_sg_settings.abs_compress;
        settings.therm_expan = data.gt_sg_settings.therm_expan;
        settings.hull_length = hullLength;
        run_flightmodel(settings)
        
        mean_w_H2O = bindata2(...
            [data.flight(dives).glide_vert_spd] - [data.flight(dives).model_vert_spd],... ( = w_H20 )
            sign([cast_dir{dives}]),[data.hydrography(dives).depth],... % (x and y positions)
            [-2 0 2],[0:4:1000]); % (x and y bins)
        mean_glide = bindata2(...
            [data.flight(dives).glide_vert_spd],...
            sign([cast_dir{dives}]),[data.hydrography(dives).depth],... % (x and y positions)
            [-2 0 2],[0:4:1000]); % (x and y bins)
        mean_model = bindata2(...
            [data.flight(dives).model_vert_spd],...
            sign([cast_dir{dives}]),[data.hydrography(dives).depth],... % (x and y positions)
            [-2 0 2],[0:4:1000]); % (x and y bins)
        
        % Plot prelim state for future comparison
        figure;
        subplot(1,2,1)
        hold on
        plot([mean_glide(:,1);mean_glide(end:-1:1,2)],[2:4:1000,1000:-4:2],'-b');
        plot([mean_model(:,1);mean_model(end:-1:1,2)],[2:4:1000,1000:-4:2],'-r');
        axis ij tight
        plot([0 0],get(gca,'YLim'),'-k');
        ylabel('Depth (m)')
        xlabel('Glider W')
        title('Blue, dP/dt | Red, buoyancy model')
        subplot(1,2,2)
        hold on
        plot([mean_w_H2O(:,1);mean_w_H2O(end:-1:1,2)],[2:4:1000,1000:-4:2],'-k');
        axis ij tight
        plot([0 0],get(gca,'YLim'),'-k');
        ylabel('Depth (m)')
        xlabel('W(H_2O)')
        title(['Run number: ',num2str(run_num),' of ',num2str(2*flight_regression_passes)]);
    end
%%%%%%%%%%%%%% END OF PLOTTING ROUTINES %%%%%%%%%%%%%%

end

function [ym,yb] = bindata2(y,x1,x2,x1rg,x2rg)
%function [ym,yb] = bindata2(y,x1,x2,x1rg,x2rg)
%Computes:
%ym(ii,jj) = mean(y(x1>=x1rg(ii) & x1 < x1rg(ii+1) & x2>=x2rg(jj) & x2 < x2rg(jj+1))
%for every ii, jj
%If a bin is empty it returns nan for that bin
%using a fast algorithm which uses no looping
%Also returns yb, the approximation of y using binning (useful for r^2
%calculations). Example:
%
%x = randn(500,2);
%y = sum(x.^2,2) + randn(500,1);
%xrg = linspace(-3,3,10)';
%[ym,yb] = bindata2(y,x(:,1),x(:,2),xrg,xrg);
%subplot(1,2,1);plot3(x(:,1),x(:,2),y,'.');
%subplot(1,2,2);h = imagesc(xrg,xrg,ym);
%set(h,'AlphaData',~isnan(ym)); box off;
%
%By Patrick Mineault
%Refs: http://xcorr.net/?p=3326
%      http://www-pord.ucsd.edu/~matlab/bin.htm
% Modified by Bastien Y. Queste to account for NaNs and to permute the matrix at the end.

good = ~isnan(x1+x2+y);

[~,whichedge1] = histc(x1(good),x1rg(:)');
[~,whichedge2] = histc(x2(good),x2rg(:)');

bins1 = min(max(whichedge1,1),length(x1rg)-1);
bins2 = min(max(whichedge2,1),length(x2rg)-1);

bins = (bins2-1)*(length(x1rg)-1)+bins1;

xpos = ones(size(bins,1),1);
ns = sparse(bins,xpos,1,(length(x1rg)-1)*(length(x2rg)-1),1);
ysum = sparse(bins,xpos,y(good),(length(x1rg)-1)*(length(x2rg)-1),1);
ym = full(ysum)./(full(ns));
yb = ym(bins);
ym = reshape(ym,length(x1rg)-1,length(x2rg)-1)';
end