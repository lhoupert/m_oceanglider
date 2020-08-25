function glider = convert_UEA_struc_to_common_glider_struc(data,missioninfo,proc_param)
%
%==========================================================================
% Conversion of matlab structure from UEA-toolbox to generic matlab structure (similar to
% the one loaded with read_seaglider_nc_basestation.m
%
% function glider = convert_UEA_struc_to_common_glider_struc(data,missioninfo,proc_param);
%
% by L. Houpert 25/01/2016
%

maxdivenb = length(data.eng);
glider.name = missioninfo.glname;
glider.mission = missioninfo.glmission;
%--------------------------------------------------------------
% Load general mission parameter from the 1st nc file
divenb = 200;
glider.version             = data.eng(divenb).version ;
glider.basestation_version = data.eng(divenb).basestation_version ;
glider.type                = 'Seaglider';
glider.snum                = data.eng(divenb).glider ; 
glider.missionnum          = data.eng(divenb).mission ; 
glider.TSdataused          = 'Processed#UEA-toolbox';

if missioninfo.rawdata.flight_regression_passes == 0
    glider.flightmodelused     = 'HDM#no_opti';
else
    glider.flightmodelused     = ['HDM#opti_' ...
                    num2str(missioninfo.rawdata.flight_regression_passes) '_regres_passes'] ;
end
%--------------------------------------------------------------
iic = 0;
for divenb = 1:maxdivenb
            iic = divenb;%iic + 1;   
            
	    if isempty(data.hydrography(divenb).time)
	                % General information about the vehicle       
            glider.dive(iic).divenum = divenb ;
            glider.dive(iic).starttime = nan;
            glider.dive(iic).start_climb_time = nan; % seconds after dive start when the 2nd apogee pump starts [s]       
            glider.dive(iic).D_TGT = nan;
            
            % Information about the location of the profile      
            glider.dive(iic).latgps1 =	nan;
            glider.dive(iic).longps1 =	nan;
            glider.dive(iic).timegps1 =	nan;
            
            glider.dive(iic).latgps2 =	nan;           
            glider.dive(iic).longps2 =	nan;
            glider.dive(iic).timegps2 =	nan;
            
            glider.dive(iic).latgps =	nan;
            glider.dive(iic).longps =	nan;
            glider.dive(iic).timegps =	nan;
                   
            glider.dive(iic).gps1_qc =  nan;
            glider.dive(iic).gps2_qc =  nan;
            glider.dive(iic).gpse_qc =  nan;

            % Results based on CTD measurements
            glider.dive(iic).latitude = [nan nan nan nan]';      
            glider.dive(iic).longitude = [nan nan nan nan]';             
            glider.dive(iic).ctdtime = [nan nan nan nan]';
            glider.dive(iic).ctddepth = [1 5 10 100]';
            glider.dive(iic).pressure = [1 5 10 100]';
            glider.dive(iic).temperature = [nan nan nan nan]';
            glider.dive(iic).temperature_qc = glider.dive(iic).temperature * 9;       
            glider.dive(iic).conductivity =[nan nan nan nan]';
            glider.dive(iic).conductivity_qc =  glider.dive(iic).temperature * 9;       
            glider.dive(iic).salinity = [nan nan nan nan]';
            glider.dive(iic).salinity_qc =  glider.dive(iic).temperature * 9;       
            glider.dive(iic).theta = [nan nan nan nan]';
            glider.dive(iic).cons_temp = [nan nan nan nan]';
            glider.dive(iic).abs_salinity = [nan nan nan nan]';
            glider.dive(iic).buoyancy = [nan nan nan nan]'; % buoyancy of vehicle, corrected for compression effects [g]           
            glider.dive(iic).ctd_qc =  glider.dive(iic).temperature * 9;       
            % Altimeter configuration and results       
            glider.dive(iic).pilotparameter.altim_pulse = data.log(divenb).ALTIM_PULSE;
            glider.dive(iic).pilotparameter.altim_sensitivity = data.log(divenb).ALTIM_SENSITIVITY;           
            if ~isempty(data.log(divenb).ALTIM_BOTTOM_PING )%isfield(data.log(divenb),'ALTIM_BOTTOM_PING')
               glider.dive(iic).altim_bottom_ping =  data.log(divenb).ALTIM_BOTTOM_PING;
            else
               glider.dive(iic).altim_bottom_ping =  [nan nan];
            end      

            % Vehicle velocities and displacements
            glider.dive(iic).horz_speed = [nan nan nan nan]'; % vehicle horizontal speed based on hdm [cm/s]
            glider.dive(iic).vert_speed = [nan nan nan nan]'; %vehicle vertical speed based on hdm [cm/s]    
            glider.dive(iic).speed_qc = glider.dive(iic).horz_speed*9;  
            glider.dive(iic).glide_angle =[nan nan nan nan]';
            
            % Positions based on displacements and computed depth-average
            % current (DAC)
            glider.dive(iic).DAC_east = nan ; 
            glider.dive(iic).DAC_north =nan;    
            glider.dive(iic).DAC_qc = 9;
            glider.dive(iic).wwater = [nan nan nan nan]';
		continue
	    else
            % General information about the vehicle       
            glider.dive(iic).divenum = data.eng(divenb).dive ;  
            glider.dive(iic).starttime = data.hydrography(divenb).time(1);
            glider.dive(iic).start_climb_time = data.hydrography(divenb).time(data.hydrography(divenb).max_pressure_index); % seconds after dive start when the 2nd apogee pump starts [s]       
            glider.dive(iic).D_TGT = data.log(divenb).D_TGT;
            
            % Information about the location of the profile      
            glider.dive(iic).latgps1 =	data.gps_postpreviousdive(divenb,1);
            glider.dive(iic).longps1 =	data.gps_postpreviousdive(divenb,2);
            glider.dive(iic).timegps1 =	data.gps_postpreviousdive(divenb,3);
            
            glider.dive(iic).latgps2 =	data.gps_predive(divenb,1);           
            glider.dive(iic).longps2 =	data.gps_predive(divenb,2);            
            glider.dive(iic).timegps2 =	data.gps_predive(divenb,3);
             
            glider.dive(iic).latgps =	data.gps_postdive(divenb,1);           
            glider.dive(iic).longps =	data.gps_postdive(divenb,2);            
            glider.dive(iic).timegps =	data.gps_postdive(divenb,3);
                   
            glider.dive(iic).gps1_qc =  data.gps_postdive(divenb,1) * 0;
            glider.dive(iic).gps2_qc =  data.gps_postdive(divenb,1) * 0;    
            glider.dive(iic).gpse_qc =  data.gps_postdive(divenb,1) * 0;

            glider.dive(iic).latitude  = data.flight(divenb).trajectory_latlon_estimated(:,1);
            glider.dive(iic).longitude = data.flight(divenb).trajectory_latlon_estimated(:,2);
            
            % Results based on CTD measurements
            glider.dive(iic).ctdtime = (data.hydrography(divenb).time)';   
            glider.dive(iic).ctddepth = (data.hydrography(divenb).depth)'; 
            glider.dive(iic).pressure = (data.hydrography(divenb).pressure)';            
            glider.dive(iic).temperature = (data.hydrography(divenb).temp)'; 
            glider.dive(iic).temperature_qc = glider.dive(iic).temperature * 0;       
            glider.dive(iic).conductivity = (data.hydrography(divenb).conductivity)'; 
            glider.dive(iic).conductivity_qc = glider.dive(iic).conductivity * 0; 
            glider.dive(iic).salinity = (data.hydrography(divenb).salinity)'; 
            glider.dive(iic).salinity_qc = glider.dive(iic).salinity * 0; 
            glider.dive(iic).theta = sw_ptmp(glider.dive(iic).salinity,glider.dive(iic).temperature, glider.dive(iic).pressure, 0); % potential temp.
            glider.dive(iic).cons_temp = (data.hydrography(divenb).cons_temp)';
            glider.dive(iic).abs_salinity = (data.hydrography(divenb).abs_salinity)';                 
            glider.dive(iic).buoyancy = (data.flight(divenb).buoyancy)'; % buoyancy of vehicle, corrected for compression effects [g]           
            glider.dive(iic).ctd_qc = glider.dive(iic).ctdtime * 0; 

            % Altimeter configuration and results       
            glider.dive(iic).pilotparameter.altim_pulse = data.log(divenb).ALTIM_PULSE;
            glider.dive(iic).pilotparameter.altim_sensitivity = data.log(divenb).ALTIM_SENSITIVITY;           
            if ~isempty(data.log(divenb).ALTIM_BOTTOM_PING )%isfield(data.log(divenb),'ALTIM_BOTTOM_PING')
               glider.dive(iic).altim_bottom_ping =  data.log(divenb).ALTIM_BOTTOM_PING;
            else
               glider.dive(iic).altim_bottom_ping =  [nan nan];
            end      

            % Vehicle velocities and displacements
            glider.dive(iic).horz_speed = (data.flight(divenb).model_horz_spd)'; % vehicle horizontal speed based on hdm [cm/s]
            glider.dive(iic).vert_speed = (data.flight(divenb).model_vert_spd)'; %vehicle vertical speed based on hdm [cm/s]    
            glider.dive(iic).speed_qc = glider.dive(iic).horz_speed *0;  
            glider.dive(iic).glide_angle = (data.flight(divenb).model_slope)';

            % Positions based on displacements and computed depth-average
            % current (DAC)
            glider.dive(iic).DAC_east = data.hydrography(divenb).DAC_u; % Eastward component of DAC based on hdm [m/s]        
            glider.dive(iic).DAC_north = data.hydrography(divenb).DAC_v; % Northward component of DAC based on hdm [m/s]       
            glider.dive(iic).DAC_qc = zeros(size(glider.dive(iic).DAC_north));
            glider.dive(iic).wwater = (data.hydrography(divenb).w_H2O)';


%            % Computed surface current
%            glider.dive(iic).surface_curr_east = results.surface_curr_east; % Eastward component of surface current [cm/s]        
%            glider.dive(iic).surface_curr_north = results.surface_curr_north; % Northward component of surface current [cm/s]   
%            glider.dive(iic).surface_curr_qc = results.surface_curr_qc;        
    	end
end

end
