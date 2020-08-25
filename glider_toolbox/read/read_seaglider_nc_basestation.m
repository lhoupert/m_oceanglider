function [data,processlog] = read_seaglider_nc_basestation(missioninfo)
%
%=========================================================================%
% Function to read seaglider data from netcdf file (using Koengsberd
% DiveData routines, included in the toolbox)
%
% L. Houpert, SAMS 22/01/2016
%=========================================================================%
%
% [data,processlog] = read_seaglider_nc_basestation(missioninfo,(variabletoload))
% Function that read seaglider data (TS, flightmodel) from the basestation
% netcdf and concatenate all the dive in a structure. 
%   Inputs: - missioninfo: [1x1] structure with mission details
%                where missioninfo.glname     : name of the glider (e.g. sg605);
%                      missioninfo.glmission : name of the mission (e.g. OSNAP1);
%                      missioninfo.rawdata.dataDir : path of the netcdf files processed by the basestation;
%                      missioninfo.rawdata.TSlevel  = 'processed' or 'raw'
%                      missioninfo.rawdata.flightmod = 'HDM' or 'GSM'
%  Output:  - data : structure with the glider data from log and eng files
%           - processlog: structure with 2 fields giving the dive number of 
%           the missing nc files or the nc files without temperature
% 
% created by L. Houpert (houpertloic@gmail.com), 22/01/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
%  
% If one wants to add others field to the output structure (biogeochemical sensor data, flight model output,...), just modify the
% script. Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote
% repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)
%



   TSprocesstouse   = missioninfo.rawdata.TSlevel;
   flightmodeltouse = missioninfo.rawdata.flightmod;


disp(['TS data to load: ''' TSprocesstouse ''''])
disp(['Flight Model data to load: ''' flightmodeltouse ''''])

glname           = missioninfo.glname;
glmission        = missioninfo.glmission;
dataDir          = missioninfo.rawdata.dataDir;
  
nclist=dir([dataDir '/p' glname(end-2:end) '*.nc']);
ncfilesname = {nclist.name};
if isempty(ncfilesname)
    disp(' ')
    disp(['Error: No nc files in ' dataDir ])
    disp(' ')
end
maxdivenb = str2num(ncfilesname{end}(end-6:end-3));

glider.name            = glname ; 
glider.mission         = glmission ; 
glider.TSdataused      = TSprocesstouse;
glider.flightmodelused = flightmodeltouse;
%--------------------------------------------------------------
% Load general mission parameter from the 1st nc file
divenb = 200;
id_dive = [glname(end-2:end) num2str(divenb,'%04d')];
[loginfo,eng,results] = dd_get_dive_data(id_dive,1,dataDir);
glider.version = eng.version ;
glider.type = 'Seaglider';
glider.basestation_version = eng.basestation_version ;
glider.snum = eng.glider ; 
glider.missionnum = eng.mission ; 

%--------------------------------------------------------------
iic = 0;
list_no_ncfile = [];
list_no_temperature_field = [];
for divenb = 1:maxdivenb
    if rem(divenb,50)==0
        id_dive = [glname(end-2:end) num2str(divenb,'%04d')];
        disp(['Dive ' num2str(id_dive) ' / '  glname(end-2:end) num2str(maxdivenb,'%04d')])
    else
        id_dive = [glname(end-2:end) num2str(divenb,'%04d')];        
    end
    ssf = strfind(ncfilesname,['p' id_dive '.nc' ]);
    if ~isempty(cell2mat(ssf))
        
            [loginfo,eng,results] = dd_get_dive_data(id_dive,1,dataDir);
            
        if ~isfield(results,'temperature')
            disp([' No temperature field for the dive ' num2str(divenb,'%04d') ])
            list_no_temperature_field = [ list_no_temperature_field divenb];            
        else
            iic = iic + 1;    
            % General information about the vehicle       
            glider.dive(iic).divenum =  eng.dive ;  
            glider.dive(iic).starttime = eng.start_ts_dn;
            glider.dive(iic).start_climb_time = eng.start_ts_dn + results.start_of_climb_time/3600/24; % seconds after dive start when the 2nd apogee pump starts [s]       
            glider.dive(iic).D_TGT = loginfo.D_TGT;
            
            % Information about the location of the profile      
            glider.dive(iic).timegps1 = loginfo.GPS1_dn;
            glider.dive(iic).latgps1 = loginfo.GPS1_lat;      
            glider.dive(iic).longps1 = loginfo.GPS1_lon;
            glider.dive(iic).timegps2 = loginfo.GPS2_dn;
            glider.dive(iic).latgps2 = loginfo.GPS2_lat;      
            glider.dive(iic).longps2 = loginfo.GPS2_lon;        
            glider.dive(iic).timegps = loginfo.GPS_dn;
            glider.dive(iic).latgps = loginfo.GPS_lat;      
            glider.dive(iic).longps = loginfo.GPS_lon;        
            glider.dive(iic).gps1_qc = str2num(results.GPS1_qc); 
            glider.dive(iic).gps2_qc = str2num(results.GPS2_qc);        
            glider.dive(iic).gpse_qc = str2num(results.GPSE_qc);        

            glider.dive(iic).latitude = results.latitude;        
            glider.dive(iic).longitude = results.longitude; 
            
            % Results based on CTD measurements           
            glider.dive(iic).ctdtime = (results.ctd_time)/24/3600 + datenum(1970,01,01,00,00,00);   
            glider.dive(iic).ctddepth = results.ctd_depth; % ctd depth below the surface corrected from the average latitude [m]
            
            if strcmp(TSprocesstouse,'processed') == 1
                glider.dive(iic).temperature = results.temperature; % temperature corrected for the thermistor 1st order lag [deg C]
                glider.dive(iic).temperature_qc = str2num(results.temperature_qc);        
                glider.dive(iic).conductivity = results.conductivity; % conductivity corrected for anomalies [ms/cm]
                glider.dive(iic).conductivity_qc = str2num(results.conductivity_qc);
                glider.dive(iic).salinity = results.salinity; % salinity corrected for thermal-inertia effects [PSU]
                glider.dive(iic).salinity_qc = str2num(results.salinity_qc);       
            elseif strcmp(TSprocesstouse,'raw')   == 1
                glider.dive(iic).temperature = results.temperature_raw; % temperature raw
                glider.dive(iic).temperature_qc = str2num(results.temperature_raw_qc);        
                glider.dive(iic).conductivity = results.conductivity_raw; % conductivity raw
                glider.dive(iic).conductivity_qc = str2num(results.conductivity_raw_qc);
                glider.dive(iic).salinity = results.salinity_raw; % salinity raw
                glider.dive(iic).salinity_qc = str2num(results.salinity_raw_qc);           
            else
                disp('Bad definition of the TS variables to load')
                disp('Possible values for variabletoload.TSdata: ''processed'' or ''raw''')
                break
            end
            glider.dive(iic).theta = results.theta; % potential temp.         
            glider.dive(iic).buoyancy = results.buoyancy; % buoyancy of vehicle, corrected for compression effects [g]           
            glider.dive(iic).ctd_qc = str2num(results.CTD_qc); 

            % Altimeter configuration and results       
            glider.dive(iic).pilotparameter.altim_pulse = loginfo.ALTIM_PULSE;
            glider.dive(iic).pilotparameter.altim_sensitivity = loginfo.ALTIM_SENSITIVITY;           
            if isfield(loginfo,'ALTIM_BOTTOM_PING')
               glider.dive(iic).altim_bottom_ping =  loginfo.ALTIM_BOTTOM_PING;
            else
               glider.dive(iic).altim_bottom_ping =  [nan nan];
            end      

            if strcmp(flightmodeltouse,'HDM') == 1
                % Vehicle velocities and displacements
                glider.dive(iic).horz_speed = results.horz_speed; % vehicle horizontal speed based on hdm [cm/s]
                glider.dive(iic).vert_speed = results.vert_speed; % vehicle vertical speed based on hdm [cm/s]    
                glider.dive(iic).speed_qc = str2num(results.speed_qc);  
                glider.dive(iic).glide_angle = results.glide_angle;
                % Positions based on displacements and computed depth-average
                % current (DAC)
                glider.dive(iic).DAC_east = results.depth_avg_curr_east; % Eastward component of DAC based on hdm [m/s]        
                glider.dive(iic).DAC_north = results.depth_avg_curr_north; % Northward component of DAC based on hdm [m/s]       
                glider.dive(iic).DAC_qc = str2num(results.depth_avg_curr_qc);
            elseif strcmp(flightmodeltouse,'GSM') == 1         
                % Vehicle velocities and displacements            
                glider.dive(iic).horz_speed = results.horz_speed_gsm;           
                glider.dive(iic).vert_speed = results.vert_speed_gsm;            
                glider.dive(iic).glide_angle = results.glide_angle_gsm;    
                % Positions based on displacements and computed depth-average
                % current (DAC)             
                glider.dive(iic).DAC_east = results.depth_avg_curr_east_gsm; % Eastward component of DAC based on gsm [m/s]        
                glider.dive(iic).DAC_north = results.depth_avg_curr_north_gsm; % Northward component of DAC based on gsm [m/s]       
                glider.dive(iic).DAC_qc = str2num(results.depth_avg_curr_qc)*0;
             else
                disp('Bad definition of the Flight Model data to load')
                disp('Possible values for variabletoload.flightmodel: ''HDM'' or ''GSM''')
                break
            end
            depthgl = glider.dive(iic).ctddepth;
            timegl  = glider.dive(iic).ctdtime;
            dPdt = [nan; (depthgl(3:end)-depthgl(1:end-2))./(3600*24*(timegl(3:end)-timegl(1:end-2)));nan];
            glider.dive(iic).wwater = -dPdt - glider.dive(iic).vert_speed/100;
            % Computed surface current
            glider.dive(iic).surface_curr_east = results.surface_curr_east; % Eastward component of surface current [cm/s]        
            glider.dive(iic).surface_curr_north = results.surface_curr_north; % Northward component of surface current [cm/s]   
            glider.dive(iic).surface_curr_qc = str2num(results.surface_curr_qc);         
            
            
            if isfield(results,'dissolved_oxygen_sat')
                glider.dive(iic).doxy       = results.dissolved_oxygen_sat;
                glider.dive(iic).doxy_qc    = glider.dive(iic).doxy*0;
            end
             
        end      
        
    else
        disp([' No NetCdf File for the dive ' num2str(divenb,'%04d') ])
        list_no_ncfile = [ list_no_ncfile divenb];
    end
end

data = glider;

processlog.list_no_ncfile = list_no_ncfile;
processlog.list_no_temperature_field = list_no_temperature_field;

end
