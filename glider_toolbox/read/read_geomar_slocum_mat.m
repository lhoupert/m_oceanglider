function [data] = read_geomar_slocum_mat(missioninfo)
%
%=========================================================================%
% Function to read slocum glider data created from the geomartoolbox
%
% L. Houpert, SAMS 21/05/2018
%=========================================================================%
%
%  [data] = read_geomar_slocum_mat(missioninfo)
% Function that read glider data from EGO netcdf realtime files
%  Inputs: - missioninfo: [1x1] structure with mission details
%                where missioninfo.glname     : name of the glider (e.g. sg605);
%                      missioninfo.glmission : name of the mission (e.g. OSNAP1);
%                      missioninfo.rawdata.dataDir : path of the netcdf file;
%                      missioninfo.rawdata.TSlevel  = 'raw' (or 'processed' if adjusted data are available in the EGO netcdf, soon?)
%                      missioninfo.rawdata.flightmod = 'NA', 'HDM' or 'GSM'
%  Output:  - data : structure with the glider data from log and eng files
%           - processlog: structure with 2 fields giving the dive number of 
%           the missing nc files or the nc files without temperature
%
% 
%
% created by L. Houpert (houpertloic@gmail.com), 21/05/2018, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
%  
%
% License: This code is licensed under the terms of the GNU General Public License v3.0 ( When distributing derived works, the source code of the work must be made available under the same license).


   TSprocesstouse   = missioninfo.rawdata.TSlevel;
   flightmodeltouse = missioninfo.rawdata.flightmod;


disp(['TS data to load: ''' TSprocesstouse ''''])
disp(['Flight Model data to load: ''' flightmodeltouse ''''])

glname           = missioninfo.glname;
glmission        = missioninfo.glmission;
dataDir          = missioninfo.rawdata.dataDir;
                 

data = load([dataDir '_final_1sec']);
data0 = load([dataDir '_1sec']);
data_der = load([dataDir '_1sec_derived']);
data_vec = load([dataDir '_1vector']);
data_yos = load([dataDir '_yos']);

%=====================================================================
% Depth-Average current checl + cleaning
% (made for Growler on OSNAP7)
DACu = nan(1,length(data_yos.n_lat));
DACv = nan(1,length(data_yos.n_lat));
nbpoint = nan(1,length(data_yos.n_lat));
line_vx = nan(1,length(data_yos.n_lat));
nbpointf = nan(1,length(data_yos.n_lat));
line_fvx = nan(1,length(data_yos.n_lat));
nbpointi = nan(1,length(data_yos.n_lat));
line_ivx = nan(1,length(data_yos.n_lat));
pfdos_id = nan(1,length(data_yos.n_lat));
m_finalwater_vx = nan(1,length(data_yos.n_lat));
m_finalwater_vy = nan(1,length(data_yos.n_lat));
% concatenate the DAC information for all yos
for ivv=1:length(data_yos.n_lat)
    pfdos_id(ivv) = unique(data_yos.dos_id{ivv});
     if length(find(~isnan(data_yos.m_water_vx{ivv})))==1
        line_vx(ivv) = find(~isnan(data_yos.m_water_vx{ivv}));      
        nbpoint(ivv)=length(data_yos.m_water_vx{ivv});       
     end     
      if length(find(~isnan(data_yos.m_final_water_vx{ivv})))>0
        line_fvx(ivv) = find(~isnan(data_yos.m_final_water_vx{ivv}));      
        nbpointf(ivv) = length(data_yos.m_final_water_vx{ivv});      
        m_finalwater_vx(ivv) = data_yos.m_final_water_vx{ivv}(line_fvx(ivv));
        m_finalwater_vy(ivv) = data_yos.m_final_water_vy{ivv}(line_fvx(ivv));       
      end        
      if length(find(~isnan(data_yos.m_initial_water_vx{ivv})))==1
        line_ivx(ivv) = find(~isnan(data_yos.m_initial_water_vx{ivv}));      
        nbpointi(ivv)=length(data_yos.m_initial_water_vx{ivv});       
      end           
end      

% identify the different dive cycles and the associated DAC
[divecycleid,IA,IC] = unique(pfdos_id);
sameid = nan(1,length(divecycleid));
nb_f_wvx = nan(1,length(divecycleid));
m_f_vx  = nan(1,length(divecycleid));
m_f_vy  = nan(1,length(divecycleid));
for ill=1:length(divecycleid)
    iinduniqid = find(pfdos_id==divecycleid(ill));
    sameid(ill)   = length(iinduniqid);
    nb_f_wvx(ill) = length(find(~isnan(line_fvx(iinduniqid))));
    timecycle = [data_yos.m_present_time{iinduniqid}]/86400 +datenum(1970,1,1,0,0,0);
    times(ill)  =   timecycle(1);
    timee(ill)  =   timecycle(end);   
    meanlat(ill)   =  nanmean([data_yos.n_lat{iinduniqid}]);
    meanlon(ill)   =  nanmean([data_yos.n_lon{iinduniqid}]);  
    % select the "good" DAC candidates
    nb_f_wvx(ill) = length(find(~isnan(line_fvx(iinduniqid))));
    if nb_f_wvx(ill)==1 & ...
            (timee(ill) - times(ill))*24<13 % only keep DAC when only one value is available for a dive cyle and for cycle < 13h 
        indok = iinduniqid(~isnan(m_finalwater_vx(iinduniqid)));
        m_f_vx(ill) = m_finalwater_vx(indok);
        m_f_vy(ill) = m_finalwater_vy(indok);    
    end    
end
meantimecycle = times+(timee-times)/2;
deltatime = (timee - times)*24; % in hour
% figure;subplot(3,1,1); plot(meantimecycle, nlon,'+'); 
% subplot(3,1,2);plot(meantimecycle,deltatime,'+') 
% subplot(3,1,3);plot(meantimecycle,sameid,'+') 

%----------------------------------
% remove outliers:
m_f_vxsmooth=smooth(m_f_vx,7,'lowess')';
m_f_vysmooth=smooth(m_f_vy,7,'lowess')';
% figure;
% subplot(2,1,1)
% plot(meantimecycle,m_f_vxsmooth);hold on
% plot(meantimecycle,m_f_vx,'xr')
% subplot(2,1,2)
% plot(meantimecycle,m_f_vysmooth);hold on
% plot(meantimecycle,m_f_vy,'xr')

ibad= union(find(abs(m_f_vx-m_f_vxsmooth)>0.2), ...
                find(abs(m_f_vy-m_f_vysmooth)>0.2));
m_f_vy(ibad)=nan;
m_f_vx(ibad)=nan;
%===================================================================
 
pindex =nan(1,length(data.pressure));
for ijk=1:length(data.pressure)
    iok=find(~isnan(data.pressure{ijk}));
    pindex(ijk)= data.pressure{ijk}(iok(1)) - data.pressure{ijk}(iok(end));
end

baddive = find(pindex(1:2:end)>0); % if climb profile are present instead of a dive
badclimb = find(pindex(2:2:end)<0); % if dive profile are present instead of a climb
if pindex(1)>0
    disp(' ')    
    disp('PROBLEM: 1st profile start by a climb!')
    disp(' ')    
    keyboard
else
    if ~isempty(baddive)
    disp(' ')    
    disp('PROBLEM: at least one of the expected dive profile is ascending')
    disp(' ')          
    keyboard
    elseif ~isempty(badclimb)
    disp(' ')    
    disp('PROBLEM: at least one of the expected climb profile is descending')
    disp(' ')    
    keyboard   
    end
end

isoxydata = 0;
o2index =zeros(1,length(data0.oxygen));
for ijk=1:length(data0.oxygen)
    iok=find(~isnan(data0.oxygen{ijk}));
    if ~isempty(iok)
        o2index(ijk)= 1;
    end
end   
maxdivenb = length(data.time_datenum)/2;

glider.name = missioninfo.glname;
glider.mission = missioninfo.glmission;

%--------------------------------------------------------------
% Load general mission parameter
glider.version             = '';
glider.basestation_version = 'NA';
glider.type                = 'Slocum';
glider.snum                = glname; 
glider.missionnum          = glmission ; 

iic = 0;
for divenb = 1:maxdivenb
            iic = divenb;%iic + 1;   
            idive = 1 + 2*(iic-1);
            iclimb = 0 + 2*iic;
	    if isempty(find(~isnan(data.pressure{divenb}))) % if no data available
	                % General information about the vehicle       
            glider.dive(iic).divenum = divenb ;
            glider.dive(iic).starttime = nan;
            glider.dive(iic).start_climb_time = nan; % seconds after dive start when the 2nd apogee pump starts [s]       
            glider.dive(iic).D_TGT = nan;
            
            % Information about the location of the profile                
            glider.dive(iic).latgps2 =	nan;           
            glider.dive(iic).longps2 =	nan;
            glider.dive(iic).timegps2 =	nan;
            
            glider.dive(iic).latgps =	nan;
            glider.dive(iic).longps =	nan;
            glider.dive(iic).timegps =	nan;
                   
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
            
            if isoxydata == 1
                glider.dive(iic).doxy       = [nan nan nan nan]';
                glider.dive(iic).doxy_qc    = glider.dive(iic).doxy * 9; 
            end               
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
            if nanmean(data.time_datenum{idive}-data.time_datenum{idive}) ~= 0
                disp(' ');
                disp('Problem science and main datenum not equal')
                disp(' ');
            end             
            glider.dive(iic).divenum = iic;
            glider.dive(iic).starttime = data.time_datenum{idive}(1);
            glider.dive(iic).start_climb_time = data.time_datenum{iclimb}(1); % seconds after dive start when the 2nd apogee pump starts [s]       
            glider.dive(iic).D_TGT = nan;
            
            % Information about the location of the profile      
            
            glider.dive(iic).latgps2 =	data.latitude{idive}(1);           
            glider.dive(iic).longps2 =	data.longitude{idive}(1);            
            glider.dive(iic).timegps2 =	data.time_datenum{idive}(1); 
             
            glider.dive(iic).latgps =	data.latitude{iclimb}(end);           
            glider.dive(iic).longps =	data.longitude{iclimb}(end);              
            glider.dive(iic).timegps =	data.time_datenum{iclimb}(end); 
                   
            glider.dive(iic).gps2_qc =  glider.dive(iic).timegps* 0;    
            glider.dive(iic).gpse_qc =  glider.dive(iic).timegps* 0;

            glider.dive(iic).latitude  = [data.latitude{idive} data.latitude{iclimb}]';
            glider.dive(iic).longitude = [data.longitude{idive} data.longitude{iclimb}]';
            
            % Results based on CTD measurements
            glider.dive(iic).ctdtime = [data.time_datenum{idive} data.time_datenum{iclimb}]';   
            
            glider.dive(iic).pressure = [data.pressure{idive} data.pressure{iclimb}]';  
            glider.dive(iic).ctddepth = sw_dpth(glider.dive(iic).pressure, glider.dive(iic).latitude);
            
            glider.dive(iic).temperature = [data.temperature{idive} data.temperature{iclimb}]';
            glider.dive(iic).temperature_qc = glider.dive(iic).temperature * 0;      
            
            glider.dive(iic).conductivity = [data0.ctd_conductivity{idive} data0.ctd_conductivity{iclimb}]';
            glider.dive(iic).conductivity_qc = glider.dive(iic).conductivity * 0; 
            
            glider.dive(iic).salinity = [data.salinity{idive} data.salinity{iclimb}]';
            glider.dive(iic).salinity_qc = glider.dive(iic).salinity * 0; 
            
          
            if isoxydata == 1
                glider.dive(iic).doxy       = [data0.oxygen{idive} data0.oxygen{iclimb}]';
                glider.dive(iic).doxy_qc    = glider.dive(iic).doxy * 0; 
            end                   
  
            glider.dive(iic).ctd_qc = glider.dive(iic).ctdtime * 0; 

            glider.dive(iic).theta = sw_ptmp(glider.dive(iic).salinity,glider.dive(iic).temperature,glider.dive(iic).pressure,0); % potential temp.         
            glider.dive(iic).buoyancy = nan(size(glider.dive(iic).theta)); % buoyancy of vehicle, corrected for compression effects [g], not in EGO netcdf           
            glider.dive(iic).ctd_qc = 0*glider.dive(iic).ctdtime; % no such variable in EGO netcdf

            % NO Altimer data     
            glider.dive(iic).pilotparameter.altim_pulse = nan;
            glider.dive(iic).pilotparameter.altim_sensitivity = nan;           
            glider.dive(iic).altim_bottom_ping =  [nan nan];         
            
            % Vehicle velocities and displacements
            glider.dive(iic).horz_speed = [nan nan nan nan]'; % vehicle horizontal speed based on hdm [cm/s]
            glider.dive(iic).vert_speed = [nan nan nan nan]'; %vehicle vertical speed based on hdm [cm/s]    
            glider.dive(iic).speed_qc = glider.dive(iic).horz_speed*9;  
            glider.dive(iic).glide_angle =[nan nan nan nan]';
            
                        
            glider.dive(iic).wwater = [data_der.measured_w{idive} data_der.measured_w{iclimb}]' - ...
                                     [data_der.modeled_w{idive} data_der.modeled_w{iclimb}]';
 
    	end
end      


iok=find(~isnan(m_f_vx));
for iic = 1:length(glider.dive)
    glider.dive(iic).DAC_east = interp1(meantimecycle(iok),m_f_vx(iok),glider.dive(iic).start_climb_time); % Eastward component of DAC [m/s]     
    glider.dive(iic).DAC_north = interp1(meantimecycle(iok),m_f_vy(iok),glider.dive(iic).start_climb_time); % Northward component of DAC [m/s]       
    glider.dive(iic).DAC_qc = zeros(size(glider.dive(iic).DAC_north));
end

% figure;
% subplot(2,1,1)
% plot([glider.dive(:).start_climb_time],[glider.dive(:).DAC_east],'-')
% hold on
% plot(meantimecycle,m_f_vx,'xr')
% subplot(2,1,2)
% plot([glider.dive(:).start_climb_time],[glider.dive(:).DAC_north],'-')
% hold on
% plot(meantimecycle,m_f_vy,'xr')

data = glider;


end
