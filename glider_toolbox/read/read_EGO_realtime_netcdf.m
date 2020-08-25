function [data] = read_EGO_realtime_netcdf(missioninfo)
%
%=========================================================================%
% Function to read glider data from EGO netcdf realtime file (concatenated 
% timeseries file) 
%
% L. Houpert, SAMS 25/05/2016
%=========================================================================%
%
%  [data] = read_EGO_realtime_netcdf(missioninfo)
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
% Details about the EGO format available here: http://www.ego-network.org/dokuwiki/doku.php?id=public:datamanagement
%
%  Notes on the EGO netcdf format (01/06/2016):
% - The position of the glider underwater is filled by NaN, this function
% perform a linear interpolation of the latitude and longitude of the
% glider underwater according to the pre-dive and post-dive GPS positioning 
% - In this version no Depth-Averaged Current is presented in the file so all
% variable related to speed are set to nan; 
% - The number of the dives are not stored in the netcdf (= the number of
% the raw dive files) so here the dive number correspond to the dive stored
% in the netcdf
% - No QC applied on the data (navigation and science): the QC variables
% should be at 0 instead of  _FillValue
% - No Flags for the TIME_GPS, LAT_GPS, LON_GPS, but a variable POSITION_QC
% which is a flag of the same size of LONGITUDE and LATITUDE (Is it useful? as the 
% underwater position are determined by the flight model ...), the position QC
% should be on the fixes GPS variable
% - No Altimeter informations in EGO netcdf 
% - No information about the horizontal and vertical current underwater but
% there are (empty) variables for the glider position underwater (LONGITUDE and
% LATITUDE), as both of them are estimated from the flight model they should be present both or not 
% - No information about the surface currents
% - if the glider didn't take measurement during the ascendent phase: no
% data in the netcdf
% 
% new remarks (27/05/2016):
%   - if the timeseries ended with a descent phase (for example: sdeep00 -
%   Canales-sep2013), how to get the position when the glider is surfacing?
%   as the phase 0 may not be present...
%
% created by L. Houpert (houpertloic@gmail.com), 25/05/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
%  
% If one wants to add others field to the output structure (biogeochemical sensor data, flight model output,...), just modify the
% script. Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote
% repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)
%
% License: This code is licensed under the terms of the GNU General Public License v3.0 ( When distributing derived works, the source code of the work must be made available under the same license).


   TSprocesstouse   = missioninfo.rawdata.TSlevel;
   flightmodeltouse = missioninfo.rawdata.flightmod;


disp(['TS data to load: ''' TSprocesstouse ''''])
disp(['Flight Model data to load: ''' flightmodeltouse ''''])

glname           = missioninfo.glname;
glmission        = missioninfo.glmission;
dataDir          = missioninfo.rawdata.dataDir;
                 
nclist=dir([dataDir '/*.nc']);
ncfilename = {nclist.name};

glider.name            = glname ; 
glider.mission         = glmission ; 
glider.TSdataused      = TSprocesstouse;
glider.flightmodelused = flightmodeltouse;

egofilepath = [dataDir '/' ncfilename{1}];

% ncdisp(egofilepath)
info=ncinfo(egofilepath);
for ikk = 1:length(info.Variables)
    egodata.(info.Variables(ikk).Name) = ncread(egofilepath,info.Variables(ikk).Name);
end

%==========================================================================
% Split the timeseries of the deployment into glider dive cycle (dive+climb)
%
% EGO variable for the phase of the glider trajectory: % 0= surface drift; 1= descent; 2= subsurface drift; 3= inflexion; 4=ascent
%[uniqphase,iaphase,icphase] = unique(egodata.PHASE); 
[uniqphasenber,iaphasenber,icphasenber] = unique(egodata.PHASE_NUMBER); 
iaphasenber(isnan(uniqphasenber))   = [];
uniqphasenber(isnan(uniqphasenber)) = []; % remove nan if there is some in PHASE_NUMBER
indphase = cell(size(iaphasenber(1:end-1)));
phase    = nan(size(iaphasenber(1:end-1)));
for ijk= 1:length(iaphasenber(1:end-1)) 
    indphase{ijk} = iaphasenber(ijk):1:(iaphasenber(ijk+1)-1);
    phase(ijk)    = egodata.PHASE(iaphasenber(ijk));
end
% remove nan line when the PHASE was not defined:
indphase(isnan(phase)) = [];
phase(isnan(phase))    = [];

% To Match EGO_NETCDF phase definition: 0= surface drift; 1= descent; 2= subsurface drift; 3= inflexion; 4=ascent
flagsurf = 0;
flagdesc = 1;
flagasc  = 4;

% need to take care of the case when the glider is measuring only in
% descent or only in ascent
if  isempty(find(phase == flagasc)) & ~isempty(find(phase == flagdesc))  
    % no ascent data: fix gps after dive cycle is the same than the one at the beginning of the next surface drift or descent dive same 
    iphase1 = find(phase == flagdesc);
    divecyclestart = nan(1,length(iphase1)-1);
    divecycleend   = nan(1,length(iphase1)-1);   
    for icc = 1:(length(iphase1)-1) % don't take the last descent as we don't know when the glider is surfacing (if there is no phase 0)
        istart = iphase1(icc);
        if ~isempty(find(phase(istart:end) == flagsurf)) 
            if find(phase(istart+1:end) == flagsurf,1) < find(phase(istart+1:end) == flagdesc,1) % if there is a surface drift phase before the next descent
                iend   = istart + find(phase(istart+1:end) == flagsurf,1);
            else
                iend   = istart + find(phase(istart+1:end)==flagdesc,1);
            end
        else % no surface drift phase available
            iend   = istart + find(phase(istart+1:end)==flagdesc,1);
        end
        divecyclestart(icc) = indphase{istart}(1);
        divecycleend(icc)   = indphase{iend}(1);        
    end
    
elseif ~isempty(find(phase==flagasc)) & isempty(find(phase==flagdesc)) % no descent data
    iphase4 = find(phase == flagasc);
    divecyclestart = nan(1,length(iphase4)-1);
    divecycleend   = nan(1,length(iphase4)-1);   
    for icc = 1:(length(iphase4)-1) 
        iend = iphase4(icc+1); % don't take the first ascent as we don't know when the glider was surfacing before (if there is no phase 0)
        divecycleend(icc)   = indphase{iend}(end);  
        indtolook = 1:(iend-1);
        if ~isempty(find(phase(indtolook) == flagsurf))
            if find(phase(indtolook) == flagsurf,1,'last') > find(phase(indtolook) == flagasc,1,'last') % if there is a surface drift phase before the last ascent
                istart = find(phase(indtolook) == flagsurf,1,'last');
                divecyclestart(icc) = indphase{istart}(end); 
            else 
                istart = find(phase(indtolook) == flagasc,1,'last');  
                divecyclestart(icc) = indphase{istart}(end);               
            end
        else
            istart = find(phase(indtolook) == flagasc,1,'last');  
            divecyclestart(icc) = indphase{istart}(end);            
        end
    end         
elseif ~isempty(find(phase==flagasc)) & ~isempty(find(phase==flagdesc)) % decent and ascent data in the EGO file but maybe not also following each other
    iphase1 = find(phase == flagdesc);
    divecyclestart = nan(1,length(iphase1)-1);
    divecycleend   = nan(1,length(iphase1)-1);   
    for icc = 1:(length(iphase1)-1) % don't take the last descent as we don't know when the glider is surfacing (if there is no phase 0)
        istart = iphase1(icc);
        if ~isempty(find(phase(istart:end) == flagsurf)) & ~isempty(find(phase(istart:end) ==flagasc)) % if ascent and surf phase available
            if find(phase(istart+1:end) == flagsurf,1) < find(phase(istart+1:end) == flagdesc,1) % if there is a surface drift phase before the next descent
                iend   = istart + find(phase(istart+1:end) == flagsurf,1);
                divecycleend(icc)   = indphase{iend}(1);    
            elseif find(phase(istart+1:end) == flagasc,1) < find(phase(istart+1:end) == flagdesc,1)
                iend   = istart + find(phase(istart+1:end)==flagasc,1);
                divecycleend(icc)   = indphase{iend}(end);  
            else
                iend   = istart + find(phase(istart+1:end)==flagdesc,1);
                divecycleend(icc)   = indphase{iend}(1);                    
            end
        elseif ~isempty(find(phase(istart:end) == flagsurf)) & isempty(find(phase (istart:end) ==flagasc)) % if no ascent phase
             if find(phase(istart+1:end) == flagsurf,1) < find(phase(istart+1:end) == flagdesc,1) % if there is a surface drift phase before the next descent
                iend   = istart + find(phase(istart+1:end) == flagsurf,1);
                divecycleend(icc)   = indphase{iend}(1);    
            else
                iend   = istart + find(phase(istart+1:end)==flagdesc,1);
                divecycleend(icc)   = indphase{iend}(1);                    
             end 
        elseif isempty(find(phase(istart:end) == flagsurf)) & ~isempty(find(phase (istart:end) ==flagasc)) % if no surface phase
            if find(phase(istart+1:end) == flagasc,1) < find(phase(istart+1:end) == flagdesc,1)
                iend   = istart + find(phase(istart+1:end)==flagasc,1);
                divecycleend(icc)   = indphase{iend}(end);  
            else
                iend   = istart + find(phase(istart+1:end)==flagdesc,1);
                divecycleend(icc)   = indphase{iend}(1);                    
            end           
        else % no surface drift and ascent phase available
            iend   = istart + find(phase(istart+1:end)==flagdesc,1);
            divecycleend(icc)   = indphase{iend}(1);            
        end
        divecyclestart(icc) = indphase{istart}(1);
    end
else
    disp(' ')
    disp('No ascent and descent data in the EGO netcdf file')
    disp(' ')
    data = [];
    return
end
%========================================================================

%-------------------------------------------------------------------------
% add bad flag position 
%-------------------------------------------------------------------------

tgps2   = egodata.TIME(divecyclestart);
latgps2 = egodata.LATITUDE(divecyclestart);
longps2 = egodata.LONGITUDE(divecyclestart);
tgpse   = egodata.TIME(divecycleend);
latgpse = egodata.LATITUDE(divecycleend);
longpse = egodata.LONGITUDE(divecycleend);

qc_gps2 = zeros(1,length(longps2));
qc_gpse = zeros(1,length(longpse));

ibad = find((tgpse-tgps2)>3600 * 8); % if cycle > 8 hours can consider as a valid dive+climb cycle
qc_gps2(ibad) = 4;
qc_gpse(ibad) = 4;

qc_temp = egodata.TEMP_QC; 
qc_temp(isnan(qc_temp)) = 0;
qc_cond = egodata.CNDC_QC; 
qc_cond(isnan(qc_cond)) = 0;
qc_sal = egodata.PSAL_QC; 
qc_sal(isnan(qc_sal)) = 0;

pres    = egodata.PRES;
idiveok = [];
for itt = 1:length(tgps2)
    inddive{itt} = find(egodata.TIME>= tgps2(itt) & egodata.TIME< tgpse(itt));
    indstartclimb{itt} = inddive{itt}(find(pres(inddive{itt}) == max(pres(inddive{itt})),1,'last'));    
    if ~isempty(find(egodata.TIME>= tgps2(itt) & egodata.TIME< tgpse(itt)))
        idiveok = [idiveok itt];
    end
end

%--------------------------------------------------------------
% Load general mission parameter from the 1st nc file
glider.version = str2num(egodata.FIRMWARE_VERSION_NAVIGATION');
glider.basestation_version = str2num(egodata.FIRMWARE_VERSION_SCIENCE');
glider.type = egodata.PLATFORM_TYPE';
glider.snum = str2num(egodata.GLIDER_SERIAL_NO'); 
glider.missionnum = [] ; 


%--------------------------------------------------------------
iic = 0;
for divenb = idiveok
    iic = iic + 1 ;
    
    % General information about the vehicle       
    glider.dive(iic).divenum =  iic ;  % dive number is not store in the netcdf
    glider.dive(iic).starttime = tgps2(divenb)/24/3600 + datenum(1970,01,01,00,00,00);
    glider.dive(iic).start_climb_time = egodata.TIME(indstartclimb{divenb})/24/3600 + datenum(1970,01,01,00,00,00);       
    glider.dive(iic).D_TGT = nan;

    % Information about the location of the profile      
    glider.dive(iic).timegps2 = tgps2(divenb)/24/3600 + datenum(1970,01,01,00,00,00);
    glider.dive(iic).latgps2  = latgps2(divenb);      
    glider.dive(iic).longps2  = longps2(divenb);        
    glider.dive(iic).timegps  = tgpse(divenb)/24/3600 + datenum(1970,01,01,00,00,00);
    glider.dive(iic).latgps   = latgpse(divenb);      
    glider.dive(iic).longps   = longpse(divenb); 
    
    glider.dive(iic).gps2_qc  = qc_gps2(divenb);        
    glider.dive(iic).gpse_qc  = qc_gpse(divenb);        
           
    glider.dive(iic).ctdtime = egodata.TIME((inddive{divenb}))/24/3600 + datenum(1970,01,01,00,00,00);   
  
    glider.dive(iic).divephase = egodata.PHASE(inddive{divenb});
    glider.dive(iic).latitude = egodata.LATITUDE(inddive{divenb});        
    glider.dive(iic).longitude= egodata.LONGITUDE(inddive{divenb});      
 
    
    glider.dive(iic).ctddepth = sw_dpth(egodata.PRES((inddive{divenb})),glider.dive(iic).latgps2);

    if strcmp(TSprocesstouse,'raw')   == 1
        glider.dive(iic).temperature    = egodata.TEMP((inddive{divenb})); % temperature raw
        glider.dive(iic).temperature_qc = qc_temp((inddive{divenb}));        
        glider.dive(iic).conductivity   = egodata.CNDC((inddive{divenb})); % conductivity raw
        glider.dive(iic).conductivity_qc = qc_cond((inddive{divenb}));
        glider.dive(iic).salinity       = egodata.PSAL((inddive{divenb})); % salinity raw
        glider.dive(iic).salinity_qc    = qc_sal((inddive{divenb})); 
        if isfield(egodata,'WPHASE_DOXY')
            glider.dive(iic).doxy       = egodata.WPHASE_DOXY((inddive{divenb}));
            glider.dive(iic).doxy_qc    = egodata.WPHASE_DOXY_QC((inddive{divenb}));
        end        
        if isfield(egodata,'CHLA_SCALED')
            glider.dive(iic).chla       = egodata.CHLA_SCALED((inddive{divenb}));
            glider.dive(iic).chla_qc    = egodata.CHLA_SCALED_QC((inddive{divenb}));
        end
        if isfield(egodata,'BB_SCALED')
            glider.dive(iic).bb       = egodata.BB_SCALED((inddive{divenb}));
            glider.dive(iic).bb_qc    = egodata.BB_SCALED_QC((inddive{divenb}));
        end   
        if isfield(egodata,'CDOM_SCALED')
            glider.dive(iic).cdom       = egodata.CDOM_SCALED((inddive{divenb}));
            glider.dive(iic).cdom_qc    = egodata.CDOM_SCALED_QC((inddive{divenb}));
        end          
    else
        disp('Bad definition of the TS variables to load')
        disp('Possible values for variabletoload.TSdata: ''processed'' or ''raw''')
        break
    end
    glider.dive(iic).theta = sw_ptmp(glider.dive(iic).salinity,glider.dive(iic).temperature,egodata.PRES((inddive{divenb})),0); % potential temp.         
    glider.dive(iic).buoyancy = nan(size(glider.dive(iic).theta)); % buoyancy of vehicle, corrected for compression effects [g], not in EGO netcdf           
    glider.dive(iic).ctd_qc = 0*qc_temp((inddive{divenb})); % no such variable in EGO netcdf

    % NO Altimer in EGO_NETCDF        
    glider.dive(iic).pilotparameter.altim_pulse = nan;
    glider.dive(iic).pilotparameter.altim_sensitivity = nan;           
    glider.dive(iic).altim_bottom_ping =  [nan nan];     

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
    elseif strcmp(flightmodeltouse,'NA') == 1 
        glider.dive(iic).horz_speed = nan*glider.dive(iic).ctdtime;           
        glider.dive(iic).vert_speed = nan*glider.dive(iic).ctdtime;            
        glider.dive(iic).glide_angle = nan*glider.dive(iic).ctdtime;    
        % Positions based on displacements and computed depth-average
        % current (DAC)             
        glider.dive(iic).DAC_east = nan*glider.dive(iic).timegps2;        
        glider.dive(iic).DAC_north = nan*glider.dive(iic).timegps2;      
        glider.dive(iic).DAC_qc = zeros(size(glider.dive(iic).timegps2));        
    else
        disp('Bad definition of the Flight Model data to load')
        disp('Possible values for variabletoload.flightmodel: ''HDM'' or ''GSM''')
        break
    end
    
    if strcmp(flightmodeltouse,'HDM') == 1 | strcmp(flightmodeltouse,'GSM') == 1  
        depthgl = glider.dive(iic).ctddepth;
        timegl  = glider.dive(iic).ctdtime;
        dPdt = [nan; (depthgl(3:end)-depthgl(1:end-2))./(3600*24*(timegl(3:end)-timegl(1:end-2)));nan];
    
        glider.dive(iic).wwater = -dPdt - glider.dive(iic).vert_speed/100;
    
        % Computed surface current
        glider.dive(iic).surface_curr_east = results.surface_curr_east; % Eastward component of surface current [cm/s]        
        glider.dive(iic).surface_curr_north = results.surface_curr_north; % Northward component of surface current [cm/s]   
        glider.dive(iic).surface_curr_qc = str2num(results.surface_curr_qc);     
        
    elseif strcmp(flightmodeltouse,'NA') == 1
        
        glider.dive(iic).wwater = nan*glider.dive(iic).ctdtime;
        % Computed surface current
        glider.dive(iic).surface_curr_east  = nan*glider.dive(iic).timegps2;
        glider.dive(iic).surface_curr_north = nan*glider.dive(iic).timegps2;  
        glider.dive(iic).surface_curr_qc    = zeros(size(glider.dive(iic).timegps2));           
        
    end
        
end      

data = glider;


end
