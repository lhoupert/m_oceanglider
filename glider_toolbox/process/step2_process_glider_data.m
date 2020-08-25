function glider = step2_process_glider_data(missioninfo,mergedpflevel,step2opt)
%
%=========================================================================%
% Function to process seaglider data: 
%
% L. Houpert, SAMS 08/05/2016
%=========================================================================%
%
% step2_process_glider_data(******)
% Function that apply final processing (detided DAC, merged or not profiles, build a structure of the data.
% Results are save in [missioninfo.step2.procdataDir filesep missioninfo.step2.matfilename]
%
%   Inputs: o missioninfo: [1x1] structure with mission details (defined in
%             users_param/loadosnapmissionparam.m)
%                where - missioninfo.step2.procdataDir : path of the
%                directory to save step2 processing outputs
%                      - missioninfo.step2.matfilename
%                      - missioninfo.section : [kx1] structure indexing the time
%                      of the beginning and end of the k glider sections as a 
%                      [1x2] matrix defined by the users in missioninfo.section(i).time
%                      - missioninfo.rawdata.dataDir : path of the eng and log files used for gps correction;
%
%           o mergedpflevel: - 0 : climb and dive are kept
%                            - 1 : climb and dive merged
%                            - 2 : only dives
%                            - 3 : only climbs
%
%           o step2opt: [1x1] structure containing:
%                     - .bathyname: Bathymetry used to get an interpolation
%                     along the glider track
%                     - .bathydir: Directory where the bathymetry file is located
%                     - .untidedDAC: 1: if one wants to untide the DAC
%                     using a lowpass filter (48h Hamming window)
%                                    2: to untide the DAC using a tidal model
%                                    0: the DAC is not untided
%                                    
%                     otherwise
%                     - .pathtmdtoolboxwithmodel: path of the tmd
%                     toolbox with (tide model files if need)
%
% created by L. Houpert (houpertloic@gmail.com), 08/05/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 


if  exist(step2opt.bathydir,'dir')==7
    usebathyfile =1;
end


%-----------------------------------------
% Create outputdirs if didn't exist
outputdir = missioninfo.step2.procdataDir;  
if exist(outputdir,'dir')~=7   mkdir(outputdir); end %mkdir([outputdir '/bad_dives_detected']);end

%----------------------------------------
% load section info:
sec_time=nan(length(missioninfo.section),2);
for ill=1:length(missioninfo.section)
    sec_time(ill,:) = missioninfo.section(ill).time;
end


%------------------------------------------------------------------------------
eval(['load(''' missioninfo.step1.procdataDir filesep missioninfo.step1.matfilename '.mat'');' ])  

optsensornames = fieldnames(optsensor);

alti_data(1:2:length(timepf)) = altimeter_data;
alti_data(2:2:length(timepf)) = altimeter_data;  
bathy_data = zbathy;   
HH(1:2:length(timepf)) = abs(mergedbathy);
HH(2:2:length(timepf)) = abs(mergedbathy);  
numberofdivecycle(1:2:length(timepf)) = divenum;
numberofdivecycle(2:2:length(timepf)) = divenum; 

%-------------------------------------------------------------------------------
% Choice of different export of vertical grid using the dive+climb matrix
% create in step1
if mergedpflevel == 0
iokpfsel=[];
for ippf = 1:length(gldivephase)-1
    if gldivephase(ippf)==1 & gldivephase(ippf+1)==4 
        iokpfsel = [iokpfsel ippf ippf+1];
    end
end
%-------------------------------------------------------------------------------
% keep dive and climb separate
TIMED = repmat(timepf(iokpfsel),length(PP),1);
PRES = repmat(PP',1,length(timepf(iokpfsel)));
TP2 = TP(:,iokpfsel);
SS2 = SS(:,iokpfsel);
TT2 = TT(:,iokpfsel);
gldivephase2 = gldivephase(iokpfsel);
for izv=1:length(optsensornames) % qc for additional var
    if optsensor.(optsensornames{izv}) == 1
          VAR = upper(optsensornames{izv});
          eval([VAR '2 = ' VAR '(:,iokpfsel);']);
    end
end
timedive = timepf(iokpfsel)';
latdive = latpf(iokpfsel)';
londive = lonpf(iokpfsel)';
timediveini = timepfini(iokpfsel)';
latdiveini = latpfini(iokpfsel)';
londiveini = lonpfini(iokpfsel)';
timediveend = timepfend(iokpfsel)';
latdiveend = latpfend(iokpfsel)';
londiveend = lonpfend(iokpfsel)';
commentdivemerged= 'QC, vertically interpolated, dive and climb are kepts independant, but DAC still calculated for a dive-climb cycle';
% elseif mergedpflevel == 1
% %-------------------------------------------------------------------------------
% % Merged dive and climb 
% TIMED = repmat(timeDAC',length(PP),1);
% PRES = repmat(PP',1,length(tgpse));
% TP2=nan(length(PP),length(tgpse));
% SS2=nan(length(PP),length(tgpse));
% TT2=nan(length(PP),length(tgpse));
% for izv=1:length(optsensornames) 
%    if optsensor.(optsensornames{izv}) == 1
%        VAR = upper(optsensornames{izv});
%        eval([VAR '2 = nan(length(PP),length(tgpse))']);
%    end
% end    
% for icc=0:length(timeDAC)-1
% 	ijk=2*icc+1;
%     gldivephase2(icc+1) = 14;
% 	TT2(:,icc+1)=mean([TT(:,ijk) TT(:,ijk+1)],2);
% 	TP2(:,icc+1)=mean([TP(:,ijk) TP(:,ijk+1)],2);
% 	SS2(:,icc+1)=mean([SS(:,ijk) SS(:,ijk+1)],2);	
%     for izv=1:length(optsensornames) % 
%         if optsensor.(optsensornames{izv}) == 1
%           VAR = upper(optsensornames{izv});
%           eval([VAR '2(:,icc+1)=mean([' VAR '(:,ijk) ' VAR '(:,ijk+1)],2);']);
%         end
%     end    
% end
% timedive = timeDAC;
% latdive = latDAC;
% londive = lonDAC;
% timediveini = timepfini(1:2:end)'; 
% latdiveini = latpfini(1:2:end)'; %#ok<*COLND>
% londiveini = lonpfini(1:2:end)';
% timediveend = timepfend(2:2:end)';
% latdiveend = latpfend(2:2:end)';
% londiveend = lonpfend(2:2:end)';
% commentdivemerged= 'QC, vertically interpolated, dive and climb merged';
% %-------------------------------------------------------------------------------
% % If only dives 
% elseif mergedpflevel == 2
% TIMED = repmat(timeDAC',length(PP),1);
% PRES = repmat(PP',1,length(tgpse));
% TP2=nan(length(PP),length(tgpse));
% SS2=nan(length(PP),length(tgpse));
% TT2=nan(length(PP),length(tgpse));
% for izv=1:length(optsensornames) 
%    if optsensor.(optsensornames{izv}) == 1
%        VAR = upper(optsensornames{izv});
%        eval([VAR '2 = nan(length(PP),length(tgpse))']);
%    end
% end    
% for icc=0:length(tgpse)-1
% 	ijk=2*icc;
% 	TT2(:,icc+1)=TT(:,ijk+1);
% 	TP2(:,icc+1)=TP(:,ijk+1);
% 	SS2(:,icc+1)=SS(:,ijk+1);	
%     gldivephase2(icc+1) = gldivephase(ijk+1);
%     for izv=1:length(optsensornames) % qc for additional var
%         if optsensor.(optsensornames{izv}) == 1
%           VAR = upper(optsensornames{izv});
%           eval([VAR '2(:,icc+1)=' VAR '(:,ijk+1);']);
%         end
%     end     
% end
% timedive = nanmean([tgps2 nanmean([tgps2 tgpse],2)],2); 
% latdive = nanmean([lat2 nanmean([lat2 late],2)],2);
% londive = nanmean([lon2 nanmean([lon2 lone],2)],2);
% timediveini = timepfini(1:2:end)';
% latdiveini = latpfini(1:2:end)';
% londiveini = lonpfini(1:2:end)';
% timediveend = timepfend(1:2:end)';
% latdiveend = latpfend(1:2:end)';
% londiveend = lonpfend(1:2:end)';
% commentdivemerged= 'QC, vertically interpolated, dives only but DAC still calculated for a dive-climb cycle';
% %-------------------------------------------------------------------------------
% % If only climbs 
% elseif mergedpflevel == 3
% TIMED = repmat(timeDAC',length(PP),1);
% PRES = repmat(PP',1,length(tgpse));
% TP2=nan(length(PP),length(tgpse));
% SS2=nan(length(PP),length(tgpse));
% TT2=nan(length(PP),length(tgpse));
% for izv=1:length(optsensornames) 
%    if optsensor.(optsensornames{izv}) == 1
%        VAR = upper(optsensornames{izv});
%        eval([VAR '2 = nan(length(PP),length(tgpse))']);
%    end
% end    
% for icc=0:length(tgpse)-1
% 	ijk=2*icc+1;
% 	TT2(:,icc+1)=TT(:,ijk+1);
% 	TP2(:,icc+1)=TP(:,ijk+1);
% 	SS2(:,icc+1)=SS(:,ijk+1);
%     gldivephase2(icc+1) = gldivephase(ijk+1);
%     for izv=1:length(optsensornames) % qc for additional var
%         if optsensor.(optsensornames{izv}) == 1
%           VAR = upper(optsensornames{izv});
%           eval([VAR '2(:,icc+1)=' VAR '(:,ijk+1);']);
%         end
%     end      
% end
% timedive = nanmean([nanmean([tgps2 tgpse],2) tgpse],2); 
% latdive = nanmean([nanmean([lat2 late],2) late],2);
% londive = nanmean([nanmean([lon2 lone],2) lat2],2);
% timediveini = timepfini(2:2:end)';
% latdiveini = latpfini(2:2:end)';
% londiveini = lonpfini(2:2:end)';
% timediveend = timepfend(2:2:end)';
% latdiveend = latpfend(2:2:end)';
% londiveend = lonpfend(2:2:end)';
% commentdivemerged= 'QC, vertically interpolated, climbs only but DAC still calculated for a dive-climb cycle';
end %fi0

% if mergedpflevel ==0 % if climb and dive kept, altimeter data vector is doubled
%     alti_data(1:2:length(timedive)) = altimeter_data;
%     alti_data(2:2:length(timedive)) = altimeter_data;  
%     bathy_data = zbathy;   
%     HH(1:2:length(timedive)) = abs(mergedbathy);
%     HH(2:2:length(timedive)) = abs(mergedbathy);  
%     numberofdivecycle(1:2:length(timedive)) = divenum;
%     numberofdivecycle(2:2:length(timedive)) = divenum; 
% 
% else
%     bathy_data = zbathydive;
%     alti_data = mergedbathy;
%     HH = abs(mergedbathy);
%     numberofdivecycle = divenum;   
% end
alti_data = alti_data(iokpfsel);
bathy_data = bathy_data(iokpfsel);
HH = HH(iokpfsel);
numberofdivecycle = numberofdivecycle(iokpfsel);

%-------------------------------------------------------------------------
% In case of nan in timedive (ex: glider only recording during dive or climb )  
TIMED(:,isnan(timedive)) = [];
PRES(:,isnan(timedive)) = [];
TP2(:,isnan(timedive)) = [];
SS2(:,isnan(timedive)) = [];
TT2(:,isnan(timedive)) = [];
latdive(isnan(timedive)) = [];
londive(isnan(timedive)) = [];
timediveini(isnan(timedive)) = [];
latdiveini(isnan(timedive)) = [];
londiveini(isnan(timedive)) = [];
timediveend(isnan(timedive)) = [];
latdiveend(isnan(timedive)) = [];
londiveend(isnan(timedive)) = [];
HH(isnan(timedive)) = [];
numberofdivecycle(isnan(timedive)) = [];
gldivephase2(isnan(timedive)) = [];
for izv=1:length(optsensornames) % qc for additional var
    if optsensor.(optsensornames{izv}) == 1
        VAR = upper(optsensornames{izv});
        eval([VAR '2(:,isnan(timedive)) = [];']);
    end
end  
timedive(isnan(timedive)) = [];
%-------------------------------------------------------------------------

PDEN2 = sw_pden(SS2,TT2,PP',0)-1000;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Buoyancy Frequency and Rossby Radius calculated on a lighty smooth
% vertical profile (6m running mean).
TT2f=nan*TT2;
SS2f=nan*TT2;
PDEN2f=nan*TT2;
TP2f=nan*TT2;
vfen=6;
for ijk=1:length(timedive)
TT2f(:,ijk)=moy_gl(TT2(:,ijk),vfen);
TP2f(:,ijk)=moy_gl(TP2(:,ijk),vfen);
SS2f(:,ijk)=moy_gl(SS2(:,ijk),vfen);
PDEN2f(:,ijk)=moy_gl(PDEN2(:,ijk),vfen);
end
% %figure;plot(TT2,-PRES,'b');hold on;plot(TT2f,-PRES,'r')

[bfrq,vort,p_ave] = sw_bfrq(SS2f(1:vfen:end,:),TT2f(1:vfen:end,:),PRES(1:vfen:end,:),latdive');
[bfrqt,vortt,p_avet] = sw_bfrq(SS2f(1:vfen:end,:)*0 + 35,TT2f(1:vfen:end,:),PRES(1:vfen:end,:),latdive');
TIMED2=TIMED(1:vfen:end,:); TIMED2 = TIMED2(2:end,:);

rosradlvl = sqrt(bfrq).* repmat(HH,length(p_ave(:,1)),1)/pi./repmat(sw_f(latdive'),length(p_ave(:,1)),1)/1000;
rosradlvlt = sqrt(bfrqt).* repmat(HH,length(p_ave(:,1)),1)/pi./repmat(sw_f(latdive'),length(p_ave(:,1)),1)/1000;
rosrad =nanmean(rosradlvl);
rosradt = nanmean(rosradlvlt);

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Detited DAC 
if step2opt.untidedDAC == 2
    addpath(genpath(step2opt.pathtmdtoolboxwithmodel))
    timesel0= tgps2;
    timesel1= timeDAC;
    timesel2= tgpse; 
    latsel0= lat2; 
    latsel1= lonDAC;
    latsel2= late; 
    lonsel0= lon2;
    lonsel1= lonDAC; 
    lonsel2= lone; 
    commentdetidedDAC = 'DAC detided using tmd_tide toolbox with the AO model';
    [u_tide_vel0,uconlist] = tmd_tide_pred('Model_AO',timesel0,latsel0,lonsel0,'u');
    [v_tide_vel0,vconlist] = tmd_tide_pred('Model_AO',timesel0,latsel0,lonsel0,'v');
    [u_tide_vel1,uconlist] = tmd_tide_pred('Model_AO',timesel1,latsel1,lonsel1,'u');
    [v_tide_vel1,vconlist] = tmd_tide_pred('Model_AO',timesel1,latsel1,lonsel1,'v');
    [u_tide_vel2,uconlist] = tmd_tide_pred('Model_AO',timesel2,latsel2,lonsel2,'u');
    [v_tide_vel2,vconlist] = tmd_tide_pred('Model_AO',timesel2,latsel2,lonsel2,'v');  
    u_tidedive=nanmean([u_tide_vel0; u_tide_vel1; u_tide_vel1; u_tide_vel2])/100;
    v_tidedive=nanmean([v_tide_vel0; v_tide_vel1; v_tide_vel1; v_tide_vel2])/100;
    DAC_east_untided = DAC_east - u_tidedive';
    DAC_north_untided = DAC_north - v_tidedive';
    
elseif step2opt.untidedDAC == 1 %remove tide and internal waves
    Fc=24*3600/(3600*48); %48h in [days-1]
    Nn=19;
    Nn2=Nn*4;
    % for DAC
    timeDACint = timeDAC(1):1/24:timeDAC(end); % one hour time vector
    iok = (~isnan(DAC_north) & sqrt(DAC_north.^2 + DAC_east.^2)<0.8);
    DAC_northinterp=interp1(timeDAC(iok),DAC_north(iok),timeDACint);
    inonan = ~isnan(DAC_northinterp);
    DAC_north_untidedlp00 = auto_filt(DAC_northinterp(inonan),1/diff(timeDACint(1:2)),Fc);
    DAC_north_untidedlp = interp1(timeDACint(inonan),DAC_north_untidedlp00,timeDAC)';
    DAC_north_untidedlp(~iok)=nan;
 
    iok = (~isnan(DAC_east) & sqrt(DAC_north.^2 + DAC_east.^2)<0.8);    
    DAC_eastinterp=interp1(timeDAC(iok),DAC_east(iok),timeDACint);
    inonan = ~isnan(DAC_eastinterp);
    DAC_east_untidedlp00 = auto_filt(DAC_eastinterp(inonan),1/diff(timeDACint(1:2)),Fc);
    DAC_east_untidedlp = interp1(timeDACint(inonan),DAC_east_untidedlp00,timeDAC)';
    DAC_east_untidedlp(~iok)=nan;
    
    DAC_north_untided = DAC_north_untidedlp;
    DAC_east_untided = DAC_east_untidedlp;
    
    % for T, S, PDEN fields:
    TT2_lp=nan*TT2;
    SS2_lp=nan*TT2;
    PDEN2_lp=nan*TT2;
    TP2_lp=nan*TT2;
    timediveint = timedive(1):1/24:timedive(end); % one hour time vector  
    for izz = 1:length(PRES(:,1))
    	iokt = ~isnan(TT2(izz,:));   
    	ioks = ~isnan(SS2(izz,:)); 
    	if length(timedive(iokt))<18;continue;end
    	tt2int=interp1(timedive(iokt),TT2(izz,iokt),timediveint);
        inonan = ~isnan(tt2int);
        tt2lp00 = auto_filt(tt2int(inonan),1/diff(timediveint(1:2)),Fc);
    	tt2lp = interp1(timediveint(inonan),tt2lp00,timedive)';
    	tt2lp(~iokt)=nan;
    	TT2_lp(izz,:)=tt2lp;
    	
    	if length(timedive(ioks))<18;continue;end    	
     	tp2int=interp1(timedive(ioks),TP2(izz,ioks),timediveint);
        inonan = ~isnan(tp2int);
        tp2lp00 = auto_filt(tp2int(inonan),1/diff(timediveint(1:2)),Fc);
    	tp2lp = interp1(timediveint(inonan),tp2lp00,timedive)';
    	tp2lp(~ioks)=nan;
    	TP2_lp(izz,:)=tp2lp;   	
    	
     	ss2int=interp1(timedive(ioks),SS2(izz,ioks),timediveint);
        inonan = ~isnan(ss2int);
        ss2lp00 = auto_filt(ss2int(inonan),1/diff(timediveint(1:2)),Fc);
    	ss2lp = interp1(timediveint(inonan),ss2lp00,timedive)';
    	ss2lp(~ioks)=nan;
    	SS2_lp(izz,:)=ss2lp;   
    	
      	pd2int=interp1(timedive(ioks),PDEN2(izz,ioks),timediveint);        
        inonan = ~isnan(pd2int);
        pd2lp00 = auto_filt(pd2int(inonan),1/diff(timediveint(1:2)),Fc);
    	pd2lp = interp1(timediveint(inonan),pd2lp00,timedive)';
    	pd2lp(~ioks)=nan;
    	PDEN2_lp(izz,:)=pd2lp;   
    	   	  	
    end
    
    commentdetidedDAC = 'DAC, T, S, PDEN detided using a lowpass hamming filter with a cut of frequency of 48h';
else
    commentdetidedDAC = 'DAC not detided (DAC_*_untided set to nan)';
    DAC_east_untided = DAC_east*nan;
    DAC_north_untided = DAC_north*nan';  
end

% figure;plot(timeDAC,DAC_north)
% hold on
% plot(timeDAC,DAC_north_untided)
% % plot(timesel1,v_tidedive)
% %  plot(timeDAC,DAC_north_untidedlp)

%-------------------------------------------------------------------------
% Formating data into a structure
clear glider
glider.name = missioninfo.glname;
glider.mission = missioninfo.glmission;
glider.date_deploy = datestr(timedive(1),'dd/mm/yy');
glider.date_lastdata = datestr(timedive(end),'dd/mm/yy');
glider.process_step = missioninfo.step2.proclevel ;
glider.process_step_comment{1} = commentdivemerged;
glider.process_step_comment{2} = commentdetidedDAC;
glider.process_step_comment{3} = ['Buoyancy frequency and Rossby Radius calculated on vertical profiles smoothed using a ' num2str(vfen) 'm running mean'];
glider.process_step_comment{4} = ['Bathymetry database used to interpolate on glider track: ' step2opt.bathyname];
glider.merged_level = mergedpflevel;
for ijk=1:length(timedive)
	glider.data(ijk).divecyclenumber = numberofdivecycle(ijk);
	glider.data(ijk).divephase = gldivephase2(ijk);    
	glider.data(ijk).time = timedive(ijk);
	glider.data(ijk).lat = latdive(ijk);	
	glider.data(ijk).lon = londive(ijk);	
	glider.data(ijk).depth_alti = alti_data(ijk);		
	glider.data(ijk).depth_bathy = bathy_data(ijk);		
	glider.data(ijk).depth_mergedaltibathy = -HH(ijk);			
	glider.data(ijk).ptmp= TP2(:,ijk);
	glider.data(ijk).temp= TT2(:,ijk);	
	glider.data(ijk).sal = SS2(:,ijk);	
	glider.data(ijk).pden = PDEN2(:,ijk);		
	glider.data(ijk).pres = PRES(:,ijk);		
	glider.data(ijk).buoyfreq = bfrq(:,ijk);
	glider.data(ijk).buoyfreq_temp = bfrqt(:,ijk);	
	glider.data(ijk).pres_buoyfreq = p_ave(:,ijk);
	glider.data(ijk).rosrad = rosradlvl(:,ijk);	
	if exist('TT2_lp','var') == 1
	glider.data(ijk).ptmplp = TP2_lp(:,ijk);
	glider.data(ijk).templp = TT2_lp(:,ijk);	
	glider.data(ijk).sallp  = SS2_lp(:,ijk);	
	glider.data(ijk).pdenlp = PDEN2_lp(:,ijk);		
	else
	glider.data(ijk).ptmplp = [];
	glider.data(ijk).templp = [];	
	glider.data(ijk).sallp  = [];	
	glider.data(ijk).pdenlp = [];		
	end
    for izv=1:length(optsensornames) 
        if optsensor.(optsensornames{izv}) == 1
            VAR = upper(optsensornames{izv});
            eval(['glider.data(ijk).' optsensornames{izv} ' = ' VAR '2(:,ijk);']);
        end
    end    
	glider.data(ijk).timestart = timediveini(ijk);
	glider.data(ijk).timeend = timediveend(ijk);
	glider.data(ijk).lonstart = londiveini(ijk);	
	glider.data(ijk).lonend = londiveend(ijk);	
	glider.data(ijk).latstart = latdiveini(ijk);	
	glider.data(ijk).latend = latdiveend(ijk);        
end
for ijk=1:length(timeDAC)
	glider.DACdata(ijk).DAC_divecyclenumber = divenum(ijk);	    
	glider.DACdata(ijk).DAC_lat = latDAC(ijk);	
	glider.DACdata(ijk).DAC_lon = lonDAC(ijk);	
	glider.DACdata(ijk).DAC_time = timeDAC(ijk);	    
	glider.DACdata(ijk).DAC_east = DAC_east(ijk);	
	glider.DACdata(ijk).DAC_north = DAC_north(ijk);					
	glider.DACdata(ijk).DAC_qc = DAC_qc(ijk);	
	glider.DACdata(ijk).DAC_eastlp = DAC_east_untided(ijk);	
	glider.DACdata(ijk).DAC_northlp = DAC_north_untided(ijk);	    
end
for ijk=1:length(sec_time(:,1))
    id1 = near(timediveini,sec_time(ijk,1));
    id2 = near(timediveend,sec_time(ijk,2));
    glider.section(ijk).timedef = sec_time(ijk,:);
    glider.section(ijk).index = [id1 id2];
end
glider.optsensor = optsensor;
glider.date_deploy = datestr(timedive(1),'dd/mm/yy');
glider.date_recup = datestr(timedive(end),'dd/mm/yy');

eval(['save ' missioninfo.step2.procdataDir filesep missioninfo.step2.matfilename '.mat glider' ])  

end
