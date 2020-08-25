function step1_process_glider_data(missioninfo,qc_param,step1opt)
%
%=========================================================================%
% Function to process seaglider: automatic quality control tests, manual
% quality control for each specific mission and creation of matrice with vertical profiles
% of pot. temp, salinity, pot. density, vertical velocity (from the flight
% model)
%
% L. Houpert, SAMS 12/04/2016
%=========================================================================%
%
% step1_process_glider_data(missioninfo,qc_param,step1opt)
% Function that apply QC tests, flag bad data, and create matrix of the vertical profiles. Results are save in
% [missioninfo.step1.procdataDir filesep missioninfo.step1.matfilename]
%
%   Inputs: o missioninfo: [1x1] structure with mission details (defined in
%             users_param/loadobservatories.m)
%                where - missioninfo.step1.procdataDir : path of the
%                directory to save step1 processing outputs
%                      - missioninfo.step1.matfilename
%                      - missioninfo.section : [kx1] structure indexing the time
%                      of the beginning and end of the k glider sections as a 
%                      [1x2] matrix defined by the users in missioninfo.section(i).time
%                      - missioninfo.rawdata.dataDir : path of the eng and log files used for gps correction;
%
%           o qc_param: [1x1] structure containing the different automatic QC
%             parameters (variable name, depth range, tolerance level) per test and per variable
%             (defined in users_param/loadobservatories.m) but also the manual QC if it exists
%           
%           o step1opt: [1x1] structure containing:
%                     - .bathyname: Bathymetry used to get an interpolation
%                     along the glider track
%                     - .bathydir: Directory where the bathymetry file is located
%                     - .qflaglvllim : Maximum value of flag kept in the final data (e.g 2)
%                     users_param/manual_qc, if = 0 don't run the script
%                     but load the saved data files
%
%-------------------------------------------------
%
%	QC flags for temperature and salinity:  
%   1=good; 
%   2=probably good; 
%   3=probably bad; 
%   4=bad; 
%   8=Seabird basestation flag for interpolated value; 
%   9=no data; 
%   40=manually flagged bad; 
%   41=spike detection; 
%   42=gradient test failed; 
%   43=density inversion; 
%   44=inferior to fix threshold;  
%   45=outliers from pressure bin mean statistics 
%
%--------------------------------------------------
%
%   If a new deployment is process and need additional QC control (manual selection of bad data), edit the script step1_process_glider_data_sub_additional_qc
%
% created by L. Houpert (houpertloic@gmail.com), 12/04/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 
% If one wants to improve the QC process, particularly the manual one (creation of specific functions 
% that let the user do an interactive selection of the bad data (on a 10 or 20 profiles window for example) 
% and save the data flagged bad in a .mat or a text file to keep track of the bad data selection 
% (by adding an option to read these file, the function here could be improve to skip the manual QC selection and 
% just load these files and flagged the bad data). Then a good way to make this toolbox 
% evolve is to create a new branch, commit the changes and push the changes on the git remote repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)
%
% 

if  exist(step1opt.bathydir,'dir')==7
    usebathyfile =1;
end

qflaglvl = step1opt.qflaglvllim;

%-----------------------------------------
% Create outputdirs if didn't exist
outputdir = missioninfo.step1.procdataDir;  
if exist(outputdir,'dir')~=7   mkdir(outputdir); mkdir([outputdir '/bad_dives_detected']);end

%----------------------------------------
% load section info:
sec_time=nan(length(missioninfo.section),2);
for ill=1:length(missioninfo.section)
    sec_time(ill,:) = missioninfo.section(ill).time;
end
%------------------------------------------
load([missioninfo.step0.procdataDir '/' missioninfo.step0.matfilename])


if isfield(glider.dive(1),'divephase')
    phaseexistinrawdata = 1 ; 
else
    phaseexistinrawdata = 0 ;    
end
%------------------------------------------
% list additional sensors:
if isfield(glider.dive(1),'doxy')
    optsensor.doxy = 1;
else
    optsensor.doxy = 0;   
end
if isfield(glider.dive(1),'chla')
    optsensor.chla = 1;
else
    optsensor.chla = 0;   
end
if isfield(glider.dive(1),'bb')
    optsensor.bb = 1;
else
    optsensor.bb = 0;   
end
if isfield(glider.dive(1),'cdom')
    optsensor.cdom = 1;
else
    optsensor.cdom = 0;   
end
optsensornames = fieldnames(optsensor);
%------------------------------------------
% Load and concatenate glider mission data
iok=find(~isnan(vertcat(glider.dive(:).timegps2)));
tgps2 = vertcat(glider.dive(iok).timegps2);
divenum = vertcat(glider.dive(iok).divenum);
longitude = vertcat(glider.dive(iok).longitude);
latitude = vertcat(glider.dive(iok).latitude);
time = vertcat(glider.dive(iok).ctdtime);
pres = sw_pres(vertcat(glider.dive(iok).ctddepth),latitude);
temp = vertcat(glider.dive(iok).temperature);
salin = vertcat(glider.dive(iok).salinity);
qc_temp = (vertcat(glider.dive(iok).temperature_qc));
qc_salin = (vertcat(glider.dive(iok).salinity_qc));
qc_ctd = (vertcat(glider.dive(iok).ctd_qc));

for ijk=1:length(optsensornames)
    if optsensor.(optsensornames{ijk}) == 1
        var = eval(['vertcat(glider.dive(iok).' optsensornames{ijk} ');']);
        eval([optsensornames{ijk} ' =  var ;']);
        eval(['qc_' optsensornames{ijk} ' = vertcat(glider.dive(iok).' optsensornames{ijk} '_qc);']);
    end
end

if phaseexistinrawdata == 1
    phasefromraw = (vertcat(glider.dive(iok).divephase));
    phasefromraw(isnan(phasefromraw)) = 9; % replace nan (if there is some) by 9
end

tgps2 = vertcat(glider.dive(iok).timegps2);
lat2 = vertcat(glider.dive(iok).latgps2);
lon2 = vertcat(glider.dive(iok).longps2);
qc_gps2 = (vertcat(glider.dive(iok).gps2_qc));

tgpse = vertcat(glider.dive(iok).timegps);
late = vertcat(glider.dive(iok).latgps);
lone = vertcat(glider.dive(iok).longps);
qc_gpse = (vertcat(glider.dive(iok).gpse_qc));

hdm_vert_speed = vertcat(glider.dive(iok).vert_speed);
hdm_glide_angle = vertcat(glider.dive(iok).glide_angle);

w_water=  vertcat(glider.dive(iok).wwater);

www =hdm_vert_speed;

DAC_east = vertcat(glider.dive(iok).DAC_east);
DAC_north = vertcat(glider.dive(iok).DAC_north);
DAC_qc = (vertcat(glider.dive(iok).DAC_qc));

tempdata = horzcat(glider.dive(iok).altim_bottom_ping); 
altimeter_data0 = - (tempdata(1:2:end) + tempdata(2:2:end));

glname = glider.name;
glmission = glider.mission;

%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% GPS correction when several dive without surfacing are made and end dive
% position are similar to predive positions
[tgps2,lon2,lat2,tgpse,lone,late,DAC_qc]=correct_nogpsdata_ndives(tgps2,lon2,lat2,tgpse,lone,late,DAC_qc);


%-------------------------------------------------------------------------
% Determination of starting and ending time of the vertical profiles 
t1 = vertcat(glider.dive(iok).starttime);
t2 = vertcat(glider.dive(iok).start_climb_time); 
t1b = t2;
t2b = tgpse;

%figure; plot(t2-t1,'+-');
%hold on; plot(t2b-t1b,'+-r');
%hold on; plot(t2b-t1,'o-g');
ipf1a = nan(length(t1),1);
ipf2a = nan(length(t1),1);
ipf1b = nan(length(t1),1);
ipf2b = nan(length(t1),1);
ipf1=[];
ipf2=[];
for ijk = 1:length(t1)
   ipf1a(ijk) = near(time,t1(ijk))+1;
   ipf2a(ijk) = near(time,t2(ijk))-1;   
   ipf1b(ijk) = near(time,t1b(ijk))+1;
   ipf2b(ijk) = near(time,t2b(ijk))-1;     
   ipf1 = [ipf1; ipf1a(ijk); ipf1b(ijk)];
   ipf2 = [ipf2; ipf2a(ijk); ipf2b(ijk)];   
end
%figure;plot(ipf2-ipf1)
if exist(missioninfo.step1.procdataDir,'dir')~=7
    mkdir(missioninfo.step1.procdataDir)
end   
eval(['save '  missioninfo.step1.procdataDir filesep 'index_dive_climb_'  glname '_' glmission '.mat time ipf1 ipf2' ])  

% -----------------------------------------------------------------------
% Flag as 9 the index without data (essentially seen in EGO netcdf)
qc_temp(isnan(temp))   = 9;
qc_salin(isnan(salin)) = 9; 
for ijk=1:length(optsensornames)
    if optsensor.(optsensornames{ijk}) == 1
        eval(['qc_' optsensornames{ijk} '(isnan(' optsensornames{ijk} ')) = 9;']);
    end
end
%-----------------------------------------------------------------------
% keep raw data
temp00=temp;
salin00=salin;
w_water00=w_water;
temp0=temp;
salin0=salin;
for ijk=1:length(optsensornames)
    if optsensor.(optsensornames{ijk}) == 1
        eval([optsensornames{ijk} '00 = ' optsensornames{ijk} ';']);
    end
end

% --------------------------------------------------
% Don't keep data flagged already flagged as bad 
temp(qc_temp>qflaglvl)=nan;
salin(qc_salin>qflaglvl)=nan; 
w_water(qc_salin>qflaglvl | qc_temp>qflaglvl)=nan;
for ijk=1:length(optsensornames)
    if optsensor.(optsensornames{ijk}) == 1
        eval([optsensornames{ijk} '(qc_' optsensornames{ijk} '>qflaglvl)=nan;']);
    end
end


%==================================================
% Manual QC for the mission manual QC file (timeflag)
%--------------------------------------------------
 if isfield(qc_param.manualsel,'temp') & ~isempty(qc_param.manualsel.temp.tindex)
    baddataindex=[];
    for iseltime=1:length(qc_param.manualsel.temp.tindex.low)
        tsel1=qc_param.manualsel.temp.tindex.low(iseltime);
        tsel2=qc_param.manualsel.temp.tindex.upp(iseltime);
        iseltt=find(time<tsel2 & time>tsel1 & qc_salin~=9);
        baddataindex = [baddataindex ; iseltt];
        baddataindex = unique(baddataindex);
    end
    qc_temp(baddataindex)=40;
    temp(qc_temp>qflaglvl)=nan;
end
if isfield(qc_param.manualsel,'salin') & ~isempty(qc_param.manualsel.salin.tindex)
    baddataindex=[];
    for iseltime=1:length(qc_param.manualsel.salin.tindex.low)
        tsel1=qc_param.manualsel.salin.tindex.low(iseltime);
        tsel2=qc_param.manualsel.salin.tindex.upp(iseltime);
        iseltt=find(time<tsel2 & time>tsel1 & qc_salin~=9);
        baddataindex = [baddataindex ; iseltt];
        baddataindex = unique(baddataindex);
    end
    qc_salin(baddataindex)=40;
    salin(qc_salin>qflaglvl)=nan;
end
% Manual qc for additional vars
for ijk=1:length(optsensornames)  
    if optsensor.(optsensornames{ijk}) == 1
        eval(['varqc = qc_' optsensornames{ijk} ';']);  
        if isfield(qc_param.manualsel,optsensornames{ijk}) & ~isempty(qc_param.manualsel.(optsensornames{ijk}).tindex)
            baddataindex=[];
            for iseltime=1:length(qc_param.manualsel.(optsensornames{ijk}).tindex.low)
                tsel1=qc_param.manualsel.(optsensornames{ijk}).tindex.low(iseltime);
                tsel2=qc_param.manualsel.(optsensornames{ijk}).tindex.upp(iseltime);
                iseltt=find(time<tsel2 & time>tsel1 & varqc~=9);
                baddataindex = [baddataindex ; iseltt];
                baddataindex = unique(baddataindex);
            end
            eval(['qc_' optsensornames{ijk} '(baddataindex) = 40;']);                  
            eval([optsensornames{ijk} '(qc_' optsensornames{ijk} '>qflaglvl)=nan;']);                     
        end
    end
end


%==================================================
% Automatic QC process
%==================================================
%-------------------------------------------------
% untruthful values
%-------------------------------------------------
if isfield(qc_param.unthrust,'temp')
    qc_temp(temp00<qc_param.unthrust.temp.low | temp00>qc_param.unthrust.temp.high) = 44;
end
if isfield(qc_param.unthrust,'salin')
    qc_salin(salin00<qc_param.unthrust.salin.low | salin00>qc_param.unthrust.salin.high) = 44;
end
for ijk=1:length(optsensornames)
    if optsensor.(optsensornames{ijk}) == 1
        if isfield(qc_param.unthrust,optsensornames{ijk})
            varn=optsensornames{ijk};
            eval(['qc_' varn '(' varn '00<qc_param.unthrust.' varn '.low | ' ...
                var '00>qc_param.unthrust.' varn '.high) = 44 ;']);
        end
    end
end


%-------------------------------------------------
% Spike detection (need at least 3 points in a 50m layer)
%-------------------------------------------------
for iqq=1:length(qc_param.spike.temp)
[ibad] = qc_spike_detection(temp,pres,qc_param.spike.temp(iqq));
qc_temp(ibad) = 41;
end
for iqq=1:length(qc_param.spike.sal)
[ibad] = qc_spike_detection(salin,pres,qc_param.spike.sal(iqq));
qc_salin(ibad) = 41;
end
temp(qc_temp>qflaglvl)=nan;
salin(qc_salin>qflaglvl)=nan;


% No gradient QC as in the Seaglider processing
% %-------------------------------------------------
% % Gradient test (need at least 3 points in a 50m layer)
% %--------------------------------------------------
% for iqq=1:length(qc_param.gradtest.temp)
% [ibad] = qc_gradient_test(temp,pres,qc_param.gradtest.temp(iqq));
% qc_temp(ibad) = 42;
% end
% for iqq=1:length(qc_param.gradtest.sal)
% [ibad] = qc_gradient_test(salin,pres,qc_param.gradtest.sal(iqq));
% qc_salin(ibad) = 42;
% end
% temp(qc_temp>qflaglvl)=nan;
% salin(qc_salin>qflaglvl)=nan;


if 0% check that the density inversion is not removing to many good data points 
%--------------------------------------------------
% Density inversion detection
%--------------------------------------------------
pden0 = sw_pden(salin,temp,pres,0)-1000;
for iqq=1:length(qc_param.invtest.pden)
[ibad] = qc_inversion_test(pden0,pres,qc_param.invtest.pden(iqq));
qc_salin(ibad) = 43;
end
end


%----------------------------------------------------------
% Remove outliers in pressure bin average
pden0 = sw_pden(salin,temp,pres,0)-1000;
for iqq=1:length(qc_param.outlierspbin.pden)
[ibad,~,~,~] = qc_presbin_outliers(pden0,pres,qc_param.outlierspbin.pden(iqq));
qc_salin(ibad) = 45;
end


%--------------------------------------------------

for ijk=1:length(optsensornames) % qc for additional var
    if optsensor.(optsensornames{ijk}) == 1
        eval([optsensornames{ijk} '(qc_' optsensornames{ijk} '>qflaglvl)=nan;']);
    end
end
temp(qc_temp>qflaglvl)=nan;
salin(qc_salin>qflaglvl)=nan;
w_water(qc_salin>qflaglvl | qc_temp>qflaglvl)=nan;

%==================================================
% Vertical interpolatation of the profiles 
%--------------------------------------------------
ptmp00=sw_ptmp(salin00,temp00,pres,0);
ptmp=sw_ptmp(salin,temp,pres,0);
pden = sw_pden(salin,temp,pres,0)-1000;
pden00 = sw_pden(salin00,temp00,pres,0)-1000;

PP = 0:1:1010;

%--------------------------------------------------
calcdivephase = ones(size(time))*9;
%--------------------------------------------------
immax = length(ipf1);

TT = nan(length(PP),immax);
TP = nan(length(PP),immax);
SS = nan(length(PP),immax);
WWW = nan(length(PP),immax);
gldivephase = nan(1,immax);

for ijk=1:length(optsensornames) % qc for additional var
    if optsensor.(optsensornames{ijk}) == 1
        VAR = upper(optsensornames{ijk});
        eval([VAR ' = nan(length(PP),immax);']);
    end
end
timepf    = nan(1,immax);
lonpf     = nan(1,immax);
latpf     = nan(1,immax);
timepfini = nan(1,immax);
lonpfini  = nan(1,immax);
latpfini  = nan(1,immax);
timepfend = nan(1,immax);
lonpfend  = nan(1,immax);
latpfend  = nan(1,immax);
%---------------------------------------------------------------------
imm=0;
noprofile =0;
for ijk = 1:length(ipf1)
    indint = ipf1(ijk):ipf2(ijk);
    indint01 = indint(~isnan(pres(indint)));   
    imm = imm +1;
    if length(indint)>4 & length(indint01)>4 & max(pres(indint))>20  & max(diff(pres(indint)))<100 % at least 4 points, deeper than 20db and no jump of more than 100db
        %ipok = 1 + find(abs(diff(pres(indint01)) - mean(diff(pres(indint01))))< abs(mean(diff(pres(indint01))))/2);
        %ipok = 1 + find(abs(diff(pres(indint01)))< 50);
        ipok0 = find(pres(indint01) < (max(pres(indint01))-2) & pres(indint01) > (min(pres(indint01))+1)); % apogee & top 
        [pp0 ipokuni wdwd]=unique(pres(indint01(ipok0)),'stable') ;
        ipok = ipok0(ipokuni);
        indint1 = indint01(ipok); % indint1 = [indint01(1) indint01(ipok)] ;    
        if mean(diff(pres(indint01)))>0   
            if mod(imm,2)==0   
                disp(' ');
                disp(['Problem in the diving phase of the glider ' glname ' mission: ' glmission ' on the  ' ...
                    datestr(time(indint01(1)),'dd/mm/yy HH:MM') ' - ' datestr(time(indint01(end)),'dd/mm/yy HH:MM') ]);
                disp('Glider should be in an ascent phase but according to the mean pressure difference he is in a descent phase.');    
                disp('Data are skipped for vertical interpolation.');                  
                disp(' ');  
                calcdivephase(indint01) = 9; % not known                   
            else
            calcdivephase(indint01) = 1; % descent phase                    
            end
        elseif mean(diff(pres(indint01)))<0

            if mod(imm,2)==1 
                disp(' ');
                disp(['Problem in the diving phase of the glider ' glname ' mission: ' glmission ' on the  ' ...
                    datestr(time(indint01(1)),'dd/mm/yy HH:MM') ' - ' datestr(time(indint01(end)),'dd/mm/yy HH:MM') ]);
                disp('Glider should be in a descent phase but according to the mean pressure difference he is in an ascent phase.');    
                disp('Data are skipped for vertical interpolation.');                  
                disp(' ');      
                calcdivephase(indint01) = 9; % not known                 
            else
                calcdivephase(indint01) = 4; % ascent phase                 
            end
        end     
        timepf(imm) = nanmean(time(indint));
        lonpf(imm)  = nanmean(longitude(indint));
        latpf(imm)  = nanmean(latitude(indint));
        timepfini(imm) = time(indint(1));
        timepfend(imm) = time(indint(end));          
        lonpfini(imm) = longitude(indint(1));
        lonpfend(imm) = longitude(indint(end));  
        latpfini(imm) = latitude(indint(1));
        latpfend(imm) = latitude(indint(end));     
        indinttemp = indint1(~isnan(temp(indint1)));
        indint2 = indint1(~isnan(salin(indint1)));
        lpf = length(pres(indint01));
        if length(indinttemp)>4 & length(indint2)<4
        	TT(:,imm) = interp1(pres(indinttemp),temp(indinttemp),PP)';    
         	TP(:,imm) = nan*PP';
        	SS(:,imm) = nan*PP';       	
         	WWW(:,imm) = nan*PP';        
        elseif length(indinttemp)>4 & length(indint2)>4
        	TT(:,imm) = interp1(pres(indinttemp),temp(indinttemp),PP)';    
        	TP(:,imm) = interp1(pres(indint2),ptmp(indint2),PP)';
        	SS(:,imm) = interp1(pres(indint2),salin(indint2),PP)'; 
        	WWW(:,imm) = interp1(pres(indint2),w_water(indint2),PP)'; 
            for izv=1:length(optsensornames) % qc for additional var
                if optsensor.(optsensornames{izv}) == 1
                    VAR = upper(optsensornames{izv});
                    eval(['indintvar = indint1(~isnan(' optsensornames{izv}  '(indint1)));']);
                    if length(indintvar)<4
                        continue
                    else
                        eval([VAR '(:,imm) = interp1(pres(indintvar),' optsensornames{izv}  '(indintvar),PP);']);
                    end
                end
            end
        else
        	TT(:,imm) = nan*PP';
        	TP(:,imm) = nan*PP';
        	SS(:,imm) = nan*PP';
        	WWW(:,imm) = nan*PP';
        end	
        gldivephase(imm) = calcdivephase(indint01(1));                
    else
        noprofile = noprofile +1;
    end
end
% % creation of a dive-phase variable for interp profiles, following EGO format: 0=surface drift; 1= descent; 2= subsurface drift; 3= inflexion; 4=ascent
% gldivephase(1:2:length(TT)) = 1;
% gldivephase(2:2:length(TT)) = 4;
% % if the dive phase exist in the raw file, we keep the orginal otherwise we
% % used the one calculated
if phaseexistinrawdata == 1    
    divephase = phasefromraw;
    divephase(divephase==9) = calcdivephase(divephase==9); 
else
    divephase = calcdivephase;
end
%
%-------------------------------------------------------------------------
% creation of a dive id variable
gldiveid = nan(size(time));
for ijk = 1:2:length(ipf1)
    indint = ipf1(ijk):ipf1(ijk+1);
    gldiveid(indint)=divenum((ijk-1)/2 +1);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% save index of the bad salinity and temperature profiles 
% in a specific matrix, so can be reuse easily by other process scripts (like flight model optimisation)

baddivesal = union(divenum(find(isnan(nanmean(SS(:,1:2:end))))),divenum(find(isnan(nanmean(SS(:,2:2:end))))));
baddivetemp = union(divenum(find(isnan(nanmean(TT(:,1:2:end))))),divenum(find(isnan(nanmean(TT(:,2:2:end))))));

eval(['save '  missioninfo.step1.procdataDir filesep 'bad_dives_detected' filesep  glname '_' glmission '_detection_baddivenbers.mat baddivesal baddivetemp' ])  
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Set nan WW when bad temperature and salinity
divenum2 = sort([divenum ; divenum]);
ibadss = [];
for ibb = 1:length(baddivesal)
    ibadss = [ ibadss find(divenum2==baddivesal(ibb))];
end
ibadtt = [];
for ibb =  1:length(baddivetemp)
    ibadtt = [ibadtt find(divenum2==baddivetemp(ibb))];
end
WWW(:,ibadss)=nan;
WWW(:,ibadtt)=nan;
for idid = union( baddivesal, baddivetemp)'
	ibad=find(gldiveid == idid);
	w_water(ibad)=nan;
end

PDEN = sw_pden(SS,TT,PP',0)-1000;

%-------------------------------------------------------------------------
% definition position for DAC:
timeDAC = nanmean([tgps2 tgpse],2);
latDAC = nanmean([lat2 late],2);
lonDAC = nanmean([lon2 lone],2);
%-------------------------------------------------------------------------
% Interpolate bathymetrie data on glider track
if usebathyfile == 1 
    zbathydive = ship_interpbathy(step1opt.bathyname,lonDAC,latDAC,step1opt.bathydir);
    zbathy = ship_interpbathy(step1opt.bathyname,longitude,latitude,step1opt.bathydir);
else
    zbathydive = nan*lonDAC;
    zbathy = nan*longitude;
end

%===================================================================
% process altimeter data
%---------------------------
% 2nd method

spanvalue =3;
londist = cumsum([0 ;sw_dist(latDAC*0+58,lonDAC,'km')]);
altimeter_data = altimeter_data0;

% % remove bad altimeter data if manual qc defined
if isfield(qc_param.manualsel,'altimeter') & ~isempty(qc_param.manualsel.altimeter.tindex)
    baddataindex=[];
    for iseltime=1:length(qc_param.manualsel.altimeter.tindex.low)
        tsel1=qc_param.manualsel.altimeter.tindex.low(iseltime);
        tsel2=qc_param.manualsel.altimeter.tindex.upp(iseltime);
        iseltt=find(timeDAC<tsel2 & timeDAC>tsel1);
        baddataindex = [baddataindex ; iseltt];
        baddataindex = unique(baddataindex);
    end
    altimeter_data(baddataindex)=nan;
end

% merge altimeter and bathymetry data
mergedbathy = altimeter_data';
if usebathyfile == 1 
    mergedbathy(isnan(altimeter_data)) = zbathydive(isnan(altimeter_data));
    mergedbathy = smooth(mergedbathy,spanvalue,'moving');
end
%-------------------------------------------------------------------------
additionalvarstr = [];
for izv=1:length(optsensornames) % qc for additional var
   if optsensor.(optsensornames{izv}) == 1
        VAR = upper(optsensornames{izv});
        additionalvarstr = [additionalvarstr ' ' optsensornames{izv} ' ' optsensornames{izv} '00 ' VAR];
   end
end
eval(['save ' missioninfo.step1.procdataDir filesep missioninfo.step1.matfilename '.mat timepf lonpf latpf timepfini lonpfini latpfini timepfend lonpfend latpfend PP TP TT SS PDEN tgps2 lat2 lon2 tgpse late lone altimeter_data zbathydive zbathy mergedbathy time longitude latitude ptmp pden temp salin pres ptmp00 temp00 pden00 salin00 w_water00 ipf1 ipf2 hdm* divenum w_water gldiveid WWW baddivetemp baddivesal lonDAC latDAC timeDAC DAC_* qc* divephase gldivephase optsensor ' additionalvarstr])  



end
