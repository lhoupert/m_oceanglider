function plot_step2_process_glider(missioninfo,graphparam,step2opt,refDB)
%=========================================================================%
% Function to plot post-processed data for different time window (sections defined by the user), as shifted profiles, time series,
% profiles, surface plot, map, ... 
%
% L. Houpert, SAMS 19/05/2016
%=========================================================================%
%
% plot_step2_process_glider(missioninfo,param_graph,step2opt,refDB)
%
% TODO: % See how merge this function with plot_step1_process_glider_data 
%
% created by L. Houpert (houpertloic@gmail.com), 19/05/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 
% Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)

load([missioninfo.step2.procdataDir filesep missioninfo.step2.matfilename '.mat'])

fnames=fieldnames(glider.data(1));
for ikk=1:length(fnames)
    eval([fnames{ikk} ' = [glider.data.(fnames{ikk})];']);
end
fnames=fieldnames(glider.DACdata(1));
for ikk=1:length(fnames)
    eval([fnames{ikk} ' = [glider.DACdata.(fnames{ikk})];']);
end
clear glider

DAC_lat(DAC_qc~=1 & DAC_qc~=0)=nan;
DAC_lon(DAC_qc~=1 & DAC_qc~=0)=nan;
DAC_east(DAC_qc~=1 & DAC_qc~=0)=nan;
DAC_north(DAC_qc~=1 & DAC_qc~=0)=nan;
DAC_east_untided(DAC_qc~=1 & DAC_qc~=0)=nan;
DAC_north_untided(DAC_qc~=1 & DAC_qc~=0)=nan;



for isecc = 1:length(missioninfo.section)
close all
    
    id1 = near(timestart,missioninfo.section(isecc).time(1));
    id2 = near(timeend,missioninfo.section(isecc).time(2));
    
    ip1 = near(time,missioninfo.section(isecc).time(1));
    ip2 = near(time,missioninfo.section(isecc).time(2));
    
    i1 = near(time,timestart(id1));
    i2 = near(time,timeend(id2));
    
    idac1 = near(DAC_time,missioninfo.section(isecc).time(1));
    idac2 = near(DAC_time,missioninfo.section(isecc).time(2));

    date1 = min(time(i1:i2));
    date2 = max(time(i1:i2));
    
    
    repfigmain = 'time_series_profiles_sections';
    graphparam.mainplot.dir  = [missioninfo.step2.plotdir filesep missioninfo.glmission filesep repfigmain filesep missioninfo.step2.proclevel];
    graphparam.mainplot.name = missioninfo.step2.plotname;
    if exist(graphparam.mainplot.dir ,'dir')~=7 ; mkdir(graphparam.mainplot.dir);end

    %----------------------------------------------------------------------
    % Plot time series 
    for ivar=1:length(graphparam.timeseriesparam) 
        if exist(graphparam.timeseriesparam(ivar).matlabvar_x,'var') & exist(graphparam.timeseriesparam(ivar).matlabvar_y,'var')         
            data =[];
            %-----        
            % step 1
            data.x = eval(graphparam.timeseriesparam(ivar).matlabvar_x);
            data.y = eval(graphparam.timeseriesparam(ivar).matlabvar_y);
            if ~isempty(graphparam.timeseriesparam(ivar).matlabvar_y00)
                data.y00 = eval(graphparam.timeseriesparam(ivar).matlabvar_y00);   
            else
                data.y00=[];
            end
            if ~isempty(graphparam.timeseriesparam(ivar).matlabvar_qc)
                data.qc = eval(graphparam.timeseriesparam(ivar).matlabvar_qc);
            else
                data.qc = data.x*0;
            end
            data.i1  = ip1; % index profiles
            data.i2  = ip2;
            data.date1  = date1;
            data.date2  = date2;        
            data.time  = time;
            data.gname = missioninfo.gnickname;
            data.sectionnb = isecc;

            plot_mission_data_time_series_qc(data,graphparam.timeseriesparam(ivar),graphparam)
                          
        end
    end
    %----------------------------------------------------------------------   
    

    %----------------------------------------------------------------------
    % Plot profiles 
    for ivar= 1:length(graphparam.profilparam) 
        if exist(graphparam.profilparam(ivar).matlabvar_x,'var') & exist(graphparam.profilparam(ivar).matlabvar_y,'var') 
            data =[];
            %-----        
            % step 1
            data.xvarname = graphparam.profilparam(ivar).matlabvar_x;
            data.yvarname = graphparam.profilparam(ivar).matlabvar_y;        
            data.x = eval(graphparam.profilparam(ivar).matlabvar_x);
            data.y = eval(graphparam.profilparam(ivar).matlabvar_y);
            data.ip1  = ip1; % index profiles
            data.ip2  = ip2;
            data0 = data.x(:,ip1:ip2);
            if sum(~isnan(data0(:)))==0
                continue
            end            
            data.date1  = date1;
            data.date2  = date2;        
            data.time  = time;
            data.gname = missioninfo.gnickname;
            data.sectionnb = isecc;
            data.gldivephase = divephase;

            plot_mission_data_profils(data,graphparam.profilparam(ivar),graphparam,refDB)
        end
    end

    %----------------------------------------------------------------------
    % Plot DAC time series and stickvec     
    if (date2-date1)>150
    	graphparam.timeseriesDAC.stickvecscale  = 0.05; %round(length(i1:i2)/5000);% 100 ;  
    else
    	graphparam.timeseriesDAC.stickvecscale  = 0.25;
    end   
    data =[];
    %-------------------        
    % step 2
    data.time =  eval(graphparam.timeseriesDAC.time)';
    data.lon =  eval(graphparam.timeseriesDAC.lon);   
    data.lat =  eval(graphparam.timeseriesDAC.lat); 
    data.UU =  eval(graphparam.timeseriesDAC.UU)';         
    data.VV =  eval(graphparam.timeseriesDAC.VV)'; 
    data.idac1  = idac1;
    data.idac2  = idac2;
    data.date1  = date1;
    data.date2  = date2;        
    data.gname = missioninfo.gnickname;
    data.gmission = missioninfo.glmission;   
    data.sectionnb = isecc;

    plot_mission_DAC_time_series(data,graphparam)
    
    %----------------------------------------------------------------------   
    % plot maps with DAC
    for ivar= 1:length(graphparam.mapDACparam) 
        if exist(graphparam.mapDACparam(ivar).matlabvar_c,'var') 
            data = [];
            data.lon =  eval(graphparam.mapDACparam(ivar).lon);   
            data.lat =  eval(graphparam.mapDACparam(ivar).lat); 
            data.UU =  eval(graphparam.mapDACparam(ivar).UU);       
            data.VV =  eval(graphparam.mapDACparam(ivar).VV);
            cdata = eval(graphparam.mapDACparam(ivar).matlabvar_c); 
            if isfield(graphparam.mapDACparam(ivar),'matlabvar_c_zlevel') % if a vertical level is defined
                if ~isempty(graphparam.mapDACparam(ivar).matlabvar_c_zlevel)
                    zlvl = graphparam.mapDACparam(ivar).matlabvar_c_zlevel;      
                    cdata = cdata(near(nanmean(pres,2),zlvl),:);
                end
            end
            cdatatime = eval(graphparam.mapDACparam(ivar).matlabvar_ctime); 
            data.cdata = interp1(cdatatime(~isnan(cdata)),cdata(~isnan(cdata)),DAC_time);
            data.idac1  = idac1;
            data.idac2  = idac2;    
            data.date1  = date1;
            data.date2  = date2;        
            data.gname = missioninfo.gnickname;
            data.gmission = missioninfo.glmission;   
            data.sectionnb = isecc;


            plot_mission_DAC_track_map(data,graphparam.mapDACparam(ivar),graphparam,step2opt)
        end
    end   
    

    %----------------------------------------------------------------------    
    % Surface plot

     for ivar= 1:length(graphparam.sectionparam) 
        if exist(graphparam.sectionparam(ivar).matlabvar_x,'var') & exist(graphparam.sectionparam(ivar).matlabvar_y,'var') & exist(graphparam.sectionparam(ivar).matlabvar_c,'var') 
            data = [];
            %-----        
            % step 2  
            graphparam.sectionparam(ivar).matlabvar_c
            data.yvarname = graphparam.sectionparam(ivar).matlabvar_y;         
            data.x = eval(graphparam.sectionparam(ivar).matlabvar_x);
            data.y = eval(graphparam.sectionparam(ivar).matlabvar_y); 
            data.c = eval(graphparam.sectionparam(ivar).matlabvar_c); 
            data.id1  = ip1;
            data.id2  = ip2;        
            data0 = data.c(:,ip1:ip2);
            if sum(~isnan(data0(:)))==0
                continue
            end                
            data.date1  = date1;
            data.date2  = date2;        
            data.time  = eval(graphparam.sectionparam(ivar).matlabvar_x); % time of each profile
            data.gname = missioninfo.gnickname;
            data.gmission = missioninfo.glmission;        
            data.sectionnb = isecc;
            data.zalti     = -depth_alti;
            usebathy       = 1;     
            data.zbathy    = -depth_bathy;

            if isfield(graphparam.section,'contourlinevar') % if one wants to add contour lines on the section plot
                data.contourvar = graphparam.section.contourlinevar;
                data.contdata   = eval(graphparam.section.contourlinevar);
                data.xcontdata  = eval(graphparam.section.contourlinex);
                data.ycontdata  = eval(graphparam.section.contourlinex);          
            end
            plot_mission_contourf_section(data,usebathy,graphparam.sectionparam(ivar),graphparam)
            
            data.timeini  = timestart;           
            data.timeend  = timeend;   
            plot_mission_surfplus_section(data,usebathy,graphparam.sectionparam(ivar),graphparam)
        end
     end   
    %----------------------------------------------------------------------    


end

end

