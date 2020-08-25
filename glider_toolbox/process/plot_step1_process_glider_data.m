function plot_step1_process_glider_data(missioninfo,graphparam,step1opt,refDB)
%=========================================================================%
% Function to plot post-processed data for different time window (sections defined by the user), as shifted profiles, time series,
% profiles, surface plot, map, ... 
%
% L. Houpert, SAMS 12/04/2016
%=========================================================================%
%
% step1_process_glider_data(glider,missioninfo,qc_param,step1opt)
% Function that apply QC tests, flag bad data, and create matrix of the vertical profiles. Results are save in
% [missioninfo.step1.procdataDir filesep missioninfo.step1.matfilename]
%
%   Inputs: o missioninfo: [1x1] structure with mission details (defined in
%             users_param/loadosnapmissionparam.m)
%              
%           o graphparam: [1x1] structure containing all the general and specific figure and
%           graph parameters (defined in users_param/graphparamgliderproc.m)
%           
%           o step1opt: [1x1] structure containing:
%                     - .bathyname: Bathymetry used to get an interpolation
%                     along the glider track
%                     - .bathydir: Directory where the bathymetry file is located
%                     - .qflaglvllim : Maximum value of flag kept in the final data (e.g 2)
%                     - .usemanualflag : Use or not the manual QC scripts
%                     or datafile (value = 1 or 0)
%                     - .manualselect_bad: if = 1 run the additional manual QC scripts in 
%                     users_param/manual_qc, if = 0 don't run the script
%                     but load the saved data files
%
%           o refDB: if not empty, this structure contains reference ctd
%           profiles to be plotted together with the glider data. Details of
%           the field of refDB can be seen in clean_data_DBref.m 
%
%
% Soon, planned to install an export function so the data with QC can be
% exported at the EGO_netdf format
% (http://www.ego-network.org/dokuwiki/doku.php?id=public:datamanagement)
% 
% created by L. Houpert (houpertloic@gmail.com), 12/04/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 
% Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)
%
load([missioninfo.step1.procdataDir filesep missioninfo.step1.matfilename '.mat'])

for isecc = 1:length(missioninfo.section)
close all
    
    id1 = near(tgps2,missioninfo.section(isecc).time(1));
    id2 = near(tgpse,missioninfo.section(isecc).time(2));
    
    ip1 = near(timepf,missioninfo.section(isecc).time(1));
    ip2 = near(timepf,missioninfo.section(isecc).time(2));

    idac1 = near(timeDAC,missioninfo.section(isecc).time(1));
    idac2 = near(timeDAC,missioninfo.section(isecc).time(2));
    
    i1 = near(time,tgps2(id1));
    i2 = near(time,tgpse(id2));

    date1 = min(time(i1:i2));
    date2 = max(time(i1:i2));
    
    repfigmain = 'time_series_profiles_sections';
    graphparam.mainplot.dir  = [missioninfo.step1.plotdir filesep missioninfo.glmission filesep repfigmain filesep missioninfo.step1.proclevel];
    graphparam.mainplot.name = missioninfo.step1.plotname;
    if exist(graphparam.mainplot.dir ,'dir')~=7 ; mkdir(graphparam.mainplot.dir);end
    
    %========================================================================%
    % Plot shifted profiles
    if isecc==1 %1 % for the 1st section that correspond to the whole deployment
        
      for ivar=1:length(graphparam.shiftplotparam) 
        data=[];
        %-----        
        % step 1
        shiftrep=['shifted_profiles_' num2str(graphparam.shiftplot.shiftpf_nb)];
        graphparam.shiftplot.dir  = [missioninfo.step1.plotdir filesep missioninfo.glmission filesep shiftrep filesep missioninfo.step1.proclevel];
        graphparam.shiftplot.name = [missioninfo.step1.plotname '_sec' num2str(isecc)];
        
        if exist(graphparam.shiftplot.dir ,'dir')~=7 ; mkdir(graphparam.shiftplot.dir);end
            %--
        if exist(graphparam.shiftplotparam(ivar).matlabvar_x,'var') & exist(graphparam.shiftplotparam(ivar).matlabvar_y,'var') 
            data.x = eval(graphparam.shiftplotparam(ivar).matlabvar_x);
            data.y = eval(graphparam.shiftplotparam(ivar).matlabvar_y);
            data.qc = eval(graphparam.shiftplotparam(ivar).matlabvar_qc);
            data.x00 = eval([graphparam.shiftplotparam(ivar).matlabvar_x '00']);
            if strcmp(graphparam.shiftplotparam(ivar).matlabvar_y,'pres')
                data.y00 = eval([graphparam.shiftplotparam(ivar).matlabvar_y ]);        
                else
                data.y00 = eval([graphparam.shiftplotparam(ivar).matlabvar_y '00']);       
            end
            data.ipf1 = ipf1;
            data.ipf2 = ipf2;        
            data.date1 = date1;
            data.date2 = date2;
            data.divephase = divephase;            
            data.time  = time;
            data.gname = missioninfo.gnickname;
            data.sectionnb = isecc;

            plot_mission_data_shifted_profil_qc(data,graphparam.shiftplotparam(ivar),graphparam)

        else
            continue
        end
        
      end   
    
    %=======================================================================%    
    end     
   
    %----------------------------------------------------------------------
    % Plot time series 
    for ivar=1:length(graphparam.timeseriesparam) 
        if exist(graphparam.timeseriesparam(ivar).matlabvar_x,'var') & exist(graphparam.timeseriesparam(ivar).matlabvar_y,'var')         
            data =[];
            %-----        
            % step 1
            data.x = eval(graphparam.timeseriesparam(ivar).matlabvar_x);
            data.y = eval(graphparam.timeseriesparam(ivar).matlabvar_y);
            data.y00 = eval(graphparam.timeseriesparam(ivar).matlabvar_y00);        
            data.qc = eval(graphparam.timeseriesparam(ivar).matlabvar_qc);
            data.i1  = i1;
            data.i2  = i2;
            data.date1  = date1;
            data.date2  = date2;        
            data.time  = time;
            data.gname = missioninfo.gnickname;
            data.sectionnb = isecc;

            plot_mission_data_time_series_qc(data,graphparam.timeseriesparam(ivar),graphparam)
        else
            continue
        end
    end
    %----------------------------------------------------------------------   
    

    %----------------------------------------------------------------------
    % Plot profiles 
    for ivar= 1:length(graphparam.profilparam) 
        if exist(graphparam.profilparam(ivar).matlabvar_x,'var') & exist(graphparam.profilparam(ivar).matlabvar_y,'var') 
            data = [];
            %-----        
            % step 1
            data.xvarname = graphparam.profilparam(ivar).matlabvar_x;
            data.yvarname = graphparam.profilparam(ivar).matlabvar_y;        
            data.x = eval(graphparam.profilparam(ivar).matlabvar_x);
            data.y = eval(graphparam.profilparam(ivar).matlabvar_y);     
            data.ip1  = ip1;
            data.ip2  = ip2;
            data0 = data.x(:,ip1:ip2);
            if sum(~isnan(data0(:)))==0
                continue
            end
            data.date1  = date1;
            data.date2  = date2;        
            data.time  = timepf;
            data.gldivephase = gldivephase;
            data.gname = missioninfo.gnickname;
            data.sectionnb = isecc;

            plot_mission_data_profils(data,graphparam.profilparam(ivar),graphparam,refDB)
        else
            continue
        end
    end

    %----------------------------------------------------------------------    
    % plot sections
    if (date2-date1)>100 % to reduce the number of dots plotted in scatter plots (if more than 100 days of data to long to process)
    	graphparam.scattersection.scatter_int  = 100; %round(length(i1:i2)/5000);% 100 ;  
    else
    	graphparam.scattersection.scatter_int  = 10;
    end
   
     for ivar= 1:length(graphparam.scatsectionparam) 
        if exist(graphparam.scatsectionparam(ivar).matlabvar_x,'var') & exist(graphparam.scatsectionparam(ivar).matlabvar_y,'var') & exist(graphparam.scatsectionparam(ivar).matlabvar_c,'var') 
            data = [];
            %-----        
            % step 1  
            data.yvarname = graphparam.scatsectionparam(ivar).matlabvar_y;         
            data.x = eval(graphparam.scatsectionparam(ivar).matlabvar_x);
            data.y = eval(graphparam.scatsectionparam(ivar).matlabvar_y); 
            data.c = eval(graphparam.scatsectionparam(ivar).matlabvar_c); 
            data.c00 = eval(graphparam.scatsectionparam(ivar).matlabvar_c00);        
            data.qc = eval(graphparam.scatsectionparam(ivar).matlabvar_qc);
            data.id1  = id1;
            data.id2  = id2;        
            data.i1  = i1;
            data.i2  = i2;
            data.date1  = date1;
            data.date2  = date2;        
            data.time  = timeDAC; % time for the cycle dive+climb
            data.gname = missioninfo.gnickname;
            data.sectionnb = isecc;
            data.zalti     = -altimeter_data;
            usebathy       = 1;     
            data.zbathy    = -zbathydive;

            plot_mission_data_section(data,usebathy,graphparam.scatsectionparam(ivar),graphparam)
        else
            continue
        end
     end   
    %----------------------------------------------------------------------      
    

    
    
    %----------------------------------------------------------------------
    % Plot DAC time series and stickvec     
    if (date2-date1)>200
    	graphparam.timeseriesDAC.stickvecscale  = 0.05; %round(length(i1:i2)/5000);% 100 ;  
    else
    	graphparam.timeseriesDAC.stickvecscale  = 0.25;
    end   
    data =[];
    %-----        
    % step 1
    data.time = timeDAC;
    data.lon = lonDAC;   
    data.lat = latDAC; 
    data.UU = DAC_east;        
    data.VV = DAC_north;
    data.idac1  = idac1;
    data.idac2  = idac2;
    data.date1  = date1;
    data.date2  = date2;        
    data.gname = missioninfo.gnickname;
    data.sectionnb = isecc;

    plot_mission_DAC_time_series(data,graphparam)
  
 
      
    %----------------------------------------------------------------------       
       

    %----------------------------------------------------------------------   
    % plot maps with DAC
    for ivar= 1:length(graphparam.mapDACparam) 
        if exist(graphparam.mapDACparam(ivar).matlabvar_c,'var') 
            data = [];
            data.lon = lonDAC;   
            data.lat = latDAC; 
            data.UU = DAC_east;        
            data.VV = DAC_north;
            cdata = eval(graphparam.mapDACparam(ivar).matlabvar_c); 
            cdatatime = eval(graphparam.mapDACparam(ivar).matlabvar_ctime);           
            data.cdata = interp1(cdatatime(~isnan(cdata)),cdata(~isnan(cdata)),timeDAC);
            data.idac1  = idac1;
            data.idac2  = idac2;  
            data.date1  = date1;
            data.date2  = date2;        
            data.gname = missioninfo.gnickname;
            data.gmission = missioninfo.glmission;
            data.sectionnb = isecc;

            plot_mission_DAC_track_map(data,graphparam.mapDACparam(ivar),graphparam,step1opt)
        else
            continue
        end
    end   
    
    
%----------------------------------------------------------------------       
end

end



