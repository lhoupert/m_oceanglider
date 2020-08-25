function offsetglider = step1b_glider_offset_from_obs(missioninfo,otherobs,tol,figparam)
%
%=========================================================================%
% Function to compare ctd and DAC from glider data to other obs (e.g. mooring), 
% plot and save the offsets 
%
% L. Houpert, SAMS 25/04/2016
%=========================================================================%
%
% offsetglider = step1b_glider_offset_from_obs(missioninfo,otherobs,tol,figparam)
% Function that compares ctd data from glider with other obs (e.g. mooring), plots and saves the offsets 
%
%   Inputs: o missioninfo: [1x1] structure with mission details (defined in
%             users_param/loadosnapmissionparam.m)
%
%           o otherobs: structure with other data to be compared to glider data
%
%           o tol: define the threshold for determining the closest point in time and space 
%          (defined by tol.time (in days), tol.dist (in km) and tol.dist2; 
%          eg: tol.time = 5, tol.dist=50, tol.time2 = 2, tol.dist2 =10)           
%
%           o figparam: [1x1] structure containing all the general and specific figure and
%           graph parameters (defined in users_param/graphparamgliderproc.m)
%   
%   Output: offsetglider: structure containing the different glider offset
%   estimated from the crossing point(s) (can be several for a mooring) between a 
%   glider mission and another dataset, the temperature and salinity
%   profiles for both dataset and for the selected time period are also
%   include
%
%
%
% created by L. Houpert (houpertloic@gmail.com), 25/04/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 
% Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)


%-------------------------------------------------------------------------
eval(['gldata=load(''' missioninfo.step1.procdataDir filesep missioninfo.step1.matfilename '.mat'');' ])  


glPP = (gldata.PP)';
glSS = gldata.SS;%-0.2;
glTP = gldata.TP; %glTP = sw_ptmp(glSS,gldata.TT,glPP,0);
glbathydive2=nan(length( gldata.zbathydive)*2,1)';
glbathydive2(1:2:length(glbathydive2)) = gldata.zbathydive;
glbathydive2(2:2:length(glbathydive2)) = gldata.zbathydive;


switch otherobs.type
    case 'mooring'
        stepmoor = 1 ;%round(0.25/(otherobs.time(1,10)-otherobs.time(1,9))); %2-day timestep so not to many comparison point
        timeotherobs2 = otherobs.time(1,1:stepmoor:end)'; 
        lonotherobs2  = repmat(otherobs.lon,1,length(timeotherobs2))';
        latotherobs2  = repmat(otherobs.lat,1,length(timeotherobs2))';
        zotherobs2    = repmat(otherobs.depth(1),1,length(timeotherobs2))';
        for ijk=1:length(gldata.lonpf)
            distglmoor(ijk)=sw_dist([gldata.latpf(ijk) otherobs.lat],[gldata.lonpf(ijk) otherobs.lon],'km');
        end
        presmean = nanmean(otherobs.pres,2);
        ipres1000 = find(presmean<1001);
        presotherobs2  = otherobs.pres(ipres1000,1:stepmoor:end); 
        salotherobs2  = otherobs.sal(ipres1000,1:stepmoor:end); 
        ptempotherobs2  = otherobs.ptemp(ipres1000,1:stepmoor:end); 

        depththreshold = 1/5; %1/5 of the depth of the mooring

        isel0 = find(gldata.lonpf>otherobs.lon-2 & gldata.lonpf<otherobs.lon+2) ; % keep inly the data close to the mooring at 1deg in longitude
        [indX00,indY00,distmat00] = findclosest(gldata.timepf(isel0),gldata.lonpf(isel0),gldata.latpf(isel0),timeotherobs2,lonotherobs2,latotherobs2,tol.time,tol.dist);

        if ~isempty(indX00)
            for ijk=1:length(indX00)
                indXnodepth{ijk}  = isel0(unique(indX00{ijk}));
                indYnodepth{ijk} = indY00{ijk}; % not true if stepmoor ~= 1
                distgl=[];       
                for ikk=1:length(indXnodepth{ijk})
                    igg = indXnodepth{ijk}(ikk);
                    imm = indYnodepth{ijk}(ikk);
                    distgl(ikk) = sw_dist([gldata.latpf(igg) latotherobs2(imm)],[gldata.lonpf(igg) lonotherobs2(imm)],'km');
                end 
                distspacenodepth{ijk} = distgl;
            end
        end
        taxmin = min(gldata.timepf(isel0)); 
        taxmax = max(gldata.timepf(isel0)); 

        % taking into account bathymetry
        isel = find(gldata.lonpf>otherobs.lon-2 & gldata.lonpf<otherobs.lon+2 & glbathydive2 < otherobs.depth(1)*(1 - depththreshold) &  glbathydive2 > otherobs.depth(1)*(1 + depththreshold) ) ; % keep inly the data close to the mooring at 1deg in longitude
        [indX0,indY0,distmat] = findclosest(gldata.timepf(isel),gldata.lonpf(isel),gldata.latpf(isel),timeotherobs2,lonotherobs2,latotherobs2,tol.time,tol.dist);

        if ~isempty(indX0)

            for ijk=1:length(indX0)
                indX{ijk} = isel(unique(indX0{ijk}));
                indY{ijk} = unique(indY0{ijk}); % not true if stepmoor ~= 1
                distgl=[];
                for ikk=1:length(indX{ijk})
                    igg = indX{ijk}(ikk);
                    imm = indY{ijk}(ikk);
                    distgl(ikk) = sw_dist([gldata.latpf(igg) latotherobs2(imm)],[gldata.lonpf(igg) lonotherobs2(imm)],'km');
                end 
                distspace{ijk} = distgl;        
            end


            dirplot  = [missioninfo.step1.plotdir filesep missioninfo.glmission filesep otherobs.figdirname];
            if exist(dirplot ,'dir')~=7 ; mkdir(dirplot);end

            figsuf   = [missioninfo.glmission '_' otherobs.name];

            %========================================================================
            % plot glider-otherobs
            %------------------------------------------------------------------------
            figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
            set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
              'Paperposition',figparam.figpos)  
            subplot(1,2,1)
            hold on
            plot(gldata.lonpf(isel0),gldata.timepf(isel0),'-.','color',[1 1 1]*0.65)     
            for ijk=1:length(indX)
                plot(gldata.lonpf(indX{ijk}),gldata.timepf(indX{ijk}),'-x')
                plot(lonotherobs2(indY{ijk}),timeotherobs2(indY{ijk}),'k--')
            end
            set(gca,'ylim',[taxmin taxmax])
            datetick('y','keeplimits')
            xlabel('Longitude')
            ylabel('Time')

            %-
            subplot(1,2,2)
            hold on
            plot(gldata.timepf(isel0),gldata.latpf(isel0),'-.','color',[1 1 1]*0.65)     
            for ijk=1:length(indX)
                plot(gldata.timepf(indX{ijk}),gldata.latpf(indX{ijk}),'-x')
                plot(timeotherobs2(indY{ijk}),latotherobs2(indY{ijk}),'k--')
            end
            set(gca,'xlim',[taxmin taxmax])
            datetick('x','keeplimits')
            xlabel('Time')
            ylabel('Latitude')

            set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
            print(gcf,'-dpng',[dirplot filesep figsuf '_' figparam.crossing.lontimeplot.name ]);

            %------------------------------------------------------------------------
            figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
            set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
              'Paperposition',figparam.figpos)  
            subplot(2,1,1)
            hold on
            plot(gldata.timepf(isel0),glbathydive2(isel0),'-.','color',[1 1 1]*0.65)        
            for ijk=1:length(indX)
                plot(timeotherobs2(indY{ijk}),zotherobs2(indY{ijk}),'k--')        
                plot(gldata.timepf(indX{ijk}), glbathydive2(indX{ijk}),'-x');
            end
            set(gca,'xlim',[taxmin taxmax])
            datetick('x','keeplimits')
            xlabel('Time')
            ylabel('Depth (m)')       
            %-
            subplot(2,1,2)
            hold on
            plot(gldata.timepf(isel0),distglmoor(isel0),'-.','color',[1 1 1]*0.65)       
            for ijk=1:length(indX)
                plot(gldata.timepf(indX{ijk}), distspace{ijk},'-x');
            end
            set(gca,'xlim',[taxmin taxmax])
            datetick('x','keeplimits') 
            xlabel('Time')
            ylabel('Distance (km)') 
            set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
            print(gcf,'-dpng',[dirplot filesep figsuf  '_' figparam.crossing.timedistplot.name ]);
            %-------------------------------------------------------------------------

            
            
            %=========================================================================
            %=========================================================================
            if otherobs.CTDcomparison == 1
            %======================================================================================
            % Creation of a structure summarizing the different offset for each level and 
            % plot glider-otherobs TS and profils
            %--------------------------------------------------------------------------------------
            for ijk=1:length(indX)
                indglider = indX{ijk};
                indother  = indY{ijk};
                % define closer profiles in space and time
                indglider2 = indglider(distspace{ijk} < tol.dist2); % glider profile less than 10km
                indglider3 = indglider(distspace{ijk} == min(distspace{ijk}));
                dmoorini = timeotherobs2(indother(1));
                dmoorfff = timeotherobs2(indother(end));
                if ~isempty(indglider2)        
                    indother2  = indother(timeotherobs2(indother) > min(gldata.timepf(indglider2)) & ...
                                     timeotherobs2(indother) < max(gldata.timepf(indglider2)));
                    indother3  = indother(near(timeotherobs2(indother), gldata.timepf(indglider3)));
                    dmoorini2 = timeotherobs2(indother2(1));
                    dmoorfff2 = timeotherobs2(indother2(end));   
                    dmoorini3 = timeotherobs2(indother3);
                    dmoorfff3 = timeotherobs2(indother3);
                    differentcase = 1:3';
                else
                    differentcase = 1;
                end    
                    %-----------------------------------------------------------------------------------
                    % Calculate difference for each level of measurement of the
                    % other dataset statistics for the different distance criteria
                    % definition of the structure summarizing the offset
                    offsetglider(ijk).glidername = missioninfo.gnickname;
                    offsetglider(ijk).glidermission = missioninfo.glmission;
                    offsetglider(ijk).otherdatasettype = otherobs.type;
                    offsetglider(ijk).otherdatasetname = otherobs.name;
                    %------------------------------------------------------------------------                     
                    for issel = differentcase
                        switch issel
                            case 1 
                                indg = indglider;
                                indo = indother;
                                offsetglider(ijk).selcrit(issel).toltime = tol.time;
                                offsetglider(ijk).selcrit(issel).toldist = tol.dist;
                                offsetglider(ijk).selcrit(issel).timerange = [dmoorini dmoorfff];
                            case 2
                                indg = indglider2;
                                indo = indother2;
                                offsetglider(ijk).selcrit(issel).toltime = tol.time2;
                                offsetglider(ijk).selcrit(issel).toldist = tol.dist2;   
                                offsetglider(ijk).selcrit(issel).timerange = [dmoorini2 dmoorfff2];                        
                            case 3
                                indg = indglider3;
                                indo = indother3;                 
                                offsetglider(ijk).selcrit(issel).toltime = 0;
                                offsetglider(ijk).selcrit(issel).toldist = min(distspace{ijk});
                                offsetglider(ijk).selcrit(issel).timerange = [dmoorini3 dmoorfff3];                           
                        end
                        meanmoorlvl = nanmean(presotherobs2(:,indo),2);
                        stdmoorlvl  = nanstd(presotherobs2(:,indo),0,2);
                        iglpres = [];
                        indotherpres = [];
                        for ipp=1:length(meanmoorlvl)
                            if abs(glPP(near(glPP,meanmoorlvl(ipp))) - meanmoorlvl(ipp))<2
                                iglpres = [iglpres near(glPP,meanmoorlvl(ipp))];
                                indotherpres = [indotherpres ipp];
                            end
                        end

                        % use of nanmean and nanstd in case no data for some
                        % specific profiles
                        glmeantpot = nanmean(glTP(iglpres,indg),2);
                        glstdtpot  = nanstd(glTP(iglpres,indg),0,2);
                        glmeansal  = nanmean(glSS(iglpres,indg),2);
                        glstdsal   = nanstd(glSS(iglpres,indg),0,2);
                        omeantpot  = nanmean(ptempotherobs2(indotherpres,indo),2);
                        ostdtpot  = nanstd(ptempotherobs2(indotherpres,indo),0,2);
                        omeansal  = nanmean(salotherobs2(indotherpres,indo),2);
                        ostdsal   = nanstd(salotherobs2(indotherpres,indo),0,2);

                        i700 = find(meanmoorlvl(indotherpres) > 700); %700m due to RTEB1

                        offsetglider(ijk).selcrit(issel).meanpreslvl = meanmoorlvl(indotherpres);
                        offsetglider(ijk).selcrit(issel).stdpreslvl = stdmoorlvl(indotherpres);
                        offsetglider(ijk).selcrit(issel).glider.tpot.mean = glmeantpot;
                        offsetglider(ijk).selcrit(issel).glider.tpot.std = glstdtpot;
                        offsetglider(ijk).selcrit(issel).glider.tpot.numpf = sum(~isnan(glTP(iglpres,indg)),2);                 
                        offsetglider(ijk).selcrit(issel).glider.sal.mean = glmeansal;
                        offsetglider(ijk).selcrit(issel).glider.sal.std = glstdsal;  
                        offsetglider(ijk).selcrit(issel).glider.sal.numpf = sum(~isnan(glSS(iglpres,indg)),2);                  
                        offsetglider(ijk).selcrit(issel).otherobs.tpot.mean = omeantpot;
                        offsetglider(ijk).selcrit(issel).otherobs.tpot.std = ostdtpot; 
                        offsetglider(ijk).selcrit(issel).otherobs.tpot.numpf = sum(~isnan(ptempotherobs2(indotherpres,indo)),2);                   
                        offsetglider(ijk).selcrit(issel).otherobs.sal.mean = omeansal;
                        offsetglider(ijk).selcrit(issel).otherobs.sal.std = ostdsal;   
                        offsetglider(ijk).selcrit(issel).otherobs.sal.numpf = sum(~isnan(salotherobs2(indotherpres,indo)),2);                   


                        % difference
                        offsetglider(ijk).selcrit(issel).difftpot = glmeantpot -omeantpot ; 
                        offsetglider(ijk).selcrit(issel).difftpotsumstd = sqrt(glstdtpot.^2 + ostdtpot.^2) ; 
                        offsetglider(ijk).selcrit(issel).diffsal = glmeansal -omeansal ;
                        offsetglider(ijk).selcrit(issel).diffsalsumstd = sqrt(glstdsal.^2 + ostdsal.^2) ; ;               

                        offsetglider(ijk).selcrit(issel).tpotdiffmean = nanmean(glmeantpot -omeantpot);
                        offsetglider(ijk).selcrit(issel).tpotdiffstd  = [];%nanstd(glmeantpot -omeantpot,0);
                        offsetglider(ijk).selcrit(issel).tpotdiffmeanstd = sqrt(nanmean(offsetglider(ijk).selcrit(issel).difftpotsumstd.^2));    
                        offsetglider(ijk).selcrit(issel).saldiffmean  = nanmean(glmeansal -omeansal);
                        offsetglider(ijk).selcrit(issel).saldiffstd   = [];%nanstd(glmeansal -omeansal,0);    
                        offsetglider(ijk).selcrit(issel).saldiffmeanstd = [];%sqrt(nanmean(offsetglider(ijk).selcrit(issel).diffsalsumstd.^2));

                        offsetglider(ijk).selcrit(issel).tpotdiffmeanbelow700    = nanmean(glmeantpot(i700) -omeantpot(i700));
                        offsetglider(ijk).selcrit(issel).tpotdiffstdbelow700     = [];%nanstd(glmeantpot(i700) -omeantpot(i700),0);
                        offsetglider(ijk).selcrit(issel).tpotdiffmeanstdbelow700 = sqrt(nanmean(offsetglider(ijk).selcrit(issel).difftpotsumstd(i700).^2));
                        offsetglider(ijk).selcrit(issel).saldiffmeanbelow700     = nanmean(glmeansal(i700) -omeansal(i700));
                        offsetglider(ijk).selcrit(issel).saldiffstdbelow700      = [];%nanstd(glmeansal(i700) -omeansal(i700),0);      
                        offsetglider(ijk).selcrit(issel).saldiffmeanstdbelow700  =  sqrt(nanmean(offsetglider(ijk).selcrit(issel).diffsalsumstd(i700).^2));             
                    end
                                
                %-----------------------------------------------------------------------------------


                lim_ax_sal = [min([min(glSS(:,indg)) min(salotherobs2(:,indo))]) ...
                              max([max(glSS(:,indg)) max(salotherobs2(:,indo))])] ; 
                lim_ax_ptemp = [min([min(glTP(:,indg)) min(ptempotherobs2(:,indo))]) ...
                              max([max(glTP(:,indg)) max(ptempotherobs2(:,indo))])] ; 

                if ~isempty(indglider2) 
                    titlestr = {[otherobs.name ' (+-) - ' missioninfo.gnickname ' (-) comparisons for ' ...
                            datestr(dmoorini,'dd/mm/yy') ' - ' datestr(dmoorfff,'dd/mm/yy') ' (grey lines and crosses), '], ...
                            [datestr(dmoorini2,'dd/mm/yy') ' - ' datestr(dmoorfff2,'dd/mm/yy') ' (color lines and crosses)' ]};
                else
                    titlestr = {[otherobs.name ' (+-) - ' missioninfo.gnickname ' (-) comparisons for ' ...
                            datestr(dmoorini,'dd/mm/yy') ' - ' datestr(dmoorfff,'dd/mm/yy') ' (grey lines and crosses)']};            
                end
                %---------------------------------------------------------------------
                % Plot Profils
                figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
                set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
                title(titlestr)
                subplot(1,2,1)
                hold on
                plot(glTP(:,indglider),glPP,'color',[1 1 1]*0.65)
                plot(ptempotherobs2(:,indother),presotherobs2(:,indother),'+','color',[1 1 1]*0.45)
                if ~isempty(indglider2)
                    plot(glTP(:,indglider2),glPP,'b')    
                    plot(glTP(:,indglider3),glPP,'color','k','linewidth',2)      
                    plot(ptempotherobs2(:,indother2),presotherobs2(:,indother2),'-.+r')  
                    plot(ptempotherobs2(:,indother3),presotherobs2(:,indother3),'--xk','linewidth',2)    
                end
                set(gca,'ydir','reverse')
                ylabel('Pressure')
                xlabel('Pot. Temp')

                subplot(1,2,2)
                hold on
                plot(glSS(:,indglider),glPP,'color',[1 1 1]*0.65)
                plot(salotherobs2(:,indother),presotherobs2(:,indother),'+','color',[1 1 1]*0.45)    
                if ~isempty(indglider2)  
                    plot(glSS(:,indglider2),glPP,'b')    
                    plot(glSS(:,indglider3),glPP,'color','k','linewidth',2)      
                    plot(salotherobs2(:,indother2),presotherobs2(:,indother2),'-.+r')  
                    plot(salotherobs2(:,indother3),presotherobs2(:,indother3),'--xk','linewidth',2)    
                end
                ylabel('Pressure')
                xlabel('Salinity')   
                set(gca,'ydir','reverse')

                annotation('textbox', [0 0.9 1 0.1], ...
                'String', titlestr, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center')

                set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,   

                set(gcf,'render','painters')
                print(gcf,'-dpng',[dirplot filesep figsuf  '_' figparam.crossing.profilsplot.name '_' num2str(ijk)]);
                if figparam.epsgraph == 1
                    print(gcf,'-depsc2',[dirplot filesep figsuf  '_' figparam.crossing.profilsplot.name '_' num2str(ijk)]);    
                end
                %---------------------------------------------------------------------
                % Plot TS    
                figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
                set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
                hold on
                contour(figparam.contTS.ssc,figparam.contTS.ttc,...
                          figparam.contTS.ddc,figparam.contTS.vcont,'k','linewidth',figparam.contTS.linewidth);       
                plot(glSS(:,indglider),glTP(:,indglider),'color',[1 1 1]*0.65)
                plot(salotherobs2(:,indother),ptempotherobs2(:,indother),'+','color',[1 1 1]*0.45)    
                if ~isempty(indglider2)
                    plot(glSS(:,indglider2),glTP(:,indglider2),'b')    
                    plot(glSS(:,indglider3),glTP(:,indglider3),'color','k','linewidth',2)      
                    plot(salotherobs2(:,indother2),ptempotherobs2(:,indother2),'-.+r')  
                    plot(salotherobs2(:,indother3),ptempotherobs2(:,indother3),'--xk','linewidth',2)    
                end       
                ylabel('Pot. Temperature')
                xlabel('Salinity')
                set(gca,'xlim',lim_ax_sal,'ylim',lim_ax_ptemp)

                title(titlestr);

                set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%, 
                set(gcf,'render','painters')
                print(gcf,'-dpng',[dirplot filesep figsuf  '_' figparam.crossing.TSplot.name '_' num2str(ijk)]);
                if figparam.epsgraph == 1
                    print(gcf,'-depsc2',[dirplot filesep figsuf  '_' figparam.crossing.TSplot.name '_' num2str(ijk)]);    
                end
            end
            %=====================================================================
            % Print statistics file
            %--------------------------------------------------------------------
            % 1. calculate difference between glider data and the other dataset at the
            % level of measurement of the other data set and the mean difference of
            % the whole profile
            outputdir = [missioninfo.step1.procdataDir filesep 'calib' filesep missioninfo.glmission]
            if exist(outputdir)~=7
                mkdir(outputdir)
            end
            statsfile = figsuf;
            for ijk=1:length(offsetglider)
                for ilm=1:length(offsetglider(ijk).selcrit)
                    offsetglider(ijk).glidername = missioninfo.gnickname;
                    offsetglider(ijk).otherdatasettype = otherobs.type;
                    offsetglider(ijk).otherdatasetname = otherobs.name;
                    selectparamstr = [ num2str(offsetglider(ijk).selcrit(ilm).toltime,'%4.0f') 'days_' num2str(offsetglider(ijk).selcrit(ilm).toldist,'%4.0f') 'km'];
                    %-------------------------------------------------------------------------
                    fidlog    = fopen([outputdir filesep statsfile '_' num2str(ijk) '_' selectparamstr '.log'],'w');
                    fprintf(fidlog,['Glider Name                  : ' missioninfo.gnickname '\n']);
                    fprintf(fidlog,['Glider Mission               : ' missioninfo.glmission '\n']);
                    numglpfs = nanmean(offsetglider(ijk).selcrit(ilm).glider.sal.numpf);
                    numotherpfs = nanmean(offsetglider(ijk).selcrit(ilm).otherobs.sal.numpf);                   
                    fprintf(fidlog,['Nber of glider profiles used : ' num2str(numglpfs,'%4.0f') '\n']);                    
                    fprintf(fidlog,['Other Dataset Type           : ' offsetglider(ijk).otherdatasettype '\n']);
                    fprintf(fidlog,['Other Dataset Name           : ' offsetglider(ijk).otherdatasetname '\n']);
                    fprintf(fidlog,['Nber of other dataset pfs.   : ' num2str(numotherpfs,'%4.0f') '\n']);                       
                    fprintf(fidlog,['Time Window for comparison   : ' num2str(offsetglider(ijk).selcrit(ilm).toltime,'%4.0f') ' (days) \n']);     
                    fprintf(fidlog,['Max. Distance for comparison : ' num2str(offsetglider(ijk).selcrit(ilm).toldist,'%4.0f') ' (km) \n']);                 
                    fprintf(fidlog,'\n');  
                    if offsetglider(ijk).selcrit(ilm).toltime ~=0 % 
                        fprintf(fidlog,[offsetglider(ijk).otherdatasettype ' mean depth;' offsetglider(ijk).otherdatasettype ' std depth; Std Mean Pot. Temp ;Mean Pot. Temp. Difference; Std Pot. Temp. Difference;' ...
                                    'Std Mean Salinity; Mean Salinity Difference; Std Salinity Difference \n']);
                        fprintf(fidlog,[' ; ; (other dataset); (glider - otherdataset); sqrt(glstdtpot^2 + ostdtpot^2) ;' ...
                                    '(other dataset); (glider - otherdataset); sqrt(glstdtpot^2 + ostdtpot^2)  \n']);                                
                        for ipp=1:length( offsetglider(ijk).selcrit(ilm).meanpreslvl)
                               fprintf(fidlog,'%4.3f;%4.3f;%4.3f;%4.3f;%4.3f;%4.3f;%4.3f;%4.3f\n',offsetglider(ijk).selcrit(ilm).meanpreslvl(ipp),...
                                   offsetglider(ijk).selcrit(ilm).stdpreslvl(ipp), offsetglider(ijk).selcrit(ilm).otherobs.tpot.std(ipp),...
                                   offsetglider(ijk).selcrit(ilm).difftpot(ipp), offsetglider(ijk).selcrit(ilm).difftpotsumstd(ipp), ...
                                   offsetglider(ijk).selcrit(ilm).otherobs.sal.std(ipp),offsetglider(ijk).selcrit(ilm).diffsal(ipp),...
                                   offsetglider(ijk).selcrit(ilm).diffsalsumstd(ipp));
                        end

                    else
                         fprintf(fidlog,[offsetglider(ijk).otherdatasettype ' depth;Pot. Temp. Difference;' ...
                                    'Salinity Difference\n']);               
                         for ipp=1:length( offsetglider(ijk).selcrit(ilm).meanpreslvl)
                               fprintf(fidlog,'%4.3f;%4.3f;%4.3f\n',offsetglider(ijk).selcrit(ilm).meanpreslvl(ipp),...
                                   offsetglider(ijk).selcrit(ilm).difftpot(ipp),offsetglider(ijk).selcrit(ilm).diffsal(ipp));
                         end               

                    end

                    fprintf(fidlog,'\n');   
                    fprintf(fidlog,'Depth levels Range;;;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;;Mean Sal. Diff.; Std Sal. Diff.\n');   
                    fprintf(fidlog,'%s;;;%4.3f;%4.3f;;%4.3f;%4.3f\n','All  Depths',...
                                   offsetglider(ijk).selcrit(ilm).tpotdiffmean,offsetglider(ijk).selcrit(ilm).tpotdiffstd,...
                                   offsetglider(ijk).selcrit(ilm).saldiffmean,offsetglider(ijk).selcrit(ilm).saldiffstd); 
                    fprintf(fidlog,'\n');   
                    fprintf(fidlog,'Depth levels Range;;;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;;Mean Sal. Diff.; Std Sal. Diff.\n');    
                    fprintf(fidlog,'%s;;;%4.3f;%4.3f;;%4.3f;%4.3f\n','700-1000m',...
                                   offsetglider(ijk).selcrit(ilm).tpotdiffmeanbelow700,offsetglider(ijk).selcrit(ilm).tpotdiffstdbelow700,...
                                   offsetglider(ijk).selcrit(ilm).saldiffmeanbelow700,offsetglider(ijk).selcrit(ilm).saldiffstdbelow700); 

                    fclose(fidlog);                

                end
            end
            % end

            % save matlab version of the calib
            save([outputdir filesep statsfile],'offsetglider')

            %-------------------------------------------------------------------------
          end      % plot glider-otherobs salinity-salinity
     else

        disp(' ')
        disp(['No crossing points for ' missioninfo.gnickname ' and ' otherobs.name '. Check time limits of the 2 datasets'])
        disp(' ')
        offsetglider = [];
    end   

    if otherobs.UVspeedcomparison == 1
    %-------------------------------------------------------------------------
    % plot glider-otherobs speed

    end        
        
        
   case 'ctd'   
       
        timediff = nan(length(gldata.timepf),length(otherobs.time));
        distdiff = nan(length(gldata.timepf),length(otherobs.time));
        for ijk=1:length(timediff(:,1))
           timediff(ijk,:) =  abs(repmat(gldata.timepf(ijk),1,length(otherobs.time)) - otherobs.time);
           for icc = 1:length(otherobs.time)
                distdiff(ijk,icc) =  sw_dist([gldata.latpf(ijk) otherobs.lat(icc)],[gldata.lonpf(ijk) otherobs.lon(icc)],'km');
           end
        end
        
  
        distlim = 20; %in km
        timelim = 3; %in days

        [iglt,ictd,distdiff,timediff] = findcross_space(gldata.timepf,gldata.lonpf,gldata.latpf,otherobs.time,otherobs.lon,otherobs.lat,timelim,distlim);     
        
  if ~isempty(ictd)
        timeotherobs2 = otherobs.time(ictd);
        lonotherobs2 = otherobs.lon(ictd);
        latotherobs2 = otherobs.lat(ictd);        
        zotherobs2   = otherobs.depth(ictd);    
        for ijk=1:length(gldata.lonpf)
            distglmoor(ijk)=sw_dist([gldata.latpf(ijk) latotherobs2],[gldata.lonpf(ijk) latotherobs2],'km');
        end 
        
        
        presmean = otherobs.pres(:,ictd);
        ipres1000 = find(presmean<1001);
        presotherobs2  = otherobs.pres(ipres1000,ictd); 
        salotherobs2  = otherobs.sal(ipres1000,ictd); 
        ptempotherobs2  = otherobs.ptemp(ipres1000,ictd); 

        depththreshold = 1/5; %1/5 of the depth of the mooring

        isel0 = find(gldata.lonpf>lonotherobs2-2 & gldata.lonpf<lonotherobs2+2) ; % keep inly the data close to the mooring at 1deg in longitude
        [indX00,indY00,distmat00] = findclosest(gldata.timepf(isel0),gldata.lonpf(isel0),gldata.latpf(isel0),timeotherobs2,lonotherobs2,latotherobs2,tol.time,tol.dist);
  
  
        if ~isempty(indX00)
            for ijk=1:length(indX00)
                indXnodepth{ijk}  = isel0(unique(indX00{ijk}));
                indYnodepth{ijk} = indY00{ijk}; % not true if stepmoor ~= 1
                distgl=[];       
                for ikk=1:length(indXnodepth{ijk})
                    igg = indXnodepth{ijk}(ikk);
                    imm = indYnodepth{ijk}(ikk);
                    distgl(ikk) = sw_dist([gldata.latpf(igg) latotherobs2(imm)],[gldata.lonpf(igg) lonotherobs2(imm)],'km');
                end 
                distspacenodepth{ijk} = distgl;
            end
        end
        taxmin = min(gldata.timepf(isel0)); 
        taxmax = max(gldata.timepf(isel0)); 

        % taking into account bathymetry
%        isel = find(gldata.lonpf>lonotherobs2-2 & gldata.lonpf<lonotherobs2+2 & glbathydive2 < otherobs.depth(1)*(1 - depththreshold) &  glbathydive2 > otherobs.depth(1)*(1 + depththreshold) ) ; % keep inly the data close to the mooring at 1deg in longitude
        isel = find(gldata.lonpf>lonotherobs2-2 & gldata.lonpf<lonotherobs2+2) ; % keep inly the data close to the mooring at 1deg in longitude
 
        [indX0,indY0,distmat] = findclosest(gldata.timepf(isel),gldata.lonpf(isel),gldata.latpf(isel),timeotherobs2,lonotherobs2,latotherobs2,tol.time,tol.dist);
  
  else
      indX0=[];
  end
  
        if ~isempty(indX0)

            for ijk=1:length(indX0)
                indX{ijk} = isel(unique(indX0{ijk}));
                indY{ijk} = unique(indY0{ijk}); % not true if stepmoor ~= 1
                distgl=[];
                for ikk=1:length(indX{ijk})
                    igg = indX{ijk}(ikk);
                    imm = indY{ijk};
                    distgl(ikk) = sw_dist([gldata.latpf(igg) latotherobs2(imm)],[gldata.lonpf(igg) lonotherobs2(imm)],'km');
                end 
                distspace{ijk} = distgl;        
            end


            dirplot  = [missioninfo.step1.plotdir filesep missioninfo.glmission filesep otherobs.figdirname];
            if exist(dirplot ,'dir')~=7 ; mkdir(dirplot);end

            figsuf   = [missioninfo.glmission '_' otherobs.name];

            %========================================================================
            % plot glider-otherobs
            %------------------------------------------------------------------------
            figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
            set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
              'Paperposition',figparam.figpos)  
            subplot(1,2,1)
            hold on
            plot(gldata.lonpf(isel0),gldata.timepf(isel0),'-.','color',[1 1 1]*0.65)     
            for ijk=1:length(indX)
                plot(gldata.lonpf(indX{ijk}),gldata.timepf(indX{ijk}),'-x')
                plot(lonotherobs2(indY{ijk}),timeotherobs2(indY{ijk}),'k--')
            end
            set(gca,'ylim',[taxmin taxmax])
            datetick('y','keeplimits')
            xlabel('Longitude')
            ylabel('Time')

            %-
            subplot(1,2,2)
            hold on
            plot(gldata.timepf(isel0),gldata.latpf(isel0),'-.','color',[1 1 1]*0.65)     
            for ijk=1:length(indX)
                plot(gldata.timepf(indX{ijk}),gldata.latpf(indX{ijk}),'-x')
                plot(timeotherobs2(indY{ijk}),latotherobs2(indY{ijk}),'k--')
            end
            set(gca,'xlim',[taxmin taxmax])
            datetick('x','keeplimits')
            xlabel('Time')
            ylabel('Latitude')

            set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
            print(gcf,'-dpng',[dirplot filesep figsuf '_' figparam.crossing.lontimeplot.name ]);

            %------------------------------------------------------------------------
            figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
            set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
              'Paperposition',figparam.figpos)  
            subplot(2,1,1)
            hold on
            plot(gldata.timepf(isel0),glbathydive2(isel0),'-.','color',[1 1 1]*0.65)        
            for ijk=1:length(indX)
                plot(timeotherobs2(indY{ijk}),zotherobs2(indY{ijk}),'k--')        
                plot(gldata.timepf(indX{ijk}), glbathydive2(indX{ijk}),'-x');
            end
            set(gca,'xlim',[taxmin taxmax])
            datetick('x','keeplimits')
            xlabel('Time')
            ylabel('Depth (m)')       
            %-
            subplot(2,1,2)
            hold on
            plot(gldata.timepf(isel0),distglmoor(isel0),'-.','color',[1 1 1]*0.65)       
            for ijk=1:length(indX)
                plot(gldata.timepf(indX{ijk}), distspace{ijk},'-x');
            end
            set(gca,'xlim',[taxmin taxmax])
            datetick('x','keeplimits') 
            xlabel('Time')
            ylabel('Distance (km)') 
            set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
            print(gcf,'-dpng',[dirplot filesep figsuf  '_' figparam.crossing.timedistplot.name ]);
            %-------------------------------------------------------------------------

            
            
            %=========================================================================
            %=========================================================================
            if otherobs.CTDcomparison == 1
            %======================================================================================
            % Creation of a structure summarizing the different offset for each level and 
            % plot glider-otherobs TS and profils
            %--------------------------------------------------------------------------------------

            for ijk=1:length(indX)
                indglider = indX{ijk};
                indother  = indY{ijk};
                % define closer profiles in space and time
                indglider2 = indglider(distspace{ijk} < tol.dist2); % glider profile less than 10km
                indglider3 = indglider(distspace{ijk} == min(distspace{ijk}));
                dmoorini = timeotherobs2(indother(1));
                dmoorfff = timeotherobs2(indother(end));
                if ~isempty(indglider2)        
                    indother2  = indother; % indother(near(timeotherobs2(indother), gldata.timepf(indglider2)));
                    indother3  = indother; % indother(near(timeotherobs2(indother), gldata.timepf(indglider3)));
                    dmoorini2 = timeotherobs2(indother2(1));
                    dmoorfff2 = timeotherobs2(indother2(end));   
                    dmoorini3 = timeotherobs2(indother3);
                    dmoorfff3 = timeotherobs2(indother3);
                    differentcase = 1:3';
                else
                    differentcase = 1;
                end    
                    %-----------------------------------------------------------------------------------
                    % Calculate difference for each level of measurement of the
                    % other dataset statistics for the different distance criteria
                    % definition of the structure summarizing the offset
                    offsetglider(ijk).glidername = missioninfo.gnickname;
                    offsetglider(ijk).glidermission = missioninfo.glmission;
                    offsetglider(ijk).otherdatasettype = otherobs.type;
                    offsetglider(ijk).otherdatasetname = otherobs.name;
                    %-------------------------------------------------------------------------
                    for issel = differentcase
                        switch issel
                            case 1 
                                indg = indglider;
                                indo = indother;
                                offsetglider(ijk).selcrit(issel).toltime = tol.time;
                                offsetglider(ijk).selcrit(issel).toldist = tol.dist;
                                offsetglider(ijk).selcrit(issel).timerange = [dmoorini dmoorfff];
                            case 2
                                indg = indglider2;
                                indo = indother2;
                                offsetglider(ijk).selcrit(issel).toltime = tol.time2;
                                offsetglider(ijk).selcrit(issel).toldist = tol.dist2;   
                                offsetglider(ijk).selcrit(issel).timerange = [dmoorini2 dmoorfff2];                        
                            case 3
                                indg = indglider3;
                                indo = indother3;                 
                                offsetglider(ijk).selcrit(issel).toltime = 0;
                                offsetglider(ijk).selcrit(issel).toldist = min(distspace{ijk});
                                offsetglider(ijk).selcrit(issel).timerange = [dmoorini3 dmoorfff3];                           
                        end
                        meanmoorlvl = nanmean(presotherobs2(:,indo),2);
                        stdmoorlvl  = nanstd(presotherobs2(:,indo),0,2);
                        iglpres = [];
                        indotherpres = [];
                        for ipp=1:length(meanmoorlvl)
                            if abs(glPP(near(glPP,meanmoorlvl(ipp))) - meanmoorlvl(ipp))<2
                                iglpres = [iglpres near(glPP,meanmoorlvl(ipp))];
                                indotherpres = [indotherpres ipp];
                            end
                        end

                        % use of nanmean and nanstd in case no data for some
                        % specific profiles
                        glmeantpot = nanmean(glTP(iglpres,indg),2);
                        glstdtpot  = nanstd(glTP(iglpres,indg),0,2);
                        glmeansal  = nanmean(glSS(iglpres,indg),2);
                        glstdsal   = nanstd(glSS(iglpres,indg),0,2);
                        omeantpot  = nanmean(ptempotherobs2(indotherpres,indo),2);
                        ostdtpot  = nanstd(ptempotherobs2(indotherpres,indo),0,2);
                        omeansal  = nanmean(salotherobs2(indotherpres,indo),2);
                        ostdsal   = nanstd(salotherobs2(indotherpres,indo),0,2);

                        i700 = find(meanmoorlvl(indotherpres) > 700); %700m due to RTEB1

                        offsetglider(ijk).selcrit(issel).meanpreslvl = meanmoorlvl(indotherpres);
                        offsetglider(ijk).selcrit(issel).stdpreslvl = stdmoorlvl(indotherpres);
                        offsetglider(ijk).selcrit(issel).glider.tpot.mean = glmeantpot;
                        offsetglider(ijk).selcrit(issel).glider.tpot.std = glstdtpot;
                        offsetglider(ijk).selcrit(issel).glider.tpot.numpf = sum(~isnan(glTP(iglpres,indg)),2);                 
                        offsetglider(ijk).selcrit(issel).glider.sal.mean = glmeansal;
                        offsetglider(ijk).selcrit(issel).glider.sal.std = glstdsal;  
                        offsetglider(ijk).selcrit(issel).glider.sal.numpf = sum(~isnan(glSS(iglpres,indg)),2);                  
                        offsetglider(ijk).selcrit(issel).otherobs.tpot.mean = omeantpot;
                        offsetglider(ijk).selcrit(issel).otherobs.tpot.std = ostdtpot; 
                        offsetglider(ijk).selcrit(issel).otherobs.tpot.numpf = sum(~isnan(ptempotherobs2(indotherpres,indo)),2);                   
                        offsetglider(ijk).selcrit(issel).otherobs.sal.mean = omeansal;
                        offsetglider(ijk).selcrit(issel).otherobs.sal.std = ostdsal;   
                        offsetglider(ijk).selcrit(issel).otherobs.sal.numpf = sum(~isnan(salotherobs2(indotherpres,indo)),2);                   


                        % difference
                        offsetglider(ijk).selcrit(issel).difftpot = glmeantpot -omeantpot ; 
                        offsetglider(ijk).selcrit(issel).difftpotsumstd = sqrt(glstdtpot.^2 + ostdtpot.^2) ; 
                        offsetglider(ijk).selcrit(issel).diffsal = glmeansal -omeansal ;
                        offsetglider(ijk).selcrit(issel).diffsalsumstd = sqrt(glstdsal.^2 + ostdsal.^2) ; ;               

                        offsetglider(ijk).selcrit(issel).tpotdiffmean = nanmean(glmeantpot -omeantpot);
                        offsetglider(ijk).selcrit(issel).tpotdiffstd  = [];%nanstd(glmeantpot -omeantpot,0);
                        offsetglider(ijk).selcrit(issel).tpotdiffmeanstd = sqrt(nanmean(offsetglider(ijk).selcrit(issel).difftpotsumstd.^2));    
                        offsetglider(ijk).selcrit(issel).saldiffmean  = nanmean(glmeansal -omeansal);
                        offsetglider(ijk).selcrit(issel).saldiffstd   = [];%nanstd(glmeansal -omeansal,0);    
                        offsetglider(ijk).selcrit(issel).saldiffmeanstd = [];%sqrt(nanmean(offsetglider(ijk).selcrit(issel).diffsalsumstd.^2));

                        offsetglider(ijk).selcrit(issel).tpotdiffmeanbelow700    = nanmean(glmeantpot(i700) -omeantpot(i700));
                        offsetglider(ijk).selcrit(issel).tpotdiffstdbelow700     = [];%nanstd(glmeantpot(i700) -omeantpot(i700),0);
                        offsetglider(ijk).selcrit(issel).tpotdiffmeanstdbelow700 = sqrt(nanmean(offsetglider(ijk).selcrit(issel).difftpotsumstd(i700).^2));
                        offsetglider(ijk).selcrit(issel).saldiffmeanbelow700     = nanmean(glmeansal(i700) -omeansal(i700));
                        offsetglider(ijk).selcrit(issel).saldiffstdbelow700      = [];%nanstd(glmeansal(i700) -omeansal(i700),0);      
                        offsetglider(ijk).selcrit(issel).saldiffmeanstdbelow700  =  sqrt(nanmean(offsetglider(ijk).selcrit(issel).diffsalsumstd(i700).^2));             
                                 end

                %-----------------------------------------------------------------------------------


                lim_ax_sal = [min([min(glSS(:,indg)) min(salotherobs2(:,indo))]) ...
                              max([max(glSS(:,indg)) max(salotherobs2(:,indo))])] ; 
                lim_ax_ptemp = [min([min(glTP(:,indg)) min(ptempotherobs2(:,indo))]) ...
                              max([max(glTP(:,indg)) max(ptempotherobs2(:,indo))])] ; 

                if ~isempty(indglider2) 
                    titlestr = {[otherobs.name ' (+-) - ' missioninfo.gnickname ' (-) comparisons for ' ...
                            datestr(dmoorini,'dd/mm/yy') ' - ' datestr(dmoorfff,'dd/mm/yy') ' (grey lines and crosses), '], ...
                            [datestr(dmoorini2,'dd/mm/yy') ' - ' datestr(dmoorfff2,'dd/mm/yy') ' (color lines and crosses)' ]};
                else
                    titlestr = {[otherobs.name ' (+-) - ' missioninfo.gnickname ' (-) comparisons for ' ...
                            datestr(dmoorini,'dd/mm/yy') ' - ' datestr(dmoorfff,'dd/mm/yy') ' (grey lines and crosses)']};            
                end
                %---------------------------------------------------------------------
                % Plot Profils
                figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
                set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
                title(titlestr)
                subplot(1,2,1)
                hold on
                plot(glTP(:,indglider),glPP,'color',[1 1 1]*0.65)
                plot(ptempotherobs2(:,indother),presotherobs2(:,indother),'+','color',[1 1 1]*0.45)
                if ~isempty(indglider2)
                    plot(glTP(:,indglider2),glPP,'b')    
                    plot(glTP(:,indglider3),glPP,'color','k','linewidth',2)      
                    plot(ptempotherobs2(:,indother2),presotherobs2(:,indother2),'-.+r')  
                    plot(ptempotherobs2(:,indother3),presotherobs2(:,indother3),'--xk','linewidth',2)    
                end
                set(gca,'ydir','reverse')
                ylabel('Pressure')
                xlabel('Pot. Temp')

                subplot(1,2,2)
                hold on
                plot(glSS(:,indglider),glPP,'color',[1 1 1]*0.65)
                plot(salotherobs2(:,indother),presotherobs2(:,indother),'+','color',[1 1 1]*0.45)    
                if ~isempty(indglider2)  
                    plot(glSS(:,indglider2),glPP,'b')    
                    plot(glSS(:,indglider3),glPP,'color','k','linewidth',2)      
                    plot(salotherobs2(:,indother2),presotherobs2(:,indother2),'-.+r')  
                    plot(salotherobs2(:,indother3),presotherobs2(:,indother3),'--xk','linewidth',2)    
                end
                ylabel('Pressure')
                xlabel('Salinity')   
                set(gca,'ydir','reverse')

                annotation('textbox', [0 0.9 1 0.1], ...
                'String', titlestr, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center')

                set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,   

                set(gcf,'render','painters')
                print(gcf,'-dpng',[dirplot filesep figsuf  '_' figparam.crossing.profilsplot.name '_' num2str(ijk)]);
                if figparam.epsgraph == 1
                    print(gcf,'-depsc2',[dirplot filesep figsuf  '_' figparam.crossing.profilsplot.name '_' num2str(ijk)]);    
                end
                %---------------------------------------------------------------------
                % Plot TS    
                figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
                set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
                hold on
                contour(figparam.contTS.ssc,figparam.contTS.ttc,...
                          figparam.contTS.ddc,figparam.contTS.vcont,'k','linewidth',figparam.contTS.linewidth);       
                plot(glSS(:,indglider),glTP(:,indglider),'color',[1 1 1]*0.65)
                plot(salotherobs2(:,indother),ptempotherobs2(:,indother),'+','color',[1 1 1]*0.45)    
                if ~isempty(indglider2)
                    plot(glSS(:,indglider2),glTP(:,indglider2),'b')    
                    plot(glSS(:,indglider3),glTP(:,indglider3),'color','k','linewidth',2)      
                    plot(salotherobs2(:,indother2),ptempotherobs2(:,indother2),'-.+r')  
                    plot(salotherobs2(:,indother3),ptempotherobs2(:,indother3),'--xk','linewidth',2)    
                end       
                ylabel('Pot. Temperature')
                xlabel('Salinity')
                set(gca,'xlim',lim_ax_sal,'ylim',lim_ax_ptemp)

                title(titlestr);

                set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%, 
                set(gcf,'render','painters')
                print(gcf,'-dpng',[dirplot filesep figsuf  '_' figparam.crossing.TSplot.name '_' num2str(ijk)]);
                if figparam.epsgraph == 1
                    print(gcf,'-depsc2',[dirplot filesep figsuf  '_' figparam.crossing.TSplot.name '_' num2str(ijk)]);    
                end
            end
            %=====================================================================
            % Print statistics file
            %--------------------------------------------------------------------
            % 1. calculate difference between glider data and the other dataset at the
            % level of measurement of the other data set and the mean difference of
            % the whole profile
            outputdir = [missioninfo.step1.procdataDir filesep 'calib' filesep missioninfo.glmission]
            if exist(outputdir)~=7
                mkdir(outputdir)
            end
            statsfile = figsuf;
            for ijk=1:length(offsetglider)
                for ilm=1:length(offsetglider(ijk).selcrit)
                    offsetglider(ijk).glidername = missioninfo.gnickname;
                    offsetglider(ijk).otherdatasettype = otherobs.type;
                    offsetglider(ijk).otherdatasetname = otherobs.name;
                    selectparamstr = [ num2str(offsetglider(ijk).selcrit(ilm).toltime,'%4.0f') 'days_' num2str(offsetglider(ijk).selcrit(ilm).toldist,'%4.0f') 'km'];
                    %-------------------------------------------------------------------------
                    fidlog    = fopen([outputdir filesep statsfile '_' num2str(ijk) '_' selectparamstr '.log'],'w');
                    fprintf(fidlog,['Glider Name                  : ' missioninfo.gnickname '\n']);
                    fprintf(fidlog,['Glider Mission               : ' missioninfo.glmission '\n']);
                    fprintf(fidlog,['Other Dataset Type           : ' offsetglider(ijk).otherdatasettype '\n']);
                    fprintf(fidlog,['Other Dataset Name           : ' offsetglider(ijk).otherdatasetname '\n']);
                    fprintf(fidlog,['Time Window for comparison   : ' num2str(offsetglider(ijk).selcrit(ilm).toltime,'%4.0f') ' (days) \n']);     
                    fprintf(fidlog,['Max. Distance for comparison : ' num2str(offsetglider(ijk).selcrit(ilm).toldist,'%4.0f') ' (km) \n']);                 
                    fprintf(fidlog,'\n');  
                    if offsetglider(ijk).selcrit(ilm).toltime ~=0 % 
                        fprintf(fidlog,[offsetglider(ijk).otherdatasettype ' mean depth;' offsetglider(ijk).otherdatasettype ' std depth;Mean Pot. Temp. Difference; Std Pot. Temp. Difference;' ...
                                    'Mean Salinity Difference; Std Salinity Difference \n']);
                        for ipp=1:length( offsetglider(ijk).selcrit(ilm).meanpreslvl)
                               fprintf(fidlog,'%4.3f;%4.3f;%4.3f;%4.3f;%4.3f;%4.3f\n',offsetglider(ijk).selcrit(ilm).meanpreslvl(ipp),...
                                   offsetglider(ijk).selcrit(ilm).stdpreslvl(ipp),offsetglider(ijk).selcrit(ilm).difftpot(ipp),...
                                   offsetglider(ijk).selcrit(ilm).difftpotsumstd(ipp), offsetglider(ijk).selcrit(ilm).diffsal(ipp),...
                                   offsetglider(ijk).selcrit(ilm).diffsalsumstd(ipp));
                        end

                    else
                         fprintf(fidlog,[offsetglider(ijk).otherdatasettype ' depth;Pot. Temp. Difference;' ...
                                    'Salinity Difference\n']);               
                         for ipp=1:length( offsetglider(ijk).selcrit(ilm).meanpreslvl)
                               fprintf(fidlog,'%4.3f;%4.3f;%4.3f\n',offsetglider(ijk).selcrit(ilm).meanpreslvl(ipp),...
                                   offsetglider(ijk).selcrit(ilm).difftpot(ipp),offsetglider(ijk).selcrit(ilm).diffsal(ipp));
                         end               

                    end

                    fprintf(fidlog,'\n');   
                    fprintf(fidlog,'Depth levels Range;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;Mean Sal. Diff.; Std Sal. Diff.\n');   
                    fprintf(fidlog,'%s;%4.3f;%4.3f;%4.3f;%4.3f\n','All  Depths',...
                                   offsetglider(ijk).selcrit(ilm).tpotdiffmean,offsetglider(ijk).selcrit(ilm).tpotdiffstd,...
                                   offsetglider(ijk).selcrit(ilm).saldiffmean,offsetglider(ijk).selcrit(ilm).saldiffstd); 
                    fprintf(fidlog,'\n');   
                    fprintf(fidlog,'Depth levels Range;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;Mean Sal. Diff.; Std Sal. Diff.\n');    
                    fprintf(fidlog,'%s;%4.3f;%4.3f;%4.3f;%4.3f\n','700-1000m',...
                                   offsetglider(ijk).selcrit(ilm).tpotdiffmeanbelow700,offsetglider(ijk).selcrit(ilm).tpotdiffstdbelow700,...
                                   offsetglider(ijk).selcrit(ilm).saldiffmeanbelow700,offsetglider(ijk).selcrit(ilm).saldiffstdbelow700); 

                    fclose(fidlog);                

                end
            end
            % end

            % save matlab version of the calib
            save([outputdir filesep statsfile],'offsetglider')

            %-------------------------------------------------------------------------
          end      % plot glider-otherobs salinity-salinity
     else

        disp(' ')
        disp(['No crossing points for ' missioninfo.gnickname ' and ' otherobs.name '. Check time limits of the 2 datasets'])
        disp(' ')
        offsetglider = [];
    end   

    if otherobs.UVspeedcomparison == 1
    %-------------------------------------------------------------------------
    % plot glider-otherobs speed

    end         
        
   
        
        
    otherwise
        disp('Dataset not mooring or CTD cruise, precise the way to handle the lat, lon and depth of the other data set')
end
    

    
    

    
