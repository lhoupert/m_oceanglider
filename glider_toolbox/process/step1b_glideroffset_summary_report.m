function step1b_glideroffset_summary_report(offsetglider,missioninfo,tol)
%
%=========================================================================%
% Function to writte a summary of all the apparent offsets for a entire glider mission,
% calculated using the output of the function step1b_glider_offset_from_obs
%
% L. Houpert, SAMS 02/06/2017
%=========================================================================%
%
%  step1b_plot_glideroffset_mission(offsetglider,missioninfo,figparam)
%  
%  - offsetglider: a [1xN] structure where N correspond to the number of  
%  different dataset for which the glider comparisons are made (each
%  dataset comparison is produced by the fonction step1b_glider_offset_from_obs)
%
%  - missioninfo: [1x1] structure with mission details (defined in
%             users_param/loadosnapmissionparam.m)
%  
%  - figparam: [1x1] structure containing all the general and specific figure and
%           graph parameters (defined in users_param/graphparamgliderproc.m)
%
% created by L. Houpert (houpertloic@gmail.com), 03/05/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 
% Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)

outputdir = [missioninfo.step1.procdataDir filesep 'calib' filesep missioninfo.glmission];
fname   = ['Summary_Offsets_' missioninfo.glmission ];
statsfile = fname;    

for ilm=1:3
           
    switch(ilm)
        case 1
            selectparamstr = [ num2str(tol.time,'%4.0f') 'days_' num2str(tol.dist,'%4.0f') 'km'];              
        case 2
            selectparamstr = [ num2str(tol.time2,'%4.0f') 'days_' num2str(tol.dist2,'%4.0f') 'km']; 
        case 3    
            selectparamstr = 'closest';            
    end

    fidlog(ilm)    = fopen([outputdir filesep statsfile '_' selectparamstr '.log'],'w');  
    fprintf(fidlog(ilm),['Glider Name                  : ' missioninfo.gnickname '\n']);                      
    fprintf(fidlog(ilm),['Glider Mission               : ' missioninfo.glmission '\n\n']);      
    
    for ijk = 1:length(offsetglider)

        % the whole profile
        outputdir = [missioninfo.step1.procdataDir filesep 'calib' filesep missioninfo.glmission];
        if exist(outputdir)~=7
             mkdir(outputdir)
        end
  
if length(offsetglider(ijk).selcrit)<ilm, continue;end

        switch(offsetglider(ijk).selcrit(ilm).toltime)
            case tol.time
                fidlogtowrite = fidlog(1);              
            case tol.time2
                fidlogtowrite = fidlog(2);  
            case 0    
                fidlogtowrite = fidlog(3);           
        end

        
        numglpfs = nanmean(offsetglider(ijk).selcrit(ilm).glider.sal.numpf);
        numotherpfs = nanmean(offsetglider(ijk).selcrit(ilm).otherobs.sal.numpf);    
        comptimerange = offsetglider(ijk).selcrit(ilm).timerange;
        timecp1 = datestr(comptimerange(1),'dd/mm/yy HH:MM');
        timecp2 = datestr(comptimerange(2),'dd/mm/yy HH:MM');
        
        fprintf(fidlogtowrite,['Other Dataset Type           : ' offsetglider(ijk).otherdatasettype '\n']);
        fprintf(fidlogtowrite,['Other Dataset Name           : ' offsetglider(ijk).otherdatasetname '\n']);
        fprintf(fidlogtowrite,['Nber of other dataset pfs.   : ' num2str(numotherpfs,'%4.0f') '\n']);   
        fprintf(fidlogtowrite,['Nber of glider profiles used : ' num2str(numglpfs,'%4.0f') '\n']);   
        fprintf(fidlogtowrite,['Time range                   : ' timecp1 ' - ' timecp2 '\n']);              
        fprintf(fidlogtowrite,['Time Window for comparison   : ' num2str(offsetglider(ijk).selcrit(ilm).toltime,'%4.0f') ' (days) \n']);     
        fprintf(fidlogtowrite,['Max. Distance for comparison : ' num2str(offsetglider(ijk).selcrit(ilm).toldist,'%4.0f') ' (km) \n']);
        %fprintf(fidlogtowrite,'\n');   
                    
        %=====================================================================
        % Print statistics file     
        %--------------------------------------------------------------------
        % 1. calculate difference between glider data and the other dataset at the
        % level of measurement of the other data set and the mean difference of
        %-------------------------------------------------------------------------

 
                        if offsetglider(ijk).selcrit(ilm).toltime ~=0 % 
                            fprintf(fidlogtowrite,[offsetglider(ijk).otherdatasettype ' mean depth;' offsetglider(ijk).otherdatasettype ' std depth; Std Mean Pot. Temp ;Mean Pot. Temp. Difference; Std Pot. Temp. Difference;' ...
                                        'Std Mean Salinity; Mean Salinity Difference; Std Salinity Difference \n']);
                            fprintf(fidlogtowrite,[' ; ; (other dataset); (glider - otherdataset); sqrt(glstdtpot^2 + ostdtpot^2) ;' ...
                                        '(other dataset); (glider - otherdataset); sqrt(glstdtpot^2 + ostdtpot^2)  \n']);   
    
                            zmaxsel = find(~isnan(offsetglider(ijk).selcrit(ilm).difftpot) & ...
                                offsetglider(ijk).selcrit(ilm).difftpotsumstd~=0 ,1,'last');
                            
                            for ipp = zmaxsel-1:zmaxsel
                                   fprintf(fidlogtowrite,'%4.3f;%4.3f;%4.3f;%4.3f;%4.3f;%4.3f;%4.3f;%4.3f\n',offsetglider(ijk).selcrit(ilm).meanpreslvl(ipp),...
                                       offsetglider(ijk).selcrit(ilm).stdpreslvl(ipp), offsetglider(ijk).selcrit(ilm).otherobs.tpot.std(ipp),...
                                       offsetglider(ijk).selcrit(ilm).difftpot(ipp), offsetglider(ijk).selcrit(ilm).difftpotsumstd(ipp), ...
                                       offsetglider(ijk).selcrit(ilm).otherobs.sal.std(ipp),offsetglider(ijk).selcrit(ilm).diffsal(ipp),...
                                       offsetglider(ijk).selcrit(ilm).diffsalsumstd(ipp));
                            end
                            
 
                            %fprintf(fidlogtowrite,'\n');   
                            %fprintf(fidlogtowrite,'Depth levels Range;;;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;;Mean Sal. Diff.; Std Sal. Diff.\n');   
                            fprintf(fidlogtowrite,'%s;;;%4.3f;%4.3f;;%4.3f;%4.3f\n',...
                                            'All  Depths', offsetglider(ijk).selcrit(ilm).tpotdiffmean,offsetglider(ijk).selcrit(ilm).tpotdiffstd,...
                                           offsetglider(ijk).selcrit(ilm).saldiffmean,offsetglider(ijk).selcrit(ilm).saldiffstd); 
                            %fprintf(fidlogtowrite,'\n');   
                            %fprintf(fidlogtowrite,'Depth levels Range;;;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;;Mean Sal. Diff.; Std Sal. Diff.\n');    
                            fprintf(fidlogtowrite,'%s;;;%4.3f;%4.3f;;%4.3f;%4.3f\n',...
                                            '700-1000m',offsetglider(ijk).selcrit(ilm).tpotdiffmeanbelow700,offsetglider(ijk).selcrit(ilm).tpotdiffstdbelow700,...
                                           offsetglider(ijk).selcrit(ilm).saldiffmeanbelow700,offsetglider(ijk).selcrit(ilm).saldiffstdbelow700); 

                           

                        else
                             fprintf(fidlogtowrite,[offsetglider(ijk).otherdatasettype ' depth;Pot. Temp. Difference;' ...
                                        'Salinity Difference\n']);   
                             zmaxsel = find(~isnan(offsetglider(ijk).selcrit(ilm).difftpot),1,'last');                                    
                             for ipp = zmaxsel-1:zmaxsel
                                   fprintf(fidlogtowrite,'%4.3f;%4.3f;%4.3f\n',offsetglider(ijk).selcrit(ilm).meanpreslvl(ipp),...
                                       offsetglider(ijk).selcrit(ilm).difftpot(ipp),offsetglider(ijk).selcrit(ilm).diffsal(ipp));
                             end     
                             
                            %fprintf(fidlogtowrite,'\n');   
                            %fprintf(fidlogtowrite,'Depth levels Range;;;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;;Mean Sal. Diff.; Std Sal. Diff.\n');   
                            fprintf(fidlogtowrite,'%s;%4.3f;%4.3f\n',...
                                            'All  Depths', offsetglider(ijk).selcrit(ilm).tpotdiffmean,...
                                           offsetglider(ijk).selcrit(ilm).saldiffmean); 
                            %fprintf(fidlogtowrite,'\n');   
                            %fprintf(fidlogtowrite,'Depth levels Range;;;Mean Pot. Temp. Diff.; Std Pot. Temp. Diff.;;Mean Sal. Diff.; Std Sal. Diff.\n');    
                            fprintf(fidlogtowrite,'%s;%4.3f;%4.3f\n',...
                                            '700-1000m',offsetglider(ijk).selcrit(ilm).tpotdiffmeanbelow700,...
                                           offsetglider(ijk).selcrit(ilm).saldiffmeanbelow700); 

                        end


                        
                fprintf(fidlogtowrite,'\n\n');         
    end
                         
end
                
fclose(fidlog(1));
fclose(fidlog(2));
fclose(fidlog(3));

