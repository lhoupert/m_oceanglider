function step1b_plot_glideroffset_mission(offsetglider,missioninfo,figparam)
%
%=========================================================================%
% Function to plot all the apparent offsets for a entire glider mission,
% calculated using the output of the function step1b_glider_offset_from_obs
%
% L. Houpert, SAMS 03/05/2016
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

for ivar = 1:2 % loop on the different variable (pot. temp or salinity)

 figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
 set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                'Paperposition',figparam.figpos)  
 hold on
%subplot(2,1,1);hold on 
 for loopplot=1:2
    for isel = 1:length(offsetglider) % loop on the different datasets used
        datasetss = regexp(offsetglider(isel).otherdatasetname,' ','split');
        otherdatasetname = datasetss{1};
        refvarmean = 0;
        refvarstd  = 0;
        refvartime = 0;
        varpfsgl   = [];
        varpfsother= [];
        for icrit = 1:length(offsetglider(isel).selcrit) % loop on the different spatial/temporal criterion
            izz = length(offsetglider(isel).selcrit(icrit).meanpreslvl);
            vardepth = offsetglider(isel).selcrit(icrit).meanpreslvl(izz);     
            %vartime0  = offsetglider(isel).selcrit(icrit).timerange;  vartime  = [mean(vartime0) - (vartime0(2) - vartime0(1))/4 ; mean(vartime0) + (vartime0(2) - vartime0(1))/4];
            vartime = offsetglider(isel).selcrit(icrit).timerange;
            switch ivar
                case 1
                    varmean     = offsetglider(isel).selcrit(icrit).difftpot(izz);
                    varstd      = offsetglider(isel).selcrit(icrit).difftpotsumstd(izz);
                    varpfsgl(icrit)    = offsetglider(isel).selcrit(icrit).glider.tpot.numpf(izz);
                    varpfsother(icrit) = offsetglider(isel).selcrit(icrit).otherobs.tpot.numpf(izz);
                    varstr      = 'PotTemp';
                    ylimdef     = figparam.missionoffset.ylimtpot;
                case 2
                    varmean     = offsetglider(isel).selcrit(icrit).diffsal(izz);
                    varstd      = offsetglider(isel).selcrit(icrit).diffsalsumstd(izz);
                    varpfsgl(icrit)    = offsetglider(isel).selcrit(icrit).glider.sal.numpf(izz);
                    varpfsother(icrit) = offsetglider(isel).selcrit(icrit).otherobs.sal.numpf(izz);
                    varstr      = 'Salinity';
                    ylimdef     = figparam.missionoffset.ylimsal;                   
            end
                       
            switch icrit
                case 1
                    refvarmean = varmean; % for text display
                    refvarstd  = varstd;
                    refvartime = vartime;
                    recfacecol = [1 1 1]*0.6;
                    texcol = 'k';
                    %rectposition = [vartime(1),varmean-varstd,vartime(2)-vartime(1),2*varstd];
                    %rectangle('Position',rectposition,'FaceColor',recfacecol,'EdgeColor',recedgecol)       
                    patchx = [vartime(1) vartime(2) vartime(2) vartime(1)];
                    patchy = [varmean-varstd varmean-varstd varmean+varstd varmean+varstd];
                    patchyb = [varmean varmean varmean varmean];
                    if loopplot == 1 
                        pl1=patch(patchx,patchy,recfacecol);
                        pl1b=patch(patchx,patchyb,'k');                        
                    end
                    toldist1 = offsetglider(isel).selcrit(icrit).toldist;
                    toltime1 = offsetglider(isel).selcrit(icrit).toltime;                     
                case 2
                    recfacecol = 'b';
                    texcol = recfacecol;
                    %rectposition = [vartime(1),varmean-varstd,vartime(2)-vartime(1),2*varstd];
                    %rectangle('Position',rectposition,'FaceColor',recfacecol,'EdgeColor',recedgecol)  
                    patchx = [vartime(1) vartime(2) vartime(2) vartime(1)];
                    patchy = [varmean-varstd varmean-varstd varmean+varstd varmean+varstd];
                    if loopplot == 2 
                        pl2 = patch(patchx,patchy,recfacecol);   
                        pl2b=plot(mean(vartime),varmean,'marker','+','markersize',12,'color','y','linewidth',3);
                    end
                    toldist2 = offsetglider(isel).selcrit(icrit).toldist;      
                    toltime2 = offsetglider(isel).selcrit(icrit).toltime;                          
                case 3 % only the closest profiles
                    recfacecol = 'r';
                    texcol = recfacecol;  
                    if loopplot ==2
                        pl3=plot(vartime(1),varmean,'marker','o','markersize',8,'color','r','linestyle','none','linewidth',2);
                    end
            end
%             if loopplot == 2
%                 ytextdepth = (refvarmean + refvarstd - icrit.*refvarstd/3);
%                 xtextdepth = refvartime(2);
%                 text(xtextdepth,ytextdepth,num2str(vardepth,'%4.0f'),'color',texcol)
%             end
        end   
        if loopplot == 2        
            xtext = offsetglider(isel).selcrit(1).timerange(1);
            ytext1=ylimdef(1) + 4*((refvarmean-refvarstd) - ylimdef(1))/8;
            text(xtext,ytext1, ...
                [otherdatasetname '(' num2str(vardepth,'%4.0f') 'm):'],...
                'horizontalalignment','left')        % refvarmean + 1.1.*refvarstd
            ytext2=ylimdef(1) + 3*((refvarmean-refvarstd) - ylimdef(1))/8;
            text(xtext,ytext2, ...
                [num2str(varpfsgl(1),'%4.0f') ' gl data'],...
                'horizontalalignment','left','color',[1 1 1]*0.6)        % refvarmean + 1.1.*refvarstd            
            if length(offsetglider(isel).selcrit)>1
                ytext3=ylimdef(1) + 2*((refvarmean-refvarstd) - ylimdef(1))/8;
                text(xtext,ytext3, ...
                    [num2str(varpfsgl(2),'%4.0f') ' gl data' ],...
                    'horizontalalignment','left','color','b')   
                ytext4=ylimdef(1) + 1*((refvarmean-refvarstd) - ylimdef(1))/8;
                text(xtext,ytext4, ...
                    ['Gl. at: ' num2str(offsetglider(isel).selcrit(3).toldist,'%4.2f') ...
                    'km ' ],...
                    'horizontalalignment','left','color','r')      
            end
        end
    end
  end % rof loop plot
    legend([pl1,pl2,pl3],['Mean(-k) +- Std,  ' num2str(toldist2,'%3.0f') 'km < dist < ' num2str(toldist1,'%3.0f') 'km ; time < ' num2str(toltime1,'%2.0f') 'days'],...
                    ['Mean(+y) +- Std, dist < ' num2str(toldist2,'%3.0f') 'km ; time < ' num2str(toltime2,'%2.0f') ' days'],...
                    'Closest glider profile','location','northwest')
                        
    xlabel('Time')
    ylabel(['Apparent ' varstr ' Offset Glider'])
    set(gca,'ylim',ylimdef)
    datetick('x','mmm')
    title(['Apparent ' varstr ' Offset for ' offsetglider(isel).glidername ' - ' offsetglider(isel).glidermission])
    
    dirplot  = [missioninfo.step1.plotdir filesep figparam.missionoffset.dirfig  ];    
    if exist(dirplot ,'dir')~=7 ; mkdir(dirplot);end
    figname   = [missioninfo.glmission '_tol_' num2str(toltime1,'%2.0f') 'd_' num2str(toldist1,'%2.0f') 'km_' ...
        num2str(toldist2,'%2.0f') 'km_offset_' varstr ];
    
    set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%, 
    set(gcf,'render','painters')
    print(gcf,'-dpng',[dirplot filesep figname]);
    
    if figparam.pdfgraph == 1
        set(gcf,'paperorientation','landscape')
        print(gcf,'-dpdf',[dirplot filesep figname]);        
    end
end % rof ivar

end
