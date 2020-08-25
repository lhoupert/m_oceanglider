function plot_mission_data_profils(data,paramplot,figparam,refDB)
% function plot_mission_data_profils(data,paramplot,figparam,refDB)
%
% refDB structure has to have the same field name than data ('TP','SS','PP') 
% or if no reference data used, just set refDB=[];
%
% created by L. Houpert (houpertloic@gmail.com), 05/02/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox


ip1 = data.ip1;
ip2 = data.ip2;

[mm nn] = size(data.x);
xdata   = data.x(:,ip1:ip2) ;

ydata0  = data.y;
if length(data.y(:)) ~= length(data.x(:)) % for case when pressure is given as a vector and not a matrix
    ydata0 = repmat(data.y(:),1,nn);
end
ydata   = ydata0(:,ip1:ip2) ;

date1 = data.date1;
date2 = data.date2;

gldivephase = data.gldivephase(ip1:ip2);

varlongname  = paramplot.varlongname;       
figrep       = figparam.mainplot.dir ;
figsuf       = figparam.mainplot.name ;
figsecnb     = ['sec' num2str(data.sectionnb)];
figtype      = figparam.profil.figname ;
figext       = paramplot.plotext;


lim_ax_x = [min(min(xdata(:,1:end)))-0.1 max(max(xdata(:,1:end)))+0.1];
lim_ax_y = [min(min(ydata(:,1:end)))-0.1 max(max(ydata(:,1:end)))+0.1];



figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
hold on

if ~isempty(refDB) && (strcmp(data.xvarname,'TP') || strcmp(data.xvarname,'SS')) 
      plot(refDB.(data.xvarname),refDB.(data.yvarname),'color',figparam.color.greycol{2});     
end

if (strcmp(data.yvarname,'TP') && strcmp(data.xvarname,'SS')) 
    contour(figparam.contTS.ssc,figparam.contTS.ttc,...
        figparam.contTS.ddc,figparam.contTS.vcont,'k','linewidth',figparam.contTS.linewidth);
end

plot(xdata(:,gldivephase==1),ydata(:,gldivephase==1),'b.','markersize',figparam.mksize1); 
plot(xdata(:,gldivephase==4),ydata(:,gldivephase==4),'m.','markersize',figparam.mksize1);   
set(gca,'xlim',lim_ax_x ,'ylim',lim_ax_y)
if (strcmp(data.yvarname,'PP')) | (strcmp(data.yvarname,'pres')) 
    set(gca,'Ydir','reverse')
end
title(['Profiles of ' varlongname ' (b: dive, m: climb)  ---  ' data.gname '  '...
    datestr(date1,'dd/mm/yyyy HH:MM') ' - ' datestr(date2,'dd/mm/yyyy HH:MM') ],'fontsize',figparam.fstitle)         


set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
print(gcf,'-dpng',[figrep filesep figsuf '_' figsecnb '_' figtype  '_' figext ]);
close

end
    
