function plot_mission_surfplus_section(data,usebathyfile,paramplot,figparam)
% function plot_mission_surfplus_section(data,usebathyfile,paramplot,figparam)
% create a surface plot of the glider data and double the data matrix to see for all data points a color patch
%     
% created by L. Houpert (houpertloic@gmail.com), 05/02/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox


id1 = data.id1;
id2 = data.id2;
iok = find(~isnan(nanmean(data.c,1)));
indselect = intersect(id1:id2,iok);

if usebathyfile == 1;
  zbathydive =   data.zbathy(indselect) ; 
end
timedive = data.time(indselect);
altimeter_data = data.zalti(indselect);

xdata   = data.x(:,indselect) ;
ydata   = data.y(:,indselect) ;  
cdata   = data.c(:,indselect) ;
[mmx nnx] = size(xdata)
if mmx==1
    xdata=repmat(xdata,size(ydata(:,1)));
end
if isfield(data,'contourvar') % if the user want to plot additional contour lines on the surface plot
    contdata  = data.contdata(:,indselect) ;
    xcontdata = data.xcontdata(:,indselect)
    ycontdata = data.ycontdata(:,indselect)
end
date1 = data.date1;
date2 = data.date2;

datestr(date1)
datestr(date2)

varlongname  = paramplot.varlongname;       
figrep       = figparam.mainplot.dir ;
figsuf       = figparam.mainplot.name ;
figtype      = figparam.section.figname;
figsecnb     = ['sec' num2str(data.sectionnb) '_surfplus'];
figext       = paramplot.plotext;
nbdatetik    = figparam.nbdatetik;


facctik   = paramplot.facctik;
c1        = paramplot.clim1;
c2        = paramplot.clim2;  
stepctick = paramplot.stepctick;
col=lines;close;cn=col(1:10,:);
if mod(c1,stepctick)==0
    c1tik=c1;
else
   c1tik=c1 + (stepctick-mod(c1,stepctick));
end

if mod(c2,stepctick)==0
   c2tik=c2;
else
   c2tik=c2 - mod(c2,stepctick);
end     
ctick=c1tik:stepctick:c2tik;
contctick=c1:(stepctick/facctik):c2;
clip=[c1-(contctick(2)-contctick(1)) contctick c2+(contctick(2)-contctick(1))];
nlct=length(contctick); 


figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
ax=axes('position',figparam.section.axpos);
colormap(paramplot.colmap)

[xxxx yyyy zzzz]=preparetodblesurf(xdata,ydata,cdata,figparam.dblesurf.ddxx,figparam.dblesurf.ddyy);
surf(xxxx,yyyy,zzzz)     
view([0 90]); shading interp
hold on         

plot(timedive,0*timedive,'.','color','k','linewidth',1)
if usebathyfile==1;
    plot(timedive,zbathydive,'color',figparam.section.colbathy ,'linewidth',figparam.section.linewidthbathy)
end
plot(timedive,altimeter_data,'color',figparam.section.colalti,'linewidth',figparam.section.linewidthalti)
set(gca,'xlim',[date1 date2],'xtick',date1:round((date2-date1)/nbdatetik):date2,'ylim',paramplot.ylim)
datetick('x',19,'keeplimits','keepticks')
caxis([c1 c2])    

if (strcmp(data.yvarname,'pres')) | (strcmp(data.yvarname,'pres_buoyfreq'))
    set(gca,'Ydir','reverse')
end

if isfield(data,'contourvar') & (strcmp(data.contourvar,'pden'))
    [cb0 hb0]=contour(xcontdata,ycontdata,contdata,figparam.contpden.lvl1,'linewidth',figparam.contpden.width);	
    set(hb0,'linestyle',figparam.contpden.style,'color',figparam.contpden.col)
    [cb hb]=contour(xcontdata,ycontdata,contdata,figparam.contpden.lvl2,'linewidth',figparam.contpden.width);	
    set(hb,'linestyle',figparam.contpden.style,'color',figparam.contpden.col)
    clabel(cb,hb,'LabelSpacing',figparam.contpden.labspace,'fontsize',figparam.contpden.fs,'color',figparam.contpden.col)
end

title([ data.gname ' - ' data.gmission ' (' datestr(date1,'dd/mm/yy') ...
        '-' datestr(date2,'dd/mm/yy') ')'])
    
clipref=clip;
dlen=0.25;
colo=clrbr2b(gca,clipref,[0.85 0.25 0.02 0.5],paramplot.cbarlabel);    
c1bis=c1-(clipref(3)-clipref(2));
c2bis=c2+(clipref(end-1)-clipref(end-2));
set(colo,'ytick',ctick,'ylim',[c1bis c2bis])
caxis([c1 c2])


set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
print(gcf,'-dpng','-r300',[figrep filesep figsuf '_' figsecnb '_' figtype  '_' figext ]);
close



    
end
