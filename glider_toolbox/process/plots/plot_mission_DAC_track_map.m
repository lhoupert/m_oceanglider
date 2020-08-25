function plot_mission_DAC_track_map(data,paramplot,figparam,step1opt)
% function plot_mission_DAC_track_map(data,paramplot,figparam,step1opt)
%
% created by L. Houpert (houpertloic@gmail.com), 05/02/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox


if  exist(step1opt.bathydir,'dir')==7
    usebathy =1;
end

ind1 = data.idac1;
ind2 = data.idac2;
itt  = ind1:ind2;
lon   = data.lon ;
lat   = data.lat ;  
UU    = data.UU ;  
VV    = data.VV ;  
cdata = data.cdata ;  
date1 = data.date1;
date2 = data.date2;

datestr(date1)
datestr(date2)

varlongname  = paramplot.varlongname;       
figrep       = figparam.mainplot.dir ;
figsuf       = figparam.mainplot.name ;
figsecnb     = ['sec' num2str(data.sectionnb)];
figext       = paramplot.plotext;
nbdatetik    = figparam.nbdatetik;
scale        = figparam.mapDAC.scale;
stikstep     = figparam.mapDAC.stickstep;
figtype      = figparam.mapDAC.figname;

x0 = 50; % km 
y0 = 50; % km
lonmin=min(data.lon);
lonmax=max(data.lon);
latmin=min(data.lat);
latmax=max(data.lat);
[platmax,plonmax]=xy2latlon(x0,y0,latmax,lonmax);
[platmin,plonmin]=xy2latlon(-x0,-y0,latmin,lonmin);
limitcmap=[plonmin-0.5  plonmax+0.5 platmin-0.5 platmax+0.5];

figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.mapDAC.figpos )  
hold on

m_proj...
('mercator','longitudes',[plonmin plonmax],'lat',[platmin  platmax]);
%m_coast('patch',[.2 .2 .2],'edgecolor','k');


if usebathy == 1
    [cs ,h]=m_drawbathy('contour',step1opt.bathyname,-3500:500:0,limitcmap,step1opt.bathydir);
end
% clab=clabel(cs,h,'manual');
% set(clab,'color',[0.5 0.5 0.5])
m_gshhs_h('patch',[.5 .5 .5],'edgecolor','k');
m_grid('xtick',10,'tickdir','in','fontsize',12);
m_track(lon(itt),lat(itt),'color','r','linew',1);
hold on
[HP, HT]=m_vec(scale,lon(itt(1:stikstep:end)),lat(itt(1:stikstep:end)),...
      UU(itt(1:stikstep:end)),VV(itt(1:stikstep:end)),cdata(itt(1:stikstep:end)));      
cba=colorbar;
if strcmp(paramplot.matlabvar_c,'timedive') | strcmp(paramplot.matlabvar_c,'DAC_time') | strcmp(paramplot.matlabvar_c,'timeDAC')
    caxis([date1 date2])
    set(cba,'yticklabel',datestr(get(cba,'ytick'),'dd mmm yy') );
else
    caxis([min(cdata) max(cdata)])
    
end
title([ data.gname ' - ' data.gmission ' (' datestr(date1,'dd/mm/yy') ...
        '-' datestr(date2,'dd/mm/yy') ') - ' paramplot.subtitle  ')'])
% Scale.
dlat=(platmax-platmin)/5;
dlon=(plonmax-plonmin)/10;
xpos1 = plonmin + dlon;
ypos1 = platmin + dlat;
[hp,ht]=m_vec(scale,xpos1,ypos1,0.2,0);
[hp2,ht2]=m_vec(scale,xpos1,ypos1,0,0.2);
m_text(xpos1,ypos1-dlat*0.15,'20 cm/s');


set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
print(gcf,'-dpng',[figrep filesep figsuf '_' figsecnb '_' figtype  '_' figext ]);
close


    
end
    
