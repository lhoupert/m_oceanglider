function plot_mission_DAC_time_series(data,figparam)
% plot_mission_DAC_time_series(data,figparam)
%
% created by L. Houpert (houpertloic@gmail.com), 05/02/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox


ind1 = data.idac1;
ind2 = data.idac2;

lon   = data.lon(ind1:ind2) ;
lat   = data.lat(ind1:ind2) ;  
time  = data.time(ind1:ind2);
UU    = data.UU(ind1:ind2);  
VV    = data.VV(ind1:ind2);  
date1 = data.date1;
date2 = data.date2;

datestr(date1)
datestr(date2)

varlongname  = figparam.timeseriesDAC.varlongname;       
figrep       = figparam.mainplot.dir ;
figsuf       = figparam.mainplot.name ;
figtype      = figparam.timeseriesDAC.figname;
figsecnb     = ['sec' num2str(data.sectionnb)];
figext       = figparam.timeseriesDAC.plotext;
nbdatetik    = figparam.nbdatetik;



figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
hold on

subplot(3,1,1)
stickvect(time,figparam.timeseriesDAC.stickvecscale,UU,VV);
title('Stick Vector of the depth Average Current (m.s-1)','fontsize',figparam.fstitle)	  
set(gca,'xlim',[date1 date2],'xtick',date1:round((date2-date1)/nbdatetik):date2)
datetick('x',19,'keeplimits','keepticks')

subplot(3,1,2)
plot(time,UU,'-+','color',figparam.color.coldef{1},'markersize',figparam.mksize1);    
title('Time series of U component (m.s-1)','fontsize',figparam.fstitle)	  
set(gca,'xlim',[date1 date2],'xtick',date1:round((date2-date1)/nbdatetik):date2)
datetick('x',19,'keeplimits','keepticks')   

subplot(3,1,3)
plot(time,VV,'-+','color',figparam.color.coldef{1},'markersize',figparam.mksize1);    
title('Time series of V component (m.s-1)','fontsize',figparam.fstitle)	  
set(gca,'xlim',[date1 date2],'xtick',date1:round((date2-date1)/nbdatetik):date2)
datetick('x',19,'keeplimits','keepticks') 

  


set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
print(gcf,'-dpng',[figrep filesep figsuf '_' figsecnb '_' figtype  '_' figext ]);
close


    
end
    
