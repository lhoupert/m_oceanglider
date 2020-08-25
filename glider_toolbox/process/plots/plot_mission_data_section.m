function plot_mission_data_section(data,usebathyfile,paramplot,figparam)
% function plot_mission_data_section(data,usebathyfile,paramplot,figparam)
%
% created by L. Houpert (houpertloic@gmail.com), 05/02/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox

id1 = data.id1;
id2 = data.id2;

if usebathyfile == 1;
  zbathydive =   data.zbathy(id1:id2) ; 
end
timedive = data.time(id1:id2);
altimeter_data = data.zalti(id1:id2);

i1 = data.i1;
i2 = data.i2;

xdata   = data.x(i1:i2) ;
ydata   = data.y(i1:i2) ;  
cdata   = data.c(i1:i2) ;
if isempty(data.qc)
    useqc = 0;
else
    useqc = 1;
    qc_flag = data.qc(i1:i2);
end


date1 = data.date1;
date2 = data.date2;

datestr(date1)
datestr(date2)

varlongname  = paramplot.varlongname;       
figrep       = figparam.mainplot.dir ;
figsuf       = figparam.mainplot.name ;
figtype      = figparam.scattersection.figname;
figsecnb     = ['sec' num2str(data.sectionnb)];
figext       = paramplot.plotext;
nbdatetik    = figparam.nbdatetik;
isc           = figparam.scattersection.scatter_int;

idataok      = find(~isnan(cdata));       
idatabad = find(isnan(cdata));

ibadqc{1} = idatabad(qc_flag(idatabad) == 3); symbqc{1} = 'o'; legstr{1}='flag 3';
ibadqc{2} = idatabad(qc_flag(idatabad) == 4); symbqc{2} = '+'; legstr{2}='flag 4';
ibadqc{3} = idatabad(qc_flag(idatabad) == 8); symbqc{3} = 'd'; legstr{3}='flag 8'; % Seabird basestation flag for interpolated value
ibadqc{4} = idatabad(qc_flag(idatabad) == 41); symbqc{4} = '^'; legstr{4}='flag 41';% spike     
ibadqc{5} = idatabad(qc_flag(idatabad) == 42); symbqc{5} = 'h'; legstr{5}='flag 42'; % gradient test failed
ibadqc{6} = idatabad(qc_flag(idatabad) == 43); symbqc{6} = 'v'; legstr{6}='flag 43'; % density inversion
ibadqc{7} = idatabad(qc_flag(idatabad) == 44); symbqc{7} = '>'; legstr{7}='flag 44'; %  inferior to fix threshold
ibadqc{8} = idatabad(qc_flag(idatabad) == 40); symbqc{8} = 'x'; legstr{8}='flag 40'; %  manually flagged bad    

indlegok =[];

figg=figure('visible',figparam.vsblfig,'position',figparam.scrsz);
set(figg,'PaperUnits','centimeters','PaperOrientation','portrait',... 
                    'Paperposition',figparam.figpos)  
hold on
if useqc == 1

    pl2=scatter(xdata(idataok(1:isc:end)),ydata(idataok(1:isc:end)),figparam.scattersection.scatsize,cdata(idataok(1:isc:end)),'filled');           
    for ipp=1:length(ibadqc)

          pl3{ipp} =plot(xdata(ibadqc{ipp}) ,ydata(ibadqc{ipp}),symbqc{ipp},...
                'markeredgecolor',figparam.scattersection.mkedgecolor,...
                'markerfacecolor',figparam.scattersection.mkfacecolor,...
                'linewidth',figparam.scattersection.mklinewidth,...
                'markersize',figparam.scattersection.mksize);      
        if ~isempty(pl3{ipp}) ; 
            indlegok = [indlegok ipp];
        end
    end   

    indqfok = unique(indlegok);
    legcol{1}='flag 0, 1 or 2';
    legend([pl2 pl3{indqfok}],[legcol legstr(indqfok)],'location','west')              

    for ipp=1:length(legstr)
        set(pl3{ipp},'zdata',get(pl3{ipp},'xdata')*0+100)            
    end

else % no raw data
       
    scatter(xdata(1:isc:end),ydata(1:isc:end),figparam.scattersection.scatsize,cdata(1:isc:end),'filled')                    
    
end


if usebathyfile==1;
    plot(timedive,zbathydive,'color',figparam.scattersection.colbathy ,'linewidth',figparam.scattersection.linewidthbathy)
end
plot(timedive,altimeter_data,'color',figparam.scattersection.colalti,'linewidth',figparam.scattersection.linewidthalti)
set(gca,'xlim',[date1 date2],'xtick',date1:round((date2-date1)/nbdatetik):date2,'ylim',figparam.yaxispreslimit)

datetick('x',19,'keeplimits','keepticks')
if (strcmp(data.yvarname,'pres'))
    set(gca,'Ydir','reverse')
end
title([varlongname ' section  ---  ' data.gname '  '...
    datestr(date1,'dd/mm/yyyy HH:MM') ' - ' datestr(date2,'dd/mm/yyyy HH:MM') ],'fontsize',figparam.fstitle)         

colormap(paramplot.colmap)
caxis([paramplot.clim1 paramplot.clim2])      
colorbar


set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
print(gcf,'-dpng',[figrep filesep figsuf '_' figsecnb '_' figtype  '_' figext ]);
close



    
end
