function plot_mission_data_time_series_qc(data,paramplot,figparam)
% function plot_mission_data_time_series_qc(data,paramplot,figparam)
%     
% created by L. Houpert (houpertloic@gmail.com), 05/02/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox

i1 = data.i1;
i2 = data.i2;

xdata   = data.x(i1:i2) ;
ydata   = data.y(i1:i2) ;  
qc_flag = data.qc(i1:i2);
if isempty(data.y00)
    useqc = 0;
else
    useqc = 1;
    ydata00 = data.y00(i1:i2) ;
end
date1 = data.date1;
date2 = data.date2;

datestr(date1)
datestr(date2)

varlongname  = paramplot.varlongname;       
figrep       = figparam.mainplot.dir ;
figsuf       = figparam.mainplot.name ;
figtype      = figparam.timeseries.figname;
figsecnb     = ['sec' num2str(data.sectionnb)];
figext       = paramplot.plotext;
nbdatetik    = figparam.nbdatetik;


idatabad  = find(isnan(ydata));
idataok  = find(~isnan(ydata));
inotnan   = intersect(idatabad,find(qc_flag~=9));          

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

    pl1=plot(xdata(inotnan),ydata00(inotnan),figparam.linestyle{2},'color',figparam.color.greycol{2},'markersize',figparam.mksize2);        
    pl2=plot(xdata(idataok) ,ydata(idataok),figparam.linestyle{3},'color',figparam.color.coldef{1},'markersize',figparam.mksize1);               
    for ipp=1:length(ibadqc)

        pl3{ipp} =plot(xdata(ibadqc{ipp}) ,ydata00(ibadqc{ipp}),symbqc{ipp},'color',figparam.color.greycol{3},'markersize',figparam.mksize2);      
        if ~isempty(pl3{ipp}) ; 
            indlegok = [indlegok ipp];
        end
    end   

    indqfok = unique(indlegok);
    legcol{1}='flag 0, 1 or 2';
    legend([pl2 pl3{indqfok}],[legcol legstr(indqfok)],'location','west')              

    set(pl1,'zdata',get(pl1,'xdata')*0-1)
    for ipp=1:length(legstr)
        set(pl3{ipp},'zdata',get(pl3{ipp},'xdata')*0+1)            
    end

else % no raw data
       
    pl2=plot(xdata ,ydata,figparam.linestyle{3},'color',figparam.color.coldef{1},'markersize',figparam.mksize1);       
end

set(gca,'xlim',[date1 date2],'xtick',date1:round((date2-date1)/nbdatetik):date2)
datetick('x',19,'keeplimits','keepticks')

title(['Times series of ' varlongname ],'fontsize',figparam.fstitle)         


set(gcf,'position',figparam.pos_eps*100,'paperunits','inches','paperposition',figparam.pos_eps)%,             
print(gcf,'-dpng',[figrep filesep figsuf '_' figsecnb '_' figtype  '_' figext ]);
close


    
end
    
