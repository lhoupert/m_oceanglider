clear all; close all; limitcmap=[-16 44 28 48];
plonmin=2.5;plonmax=6;platmin=42;platmax=44; 
fs_map=12*2;
%limitcmap=[-0 10 38 44.5];
%plonmin=0;plonmax=10;platmin=38;platmax=44.5; fs_map=20;

minc=-3600; maxc=0;

fig=figure;
set(fig,'PaperUnits','centimeters','PaperOrientation','portrait',... 
            'Paperposition', [0 0 29.7 21],'Papertype','A4')

m_proj...
('mercator','lon',[plonmin plonmax],'lat',[platmin  platmax]);
%[cs ,h]=m_drawbathy_col('contourf','gebco',[minc:200:maxc],limitcmap);
[cs ,h]=m_drawbathy2('contourf','etopo6min',[minc:200:maxc],limitcmap);
%m_gshhs('fc','color','k');
%m_gshhs('fr','color','b');

m_grid('box','fancy','tickdir','in','fontsize',fs_map)
set(gca,'fontsize',fs_map)

% caxis([minc maxc])
% tikc=(maxc-minc)/100;
% ccl=minc:tikc:maxc;
% %ac=clrbr(ax2,ccl,[sx+0.07 .07 .01 sy]);
% %ac=clrbrh(gca,ccl,[0.3 0.16 .4 0.02]);
% ac=clrbrh(gca,ccl,[0.17 0.28 .4 0.016]);
% set(ac,'fontsize',fs_map,'fontweight','bold')
% %shading interp


print(gcf,'-dpng','map_lion.png')
