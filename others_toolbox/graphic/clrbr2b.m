function [a]=clrbr2b(ax,clip,pos,labelstr)
%[a]=clrbr2b(ax,clip,pos,labelstr)
%  for contoourf with extreme data values outside clip
%
clip=clip(:);

cliplim1=clip(2)-(clip(3)-clip(2));
cliplim2=clip(end-1);%+ (clip(end-1)-clip(end-2)); 
a = axes('position',pos);
sc=contourf(repmat([0 1],length(clip),1),repmat(clip,1,2),repmat(clip,1,2),...
    clip,'linestyle','none');
    %rgb_colors = get(sc, 'CData');
    %[i, map] = rgb2ind(rgb_colors, length(clip));
    %set(s, 'CData', double(i)); 
% shading flat;
% view(0,-90)
caxis([cliplim1 cliplim2])
%dclip=( clip(end)-clip(1) )/length(clip);
set(a,'xlim',[0.1 0.9],'ylim',[cliplim1 cliplim2],...
      'xtick',[0 1],'xticklabel',[],'yaxislocation','right','box','on','color','none')
ylabel(labelstr)           
%map = get(gcf,'colormap');
grid off
%colormap(map)
%freezeColors
% b = axes('position',get(a,'position'),...
%          'xtick',get(a,'xtick'),...
%          'xlim',get(a,'xlim'),...
%          'xticklabel',[],...
%          'ytick',get(a,'ytick'),...
%          'ylim',get(a,'ylim'),...
%          'yticklabel',[],...
%          'color','none','box','on');
% set(gca,'yaxisLocation','right')     

set(gcf,'currentaxes',ax);
										   
