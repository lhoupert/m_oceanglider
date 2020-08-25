function [a]=clrbr3(ax,clip,pos,labelstr)
%
% a = clrbr(ax,clip,pos);
%
clip = clip(:);
a = axes('position',pos);
sc=pcolor(repmat([0 1],length(clip),1),repmat(clip,1,2),repmat(clip,1,2));
%     rgb_colors = get(sc, 'CData');
%     [i, map] = rgb2ind(rgb_colors, length(clip));
%     set(s, 'CData', double(i)); 
shading flat;
view(0,-90)
caxis([clip(1) clip(end)])
set(a,'xlim',[0.1 0.9],'ylim',[clip(1) clip(end)],...
      'xtick',[0 1],'xticklabel',[],'yaxislocation','right','box','on','ydir','rev')
ylabel(labelstr)           
map = get(gcf,'colormap');
grid off
colormap(map)
%freezeColors
b = axes('position',get(a,'position'),...
         'xtick',get(a,'xtick'),...
         'xlim',get(a,'xlim'),...
         'xticklabel',[],...
         'ytick',get(a,'ytick'),...
         'ylim',get(a,'ylim'),...
         'yticklabel',[],...
         'color','none','box','on');
set(gca,'yaxisLocation','right')     

set(gcf,'currentaxes',ax);
										   
