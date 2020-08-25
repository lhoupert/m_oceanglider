function addVersionNumberToPlot()
    h = uicontrol('Style', 'text', 'String', 'Ver. Trunk (post 3.06)',...
                  'Units', 'normalized', 'horizontalalignment', 'left',...
                  'FontName', 'FixedWidth', 'FontSize', 10, 'Position', [0.0025 0.975 0.075 0.025]); 

    parentColor = get(get(h, 'parent'), 'color');
    
    set(h, 'foregroundcolor', [0 0 0], 'backgroundcolor', parentColor);
end
