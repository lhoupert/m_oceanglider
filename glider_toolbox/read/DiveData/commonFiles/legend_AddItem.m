function legend_AddItem(varargin)
    global    gstrLegend;
    
    for i=1:length(varargin)
        gstrLegend{end+1} = char(varargin(i));
    end
end
