function legend_Show(handleVectors)
    global    gstrLegend;

    if ~isempty(gstrLegend)
        if exist('handleVectors', 'var')
        	hLgnd = legend(handleVectors, gstrLegend, 'Location', 'SouthEast', 'Color', 'None');
        else
        	hLgnd = legend(gstrLegend, 'Location', 'SouthEast', 'Color', 'None');
        end
        set(hLgnd, 'FontSize', 8);
    end
end