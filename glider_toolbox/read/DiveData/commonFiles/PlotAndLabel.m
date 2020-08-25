% ****************************************************************************************************
% The purpose of this function is to put labels on the lines being
% plotted so that the graph can be read and interpretted easily.
% ****************************************************************************************************
function PlotAndLabel(x, y, strColor, strLabel)
    h = plot(x, y, strColor);

    if (getAnnotatePlots() == 1)
        ch = clabel([length(x) x'; length(y) (y)'], h);
        grayValue = 240.0/255.0;
        set(ch, 'string', strLabel, 'backgroundcolor', [grayValue grayValue grayValue]);
    end
end
