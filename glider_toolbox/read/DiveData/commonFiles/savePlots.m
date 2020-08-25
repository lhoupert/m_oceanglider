% ****************************************************************************************************
% Saves the specified figure to the path specified.
% ****************************************************************************************************
function savePlots(hFig, pathname)
    global gSavePlotsInProgress;

    gSavePlotsInProgress = 1;

    imgExt = 'jpg';

    bOverwrite = 1;

    if exist(pathname, 'dir') == 0
        mkdir(pathname);
    end

    if isempty(hFig)
        [plotFigureHandles, plotFigureNames] = getPlotFigureHandles();
        if (isempty(plotFigureHandles) == 1)
            vizMsgBox('Please generate plots before trying to save them !', 'ERROR: No Plots', 'modal');
            return;
        end

        hWaitbar = waitbar(0, 'Saving plots as image files...');
        
        for i=1:length(plotFigureHandles)
            name = strtrim(plotFigureNames(i,:));
            ns_name = strrep(name, ' ', '_');

            if (bOverwrite == 1)
               image_filename = fullfile(pathname, sprintf('%s.%s', ns_name, imgExt));
            else
                index = 1;
                image_filename = fullfile(pathname, sprintf('%s %3.3d.%s', ns_name, index, imgExt));
                while exist(image_filename)
                    index = index + 1;
                    image_filename = fullfile(pathname, sprintf('%s %3.3d.%s', ns_name, index, imgExt));
                end
            end

            strWaitMsg = name;
            waitbar(i/length(plotFigureHandles), hWaitbar, sprintf('Saving file %s...', strWaitMsg));

            % PUTTING THIS TRY/CATCH IN PLACE JUST IN CASE THE END-USER IGNORES THE "SAVING file Plot..." PROGRESS DIALOG AND
            % CLOSES PLOTS BEFORE THEY ARE SAVED
            try
                saveas(plotFigureHandles(i), image_filename, imgExt);
            catch
                vizMsgBox('savePlots unhappiness !', 'ERROR', 'modal');
            end
        end %FOR

        close(hWaitbar);
    else
        index = 1;
        image_filename = fullfile(pathname, sprintf('%s %3.3d.%s', strtrim(get(hFig, 'Name')), index, imgExt));
        while exist(image_filename)
            index = index + 1;
            image_filename = fullfile(pathname, sprintf('%s %3.3d.%s', strtrim(get(hFig, 'Name')), index, imgExt));
        end
%          print(hFig, '-djpeg', '-loose', image_filename);
         saveas(hFig, image_filename, imgExt);
    end
    
    gSavePlotsInProgress = 0;
end
