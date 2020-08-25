% **********************************************************************************************************************************
% The purpose of this function is to plot fluorscence.
% **********************************************************************************************************************************
function Plot_WL_Fluor( titlestring, WETLabsData, dmin, dmax, id_and_cal_info)

    nMarkerSize = 4;
    mrk = {'o', 's', 'd', '*', '+', 'v', '^', '<', '>'}.';
    iSubPlot = 0; % incrementer to track where we are as we build NumSubPlotsNeeded SubPlots
    iPlot = 0;    % incrementer to track the number of plots produced; this facilitates color changing and marker changing
  
    setAxesFontSize(8);
    fluor_title = sprintf('Fluorescence\n\n%s        %s', titlestring, id_and_cal_info);
    title(fluor_title, 'FontSize', 8);

    if dmin >= dmax
        dmax = dmin + 1.0;
    end

    bmin = double(intmax('int32'));
    bmax = -double(intmax('int32'));
    bdel = 0;

    % ------------------------------------------------------------------------------------------------------------------------------
    % tally if different types of units are in the data collection: ppb, ug/l, and counts are all possibilities
    % ------------------------------------------------------------------------------------------------------------------------------
    isPpb = false;
    isUgl = false;
    isCounts = false;
    isPpbPhyco = false;
    isCountsPhyco = false;
  
    for sensorIndex = 1 : length(WETLabsData)
        if isempty(WETLabsData)
            return;
        end

        % --------------------------------------------------------------------------------------------------------------------------
        %
        % --------------------------------------------------------------------------------------------------------------------------
        if isfield(WETLabsData, 'CDOM')
            for i= 1 : length(WETLabsData.CDOM)
                CDOM = WETLabsData.CDOM{i};
                if isempty(CDOM.scaled) && isempty(CDOM.counts)
                    continue;
                end

                if (strcmp(CDOM.units, 'ppb'))
                    isPpb = true;
                else
                    isCounts = true;
                end
            end
        end

        % --------------------------------------------------------------------------------------------------------------------------
        %
        % --------------------------------------------------------------------------------------------------------------------------
        if isfield(WETLabsData, 'Chlorophyll')
            for i= 1 : length(WETLabsData.Chlorophyll)
                Chlorophyll = WETLabsData.Chlorophyll{i};
                if isempty(Chlorophyll) || (isempty(Chlorophyll.scaled) && isempty(Chlorophyll.counts))
                    continue;
                end

                if (strcmp(Chlorophyll.units, '{\mu}g/l'))
                    isUgl = true;
                else
                    isCounts = true;
                end
            end
        end
        
        % --------------------------------------------------------------------------------------------------------------------------
        %
        % --------------------------------------------------------------------------------------------------------------------------
        if isfield(WETLabsData, 'Phycoerythrin')
            for i= 1 : length(WETLabsData.Phycoerythrin)
                Phycoerythrin = WETLabsData.Phycoerythrin{i};
                if isempty(Phycoerythrin) || (isempty(Phycoerythrin.scaled) && isempty(Phycoerythrin.counts))
                    continue;
                end

                if (strcmp(Phycoerythrin.units, 'ppb'))
                    isPpbPhyco = true;
                else
                    isCountsPhyco = true;
                end
            end
        end
    end
  
    % Determine the number of sub plots needed based on the number of distinct unit types found
    numSubPlotsNeeded = calcSubPlotsRequired(isPpb, isUgl, isCounts, isPpbPhyco, isCountsPhyco);
            
    % ------------------------------------------------------------------------------------------------------------------------------
    % PLOT CDOM ppb vs. DEPTH ON ITS OWN AXIS
    % ------------------------------------------------------------------------------------------------------------------------------
    if isPpb
        iSubPlot = iSubPlot + 1;
        
        ax = subplot(numSubPlotsNeeded, 1, iSubPlot);

        adjustOuterPosition(ax, numSubPlotsNeeded);
        if (iSubPlot==1)
            setAxesFontSize(8);
            title(fluor_title, 'FontSize', 8);
        end
        
        xlabel('ppb');
        ylabel('Depth (m)');
        set(ax, 'YDir', 'reverse');

        grid on;
        hold on;

        legend_Reset();

        bmin = double(intmax('int32'));
        bmax = -double(intmax('int32'));
        bdel = 0;
    
        for sensorIndex = 1 : length(WETLabsData)
            if isempty(WETLabsData)
                continue;
            end
            
            if ~isfield(WETLabsData, 'CDOM')
                continue;
            end

            for i= 1 : length(WETLabsData.CDOM)
                CDOM = WETLabsData.CDOM{i};
                if isempty(CDOM) || ~strcmp(CDOM.units, 'ppb')
                    continue;
                end

                bmin = min(bmin, min(CDOM.scaled)); %#ok<AGROW>
                bmax = max(bmax, max(CDOM.scaled)); %#ok<AGROW>

                % Descent
                if ~isempty(CDOM.iwn)
                    iPlot = iPlot + 1;
                    [R G B] = Wavelength2RGB(360 + (iPlot*50));
                    line(CDOM.scaled(CDOM.iwn) ,CDOM.depth(CDOM.iwn), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                    legend_AddItem(sprintf('CDOM #%u - Descent', sensorIndex)); %#ok<AGROW>
                end

                % Ascent
                if ~isempty(CDOM.iwp)
                    iPlot = iPlot + 1;
                    [R G B] = Wavelength2RGB(360 + (iPlot*50)+25);
                    line(CDOM.scaled(CDOM.iwp), CDOM.depth(CDOM.iwp),'Marker', mrk{iPlot},'Color', [R G B],'MarkerSize', nMarkerSize, 'LineStyle', '-');
                    legend_AddItem(sprintf('CDOM #%u - Ascent', sensorIndex)); %#ok<AGROW>
                end
            end
            
        end

        [tmpMin, tmpMax, bdel] = calcPlotScale(bmin, bmax);

        if isfield(CDOM, 'depth') && ~isempty(CDOM.depth)
            [dmin, dmax, ddel] = calcPlotScale(min(CDOM.depth), max(CDOM.depth));
        end

        set(ax, 'XAxisLocation', 'bottom', 'XLim', [tmpMin tmpMax], 'XTick', tmpMin : bdel : tmpMax, 'YLim', [dmin dmax], 'YDir', 'reverse');

        legend_Show();
    end
    
    % ------------------------------------------------------------------------------------------------------------------------------
    % PLOT Phycoerythrin ppb vs. DEPTH ON ITS OWN AXIS
    % ------------------------------------------------------------------------------------------------------------------------------
    if isPpbPhyco
        iSubPlot = iSubPlot + 1;
        ax = subplot(numSubPlotsNeeded, 1, iSubPlot);

        adjustOuterPosition(ax, numSubPlotsNeeded);

        if (iSubPlot==1)
            setAxesFontSize(8);
            title(fluor_title, 'FontSize', 8);
        end

        xlabel('ppb');
        ylabel('Depth (m)');
        set(ax, 'YDir', 'reverse');
        
        grid on;
        hold on;

        legend_Reset();

        bmin = double(intmax('int32'));
        bmax = -double(intmax('int32'));
        bdel = 0;
    
        for sensorIndex = 1 : length(WETLabsData)
            if isempty(WETLabsData)
                continue;
            end
            
            if ~isfield(WETLabsData, 'Phycoerythrin')
                continue;
            end
            
            for i= 1 : length(WETLabsData.Phycoerythrin)
                Phycoerythrin = WETLabsData.Phycoerythrin{i};
                if isempty(Phycoerythrin) || ~strcmp(Phycoerythrin.units, 'ppb')
                    continue;
                end

                bmin = min(bmin, min(Phycoerythrin.scaled));
                bmax = max(bmax, max(Phycoerythrin.scaled));

                % Descent
                if ~isempty(Phycoerythrin.iwn)
                    iPlot = iPlot + 1;
                    [R G B] = Wavelength2RGB(360 + (iPlot*50));
                    line(Phycoerythrin.scaled(Phycoerythrin.iwn) ,Phycoerythrin.depth(Phycoerythrin.iwn), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                    legend_AddItem(sprintf('Phycoerythrin #%u - Descent', sensorIndex));
                end

                % Ascent
                if ~isempty(Phycoerythrin.iwp)
                    iPlot = iPlot + 1;
                    [R G B] = Wavelength2RGB(360 + (iPlot*50)+25);
                    line(Phycoerythrin.scaled(Phycoerythrin.iwp), Phycoerythrin.depth(Phycoerythrin.iwp),'Marker', mrk{iPlot},'Color', [R G B],'MarkerSize', nMarkerSize, 'LineStyle', '-');
                    legend_AddItem(sprintf('Phycoerythrin #%u - Ascent', sensorIndex));
                end
            end
        end 
        
        [tmpMin, tmpMax, bdel] = calcPlotScale(bmin, bmax);
        
        if isfield(Phycoerythrin, 'depth') && ~isempty(Phycoerythrin.depth)
            [dmin, dmax, ddel] = calcPlotScale(min(Phycoerythrin.depth), max(Phycoerythrin.depth));
        end

        set(ax, 'XAxisLocation', 'bottom', 'XLim', [tmpMin tmpMax], 'XTick', tmpMin : bdel : tmpMax, 'YLim', [dmin dmax], 'YDir', 'reverse');

        legend_Show();
    end
    
    % ----------------------------------------------------------------------------------------------------
    % PLOT Chlorophyll ug/L vs. DEPTH ON its own axis
    % ----------------------------------------------------------------------------------------------------
    if isUgl
        iSubPlot = iSubPlot + 1;
        ax = subplot(numSubPlotsNeeded, 1, iSubPlot);

        adjustOuterPosition(ax, numSubPlotsNeeded);

        if (iSubPlot==1)
            setAxesFontSize(8);
            title(fluor_title, 'FontSize', 8);
        end
    
        xlabel('{\mu}g/l');
        ylabel('Depth (m)');

        grid on;
        hold on;

        legend_Reset();

        bmin = double(intmax('int32'));
        bmax = -double(intmax('int32'));
        bdel = 0;
    
        for sensorIndex = 1 : length(WETLabsData)
            if isempty(WETLabsData) || ~isfield(WETLabsData, 'Chlorophyll')
                continue;
            end
      
            for i= 1 : length(WETLabsData.Chlorophyll)
                Chlorophyll = WETLabsData.Chlorophyll{i};
                if isempty(Chlorophyll) || ~strcmp(Chlorophyll.units, '{\mu}g/l')
                    continue;
                end
        
                bmin = min(bmin, min(Chlorophyll.scaled));
                bmax = max(bmax, max(Chlorophyll.scaled));
        
                % Descent
                if ~isempty(Chlorophyll.iwn)
                    iPlot = iPlot + 1;
                    [R G B] = Wavelength2RGB(360 + (iPlot*50));
                    line(Chlorophyll.scaled(Chlorophyll.iwn), Chlorophyll.depth(Chlorophyll.iwn), 'Marker', mrk{iPlot},'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                    legend_AddItem(sprintf('Chlorophyll #%u - Descent', sensorIndex));
                end
        
                % Ascent
                if ~isempty(Chlorophyll.iwp)
                    iPlot = iPlot + 1;
                    [R G B] = Wavelength2RGB(360 + (iPlot*50)+25);
                    line(Chlorophyll.scaled(Chlorophyll.iwp), Chlorophyll.depth(Chlorophyll.iwp), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                    legend_AddItem(sprintf('Chlorophyll #%u - Ascent', sensorIndex));
                end
            end
        end
    
        [tmpMin, tmpMax, bdel] = calcPlotScale(bmin, bmax);

        if isfield(Chlorophyll, 'depth') && ~isempty(Chlorophyll.depth)
            [dmin, dmax, ddel] = calcPlotScale(min(Chlorophyll.depth), max(Chlorophyll.depth));
        end
        
        set(ax, 'XAxisLocation', 'bottom', 'XLim', [tmpMin tmpMax], 'XTick', tmpMin : bdel : tmpMax, 'YLim', [dmin dmax], 'YDir', 'reverse');

        legend_Show();
    end
  
    % ----------------------------------------------------------------------------------------------------
    % PLOT Fluorescence counts vs. DEPTH ON its own axis
    % ----------------------------------------------------------------------------------------------------
    if isCounts
        iSubPlot = iSubPlot + 1;
        ax = subplot(numSubPlotsNeeded, 1, iSubPlot);

        adjustOuterPosition(ax, numSubPlotsNeeded);

        if (iSubPlot==1)
            setAxesFontSize(8);
            title(fluor_title, 'FontSize', 8);
        end

        xlabel('A/D Counts');
        ylabel('Depth (m)');

        grid on;
        hold on;

        legend_Reset();

        for sensorIndex = 1 : length(WETLabsData)
            if isempty(WETLabsData)
                continue;
            end

            if isfield(WETLabsData, 'CDOM')
                for i = 1 : length(WETLabsData.CDOM)
                    CDOM = WETLabsData.CDOM{i};
                    if isempty(CDOM) || ~strcmp(CDOM.units, '(A/D Counts)')
                        continue;
                    end

                    % Descent
                    if ~isempty(CDOM.iwn)
                        iPlot = iPlot + 1;
                        [R G B] = Wavelength2RGB(360 + (iPlot*50));
                        line(CDOM.counts(CDOM.iwn), CDOM.depth(CDOM.iwn), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                        legend_AddItem(sprintf('CDOM #%u - Descent', sensorIndex));
                    end
          
                    % Ascent
                    if ~isempty(CDOM.iwp)
                        iPlot = iPlot + 1;
                        [R G B] = Wavelength2RGB(360 + (iPlot*50)+25);
                        line(CDOM.counts(CDOM.iwp), CDOM.depth(CDOM.iwp), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                        legend_AddItem(sprintf('CDOM #%u - Ascent', sensorIndex));
                    end
                end
            end
      
            if isfield(WETLabsData, 'Phycoerythrin')
                for i = 1 : length(WETLabsData.Phycoerythrin)
                    Phyco = WETLabsData.Phycoerythrin{i};
                    if isempty(Phyco) || ~strcmp(Phyco.units, '(A/D Counts)')
                        continue;
                    end

                    % Descent
                    if ~isempty(Phyco.iwn)
                        iPlot = iPlot + 1;
                        [R G B] = Wavelength2RGB(360 + (iPlot*50));
                        line(Phyco.counts(Phyco.iwn), Phyco.depth(Phyco.iwn), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                        legend_AddItem(sprintf('Phycoerythrin #%u - Descent', sensorIndex));
                    end
          
                    % Ascent
                    if ~isempty(Phyco.iwp)
                        iPlot = iPlot + 1;
                        [R G B] = Wavelength2RGB(360 + (iPlot*50)+25);
                        line(Phyco.counts(Phyco.iwp), Phyco.depth(Phyco.iwp), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                        legend_AddItem(sprintf('Phycoerythrin #%u - Ascent', sensorIndex));
                    end
                end
            end
      
            if isfield(WETLabsData, 'Chlorophyll')
                for i = 1 : length(WETLabsData.Chlorophyll)
                    Chlorophyll = WETLabsData.Chlorophyll{i};
                    if isempty(Chlorophyll) || ~strcmp(Chlorophyll.units, '(A/D Counts)')
                        continue;
                    end
          
                    % Descent
                    if ~isempty(Chlorophyll.iwn)
                        iPlot = iPlot + 1;
                        [R G B] = Wavelength2RGB(360 + (iPlot*50));
                        line(Chlorophyll.counts(Chlorophyll.iwn), Chlorophyll.depth(Chlorophyll.iwn), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                        legend_AddItem(sprintf('Chlorophyll #%u - Descent', sensorIndex));
                    end
          
                    % Ascent
                    if ~isempty(Chlorophyll.iwp)
                        iPlot = iPlot + 1;
                        [R G B] = Wavelength2RGB(360 + (iPlot*50)+25);
                        line(Chlorophyll.counts(Chlorophyll.iwp), Chlorophyll.depth(Chlorophyll.iwp), 'Marker', mrk{iPlot}, 'Color', [R G B], 'MarkerSize', nMarkerSize, 'LineStyle', '-');
                        legend_AddItem(sprintf('Chlorophyll #%u - Ascent', sensorIndex));
                    end
                end
            end
        end

        set(ax, 'XAxisLocation', 'bottom', 'YLim', [dmin dmax], 'YDir', 'reverse');

        hTxt = text(0, 0, '(Missing or corrupt sg\_calib\_constants file)', 'HorizontalAlignment', 'center', 'FontSize', 8);
        xLabelPos = get(get(gca, 'XLabel'), 'Position');
        set(hTxt, 'Position', xLabelPos);

        legend_Show();
    end
  
    if ~isPpb && ~isUgl && ~isCounts && ~isPpbPhyco
        grid off;

        set(gcf, 'ActivePositionProperty', 'OuterPosition');
        height = 0.95;
        set(gca, 'OuterPosition', [0.0 0.0 1.0 height]);

        %title(sprintf('Fluorescence\n\n%s\n%s\n', titlestring, id_and_cal_info), 'FontSize', 8);
        title(fluor_title, 'FontSize', 8);

        displayNoDataMessage('Fluorescence data does not exist for this dive');
    end
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function adjustOuterPosition(ax, numSubPlotsNeeded)
    tmp1 = get(ax, 'OuterPosition');

    if (numSubPlotsNeeded > 1)
        tmp1(4) = (1.0 / numSubPlotsNeeded) - 0.05;
    else
        tmp1(4) = tmp1(4) - 0.075;
    end
    
    set(ax, 'OuterPosition', tmp1);
end

% **********************************************************************************************************************************
%
% **********************************************************************************************************************************
function [numSubPlotsNeeded] = calcSubPlotsRequired(varargin)
    numSubPlotsNeeded = sum(cell2mat(varargin));
end
