% **********************************************************************************************************************************
% Author: John E. Stranzl Jr.
% Date:   ##-???-2011
%
% Function: Plot_PAR
%
% Purpose:  The purpose of this function is to plot the data from the PAR sensor
%
% @param titleString    char - the title to be displayed in the title bar atop the plot window
% 
% @param PAR_CalData    struct - this is the calibration data for the PAR sensor that was read from the sg_calib_constants file
%                                which is to be applied to the acquired PAR sensor data before it is plotted
% 
% @returns              there are no returns from this function
% 
% 
% Author           Date             Remarks
% ----------------------------------------------------------------------------------------------------------------------------------
% ???              ??-???-2010      INITIAL INCORPORATION
% 
% John E. Stranzl	09-Aug-2012     BUG:      Enhnacement
%                                   PHASE:    Unit Testing
%                                   Problem:  The data structures were re-architected but the corresponding changes were not made to
%                                             this function.
%                                   Solution: Modified the code to utilize the re-architected data structures.
% 
% **********************************************************************************************************************************
function Plot_PAR( titleString, PAR_CalData, eng )
    setFigureProperties();
    
%     displayTitle(titleString);

    set(gca, 'OuterPosition', [0.0 0.0 1.0 0.975]);
    setAxesFontSize(8);
    
    % Calibration info.
    if isfield(PAR_CalData, 'serialNumber') && isfield(PAR_CalData, 'calDate')
        par_id_and_cal_info = sprintf(' SN: %d, CAL: %s ', PAR_CalData.serialNumber, PAR_CalData.calDate);
    else
        par_id_and_cal_info = 'Unknown calibration';
    end
    
    % Plot title, including calibration info.
    par_title = sprintf('Photosynthetic Active Radiation (PAR)\n\n%s        %s', titleString, par_id_and_cal_info);
    title(par_title, 'FontSize', 8);
    
    %addVersionNumberToPlot();

    setAxisProperties();
    
    try
        %[~, PAR_Data, depth, pitchCtl] = dataStore.queryDoubleClean('PARuV', 'depth', 'pitchCtl');
        PAR_Data = eng.qsp_PARuV;
        depth = eng.depth;
        pitchCtl = eng.pitchCtl;

        [~,n] = size(PAR_Data);
        for i = 1:n
            set(gca, 'YDir', 'reverse');

            % FIND INDEXES FOR DOWNCAST AND DESECENT DATA POINTS
            [iwn, iwp] = findUpcastDowncastIndices(depth, pitchCtl);

            % APPLY CALIBRATION DATA TO ACQUIRED DATA
            darkOffset = [];
            scaleFactor = [];
            if ~isempty(PAR_CalData) && (length(PAR_CalData) >= i)
                if isfield(PAR_CalData(i), 'darkOffset')
                    darkOffset = PAR_CalData(i).darkOffset;
                elseif isfield(PAR_CalData(i), 'sensorDark')
                    darkOffset = PAR_CalData(i).sensorDark;
                end
                
                if isfield(PAR_CalData(i), 'scaleFactor')
                    scaleFactor = PAR_CalData(i).scaleFactor;
                elseif isfield(PAR_CalData(i), 'scalingFactor')
                    scaleFactor = PAR_CalData(i).scalingFactor;
                end
            end

            if ~isempty(darkOffset) && ~isempty(scaleFactor)
                units = '(\mu mol photons\bulletm^-^2\bullets^-^1)';
                PAR_Data = (PAR_Data - (darkOffset / 1000)) / scaleFactor;
            else
                units = '(millivolts)';
            end

            hold on;
            grid on;

            xMax = 100;

            if ~isempty(iwn)
                semilogx(PAR_Data(iwn), depth(iwn), 'Marker', 'v', 'MarkerSize', 4, 'Color', 'b', 'LineStyle', '-');
            end
            if ~isempty(iwp)
                semilogx(PAR_Data(iwp), depth(iwp), 'Marker', '^', 'MarkerSize', 4, 'Color', 'm', 'LineStyle', '-');
            end

            % SCALE AND SET THE LIMITS OF THE X- and Y-AXIS
            [tmpMin, tmpMax, tdel] = calcPlotScale(min(PAR_Data), max(PAR_Data), 6);
            set(gca, 'XLim', [tmpMin tmpMax], 'XTick', (tmpMin : tdel : tmpMax));

            labelAxis(units);

            if (exist('PAR_CalData', 'var') == 1)
                if ~isempty(PAR_CalData) && (length(PAR_CalData) >= i) && (isempty(PAR_CalData(i).serialNumber) == 0)
                    displayLegend(PAR_CalData(i).serialNumber);
                else
                    displayLegend([]);
                end
            else
                displayLegend([]);
            end
        end
    catch err  % NO NEED TO RETHROW err NOR REPORT ADDITIONAL ERROR DETAILS THAT ARE AVAILABLE FROM err
        displayNoDataAnnotation();
    end    % end try/catch
end

% ****************************************************************************************************
% The purpose of this function is to label the plot axiis.
% ****************************************************************************************************
function displayTitle(titleString)
    set(gca, 'OuterPosition', [0.0 0.0 1.0 0.975]);
    setAxesFontSize(8);
    title(sprintf('Photosynthetic Active Radiation (PAR)\n\n%s', titleString), 'FontSize', 8);
    addVersionNumberToPlot();
end

% ****************************************************************************************************
% The purpose of this function is to label the plot axiis.
% ****************************************************************************************************
function setFigureProperties()
    set(gcf, 'ActivePositionProperty', 'OuterPosition');
end

% ****************************************************************************************************
% The purpose of this function is to label the plot axiis.
% ****************************************************************************************************
function setAxisProperties()
    hold off;
    grid on;

    set(gca, 'OuterPosition', [0.0 0.0 1.0 1.0]);
    setAxesFontSize(8);
end

% ****************************************************************************************************
% The purpose of this function is to label the plot axiis.
% ****************************************************************************************************
function labelAxis(units)
    xlabel(['PAR ' units], 'FontSize', 8); 
    ylabel('Depth (m)', 'FontSize', 8);
end

% ****************************************************************************************************
% The purpose of this function is to display a legend inside the plot area.
% ****************************************************************************************************
function displayLegend(serialNumber)
    
    if isempty(serialNumber)
        legendPAR = sprintf('PAR');
    else
        legendPAR = sprintf('PAR (s/n %d)', serialNumber);
    end
    
    hLgnd = legend([legendPAR ' - Descent'], [legendPAR ' - Ascent'], 'Location', 'SouthEast');
    
    set(hLgnd, 'FontSize', 8);
end

% ****************************************************************************************************
% The purpose of this function is to display an error message in the center of the plot when
% calibration data is not available for the sensor.
% ****************************************************************************************************
function displayNoDataAnnotation()
    grid off;

    v = axis;
    xtxt = (v(2) - v(1)) / 2.0;
    ytxt = (v(4) - v(3)) / 2.0;
    comment = 'PAR data does not exist for this dive';
    text(xtxt, ytxt, comment, 'HorizontalAlignment', 'center', 'FontSize', 12);

    set(gca, 'ytick', []);
    set(gca, 'xtick', []);
end
