% *************************************************************************
% The purpose of this function is to plot backscatter.
% *************************************************************************
function Plot_WL_Scatter (titlestring, WETLabsData, id_and_cal_info)

    nMarkerSize = 4;
    mrk = {'o', 's', 'd', '*', '+', 'v', '^'}.';

    set(gca, 'OuterPosition', [0.0 0.0 1.0 0.975]);

    setAxesFontSize(8);
    title(sprintf('Scattering\n\n%s        %s', titlestring, id_and_cal_info), 'FontSize', 8);

    legend_Reset();

    cmap = rand(20,3);
    cmap(1,1:3) = [1 0 0];
    cmap(2,1:3) = [0 1 0];
    cmap(3,1:3) = [0 0 1];
    cmap(4,1:3) = [1 1 0];
    cmap(5,1:3) = [0 1 1];

    numCounts = 0;
    numScaled = 0;

    hold on;
    grid on;
    set(gca,'YDir','reverse');    

    bmin = double(intmax('int32'));
    bmax = -double(intmax('int32'));
    bdel = 0;

    dmin = 0;
    dmax = 100;

    Scattering = [];

    if isempty(WETLabsData)
        return;
    end

    if ~isfield(WETLabsData, 'Scattering')
        return;
    end

    % There is one WETLabsData.Scattering structure for
    % each backscatter channel supported by the puck.
    for i = 1 : length(WETLabsData.Scattering)
        Scattering = WETLabsData.Scattering{i};
        if isempty(Scattering)
            continue;
        end

        if ~isempty(Scattering.scaled)
            bmin = min(bmin, min(Scattering.scaled));
            bmax = max(bmax, max(Scattering.scaled));
            numScaled = numScaled + 1;
        elseif ~isempty(Scattering.counts)
            bmin = min(bmin, min(Scattering.counts));
            bmax = max(bmax, max(Scattering.counts));
            numCounts = numCounts + 1;
        end

        % DOWN CAST
        nWavelength = Scattering.refData(1);

        if ~isempty(Scattering.iwn)
            if isempty(Scattering.scaled)
                line(Scattering.counts(Scattering.iwn), Scattering.depth(Scattering.iwn), 'Marker', mrk{i}, 'Color', cmap(i,:), 'MarkerSize', nMarkerSize, 'LineStyle', '-');
            else
                line(Scattering.scaled(Scattering.iwn), Scattering.depth(Scattering.iwn), 'Marker', mrk{i}, 'Color', cmap(i,:), 'MarkerSize', nMarkerSize, 'LineStyle', '-');
            end

            legend_AddItem(sprintf('%u nm - Descent', nWavelength));
        end  

        % UP CAST
        if ~isempty(Scattering.iwp)
            if isempty(Scattering.scaled)
                line(Scattering.counts(Scattering.iwp), Scattering.depth(Scattering.iwp), 'Marker', mrk{i}, 'Color', cmap(i,:)/2.0, 'MarkerSize', nMarkerSize, 'LineStyle', '-');
            else
                line(Scattering.scaled(Scattering.iwp), Scattering.depth(Scattering.iwp), 'Marker', mrk{i}, 'Color', cmap(i,:)/2.0, 'MarkerSize', nMarkerSize, 'LineStyle', '-');
            end

            legend_AddItem(sprintf('%u nm - Ascent', nWavelength ));
        end
    end

    % ------------------------------------------------------------------------------------------------------------------------------
    % DISPLAY ADVISORIES AS NEEDED
    % ------------------------------------------------------------------------------------------------------------------------------
    if (numCounts <= 0) && (numScaled <= 0)
        displayNoDataMessage('Scattering data does not exist for this dive');
    else
        % ------------------------------------------------------------------------------------------------------------------------------
        % SCALE, FORMAT AND LABEL AXIIS
        % ------------------------------------------------------------------------------------------------------------------------------
        [tmpMin, tmpMax, bdel] = calcPlotScale(bmin, bmax);
        if isfield(Scattering, 'depth') && ~isempty(Scattering.depth)
            [dmin, dmax, ddel] = calcPlotScale(min(Scattering.depth), max(Scattering.depth));
        end
        set(gca, 'XAxisLocation', 'bottom', 'XLim', [tmpMin tmpMax], 'XTick', tmpMin : bdel : tmpMax, 'YLim', [dmin dmax]);
        if numCounts <= 0
            xlabel('\beta(\theta_c) (m^-^1 sr^-^1)');
        else
            xlabel('A/D Counts');
        end
        ylabel( 'Depth (m)');

        legend_Show();

        if (numScaled < length(WETLabsData.Scattering))
            hTxt = text(0, 0, '(Missing or corrupt sg\_calib\_constants file)', 'HorizontalAlignment', 'center', 'FontSize', 8);
            xLabelPos = get(get(gca, 'XLabel'), 'Position');
            set(hTxt, 'Position', xLabelPos+[0 0.5 0]);
        end
    end
end
