% ******************************************************************************
% Function: calcPlotScale
% 
% Purpose: the purpose of this function is to calculate the scale and
%          limits of the plot based upon the data being plotted.
%
% ******************************************************************************
function [tmpMin, tmpMax, bdel] = calcPlotScale(bmin, bmax, nPrecision, bPad)

    if exist('nPrecision', 'var') == 0
        nPrecision = 2;   % DEFAULT TO TWO (2) DECIMALS OF PRECISION
    end

    bmin = roundNumber(min(bmin), nPrecision);
    bmax = roundNumber(max(bmax), nPrecision);
    if ~isfinite(bmin) || ~isfinite(bmax) || bmin >= bmax
        bmin = bmin - (bmin / 2);
        bmax = bmax + (bmax / 2);
    end

    % IF WE ARE PLOTTING REALLY SMALL NUMBERS AND THE PRECISION SPECIFIED DOES NOT HAVE ENOUGH RESOLUTION, AUTOMATICALLY DETERMINE
    % THE OPTIMAL PRECISION
    nExponent = nPrecision;
    if (bmax-bmin) > 0
        while ((bmax-bmin) / 10) < 1/10^nExponent
            nExponent = nExponent + 1;
        end
    end
    bdel = roundNumber(((bmax-bmin) / 10), nExponent);

    % THIS WILL ADD SOME PADDING TO THE LEFT AND RIGHT OF THE PLOT SO THAT THE DATA IS NOT SMASHED UP AGAINST THE SIDE OF THE FIGURE
    if exist('bPad')
        if strcmpi(bPad, 'noPad') == 1
            tmpMin = roundNumber(bmin, nPrecision);
            tmpMax = roundNumber(bmax, nPrecision);
        else
            tmpMin = roundNumber(bmin - (1*bdel), nPrecision);
            tmpMax = roundNumber(bmax + (1*bdel), nPrecision);
        end
    else
        tmpMin = roundNumber(bmin - (1*bdel), nPrecision);
        tmpMax = roundNumber(bmax + (1*bdel), nPrecision);
    end

    % Simulated dives are not truly simulated because
    % the simulation code actually generates bad data.
    % We must "band-aid" the code to account for this
    % since the visualization scripts are also being
    % utilized for plotting simulated dives.
    if tmpMin >= tmpMax
        tmpMin = 0;
    end
end
