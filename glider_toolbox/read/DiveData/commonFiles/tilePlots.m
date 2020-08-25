% ****************************************************************************************************
%
% ****************************************************************************************************
function tilePlots(figureHandles, tile, border)
    maxpos  = get (0,'screensize');     % determine terminal size in pixels

    figureHandles = sort(figureHandles);
    numfigs       = length(figureHandles);
    
    maxfigs = 100;

    if (numfigs > maxfigs)
        disp(['Unable to to display more then  ' num2str(maxfigs) ' figures.']);
    else
        if nargin == 0
            uiwait(msgbox('ERR: Handles to figures was not passed to function!', 'CreateMode', 'modal'));
            return;
        elseif (nargin == 1) || ((nargin >= 2) && (isempty(tile) == 1))
            maxfactor = sqrt(maxfigs);       % max number of figures per row or column
            sq = [1:maxfactor].^2;           % vector of integer squares
            sq = sq(find(sq>=numfigs));      % determine square grid size
            gridsize = sq(1);                % best grid size
            nrows = sqrt(gridsize);          % figure size screen scale factor
            ncols = ceil(numfigs/nrows);                   % figure size screen scale factor
        elseif nargin >= 2
            nrows = tile(1);
            ncols = tile(2);
            if numfigs > nrows*ncols
                disp (['Requested tile size too small for ' num2str(numfigs) ' open figures ']);
                return;
            end
        end

        if nargin < 2
            border = 0;
        else
            maxpos(3) = maxpos(3) - 2*border;
            maxpos(4) = maxpos(4) - 2*border;
        end
        xlen = fix(maxpos(3)/ncols) - 30; % new tiled figure width
        ylen = fix(maxpos(4)/nrows) - 90; % subtract off window header

        % Location (1,1) is at bottom left corner
        pnum=0;
        for iy = nrows-1:-1:0
            ypos = fix((iy)*(maxpos(4)/nrows)) + border;
            for ix = 0:ncols-1
                xpos = fix(ix*(maxpos(3)/ncols) + 1) + border+7;     % figure location (column)
                pnum = pnum+1;
                if (pnum>numfigs)
                    break
                else
                    figure(figureHandles(pnum))
                    set(figureHandles(pnum),'Position',[ xpos ypos xlen ylen ]); % move figure
                end
            end
        end
    end
end
