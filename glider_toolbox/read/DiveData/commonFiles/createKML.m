function createKML(handles, filename, name, sg_id, diveInfoList, iconColor)

    showSurfaceDrift = 1;
    
    % Color of pushpin icon and track line.
    % strColor format is AAbbggrr, where
    %  AA = alpha blending
    %  bb = blue value
    %  gg = green value
    %  rr = red value
    iconRoot = 'http://maps.google.com/mapfiles/kml/paddle/';
    
    iconPathDn = strcat(iconRoot, 'wht-blank.png');
    iconPathUp = strcat(iconRoot, 'wht-stars.png');

    % color will be applied to the white paddle icons
    % as well as to the track lines, targets, etc.
    if (iconColor == 1)
        %           AAbbggrr
        strColor = 'ff0000ff'; % red
    elseif (iconColor == 2)
        strColor = 'ff00ffff'; % yellow
    elseif (iconColor == 3)
        strColor = 'ff00aa00'; % green
    elseif (iconColor == 4)
        strColor = 'ffffffff'; % white        
    elseif (iconColor == 5)
        strColor = 'ffffff00'; % lt blue
    elseif (iconColor == 6)
        strColor = 'ffff00aa'; % purple
    elseif (iconColor == 7)
        strColor = 'ffff2222'; % blue
    elseif (iconColor == 8)
        strColor = 'ffffaaff'; % pink
    else
        strColor = 'ff0000ff'; % red
    end

    fid = fopen(filename, 'wt');
    if fid == -1
        vizMsgBox(sprintf('Unable to open output .kml file\n(%s)\n', filename), 'ERROR', 'modal');
        return;
    end

    % ------------------------------------------------------------------------------------------------------------------------
    % xml header info
    % ------------------------------------------------------------------------------------------------------------------------
    header = [];
    header{end+1}  = ['<?xml version="1.0" encoding="UTF-8"?>' 10];
    header{end+1}  = ['<kml xmlns="http://earth.google.com/kml/2.0"' 10];
    header{end+1}  = ['xmlns:gx="http://www.google.com/kml/ext/2.2">' 10];
    header{end+1}  = ['<Document>' 10];
    header{end+1}  = ['<name>' name '</name>' 10];

    insertXML(fid, header);


    % ------------------------------------------------------------------------------------------------------------------------
    % Style for showing icons where SeaGlider returns to surface.
    % ------------------------------------------------------------------------------------------------------------------------
    header = [];
    
    header{end+1} = ['<Style id="sg_up">' 10];
    header{end+1} = ['    <IconStyle>' 10];
    header{end+1} = ['        <colorMode>normal</colorMode>' 10];
    header{end+1} = ['        <color>' strColor '</color>' 10];
    header{end+1} = ['        <scale>1.0</scale>' 10];
    header{end+1} = ['        <Icon>' 10];
    header{end+1} = ['            <href>' iconPathUp '</href>' 10];
    header{end+1} = ['        </Icon>' 10];
    header{end+1} = ['    </IconStyle>' 10];
    header{end+1} = ['    <BalloonStyle>' 10];
    header{end+1} = ['        <text>$[description]</text>' 10];
    header{end+1} = ['    </BalloonStyle>' 10];
    %header{end+1} = ['    <LineStyle>' 10];
    %header{end+1} = ['        <color>' strColor '</color>' 10];
    %header{end+1} = ['        <colorMode>normal</colorMode>' 10];
    %header{end+1} = ['        <width>3</width>' 10];
    %header{end+1} = ['        <gx:labelVisibility>0</gx:labelVisibility>' 10];
    %header{end+1} = ['    </LineStyle>' 10];
    header{end+1} = ['</Style>' 10];

    insertXML(fid, header);
    fprintf(fid, '\n\n');

    % ------------------------------------------------------------------------------------------------------------------------
    % Style for showing icons where SeaGlider dives from surface.
    % ------------------------------------------------------------------------------------------------------------------------
    header = [];
    
    header{end+1} = ['<Style id="sg_down">' 10];
    header{end+1} = ['    <IconStyle>' 10];
    header{end+1} = ['        <colorMode>normal</colorMode>' 10];
    header{end+1} = ['        <color>' strColor '</color>' 10];
    header{end+1} = ['        <scale>1.0</scale>' 10];
    header{end+1} = ['        <Icon>' 10];
    header{end+1} = ['            <href>' iconPathDn '</href>' 10];
    header{end+1} = ['        </Icon>' 10];
    header{end+1} = ['    </IconStyle>' 10];
    header{end+1} = ['    <BalloonStyle>' 10];
    header{end+1} = ['        <text>$[description]</text>' 10];
    header{end+1} = ['    </BalloonStyle>' 10];
    header{end+1} = ['</Style>' 10];

    insertXML(fid, header);
    fprintf(fid, '\n\n');

    % ------------------------------------------------------------------------------------------------------------------------
    % Style for showing SeaGlider track.
    % ------------------------------------------------------------------------------------------------------------------------
    header = [];
    
    header{end+1} = ['<Style id="sg_track">' 10];
    header{end+1} = ['    <LineStyle>' 10];
    header{end+1} = ['        <color>' strColor '</color>' 10];
    header{end+1} = ['        <colorMode>normal</colorMode>' 10];
    header{end+1} = ['        <width>3</width>' 10];
    header{end+1} = ['        <gx:labelVisibility>0</gx:labelVisibility>' 10];
    header{end+1} = ['    </LineStyle>' 10];
    header{end+1} = ['</Style>' 10];

    insertXML(fid, header);
    fprintf(fid, '\n\n');

    % ------------------------------------------------------------------------------------------------------------------------
    % Style for showing location of 'targets'.
    % ------------------------------------------------------------------------------------------------------------------------
    header = [];
    
    header{end+1} = ['<Style id="sg_target">' 10];
    header{end+1} = ['    <IconStyle>' 10];
    header{end+1} = ['        <colorMode>normal</colorMode>' 10];
    header{end+1} = ['        <color>' strColor '</color>' 10];
    header{end+1} = ['        <scale>2.0</scale>' 10];
    header{end+1} = ['        <Icon>' 10];
    header{end+1} = ['            <href>' 'http://maps.google.com/mapfiles/kml/shapes/target.png' '</href>' 10];
    header{end+1} = ['        </Icon>' 10];
    header{end+1} = ['    </IconStyle>' 10];
    header{end+1} = ['    <BalloonStyle>' 10];
    header{end+1} = ['        <text>$[description]</text>' 10];
    header{end+1} = ['    </BalloonStyle>' 10];
    header{end+1} = ['</Style>' 10];

    insertXML(fid, header);
    fprintf(fid, '\n\n');


    % ------------------------------------------------------------------------------------------------------------------------
    % Place push pin at position where seaglider left surface and returned to surface for each dive.
    % ------------------------------------------------------------------------------------------------------------------------
    for i=1:length(diveInfoList)
        
        % Add pushpin showing where we left the surface for this dive.
        if (diveInfoList(i).dove.pos.valid ~= 1)
            errStr = sprintf('Unable to plot position where dive %d left surface - no valid GPS fix', diveInfoList(i).diveNum);
            disp(errStr);
            vizMsgBox(errStr, 'ERROR', 'modal');
            updateResults(handles, errStr);
        else
            [dat, tim] = unixTimeToString(diveInfoList(i).dove.time);

            %titleStr = sprintf('<h2>SG%03d - Dive %d</h2>', sg_id, diveInfoList(i).diveNum);
            titleStr = sprintf('<h2>Dive %d - Start</h2>', diveInfoList(i).diveNum);
            %titleStr = sprintf('<h2>Start</h2>');
            idStr = sprintf('<h3>SG%03d</h3>', sg_id);
            dateStr = sprintf('<b>Date:</b> %s', dat);
            timeStr = sprintf('<br><b>Time:</b> %s', tim);
            latStr = sprintf('<p><b>Lat:</b> %.6f', diveInfoList(i).dove.pos.lat);
            lonStr = sprintf('<br><b>Lon:</b> %.6f', diveInfoList(i).dove.pos.lon);
            tgtStr = sprintf('<p><b>Target:</b> %s', char(diveInfoList(i).target.name));
            descrip = [titleStr idStr dateStr timeStr latStr lonStr tgtStr];
            
            %[errStr] = errorListToString(diveInfoList(i).errors);
            %descrip = strcat(descrip, errStr);
            
            header = [];
            header{end+1} = ['<Placemark>' 10];
            header{end+1} = ['    <name>' sprintf('%d', diveInfoList(i).diveNum) '</name>' 10];
            header{end+1} = ['    <description> <![CDATA[' descrip ']]></description>' 10];
            header{end+1} = ['    <styleUrl>#sg_down</styleUrl>' 10];
            header{end+1} = ['    <Point>' 10];
            header{end+1} = sprintf('        <coordinates>%.6f,%.6f,0.0</coordinates>', diveInfoList(i).dove.pos.lon, diveInfoList(i).dove.pos.lat);
            header{end+1} = [' ' 10];
            header{end+1} = ['    </Point>' 10];
            header{end+1} = ['</Placemark>' 10];
            header{end+1} = [' ' 10];

            insertXML(fid, header);
            fprintf(fid, '\n\n');
        end

        % Add pushpin showing where we surfaced after the dive.
        if (diveInfoList(i).surfaced.pos.valid ~= 1)
            errStr = sprintf('Unable to plot position where dive %d surfaced - no valid GPS fix', diveInfoList(i).diveNum);
            disp(errStr);
            vizMsgBox(errStr, 'ERROR', 'modal');
            updateResults(handles, errStr);
            continue;
        else        
            [dat, tim] = unixTimeToString(diveInfoList(i).surfaced.time);

            %titleStr = sprintf('<h2>SG%03d - Dive %d</h2>', sg_id, diveInfoList(i).diveNum);
            titleStr = sprintf('<h2>Dive %d - End</h2>', diveInfoList(i).diveNum);
            %titleStr = sprintf('<h2>End</h2>');
            idStr = sprintf('<h3>SG%03d</h3>', sg_id);
            dateStr = sprintf('<b>Date:</b> %s', dat);
            timeStr = sprintf('<br><b>Time:</b> %s', tim);
            latStr = sprintf('<p><b>Lat:</b> %.6f', diveInfoList(i).surfaced.pos.lat);
            lonStr = sprintf('<br><b>Lon:</b> %.6f', diveInfoList(i).surfaced.pos.lon);
            tgtStr = sprintf('<p><b>Target:</b> %s', char(diveInfoList(i).target.name));
            descrip = [titleStr idStr dateStr timeStr latStr lonStr tgtStr];

            [errStr] = errorListToString(diveInfoList(i).errors);
            descrip = strcat(descrip, errStr);

            header = [];
            header{end+1} = ['<Placemark>' 10];
            header{end+1} = ['    <name>' sprintf('%d', diveInfoList(i).diveNum) '</name>' 10];
            header{end+1} = ['    <description> <![CDATA[' descrip ']]></description>' 10];
            header{end+1} = ['    <styleUrl>#sg_up</styleUrl>' 10];
            header{end+1} = ['    <Point>' 10];
            header{end+1} = sprintf('        <coordinates>%.6f,%.6f,0.0</coordinates>', diveInfoList(i).surfaced.pos.lon, diveInfoList(i).surfaced.pos.lat);
            header{end+1} = [' ' 10];
            header{end+1} = ['    </Point>' 10];
            header{end+1} = ['</Placemark>' 10];
            header{end+1} = [' ' 10];

            insertXML(fid, header);
            fprintf(fid, '\n\n');
        end
    end
    
    % ------------------------------------------------------------------------------------------------------------------------
    % Draw lines connecting the push pins.
    % ------------------------------------------------------------------------------------------------------------------------
    header = [];
    header{end+1} = ['<Placemark>' 10];
    header{end+1} = ['    <name>Track</name>' 10];
    header{end+1} = ['    <styleUrl>#sg_track</styleUrl>' 10];
    header{end+1} = ['    <LineString>' 10];
    header{end+1} = ['        <tessellate>1</tessellate>' 10];
    header{end+1} = ['        <altitudeMode>relativeToGround</altitudeMode>' 10];
    header{end+1} = ['        <coordinates>' 10];
    
    for i=1:length(diveInfoList)
        % point where glider dove for this segment
        if (diveInfoList(i).dove.pos.valid == 1)
            header{end+1}  = [sprintf('            %.6f,%.6f,0', diveInfoList(i).dove.pos.lon, diveInfoList(i).dove.pos.lat) 10];
        end

        % point where glider surfaced for this segment
        if (diveInfoList(i).surfaced.pos.valid == 1)
            % Google Earth does not like spaces between lon,lat,altitude.
            header{end+1}  = [sprintf('            %.6f,%.6f,0', diveInfoList(i).surfaced.pos.lon, diveInfoList(i).surfaced.pos.lat) 10];
        end
    end
    
    header{end+1}=['        </coordinates>' 10];
    header{end+1}=['    </LineString>' 10];
    header{end+1}=['</Placemark>' 10];
    
    insertXML(fid, header);
    fprintf(fid, '\n\n');
    
    % ------------------------------------------------------------------------------------------------------------------------
    % Draw targets (waypoints).
    % ------------------------------------------------------------------------------------------------------------------------
    targetsDrawn = [];
    for i=1:length(diveInfoList)
        tgt = char(diveInfoList(i).target.name);
        % Only plot each target once.
        found = 0;
        for j=1:length(targetsDrawn)
            %if (char(targetsDrawn(j)) == char(tgt))
            if (length(tgt) == length(targetsDrawn(j)) && char(tgt) == char(targetsDrawn(j)))
                found = 1;
                break;
            end
        end
        
        if (found == 0)
            tgtName = char(diveInfoList(i).target.name);
            coords = diveInfoList(i).target.location;
            radius = diveInfoList(i).target.radius;
            
            % lat
            degrees = fix(coords(1)/100);
            minutes = rem(coords(1),100);
            lat = degrees + minutes / 60.0;

            % lon
            degrees = fix(coords(2)/100);
            minutes = rem(coords(2),100);
            lon = degrees + minutes / 60.0;

            descrip = sprintf('<h2>Target - %s</h2><b>Lat:</b> %.6f<br><b>Lon:</b> %.6f<p><b>Radius:</b> %.0f m<p>', tgtName, lat, lon, radius);

            header = [];
            header{end+1} = ['<Placemark>' 10];
            header{end+1} = ['    <name>' sprintf('%s', tgtName) '</name>' 10];
            header{end+1} = ['    <description> <![CDATA[' descrip ']]></description>' 10];
            header{end+1} = ['    <styleUrl>#sg_target</styleUrl>' 10];
            header{end+1} = ['    <Point>' 10];
            
            header{end+1} = sprintf('        <coordinates>%.6f,%.6f,0.0</coordinates>', lon, lat);
            header{end+1} = [' ' 10];
            header{end+1} = ['    </Point>' 10];
            header{end+1} = ['</Placemark>' 10];
            header{end+1} = [' ' 10];

            insertXML(fid, header);
            fprintf(fid, '\n\n');
            
            targetsDrawn{end+1} = [tgtName];
        end
    end

    % ------------------------------------------------------------------------------------------------------------------------
    %
    % ------------------------------------------------------------------------------------------------------------------------
    footer='</Document>\n</kml>';
    fprintf(fid, '%s', footer);

    fclose(fid)
end

% ************************************************************************************************************************
% Create date/time string from unix time.
% ************************************************************************************************************************
function [errStr2] = errorListToString(errorList)
    errHdr = 0;
	errStr = sprintf('<p><b>Errors:</b>');
    for errNum=1:length(errorList)
        %if (errorList(errNum) > 0)
        % If a single GPS error occurred no need to report it
        if (errNum == 15 && errorList(errNum) > 1) || (errNum ~= 15 && errorList(errNum) > 0)
            if errHdr == 0
                errStr = strcat(errStr, sprintf('<font color="red">'));
                errHdr = 1;
            end
            if errNum == 16
                errStr = strcat(errStr, sprintf('<br>%d GPS PPS errors', errorList(errNum)));
            elseif errNum == 15
                errStr = strcat(errStr, sprintf('<br>%d GPS fix errors', errorList(errNum)));
            elseif errNum == 14
                errStr = strcat(errStr, sprintf('<br>%d VBD retries', errorList(errNum)));
            elseif errNum == 13
                errStr = strcat(errStr, sprintf('<br>%d roll retries', errorList(errNum)));
            elseif errNum == 12
                errStr = strcat(errStr, sprintf('<br>%d pitch retries', errorList(errNum)));
            elseif errNum == 11
                errStr = strcat(errStr, sprintf('<br>%d VBD errors', errorList(errNum)));
            elseif errNum == 10
                errStr = strcat(errStr, sprintf('<br>%d roll errors', errorList(errNum)));
            elseif errNum == 9
                errStr = strcat(errStr, sprintf('<br>%d pitch errors', errorList(errNum)));
            elseif errNum == 8
                errStr = strcat(errStr, sprintf('<br>%d CF retries on close', errorList(errNum)));
            elseif errNum == 7
                errStr = strcat(errStr, sprintf('<br>%d CF retries on write', errorList(errNum)));
            elseif errNum == 6
                errStr = strcat(errStr, sprintf('<br>%d CF retries on open', errorList(errNum)));
            elseif errNum == 5
                errStr = strcat(errStr, sprintf('<br>%d CF errors on close', errorList(errNum)));
            elseif errNum == 4
                errStr = strcat(errStr, sprintf('<br>%d CF errors on write', errorList(errNum)));
            elseif errNum == 3
                errStr = strcat(errStr, sprintf('<br>%d CF errors on open', errorList(errNum)));
            elseif errNum == 2
                errStr = strcat(errStr, sprintf('<br>%d spurious interrupts', errorList(errNum)));
            elseif errNum == 1
                errStr = strcat(errStr, sprintf('<br>%d log buffer overruns', errorList(errNum)));
            else
                errStr = strcat(errStr, sprintf('<br>Unexpected error index %d', errorList(errNum)));
            end
        end
    end
    if errHdr == 0
        [errStr2] = strcat(errStr, sprintf('<font color="green"> None!</font>'));
    else
        [errStr2] = strcat(errStr, '</font>');
    end
end

% ************************************************************************************************************************
% Create date/time string from unix time.
% ************************************************************************************************************************
function [dateStr timeStr] = unixTimeToString(startPosTime)
    dateval = unix_to_datenum(startPosTime);
	start_vec = datevec(dateval);
    dateStr = sprintf('%02d/%02d/%d', start_vec(2), start_vec(3), start_vec(1));   % mon,d,y
    timeStr = sprintf('%02d:%02d:%02d', start_vec(4), start_vec(5), start_vec(6)); % h,min,s
end

% ************************************************************************************************************************
%
% ************************************************************************************************************************
function insertXML(fid, header)
    newHeader = [];

    [m, n] = size(header);

    for i=1:m
        newHeader = sprintf('%s\n%s', newHeader, char(cell2mat(header(i,:))));
    end

    fprintf(fid, '%s', strtrim(newHeader));
end