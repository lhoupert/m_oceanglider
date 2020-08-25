function [gps_prepredive, gps_predive, gps_pastdive] = gt_sg_sub_parse_GPS(log)
%
% [gps_prepredive, gps_predive, gps_pastdive] = gt_sg_sub_parse_GPS(log)
% Parses Seaglider GPS coordinates present in log files.
%
% Inputs:
% log = the .log substructure.
%
% Outputs:
% gps_prepredive, gps_predive, gps_pastdive =
%    n_dives x 3 matrices containing lat, long and time.
%
% pastdive    :  GPS   :  coordinates upon surfacing at the end of the dive (most recent)
% prepredive  :  GPS1  :  corresponds to gps_pastdive of the previous dive (redundant)
% predive     :  GPS2  :  last coordinates before diving
%
% B.Y.QUESTE Feb 2015

arrayLength = length(log);
gps_predive = nan(arrayLength,3);
gps_pastdive = nan(arrayLength,3);
gps_prepredive = nan(arrayLength,3);

% Decoding functions
decode_GPS = @(x) sign(x).*(fix(abs(x)/100)+ mod(abs(x),100)./60);
decode_time = @(x,y) datenum(mod(x,100) + 2000,floor(mod(x,10000)/100),floor(x/10000),...
    floor(y./10000),floor(mod(y,10000)/100),mod(y,100));

% Decode standard formatting for each dive
gt_sg_sub_echo({'Parsing GPS coordinates from log data.',...
    'gps_predive: last coordinates before diving',...
    'gps_postdive: coordinates upon surfacing at the end of the dive (most recent)',...
    'gps_postpreviousdive: corresponds to gps_pastdive of the previous dive (redundant)'});

% Occasionally, empy GPS strings are returned for one of them, so ignore
% those but do the rest.
processedDives  = intersect([log.dive],find(~cellfun(@isempty,{log.GPS})));
processedDives1 = intersect([log.dive],find(~cellfun(@isempty,{log.GPS1})));
processedDives2 = intersect([log.dive],find(~cellfun(@isempty,{log.GPS2})));

% 1. Date
% 2. Time
% 3. Lat
% 4. Lon
% 5. Time to first fix in seconds ?
% 6. Horizontal dilution of precision ?
% 7. Total time to aquire fix ?
% 8. ?
gps_pastdive(processedDives,:) = [arrayfun(@(x) decode_GPS(x.GPS(3)),log(processedDives)),...
    arrayfun(@(x) decode_GPS(x.GPS(4)),log(processedDives)),...
    arrayfun(@(x) decode_time(x.GPS(1),x.GPS(2)),log(processedDives))];

gps_prepredive(processedDives1,:) = [arrayfun(@(x) decode_GPS(x.GPS1(3)),log(processedDives1)),...
    arrayfun(@(x) decode_GPS(x.GPS1(4)),log(processedDives1)),...
    arrayfun(@(x) decode_time(x.GPS1(1),x.GPS1(2)),log(processedDives1))];

gps_predive(processedDives2,:) = [arrayfun(@(x) decode_GPS(x.GPS2(3)),log(processedDives2)),...
    arrayfun(@(x) decode_GPS(x.GPS2(4)),log(processedDives2)),...
    arrayfun(@(x) decode_time(x.GPS2(1),x.GPS2(2)),log(processedDives2))];

% In one version (where iRobot messed up), they forgot to include the date,
% in this case everything is shifted up by one - BUT ONLY FOR GPS1 and GPS2
% Test to see if latititudes/longitudes are valid, if not, assume it's this
% bad version.
if any(gps_predive(:,1) < -80 | gps_predive(:,1) > 90) || mode([log.version]) == 6606
    gt_sg_sub_echo({'Parsing old-style (v.66.06) GPS string format','Please verify gps_ arrays for errors.'});
    
    extractDate = @(x) x(2)*10000 + x(1)*100 + mod(x(3)+1900,100);
    
    gps_prepredive(processedDives,:) = [arrayfun(@(x) decode_GPS(x.GPS1(2)),log(processedDives)),...
        arrayfun(@(x) decode_GPS(x.GPS1(3)),log(processedDives)),...
        arrayfun(@(x) decode_time(nan,x.GPS1(1)),log(processedDives))];
    
    gps_predive(processedDives,:) = [arrayfun(@(x) decode_GPS(x.GPS2(2)),log(processedDives)),...
        arrayfun(@(x) decode_GPS(x.GPS2(3)),log(processedDives)),...
        arrayfun(@(x) decode_time(extractDate(x.start),x.GPS2(1)),log(processedDives))];
    
    for istep=setxor(1,processedDives)
        gps_prepredive(istep,3) = gps_pastdive(istep-1,3);
    end
    
    for istep=1:numel(log)
        if gps_predive(istep,3) > gps_pastdive(istep,3)
            gt_sg_sub_echo({'Incorrect timestamp identified for predive GPS signal due to v.66.06 format issue.',['Attempting to autocorrect for dive ',num2str(istep),'.']});
            gps_predive(istep,3) = gps_predive(istep,3)-1;
        end
    end
end
end