function [indX,indY,distXY,timediffXY] = findcross_space(timeX,lonX,latX,timeY,lonY,latY,timelim,distlim)
% [indX,indY,distXY] = findcross_space(timeX,lonX,latX,timeY,lonY,latY,timelim,distlim)
% Function that find the closest point in space between two data set. If
% the closest point is not close enough according to a time and space threshold
% (defined by timelim (in days) and distlim (in km)), indX and indY are empty. 

% X = [ [54    55    56] ; [10 9.5 8]]' ;
% Y = [ [21    14    51   223   335    54    53    22]; [11 11 12 13 14 9 0 -4]]';

X = [ latX(:)' ; lonX(:)' ]';
Y = [ latY(:)' ; lonY(:)' ]';


timeX = timeX(:);
timeY = timeY(:);

% figure;subplot(2,1,1);plot(X,'ro');hold on; plot(Y,'bs');subplot(2,1,2);plot(timeX,'ro');hold on; plot(timeY,'bs')
 
DD = pdist2(X,Y,'euclidean');
DDtime = pdist2(timeX,timeY,'euclidean');
    
indbad = find(DDtime > timelim);
DD(indbad) = nan;

[indX0,indY0] = find(DD==min(DD(:)));

if isempty(indX0)
    indX = [];
    indY = [];   
    distXY = [];
    timediffXY = [];    
else
    indX0=indX0(1);
    indY0=indY0(1);
    [distXY,phaseangleXY] = sw_dist([X(indX0,1) Y(indY0,1)],[X(indX0,2) Y(indY0,2)],'km');
    timediffXY = abs(timeX(indX0)-timeY(indY0));
    if distXY < distlim
        indX = indX0;
        indY = indY0;
    else
        indX = [];
        indY = [];
        distXY = [];
        timediffXY = [];
    end

%    figure;plot(X,'ro');hold on; plot(Y,'bs');plot(X(indX,:),'m+');plot(Y(indY,:),'kx')

end

end
