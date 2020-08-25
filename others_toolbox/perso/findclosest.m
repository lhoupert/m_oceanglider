function [ indX,indY,distmat] = findclosest(timeX,lonX,latX,timeY,lonY,latY,timelim,distlim)
% [indX,indY,distmat] = findclosest(timeX,lonX,latX,timeY,lonY,latY,timelim,distlim)
% Function that find the closest points in space between two data set. If
% the closest points are not close enough according to a time and space threshold
% (defined by timelim (in days) and distlim (in km)), indX and indY are empty. 
%
% created by L. Houpert (houpertloic@gmail.com), 10/04/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox

% X = [ [54    55    56] ; [10 9.5 8]]' ;
% Y = [ [21    14    51   223   335    54    53    22]; [11 11 12 13 14 9 0 -4]]';
X0 = [ latX(:)' ; lonX(:)' ]';
Y0 = [ latY(:)' ; lonY(:)' ]';
timeX = timeX(:);
timeY = timeY(:);

% convert degree lat lon into distance
% def ref point:
reflat = nanmean([X0(:,1); Y0(:,1)]);
reflon = nanmean([X0(:,2); Y0(:,2)]);
if isempty(X0), indX=[];indY=[];distmat=[];return;end
X=nan(size(X0(:,1)));
for ijk=1:length(X0(:,1))
    [X(ijk),angleX(ijk)] = sw_dist([X0(ijk,1) reflat],[X0(ijk,2) reflon],'km');
end
X(angleX<-90 | angleX>90) = -X(angleX<-90 | angleX>90);
Y=nan(size(Y0(:,1)));
for ijk=1:length(Y0(:,1))
    [Y(ijk),angleY(ijk)] = sw_dist([Y0(ijk,1) reflat],[Y0(ijk,2) reflon],'km');
end
Y(angleY<-90 | angleY>90) = -Y(angleY<-90 | angleY>90);
% figure;subplot(2,1,1);plot(X,'ro');hold on; plot(Y,'bs');subplot(2,1,2);plot(timeX,'ro');hold on; plot(timeY,'bs')
 

DD = pdist2(X,Y,'euclidean');
DDtime = pdist2(timeX,timeY,'euclidean');
    
[ixy] = find(DDtime > timelim);
DD(ixy) = nan;

[indX0,indY0] = find(DD < distlim);

distmat.time  = DDtime;
distmat.space = DD;
%figure; surf(DD)
%figure; plot(nanmean(DD,2))

if ~isempty(indX0)

    uniqiX0 = unique(indX0);
    indlimmin = [uniqiX0(1) ; uniqiX0(find(diff(timeX(uniqiX0))>1)+1)];
    indlimmax = [uniqiX0(find(diff(timeX(uniqiX0))>1)); uniqiX0(end)];

    indX=cell(1,length(indlimmin));
    indY=cell(1,length(indlimmin));
    % find way to extract all the index in subgroup
    for ill=1:length(indX)
        iXsel = (indX0 >= indlimmin(ill) &  indX0 <= indlimmax(ill));
        indX{ill} = indX0(iXsel);
        indY{ill} = indY0(iXsel);
    end
    
else
    indX = [];
    indY = [];
end
