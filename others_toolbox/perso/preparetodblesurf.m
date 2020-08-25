function [xxxx yyyy zzzz]=preparetodblesurf(xini,yini,zz,ddx,ddy)
%function [xxxx yyyy zzzz]=preparetodblesurf(xini,yini,zz,ddx,ddy)
% double the two dimensions of a matrix and duplicate the values in order
% to not loose data point in the interpolation of surf
%
% created by L. Houpert (houpertloic@gmail.com), 12/02/2015, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox
% 
% Then a good way to make this toolbox evolve is to create a new branch, commit the changes and push the changes on the git remote repository
% (https://bitbucket.org/Lhoupert/oceano_data_toolbox)

nx=length(nanmean(xini,1));
xxxx1(:,1:2:2*nx) = xini-ddx;
xxxx1(:,2:2:2*nx) = xini+ddx;
zzzz1(:,1:2:2*nx) = zz;
zzzz1(:,2:2:2*nx) = zz;
yyyy1(:,1:2:2*nx) = yini;
yyyy1(:,2:2:2*nx) = yini;

ny=length(nanmean(yini,2));
xxxx(1:2:2*ny,:) = xxxx1;
xxxx(2:2:2*ny,:) = xxxx1;
zzzz(1:2:2*ny,:) = zzzz1;
zzzz(2:2:2*ny,:) = zzzz1;
yyyy(1:2:2*ny,:) = yyyy1-ddy;
yyyy(2:2:2*ny,:) = yyyy1+ddy;   
end

