function [dcl,dof] = decorrelation3(x,y,diag) ;
%% Compute auto-correlation of a variable; determine the weighted integral scale; 
% compute degrees of freedom in data series
% Compute the integral time scale and degrees of freedom in a timeseries
% Emery, W. J. and R. E. Thomson, 1997: Data analysis methods in physical
% oceanography. Pergamon, 634 pp. (see p263,265,379-380) define the
% integral time scale from C(0) out to the first zero crossing. Here we
% define the integral time scale as twice this value.
% Integral time scale = 2xsum(C(tau)dt) from C(0) to first zero crossing.
% If the autocorrelation decays linearly then the area under this triangle
% is 0.5 * distance to zero crossing. So twice this integral is equal to the 
% time to the zero crossing.
% If the correlation decays fast initially but more slowly later the zero
% crossing suggests a longer time than the sum which is really just a
% weighted estimate and in this case gives less weight to a long tail with
% low correlation.
%
%
% USAGE : function [dcl,dof] = decorrelation3(x,y,diag)
% x is the temporal/spatial variable ; 
% y is the variable and is normalised prior to computing the autocorrelation
% Diagnostic diag = 1/0 for diagnostic plots.
% dcl is the decorrelation length scale in the units of x
% dof is the number of degrees of freedom in x. Calculated by length(x) /
% dcl. Stuart Cunningham, July 2017
% Aug 2018: LOH fix some bugs in the code, define minimum of the
% autocorrelation fonction, decorrelation length is express as a function
% of timestep.
% July 2019: LOH add significant level for autocorrelation


doplot = diag;

x(isnan(y))=[];
y(isnan(y))=[];
% First normalise the variable
ynorm = (y - nanmean(y)) / nanstd(y);
maxlag = length(ynorm); % Set maxlag = length of the variable
[C,lags] = xcorr(ynorm,maxlag,'coeff'); % compute normalised correlation coefficient

%------------------------------------------------------------------------------------------------
% upper rejection limit for testing (under white noise assumption) each autocorrelation = 0.
% Degree of freedom for each lag
lags_dof = length(ynorm)-abs(lags);
% Compute the inverse of Student's T cumulative distribution function for a pvalue of 0.01 or 0.001 if large dataset
if length(ynorm)>100
    alpha = 0.001;
    
else
    alpha = 0.01;
    
end
tcoef = tinv(1-alpha,lags_dof); 
rejuplim = tcoef.*(1./sqrt(lags_dof));
%rejuplim = 3.*(1./sqrt(lags_dof));
% figure; plot(lags,rejuplim);hold on; plot(lags,-rejuplim)
% %rejuplim =  (1.96)*(1/sqrt(length(y)))*ones(1,2); % lower limit is x(-1)

[X0,Y0] = intersections(lags,C,lags,rejuplim); % use functions to find the intersection of correlation with line y=0;

X0lag = (X0(((length(X0)/2)+1))); % find lag value of first negative crossing
[~,Imax]=min(abs(lags - X0lag)) ; % find index of first negative crossing
[~,Imin]=min(abs(lags - -X0lag));% find index of first positive crossing

dcl = trapz(lags(Imin:Imax),C(Imin:Imax)); % Integrate correlation between first neg and first pos crossing. This is the decorrelation length
dof= length(ynorm) ./ dcl ; % Degrees of freedom = length of y / dcl

% disp(['x has',sprintf('%5.0f',length(x)),' data cycles']);
% disp(['Integral time scale = ',sprintf('%2.1f',dcl) ' * unit of time',...
%     ' : Degrees of freedom = ',sprintf('%2.0f',dof)])

if doplot;
    figure;clf;
    subplot(211);hold on ;grid on;
    plot(x,ynorm,'k');
    title(['x has ',sprintf('%5.1f',length(x)),' data cycles']);
    xlabel('x variable');ylabel('Normalised y variable');
    
    subplot(212);hold on;grid on;
p1= plot(lags,C);
%     plot([lags(1) lags(end)], rejuplim,'--b')
%     plot([lags(1) lags(end)], -rejuplim,'--b')    
p2= plot(lags, rejuplim,'--r')
    plot(lags, -rejuplim,'--r')      
p3= plot(lags([Imin Imax]),C([Imin Imax]),'ko','MarkerSize',10,'MarkerFaceColor','k','MarkerSize',5);
    xlim([lags(Imin-round(length(ynorm)/4)) lags(Imax+round(length(ynorm)/4))]);
    xlabel('Lags');title('Normalised auto-correlation of ynorm');ylabel('Correlation coefficient');
    legend([p1,p2, p3],{'autocorrelation function',[num2str((1-alpha)*100,'%2.1f') '% confidence interval'],'decorrelation length scale'})
end

end
