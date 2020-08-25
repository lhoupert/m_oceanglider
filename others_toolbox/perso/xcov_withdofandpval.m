function [xcf,lags,pval,dofy1,dofy2] = xcov_withdofandpval(y1,y2,maxlag)
% function [xcf,lags,pval,dofy1,dofy2] = xcov_withdofandpval(y1,y2,maxlag)

lags=-maxlag:1:maxlag;
xcf = nan(size(lags));
dofy1 = nan(size(lags));
dofy2 = nan(size(lags));
pval = nan(size(lags));

for cc=1:length(lags)
    ij = lags(cc);
    y1b = y1(:);
    y2b = y2(:);
        if ij<0
%             y2b(1:-ij)=[];
%             y1b((end+ij)+1:end)=[]; 
            y2b = [ y2b' nan(1,-ij) ];
            y1b = [ nan(1,-ij) y1b' ];    
       
        else 
%             y1b(1:ij)=[];
%             y2b((end-ij)+1:end)=[];       
            y2b = [ nan(1,ij) y2b'];
            y1b = [ y1b' nan(1,ij)];                 
        end
        iok = find(~isnan(y1b) & ~isnan(y2b));
        y1b = y1b(iok);
        y2b = y2b(iok);   
        [~,dof1]=decorrelation3(1:length(y1b),y1b,0);
        [~,dof2]=decorrelation3(1:length(y2b),y2b,0);    
        dofy1(cc)= dof1;
        dofy2(cc)= dof2;        
        
        % calculate crosscovarPval
        [R,P] = corrcoef_timeseries(y1b,y2b);
     
        xcf(cc) = R(2);
        pval(cc) = P(2);
end