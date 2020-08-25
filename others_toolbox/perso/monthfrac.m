function[mf]=monthfrac(num)
%MONTHFRAC  Converts a DATENUM into 'month.fraction'.
%  
%   [MF]=monthfrac(NUM) returns MF, the fraction of the current 
%   month at each date.  FLOOR(MF) is the standard month number. 
%   
%   Recently the function monthfrac disappear from jlab. This function was 
%   created to insure a continuity between the older versions of the jlab toolbox
%   and the new yearfrac function.
%
%   Usage: mf=monthfrac(num);
%
%   L. Houpert, August 2016

      
if iscell(num)
    for i=1:length(num)
        if ~isempty(num{i})
            [yf{i,1},mf{i,1}]=yearfrac(num{i});
        else
            yf{i}=[];
            mf{i}=[];
        end
    end
else
    [yf,mf]=yearfrac(num);
end

  
