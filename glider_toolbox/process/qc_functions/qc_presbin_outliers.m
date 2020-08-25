 function [indexbad,datavecmean,datavecstd,binpres] = qc_presbin_outliers(data,pres,qc_param)
% function [indexbad,datavecmean,datavecstd,binpres] = qc_presbin_outliers(data,pres,qc_param)
% Function to detect outliers in pressure bin 
%   input: data : data to check
%          pres : pressure associated to the data to test
%          qc_param.varname: name of the QC variable (for display purpose) 
%          qc_param.binsize : size of the pressure bin (m) for the mean and std calculation 
%          qc_param.tollvl : nber of standard deviation to define outliers
%   output: indexbad : indexbad in the data vector of the outliers
%                      datavecmean,datavecstd, binpres: mean, std, and
%                      pressure lvl used to determine the outliers 
%                      (larger than datamean +/- qc_param.tollvl*datastd)
%                      
%
% created by L. Houpert (houpertloic@gmail.com), 08/08/2017, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox

%-------------------------------------------------------
% detect outliers in compare to mean pressure bin data
%-------------------------------------------------------

% Determination of the number of spikes
%--- affectation of the index
pbinvec = 0:qc_param.binsize:max(pres);
cstop=0;
cloopindex = 0;
tot_outliers_nb=0;
indexbad = [];
while cstop==0
cloopindex = cloopindex + 1;
data(indexbad)=nan;     
indok=find(~isnan(data));
nb_outliers = 0;
datavecmean = nan(length(pbinvec)-1,1);
datavecstd  = nan(length(pbinvec)-1,1) ;
binpres     = nan(length(pbinvec)-1,1);
for izz=1:length(pbinvec)-1
    isel = find(pres>pbinvec(izz) & pres<pbinvec(izz+1));
    datamean = nanmean(data(isel));
    datastd  = nanstd(data(isel),0);
    ibad = find(data(isel) > datamean + qc_param.tollvl*datastd | ...
                data(isel) < datamean - qc_param.tollvl*datastd);
    nb_outliers = nb_outliers + length(ibad);
    indexbad = [indexbad isel(ibad)'];
    datavecmean(izz) = datamean;
    datavecstd(izz)  = datastd;   
    binpres(izz)     = mean(pbinvec(izz:izz+1));
end

tot_outliers_nb = tot_outliers_nb + nb_outliers ;
if nb_outliers==0, cstop = 1; end
end

disp([ qc_param.varname ' outliers(s) detected in  ' num2str(qc_param.binsize) ...
    'm pressure bin (outside ' num2str(qc_param.tollvl)  ' standard deviation): ' num2str(tot_outliers_nb)])
disp([ qc_param.varname ' nber of loop for spike detection : ' num2str(cloopindex)])

end
