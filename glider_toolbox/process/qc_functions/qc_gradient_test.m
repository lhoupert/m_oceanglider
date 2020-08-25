function [indexbad] = qc_gradient_test(data,pres,qc_param)
% function [indexbad] = qc_gradient_test(data,pres,qc_param)
% Function to detect realistic gradient 
%   input: data : data to check
%          pres : pressure associated to the data to test
%          qc_param.varname: name of the QC variable (for display purpose) 
%          qc_param.presmin: pressure min of the layer tested
%          qc_param.presmax: pressure max of the layer tested
%          qc_param.tollvl : level of tolerance for gradient value
%   output: indexbad : index in data of the bad gradient
%
% created by L. Houpert (houpertloic@gmail.com), 28/01/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox


% --------------------------
% Gradient test (need at least 3 points in a 50m layer)
% --------------------------
% Determination of the number of failed gradient tests
%--- affectation of the index
cstop=0;
cloopindex = 0;
tot_spike_nb=0;
indexbad = [];
while cstop==0
cloopindex = cloopindex + 1;
data(indexbad)=nan; 
indok=find(~isnan(data));
nb_spike = 0;
for ijk=2:(length(indok)-1)
  data2 = data(indok(ijk));
  data1 = data(indok(ijk-1));
  data3 = data(indok(ijk+1)); 
  if abs(pres(indok(ijk+1)) - pres(indok(ijk-1))) < 50 & ...	
  	((pres(indok(ijk)) > qc_param.presmin  && pres(indok(ijk))< qc_param.presmax && ...
        (abs(data2-(data3+data1)/2)/abs(pres(indok(ijk-1)) -pres(indok(ijk+1)))/2   > qc_param.tollvl )))
     indexbad = [indexbad indok(ijk)];    
     nb_spike = nb_spike + 1;   
  end
end
tot_spike_nb = tot_spike_nb + nb_spike ;
if nb_spike==0, cstop = 1; end
end

disp([ qc_param.varname ' gradient test(s) failed in the ' num2str(qc_param.presmin) ...
    '-' num2str(qc_param.presmax)  'db layer: ' num2str(tot_spike_nb)])
disp([ qc_param.varname ' nber of loop for gradient test detection : ' num2str(cloopindex)])

