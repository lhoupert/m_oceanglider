function [mt] = mean_av_gauss3(d_mean,t,km)
% t = tableau à smoother
% d_mean position des points
% km largeur de la gaussienne (le sigma)

for i = 1:length(t)
	weight_tab = gauss_distribution(d_mean,d_mean(i),km);
	iqnan = find(weight_tab<0.01);weight_tab(iqnan) = NaN;
	mt(i,:) = wmean(t,weight_tab);    
end

end
