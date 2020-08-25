function [data] =clean_spike_DB_median(data,pres,bsize,nbstd)
%function [data] =clean_spike_DB_median(data,pres,bsize,nbstd)

nbstdgl = nbstd;
bsizegl = bsize;

data00=data;
pbin=0:bsize:max(pres);
pbinref=(pbin(2)-pbin(1))/2:bsize:pbin(end);
pbingl=0:bsizegl:200; %only for the first 200m
pbinrefgl=(pbingl(2)-pbingl(1))/2:bsizegl:pbingl(end);

[mm nn]=size(data);


total_ibad=999;
climm=0;
indpff = 1:nn;
while total_ibad~=0 & climm<6
data0=data;
dat1medbin = nan(length(pbinref),nn);
dat1meddevbin = nan*dat1medbin;
dat1medbingl = nan(length(pbinrefgl),nn);;
dat1glmeddevbin = nan(length(pbinrefgl),nn);;
dat1bad =cell(1,nn);
total_ibad=0;
length(indpff);
for icc=indpff
	dat1 = data(:,icc);
	indint = find(~isnan(dat1));
	ibaddata=[];
%	for ipp=1:length(pbinrefgl) % for the first step we used the deviation to a "mean profile" calculated with all the median profiles
%		izz0 = find(pres(indint)>=pbingl(ipp) & pres(indint)<=pbingl(ipp+1));
%		izz=indint(izz0);
%		ldata=length(izz);
%		dat1medbingl(ipp,icc) = nanmedian(dat1(izz));	
%	end
	for ipp=1:length(pbinref)
		izz0 = find(pres(indint)>=pbin(ipp) & pres(indint)<=pbin(ipp+1));
		izz=indint(izz0);
		ldata=length(izz);
		dat1medbin(ipp,icc) = nanmedian(dat1(izz));
		dat1meddevbin(ipp,icc) = sqrt(nansum( (dat1(izz) - dat1medbin(ipp,icc)).^2)/(length(izz)-1));		
		ibaddata = [ibaddata izz(find(abs(dat1(izz) - dat1medbin(ipp,icc)) >  nbstd* dat1meddevbin(ipp,icc)))' ]; 
		
	end
	dat1bad{icc} = unique(ibaddata);
	
	total_ibad = total_ibad + length(dat1bad{icc});
	
	if climm > 0 % for the first step we used the deviation to a "mean profile" calculated with all the median profiles
	data(dat1bad{icc},icc) = nan ;
	end
	%figure;
	%hold on
	%plot(dat1,-pres,'-+','color',[0.6 0.6 0.6])
	%plot(dat1(dat1bad{icc}),-pres(dat1bad{icc}),'-+r')	
	%plot(dat1medbin(:,icc) - nbstd*dat1meddevbin(:,icc),-pbinref,'k','linewidth',1)
	%plot(dat1medbin(:,icc) + nbstd*dat1meddevbin(:,icc),-pbinref,'k','linewidth',1)
	
end

if climm <1 % for the first step we used the deviation to a "mean profile" calculated with all the median profiles
gl_med_pf=nanmean(dat1medbin,2);%nanmean(dat1medbingl,2);
for icc=indpff
	dat1 = data(:,icc);
	indint = find(~isnan(dat1));
	ibaddata=[];
	for ipp=1:length(pbinrefgl)
		izz0 = find(pres(indint)>=pbingl(ipp) & pres(indint)<=pbingl(ipp+1));
		izz=indint(izz0);
		ldata=length(izz);
		dat1glmeddevbin(ipp,icc) = sqrt(nansum( (dat1(izz) - gl_med_pf(ipp)).^2)/(length(izz)-1));
	end	
end
gl_med_dev=nanmean(dat1glmeddevbin,2);
for icc=indpff
	dat1 = data(:,icc);
	indint = find(~isnan(dat1));
	ibaddata=[];
	for ipp=1:length(pbinrefgl)
		izz0 = find(pres(indint)>=pbingl(ipp) & pres(indint)<=pbingl(ipp+1));
		izz=indint(izz0);
		ldata=length(izz);
		ibaddata = [ibaddata izz(find(abs(dat1(izz) - gl_med_pf(ipp)) >  nbstdgl* gl_med_dev(ipp)))' ]; 
	end
	dat1bad{icc} = unique(ibaddata);
	data(dat1bad{icc},icc) = nan ;		
end
 
else

indpff=[];
for ijk=1:length(dat1bad)
	if ~isempty(dat1bad{ijk})
	indpff = [indpff ijk];
	end
end

end %fi climm



climm = climm+1;

if 0
figure;
hold on
plot(data,-PP,'+-','color',[0.6 0.6 0.6])
plot(data0(isnan(data)),-PP(isnan(data)),'+r')

figure;
hold on
plot(data,TP,'+-','color',[0.6 0.6 0.6])
plot(data0(isnan(data)),TP(isnan(data)),'+r')
end %fi 

end

end
