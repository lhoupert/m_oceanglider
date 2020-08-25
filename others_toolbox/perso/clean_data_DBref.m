function [refDB] = clean_data_DBref(dataDB)
%
% refDB = clean_data_DBref(dataDB)
%
%   reorganise the structure datDBref and apply some final qc to remove
%   baddata
%
reftime = [dataDB.time];
reflat  = [dataDB.lat];
reflon  = [dataDB.lon];
refTP   = [dataDB.ptemp];
refSS   = [dataDB.sal]; 
refPP   = [dataDB.ppp];
[ii,jj]=find(refSS<25 | refSS>40 | refTP<-3 | refTP>30 | (refTP<4 & refSS <26) ); 
jbad=unique(jj);jok=setdiff(1:length(refSS(1,:)),jbad);
refTP(:,jbad)=[];
refSS(:,jbad)=[];
refPP(:,jbad)=[];
if 0
refSS0 = refSS; clear refSS
refTP0 = refTP; clear refTP
nbstd = 6;
bsize = 50; % binsize for the median
%---------------------------------------------%
refSS =clean_spike_DB_median(refSS0,pref,bsize,nbstd) ;
refTP =clean_spike_DB_median(refTP0,pref,bsize,nbstd) ;
%----------------------------------------------
end

for ijk = 1:length(dataDB)
    refDB.name{ijk} = [dataDB(ijk).typinstr ' - ' datestr(dataDB(ijk).time,'dd/mm/yyyy') ' # ' dataDB(ijk).datasce];
end
refDB.time = reftime;
refDB.lat  = reflat;
refDB.lon  = reflon;
refDB.TP   = refTP;
refDB.SS   = refSS;
refDB.PP   = refPP;
end