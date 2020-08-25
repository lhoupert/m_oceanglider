function [tgps2,lon2,lat2,tgpse,lone,late,DAC_qc]=correct_nogpsdata_ndives(tgps2,lon2,lat2,tgpse,lone,late,DAC_qc)
% function [tgps2,lon2,lat2,tgpse,lone,late,DAC_qc]=correct_nogpsdata_ndives(tgps2,lon2,lat2,tgpse,lone,late,DAC_qc)
%
% GPS correction when several dive without surfacing are made, 
%
% created by L. Houpert (houpertloic@gmail.com), 10/04/2016, as part of the
% oceano_data_toolbox project: https://bitbucket.org/Lhoupert/oceano_data_toolbox


ibadt2=find(diff(tgps2)<0);
%ibadte=find(diff(tgpse)<0);
tgps2(ibadt2)=tgps2(ibadt2+1);
lat2(ibadt2)=lat2(ibadt2+1);
lon2(ibadt2)=lon2(ibadt2+1);
%------------------------------------------------
[ctdive,iat,ict] = unique(tgps2);%;unique(round(tgps2*200)/200); 
tgps2_0=tgps2;
lat2_0 =lat2;
lon2_0 =lon2;
tgpse_0=tgpse;
late_0 =late;
lone_0 =lone;
itnun = [];
inun = [];
compt=0;
for itt=1:length(ctdive)-1
	if length(find(tgps2==ctdive(itt)))>1
		compt=compt+1;
		indnonuniq = find(tgps2==ctdive(itt))';
		itnun{compt} = indnonuniq;
		inun = [ inun indnonuniq];
		deltatime2 = (tgps2(indnonuniq(end)+1) - tgps2(indnonuniq(1)))/(length(indnonuniq));
		deltalon2 = (lon2(indnonuniq(end)+1) - lon2(indnonuniq(1)))/(length(indnonuniq));	
		deltalat2 = (lat2(indnonuniq(end)+1) - lat2(indnonuniq(1)))/(length(indnonuniq));			
		deltatimee = (tgpse(indnonuniq(end)) - tgpse(indnonuniq(1)-1))/(length(indnonuniq));
		deltalone = (lone(indnonuniq(end)) - lone(indnonuniq(1)-1))/(length(indnonuniq));	
		deltalate = (late(indnonuniq(end)) - late(indnonuniq(1)-1))/(length(indnonuniq));			
		for innq = 2:length(indnonuniq)
			tgps2(indnonuniq(innq)) = tgps2(indnonuniq(1)) +  (innq-1)*deltatime2;
			lon2(indnonuniq(innq)) = lon2(indnonuniq(1)) +  (innq-1)*deltalon2;	
			lat2(indnonuniq(innq)) = lat2(indnonuniq(1)) +  (innq-1)*deltalat2;	
		end	
		for innq = 1:(length(indnonuniq)-1)					
			tgpse(indnonuniq(innq)) = tgpse(indnonuniq(1)-1) +  (innq)*deltatimee;
			lone(indnonuniq(innq)) = lone(indnonuniq(1)-1) +  (innq)*deltalone;	
			late(indnonuniq(innq)) = late(indnonuniq(1)-1) +  (innq)*deltalate;		
		end
		DAC_qc(indnonuniq) = 4; % DAC certainly bad TODO: to verify
	end
end
end
