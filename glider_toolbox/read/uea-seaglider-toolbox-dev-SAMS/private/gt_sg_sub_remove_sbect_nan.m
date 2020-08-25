function [data0,indexbaddepth,indexbadct] = gt_sg_sub_remove_sbect_nan(data0,outputinlogfile)
% function [data0,indexbaddepth,indexbadct] = gt_sg_sub_remove_sbect_nan(data0,outputinlogfile)
%                                outputinlogfile = 0 (no output file) or 1
% Remove data with nan for condFreq or tempFreq or negative depth % L Houpert 18/01/2016
if nargin<2
    outputinlogfile=0;
end

indexbaddepth = [];
indexbadct = [];
for ijk=1:length(data0.eng)
	condfr = data0.eng(ijk).sbect.condFreq;
	tempfr = data0.eng(ijk).sbect.tempFreq;
	ibadcond = find(isnan(condfr));
	ibadtemp = find(isnan(tempfr));
	ibaddepth = find(data0.eng(ijk).depth < 0);	
	if ~isempty(ibaddepth);
        indexbaddepth = [indexbaddepth ijk];
	end
	if ~isempty(ibadtemp) | ~isempty(ibadcond) ;
        indexbadct = [indexbadct ijk];
	end	 
	 	
	ibad = union(ibaddepth, union(ibadcond,ibadtemp));
	if ~isempty(ibad);
		%data0.eng(ijk).sbect.condFreq(ibad)=[];
		%data0.eng(ijk).sbect.tempFreq(ibad)=[];	
		data0.eng(ijk).vbdCC(ibad)=[];		
		data0.eng(ijk).rollCtl(ibad)=[];				
		data0.eng(ijk).rollAng(ibad)=[];		
		data0.eng(ijk).rec(ibad)=[];		
		data0.eng(ijk).pitchCtl(ibad)=[];		
		data0.eng(ijk).pitchAng(ibad)=[];				
		data0.eng(ijk).head(ibad)=[];		
		data0.eng(ijk).elaps_t_0000(ibad)=[];			
		data0.eng(ijk).elaps_t(ibad)=[];		
		data0.eng(ijk).depth(ibad)=[];				
		data0.eng(ijk).GC_phase(ibad)=[];	
		% for the data0 associated to the different sensors:
		ffl=struct2cell(data0.eng(ijk));  	
		fflname = fieldnames(data0.eng(ijk));
		sfield=[];
		for iff=1:length(ffl)
	  		if isstruct(ffl{iff})
				fflname2 = fieldnames(data0.eng(ijk).(fflname{iff}));
				for ina=1:length(fflname2)
                    if ~isempty(data0.eng(ijk).(fflname{iff}).(fflname2{ina}))
                        data0.eng(ijk).(fflname{iff}).(fflname2{ina})(ibad)=[];
                    end
				end
			end
	  	end	
	  		  	
	end	  	
end

if ~isempty(indexbaddepth)
	disp('PROBLEM:  negative PRESSURE values in dives: ')
    disp(num2str([data0.eng(indexbaddepth).dive]))
	disp('Dives processed and negative values removed')
    if outputinlogfile==1
        file1 = fopen('dives_with_negative_depth_values.log','w');
        fprintf(file1,'PROBLEM:  negative PRESSURE values in dives:\n ');
        fprintf(file1,'%d\n',[data0.eng(indexbaddepth).dive]);
        fclose(file1);
    end
end

if ~isempty(indexbadct)
    disp('PROBLEM: no TEMPERATURE/CONDUCTIVITY (NaN) for some levels in dives:')
    disp( num2str([data0.eng(indexbadct).dive]))
    disp('Dives processed and the levels without temperature/conductivity data are removed')
    if outputinlogfile==1   
        file2 = fopen('dives_with_nan_CT.log','w');
        fprintf(file2,'PROBLEM: no TEMPERATURE/CONDUCTIVITY (NaN) for some levels in dives:\n');
        fprintf(file2,'%d\n',[data0.eng(indexbadct).dive]);
        fclose(file2);
    end
end
end
