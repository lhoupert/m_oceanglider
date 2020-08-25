function [dataDB]=DBhisto_loaddata6(path_db,arealimit,timeselect,minp,instrtyp,qcselec,qlvltp,qlvlss)
%[dataDB]=DBhisto_loaddata6(path_db,arealimit,timeselect,minp,instrtyp,qcselec,qlvltp,qlvlss);
% with	    path_db= path of database
%			arealimit= limt polygon selection
% 				ex: arealimit.lat=[latmin latmin latmax latmax]
%				    arealimit.lon=[lonmin lonmax lonmax lonmin]
%			timeselect.min= first date
%			timeselect.max= last date
%           minp= depth of selected profiles (ex:1200 will keep only profiles going at least to (or deeper) 1200db)
%           instrtyp= Kind of instrument ('ALL', 'CTD', 'ARGO', 'XBT',
%                   'GLIDERS', etc...)
% quality check selection only on temperature profiles
%           qcselec= need qualtiy check data? (1 for yes, 0 for no)
%           qlvltp= quality flag level if qcselec==1 
%                 QC flag: 1: GOOD ,2: PROBABLY GOOD, 3: PROBABLY BAD, 
%                 4: BAD, 0: NO QC, 8: NO TIME(<4 pts)
%                 7: probleme localisation temps ou position
%                 8: pas de données
%                10: bonnes données brutes
%                11: bonnes données brutes ajustées (controle PI)
%                20: probablement bonne données(QC de la "forme" des profils ctd profond)
%                       ex: qlvltp=[1 2 0];
%		 41: No data in the upper 15m of the profile
%           qlvlss = same as qlvltp
%%test:
%path_db='/home/lhoupert/Data/Boulot/DB_HISTO/BOXES_BISbi/';
%arealimit.lat=[42 42 41.5000 41.5000 42]
%arealimit.lon=[4.7500 5.2500 5.2500 4.7500 4.7500]
%timeselect.min=datenum(1969,01,01);
%timeselect.max=datenum(2013,01,01);
%ktest=[[3 4 4 3];[3.5 4.8 4.8 3.5];[3 6 6 4]]
%ltest=[[40 40 41 41];[41 41 41.75 41.75];[40 40 43 43 ]]
%plot(arealimit.lon,arealimit.lat) 
%hold on
%plot(ktest(1,:),ltest(1,:),'r')
%plot(ktest(2,:),ltest(2,:),'r')
%plot(ktest(3,:),ltest(3,:),'r')
gridstp=0.5;% grid width
%=================================================================================
% READ DATABASE FILES k and l correspond to the left-bottom corner of the 0.5°box
boites=[path_db '*.mat'];
boites=dir(boites);
l=nan(length(boites),4);
k=nan(length(boites),4);
dataDB=[];
for v=1:length(boites)
    point=boites(v).name;
    l(v,1)=str2num(point(6:8))/10;
    if strcmp(point(1),'E')
        aa=1;
    elseif strcmp(point(1),'W')
        aa=-1;
    end
    k(v,1)=aa*str2num(point(2:4))/10;
end
l(:,2)=l(:,1);				k(:,2)=k(:,1) + gridstp;
l(:,3)=l(:,1) + gridstp;		k(:,3)=k(:,1) + gridstp;
l(:,4)=l(:,1) + gridstp;		k(:,4)=k(:,1);

% Check if at least one box corner is in the polygon (def in arealimit)
ii_inpoly=inpolygon(k,l,arealimit.lon,arealimit.lat);
ibox=find(sum(ii_inpoly,2)>0);

% if no corner is in the polygon but the polygon cross the grid
boxnear=find( k(:,1)<max(arealimit.lon) &  k(:,2)>min(arealimit.lon) & ...
		l(:,1)<max(arealimit.lat) &  l(:,4)>min(arealimit.lat));
iv=setdiff(boxnear,ibox);

if ~isempty(iv)
ibox2=[];
for ij=iv'
	ii_inpoly2=inpolygon(arealimit.lon,arealimit.lat,k(ij,:),l(ij,:));
	if  sum(ii_inpoly2)>0
		ibox2=[ibox2 ij];
	end
end

ibox=[ibox ibox2];
end
%============================================================================
% Read database in the polygon and selection of profiles
%============================================================================

% initialisation
pmax=0;
nn=0;
for v=ibox'
    point=boites(v).name;
    l=point(6:8);
    k=point(1:4);
%    disp(['## ' num2str(v) ' --- ' k '_N' l  ' / ' datestr(now)]); 
    eval(['load ' path_db  k '_' l ])   
    
    if isempty(TIME)
        clear TP Pref SS LAT LON TIME CAMPAIGN IDINSTR TYPINSTR IDPF OP DATASOURCE TP_QC SS_QC 
        continue
    else
    
    depth_cast=nan(1,length(TIME));
    for zi=1:length(TIME)
        ind1=find((TP(:,zi)~=0), 1, 'last' );
        ind2=find(~isnan(TP(:,zi)), 1, 'last' );
        depth_cast(zi)=min([ind1 ind2]);
    end
    
        if ~isempty(strfind(instrtyp,'ALL')) & qcselec==0
            [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
                | TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp);

        elseif qcselec==0;
            instrbol=nan(1,length(TIME));
            for zii2=1:length(TIME)
                if ~isempty(strfind(TYPINSTR{zii2},instrtyp))
                    instrbol(zii2)=1;
                else
                    instrbol(zii2)=0;
                end
            end
            instrbol=logical(instrbol);
            [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
                | TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp | ~instrbol);     

        elseif ~isempty(strfind(instrtyp,'ALL'))        
            mqc1=length(qlvltp);
            qcbol1=nan(mqc1,length(TIME));
            for zii2=1:mqc1;
                qcbol1(zii2,:)= (TP_QC==qlvltp(zii2));
            end
            qc_sel1=any(qcbol1);
            mqc2=length(qlvlss);
            qcbol2=nan(mqc2,length(TIME));
            for zii2=1:mqc2;
                qcbol2(zii2,:)= (SS_QC==qlvlss(zii2));
            end
            qc_sel2=any(qcbol2);        
            
            qc_sel=qc_sel1 & qc_sel2;
            
            [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
                | TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp | ~qc_sel);          

        else % <==> Not for all data with quality check selection
            instrbol=nan(1,length(TIME));
            for zii2=1:length(TIME)
                if ~isempty(strfind(TYPINSTR{zii2},instrtyp))
                    instrbol(zii2)=1;
                else
                    instrbol(zii2)=0;
                end
            end
            instrbol=logical(instrbol);        

            mqc1=length(qlvltp);
            qcbol1=nan(mqc1,length(TIME));
            for zii2=1:mqc1;
                qcbol1(zii2,:)= (TP_QC==qlvltp(zii2));
            end
            qc_sel1=any(qcbol1);
            mqc2=length(qlvlss);
            qcbol2=nan(mqc2,length(TIME));
            for zii2=1:mqc2;
                qcbol2(zii2,:)= (SS_QC==qlvlss(zii2));
            end
            qc_sel2=any(qcbol2);        
             
            qc_sel=qc_sel1 & qc_sel2;

            [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
                | TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp ...
                | ~instrbol| ~qc_sel);           

        end
    
    indpf=unique(j);
    TP(:,indpf)=[];
    SS(:,indpf)=[];
    LAT(indpf)=[];
    LON(indpf)=[];
    TIME(indpf)=[];
    CAMPAIGN(indpf)=[];
    IDINSTR(indpf)=[];
    TYPINSTR(indpf)=[];
    IDPF(indpf)=[];
    OP(indpf)=[];
    DATASOURCE(indpf)=[];
    TP_QC(indpf)=[];	 
    SS_QC(indpf)=[]; 
    depth_cast(indpf)=[];      
     
    if pmax<max(depth_cast)
        pmax=max(depth_cast);
    end
    
    nn=nn+length(TIME);
    end

    clear TP Pref SS LAT LON TIME CAMPAIGN IDINSTR TYPINSTR IDPF OP DATASOURCE TP_QC SS_QC 
end

AX=(1:1:pmax)';
ptempdb=nan(length(AX),nn);
saldb=nan(length(AX),nn);
latdb =nan(1,nn);
londb=nan(1,nn);
timedb=nan(1,nn);
namedb=cell(1,nn);
campaigndb=cell(1,nn);
typinstrdb=cell(1,nn);
idinstrdb=cell(1,nn);
datascedb=cell(1,nn);
operdb=cell(1,nn);
tp_qcdb=nan(1,nn);
ss_qcdb=nan(1,nn);
maxdepth=nan(1,nn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD TIME
cc=0;
pind=0;
for v=ibox'

    point=boites(v).name;
    l=point(6:8);
    k=point(1:4);
%    disp(['## ' num2str(v) ' --- ' k '_N' l  ' / ' datestr(now)]); 
    eval(['load ' path_db  k '_' l ])   

    
    if isempty(TIME)
        clear TP Pref SS LAT LON TIME CAMPAIGN IDINSTR TYPINSTR IDPF OP DATASOURCE TP_QC SS_QC 
        continue
    else
    
    depth_cast=nan(1,length(TIME));
    for zi=1:length(TIME)
        ind1=find((TP(:,zi)~=0), 1, 'last' );
        ind2=find(~isnan(TP(:,zi)), 1, 'last' );
        depth_cast(zi)=min([ind1 ind2]);
    end
%     [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
%     		| TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp);
    if ~isempty(strfind(instrtyp,'ALL')) & qcselec==0
        [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
    		| TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp);
        
    elseif qcselec==0;
        instrbol=nan(1,length(TIME));
        for zii2=1:length(TIME)
            if ~isempty(strfind(TYPINSTR{zii2},instrtyp))
                instrbol(zii2)=1;
            else
                instrbol(zii2)=0;
            end
        end
        instrbol=logical(instrbol);
        [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
    		| TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp | ~instrbol);     
        
    elseif ~isempty(strfind(instrtyp,'ALL'))        
            mqc1=length(qlvltp);
            qcbol1=nan(mqc1,length(TIME));
            for zii2=1:mqc1;
                qcbol1(zii2,:)= (TP_QC==qlvltp(zii2));
            end
            qc_sel1=any(qcbol1);
            mqc2=length(qlvlss);
            qcbol2=nan(mqc2,length(TIME));
            for zii2=1:mqc2;
                qcbol2(zii2,:)= (SS_QC==qlvlss(zii2));
            end
            qc_sel2=any(qcbol2);        
             
            qc_sel=qc_sel1 & qc_sel2;
            
        [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
    		| TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp | ~qc_sel);          
        
    else % <==> Not for all data with quality check selection
        instrbol=nan(1,length(TIME));
        for zii2=1:length(TIME)
            if ~isempty(strfind(TYPINSTR{zii2},instrtyp))
                instrbol(zii2)=1;
            else
                instrbol(zii2)=0;
            end
        end
        instrbol=logical(instrbol);        
        
            mqc1=length(qlvltp);
            qcbol1=nan(mqc1,length(TIME));
            for zii2=1:mqc1;
                qcbol1(zii2,:)= (TP_QC==qlvltp(zii2));
            end
            qc_sel1=any(qcbol1);
            mqc2=length(qlvlss);
            qcbol2=nan(mqc2,length(TIME));
            for zii2=1:mqc2;
                qcbol2(zii2,:)= (SS_QC==qlvlss(zii2));
            end
            qc_sel2=any(qcbol2);        
             
            qc_sel=qc_sel1 & qc_sel2;
        
        [i j]=find( ~inpolygon(LON,LAT,arealimit.lon,arealimit.lat) ...
    		| TIME<timeselect.min | TIME>timeselect.max | depth_cast<minp ...
            | ~instrbol| ~qc_sel);           
        
    end

        
    indpf=unique(j);
    TP(:,indpf)=[];
    SS(:,indpf)=[];
    LAT(indpf)=[];
    LON(indpf)=[];
    TIME(indpf)=[];
    CAMPAIGN(indpf)=[];
    IDINSTR(indpf)=[];
    TYPINSTR(indpf)=[];
    IDPF(indpf)=[];
    OP(indpf)=[];
    DATASOURCE(indpf)=[];
    TP_QC(indpf)=[];	 
    SS_QC(indpf)=[]; 
    depth_cast(indpf)=[];    
       
    nn=nn+length(TIME);
    
    [m n]=size(TP);
    pind=max(depth_cast);

    timedb((cc+1):cc+n)      = TIME;
    londb((cc+1):cc+n)       = LON;
    latdb((cc+1):cc+n)       = LAT;
    ptempdb(1:pind,(cc+1):cc+n) = TP(1:pind,:);
    saldb(1:pind,(cc+1):cc+n)   = SS(1:pind,:);
    namedb((cc+1):cc+n)      = IDPF;    
    campaigndb((cc+1):cc+n)  = CAMPAIGN;
    typinstrdb((cc+1):cc+n)  = TYPINSTR;
    idinstrdb((cc+1):cc+n)   = IDPF;
    datascedb((cc+1):cc+n)   = DATASOURCE;
    operdb((cc+1):cc+n)      = OP;
    tp_qcdb((cc+1):cc+n)     = TP_QC;
    ss_qcdb((cc+1):cc+n)     = SS_QC;
    maxdepth((cc+1):cc+n)    = pind;

    end

    clear TP Pref SS LAT LON TIME CAMPAIGN IDINSTR TYPINSTR IDPF OP DATASOURCE TP_QC SS_QC 

    cc=cc+n;
end    


ptempdb(ptempdb==0)=nan;
saldb(saldb==0)=nan;

for idd=1:length(timedb)
dataDB(idd).time     = timedb(idd);
dataDB(idd).lon      = londb(idd);
dataDB(idd).lat      = latdb(idd);
dataDB(idd).name     = namedb{idd};
dataDB(idd).ppp      = AX;
dataDB(idd).ptemp    = ptempdb(:,idd);
dataDB(idd).sal      = saldb(:,idd);
dataDB(idd).campaign = campaigndb{idd}; 
dataDB(idd).oper     = operdb{idd};     
dataDB(idd).idinstr	= idinstrdb{idd};
dataDB(idd).typinstr = typinstrdb{idd};
dataDB(idd).datasce  = datascedb{idd};    
dataDB(idd).qlvl_tp  = tp_qcdb(idd);    
dataDB(idd).qlvl_ss  = ss_qcdb(idd);    
dataDB(idd).maxdpth   =maxdepth(idd);
end





















