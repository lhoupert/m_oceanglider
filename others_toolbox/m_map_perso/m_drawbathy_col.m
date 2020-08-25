function [cs0,h0]=m_drawbathy_col(optn,bathyname,zlvls, limitcmap)
% function [cs,h]=m_drawbathy(optn,bathyname, zlvls, limitcmap)
% optn = 'contourf' or 'contour'
% bathyname='canyon' or 'medwest' 
% zlvls (ex: -3000:50:0)
% limitcmap = limit lat, lon pour fixer bathy : limitcmap=[minlon maxlon
% minlat maxlat]

global clong

%load('MyColormaps')
load('cmap_map')
medcmap=cmap;

if strfind(bathyname,'etopo6min')
    lgmin=limitcmap(1);
    lgmax=limitcmap(2);
    ltmin=limitcmap(3) ;
    ltmax=limitcmap(4);

    fic_bathy = '/home/lhoupert/Data/Boulot/function_MATLAB/current/CASCADE/bathymetrie/bathy6min.nc';
    if (exist(fic_bathy,'file') ~=0)
%     hbar = waitbar(0,'Please, Wait');
    nc_bathy = netcdf.open(fic_bathy,'NOWRITE');
    idlat_bathy=netcdf.inqVarID(nc_bathy,'latitude');
    lat_bathy=netcdf.getVar(nc_bathy,idlat_bathy,'double');
    idlon_bathy=netcdf.inqVarID(nc_bathy,'longitude');
    lon_bathy=netcdf.getVar(nc_bathy,idlon_bathy,'double');
    lat1 = find(lat_bathy>=ltmin);
    lat1 = lat1(1);
    lat2 = find(lat_bathy<=ltmax);
    lat2 = lat2(end);
    lon1 = find(lon_bathy>=lgmin);
    lon1 = lon1(1);
    lon2 = find(lon_bathy<=lgmax);
    lon2 = lon2(end); 

    z_id=netcdf.inqVarID(nc_bathy,'z');
    z_bathy=netcdf.getVar(nc_bathy,z_id,[lon1-1  lat1-1],[(lon2-lon1+1) (lat2-lat1+1)]); % Attention, les indices commencent a 0 dans les NetCDF (d'ou le -1).
    zbathy2=double(z_bathy');
    netcdf.close(nc_bathy);
    lat_bathy=lat_bathy(lat1:lat2);
    lon_bathy=lon_bathy(lon1:lon2);
    
    [lgbathy2 ltbathy2]=meshgrid(lon_bathy,lat_bathy);
    
    if strfind(optn,'contourf')
        if clong==0
         [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap);
         hold on
         set(h0,'Linestyle','none');
         
        else
          [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap);
          hold on
          set(h0,'Linestyle','none');           
          [cs02,h02]=m_contourf(lgbathy2+360,ltbathy2,zbathy2,zlvls);colormap(medcmap);
          hold on
          set(h02,'Linestyle','none');       
        end
%         freezeColors       
        
%          [cs00,h00]=m_contourf(lgbathy2,ltbathy2,zbathy2,[-11000 -5000]);
%          colormap([[0.4784    0.0627    0.8941];[0.4784    0.0627    0.8941]]);
%          hold on
%          set(h00,'Linestyle','none');
%          freezeColors

    elseif strfind(optn,'contour')
         [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap*0+1);
         hold on
         set(h0,'Linestyle','none');
         [cs,h]=m_contour(lgbathy2,ltbathy2,zbathy2,zlvls);
         set(h,'linecolor',[0.65 0.65 0.65]) ;
         freezeColors
    else
          disp('Preciser "optn" (contour ou contourf)');
    end
    
    end

elseif strfind(bathyname,'gebco')
    lgmin=limitcmap(1);
    lgmax=limitcmap(2);
    ltmin=limitcmap(3) ;
    ltmax=limitcmap(4);

    fic_bathy = '/home/lhoupert/Data/Boulot/function_MATLAB/current/CASCADE/bathymetrie/bathy_gebco.nc';
    if (exist(fic_bathy,'file') ~=0)
%     hbar = waitbar(0,'Please, Wait');
    nc_bathy = netcdf.open(fic_bathy,'NOWRITE');
    idlat_bathy=netcdf.inqVarID(nc_bathy,'latitude');
    lat_bathy=netcdf.getVar(nc_bathy,idlat_bathy,'double');
    idlon_bathy=netcdf.inqVarID(nc_bathy,'longitude');
    lon_bathy=netcdf.getVar(nc_bathy,idlon_bathy,'double');
    lat1 = find(lat_bathy>ltmin);
    lat1 = lat1(1)-1;
    lat2 = find(lat_bathy<ltmax);
    lat2 = lat2(end)+1;
    lon1 = find(lon_bathy>lgmin);
    lon1 = lon1(1)-1;
    lon2 = find(lon_bathy<lgmax);
    lon2 = lon2(end)+1; 

    z_id=netcdf.inqVarID(nc_bathy,'z');
    z_bathy=netcdf.getVar(nc_bathy,z_id,[lon1-1 lat1-1],[(lon2-lon1+1) (lat2-lat1+1)]); % Attention, les indices commencent a 0 dans les NetCDF (d'ou le -1).
    zbathy2=double(z_bathy');
    netcdf.close(nc_bathy);
    lat_bathy=lat_bathy(lat1:lat2);
    lon_bathy=lon_bathy(lon1:lon2);
    
    [lgbathy2 ltbathy2]=meshgrid(lon_bathy,lat_bathy);
    
    if strfind(optn,'contourf')
        if clong==0
         [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap);
         hold on
         set(h0,'Linestyle','none');
%         freezeColors
        else
          [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap);
          hold on
          set(h0,'Linestyle','none');           
          [cs02,h02]=m_contourf(lgbathy2-360,ltbathy2,zbathy2,zlvls);colormap(medcmap);
          hold on
          set(h02,'Linestyle','none');       
        end
%          [cs00,h00]=m_contourf(lgbathy2,ltbathy2,zbathy2,[-11000 -5000]);
%          colormap([[0.4784    0.0627    0.8941];[0.4784    0.0627    0.8941]]);
%          hold on
%          set(h00,'Linestyle','none');
%          freezeColors

    elseif strfind(optn,'contour')
         [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap*0+1);
         hold on
         set(h0,'Linestyle','none');
         [cs,h]=m_contour(lgbathy2,ltbathy2,zbathy2,zlvls);
         set(h,'linecolor',[0.65 0.65 0.65]) ;
         freezeColors
    else
          disp('Preciser "optn" (contour ou contourf)');
    end
    
    end
elseif strfind(bathyname,'medwest')
load '~/Data/Boulot/function_MATLAB/bathy/bathy_ifremer4.mat'
lgmin=limitcmap(1);
lgmax=limitcmap(2);
ltmin=limitcmap(3) ;
ltmax=limitcmap(4);
%lgmin=0; lgmax=10; ltmin=38 ; ltmax=46;  % juste pour fixer la colormap...
mb=find(lgbathy4(1,:)<lgmin | lgbathy4(1,:)>lgmax);
lgbathy4(:,mb)=[];zbathy3(:,mb)=[];ltbathy4(:,mb)=[];
nb=find(ltbathy4(:,1)<ltmin | ltbathy4(:,1)>ltmax);
ltbathy4(nb,:)=[];zbathy3(nb,:)=[];lgbathy4(nb,:)=[];

     if strfind(optn,'contourf')
     [cs0,h0]=m_contourf(lgbathy4,ltbathy4,zbathy3,zlvls);colormap(medcmap);
     set(h0,'Linestyle','none');
     hold on
%      [cs,h]=m_contour(lgbathy4,ltbathy4,zbathy3,zlvls(1:end-1));
%      set(h,'linecolor',[0.65 0.65 0.65])
%     freezeColors
     elseif strfind(optn,'contour')
     [cs,h]=m_contour(lgbathy4,ltbathy4,zbathy3,zlvls);
     set(h,'linecolor',[0.65 0.65 0.65])
     else
          disp('Preciser "optn" (contour ou contourf)');
     end

elseif strfind(bathyname,'canyon')
 

data=load('/home/lhoupert/Data/Boulot/function_MATLAB/bathy/GOL_base.dat');

lgbathy=data(:,1);
ltbathy=data(:,2);
zbathy=-data(:,3);
clear data;
lgmin=limitcmap(1);
lgmax=limitcmap(2);
ltmin=limitcmap(3) ;
ltmax=limitcmap(4);

%lgmin=0; lgmax=10; ltmin=38 ; ltmax=46;  % juste pour fixer la colormap...

mb=find(lgbathy<lgmin | lgbathy>lgmax | ...
     ltbathy<ltmin | ltbathy>ltmax);
lgbathy(mb)=[];zbathy(mb)=[];ltbathy(mb)=[];

%xi1=[lgmin:0.01:lgmax]; yi1=[ltmin:0.01:ltmax];
xi1=[min(lgbathy):0.002:max(lgbathy)]; yi1=[min(ltbathy):0.002:max(ltbathy)];
[lgbathy2 ltbathy2]=meshgrid(xi1,yi1);
%[lgbathy2 ltbathy2 zbathy2]=griddata(lgbathy,ltbathy,zbathy,xi1,yi1');
F = TriScatteredInterp(lgbathy,ltbathy,zbathy);
zbathy2=F(lgbathy2,ltbathy2);

%zbathy2(isnan(zbathy2))=0;%999;
%[cs,h]=m_contourf(lgbathy4,ltbathy4,zbathy3,zlvls);colormap(medcmap);
     if strfind(optn,'contourf')
     [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap);
     hold on
     set(h0,'Linestyle','none');
%       [cs,h]=m_contour(lgbathy2,ltbathy2,zbathy2,zlvls);
%       set(h,'linecolor',[0.65 0.65 0.65]) ;
%     freezeColors
     elseif strfind(optn,'contour')
     [cs0,h0]=m_contourf(lgbathy2,ltbathy2,zbathy2,zlvls);colormap(medcmap*0+1);
     hold on
     set(h0,'Linestyle','none');
     [cs,h]=m_contour(lgbathy2,ltbathy2,zbathy2,zlvls);
     set(h,'linecolor',[0.65 0.65 0.65]) ;
     freezeColors
     else
          disp('Preciser "optn" (contour ou contourf)');
     end
else
     disp('precisez le nom de la bathy');
end
