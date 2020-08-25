function zi=ship_interpbathy(bathyname,lon,lat,bathydir)
%function [zi]=ship_interpbathy(bathyname,lon,lat,bathydir)
%               example:
%               [zi]=ship_interpbathy('gebco',lon,lat,'~/Dropbox/Work/function_MATLAB/bathymetrie')
%

lgmin=min(lon)-1;
lgmax=max(lon)+1;
ltmin=min(lat)-1;
ltmax=max(lat)+1;

if strfind(bathyname,'etopo6min')


    fic_bathy = [bathydir filesep 'bathy6min.nc'];

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
    
    end
    
elseif strfind(bathyname,'gebco')

    fic_bathy = [bathydir filesep 'bathy_gebco.nc'];  
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

    end
    
elseif strfind(bathyname,'medwest')
load([bathydir filesep 'bathy_ifremer4.mat'])

%lgmin=0; lgmax=10; ltmin=38 ; ltmax=46;  % juste pour fixer la colormap...
mb=find(lgbathy4(1,:)<lgmin | lgbathy4(1,:)>lgmax);
lgbathy4(:,mb)=[];zbathy3(:,mb)=[];ltbathy4(:,mb)=[];
nb=find(ltbathy4(:,1)<ltmin | ltbathy4(:,1)>ltmax);
ltbathy4(nb,:)=[];zbathy3(nb,:)=[];lgbathy4(nb,:)=[];

lgbathy2 = lgbathy4;
ltbathy2 = ltbathy4;
zbathy2  = zbathy3;

elseif strfind(bathyname,'canyon')
 

data=load([bathydir filesep 'GOL_base.dat']);

lgbathy=data(:,1);
ltbathy=data(:,2);
zbathy=-data(:,3);
clear data;

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

end
zi=interp2(lgbathy2,ltbathy2,zbathy2,lon,lat);

