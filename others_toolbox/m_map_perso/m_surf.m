function m_surf(long,lat,data,varargin);

global MAP_PROJECTION 

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if min(size(long))==1 & min(size(lat))==1,
 [long,lat]=meshgrid(long,lat);
end;

[X,Y]=m_ll2xy(long,lat,'clip','on'); % default: 'on' %First find the points outside

i=isnan(X);      % For these we set the *data* to NaN...
data(i)=NaN;

                 % And then recompute positions without clipping. THis
                 % is necessary otherwise contouring fails (X/Y with NaN
                 % is a no-no. Note that this only clips properly down
                 % columns of long/lat - not across rows. In general this
                 % means patches may nto line up properly a right/left edges.
if any(i(:)), [X,Y]=m_ll2xy(long,lat,'clip','patch'); end;  

if any(~i(:)),

 % Bug in contourf call - Solution Number: 1-1W36E8
 % (that involved whether or not the X/Y matrices were handled
 % correctly). In later versions a problem comes up with the renderer
 % that could be solved here, but is better handled in m_grid (see comments
 % in code there).
  if any(strfind(version,'(R14) Service Pack 3')) | ...
     any(strfind(version,'7.2.0.294 (R2006a)')) | ....
     any(strfind(version,'7.2.0.294 (R2006a)')) %%| ....
   %%  any(strfind(version,'7.4.0.336 (R2007a)')),
   surf('v6',X,Y,data,varargin{:});
  else
   surf(X,Y,data,varargin{:});
  end;
   
 
else
  cs=[];h=[];
end;

if nargout==0,
 clear cs h
end;
