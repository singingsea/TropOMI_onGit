function[lat_r,lon_r,x0,y0,x1,y1,sites_distance] = wind_rotation_2sites(u,v,stlat1,stlon1,stlat2,stlon2,lat,lon)

ws=sqrt(u*u + v*v);
 
if ( v>0 )        
    wd=180 + atan(u/v) * 180/ 3.1415926;
elseif ( v<0 && u>0 )
    wd=360 + atan(u/v)* 180/ 3.1415926;
elseif ( v<0 && u<0 ) 
    wd=atan(u/v)* 180/ 3.1415926;
elseif ( u==0 && v<0 ) 
    wd=0;
elseif ( v==0 && u>0 )
    wd=270;
elseif ( v==0 && u<0 ) 
    wd=90;
elseif ( v==0 && u==0 )
    wd=0;
end
 
ws=ws*3.6;  % WS in km per h  ;
 

[arclen,az0] = distance(stlat1,stlon1,stlat2,stlon2);% get archlength and az between two sites
sites_distance = earthRadius('km').*deg2rad(arclen); % get distance between two sites note --> here sites_distance is in [km]

% since the wind is always from station1 towards station2, the site
% distance is always negative!!!
sites_distance = -sites_distance;

% the following lines make sure the sites distance follow the coordination
% --> negative y corresponding to down wind pixel
% if stlat1 > stlat2
%     sites_distance = -sites_distance;
% elseif stlat1 < stlat2
%     sites_distance = sites_distance;
% elseif stlat1 == stlat2
%     if stlon1 < stlon2
%         sites_distance = -sites_distance;
%     elseif stlon1 > stlon2
%         sites_distance = sites_distance;
%     end
% end
    
    
az1 = rem(az0 + 180,360);% convert azimuth of two sites to Vitali's North-South coordinate (say, point south is positive, not the meterological defination!)
%wd = rem((wd + (360-az1) + 360),360);% here the new "wind direction" is compensated by az1
wd = rem((wd + (360+az1) + 360),360);% here the new "wind direction" is compensated by az1

x0=111.3*(lon-stlon1)*cos(stlat1/180*3.1415926);
y0=111.3*(lat-stlat1);
 
coswd=cos(-wd/180*3.1415926);
sinwd=sin(-wd/180*3.1415926);
% coswd=cos(wd/180*3.1415926);
% sinwd=sin(wd/180*3.1415926);
x1= x0*coswd + y0*sinwd;
y1=-x0*sinwd + y0*coswd;
 
lat_r=y1/111.3 + stlat1;
lon_r=x1/111.3/cos(stlat1/180*3.1415926) + stlon1;
