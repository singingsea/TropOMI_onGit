function[lat_r,lon_r,x0,y0,x1,y1] = wind_rotation(u,v,stlat,stlon,lat,lon)

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
 

x0=111.3*(lon-stlon)*cos(stlat/180*3.1415926);
y0=111.3*(lat-stlat);
 
coswd=cos(-wd/180*3.1415926);
sinwd=sin(-wd/180*3.1415926);
%coswd=cos(wd/180*3.1415926);
%sinwd=sin(wd/180*3.1415926);
x1= x0*coswd + y0*sinwd;
y1=-x0*sinwd + y0*coswd;
 
lat_r=y1/111.3 + stlat;
lon_r=x1/111.3/cos(stlat/180*3.1415926) + stlon;
