function[lat_r,lon_r,x0,y0,x1,y1] = wind_rotation_2sites_v2(u,v,stlat,stlon,stlat2,stlon2,lat,lon)
% this function is trying to repeat Vitali's wind_rotation function, which
% has some different defination of wind direction
% please note, in this new wind_rotation function, the defication of wind
% direction is following meterological defination (true north is 0 degree, west wind (wind from west, towards east) is 90 degree)


wd = atan2d(u,v);% get wind direction in degree
wd = rem((wd+360),360);% convert to [0 360] range

[arclen,az0] = distance(stlat,stlon,stlat2,stlon2);% get archlength and az between two sites
sites_distance = earthRadius('km').*deg2rad(arclen); % get distance between two sites note --> here sites_distance is in [km]

% since the wind is always from station1 towards station2, the site
% distance is always negative!!!
sites_distance = -sites_distance;

% this function need mapping tool ... 
[arclen,az] = distance(stlat,stlon,lat,lon);% get the arclength and azimuth angle between station and pixel
az = az-wd-az0+180; % based on pixel rotation method, after rotation, all wind is from the north (north wind) --> i.e. point toward south is positive!


% this is from mapping tool again.
[lat_r,lon_r] = reckon(stlat,stlon,arclen,az); %get the new location


% calculate y0 --> new north-south distance between site and original pixel
[arclen_y0,az_y0] = distance(stlat,stlon,lat,stlon);% note --> here arclen is in degree
y0 = earthRadius('km').*deg2rad(arclen_y0); % note --> here y1 is in [km]
if az_y0 == 0
    y0 = -y0;
elseif az_y0 == 180
    y0 = y0;
else
    disp('Warnning ... the sign of y0 is not determined! Please check!');
end

% calculate x0 --> new east-west distance between site and original pixel
[arclen_x0,az_x0] = distance(stlat,stlon,stlat,lon);% note --> here arclen is in degree
x0 = earthRadius('km').*deg2rad(arclen_x0); % note --> here x1 is in [km]
if round(az_x0) == 90
    x0 = -x0;
elseif round(az_x0) == 270
    x0 = x0;
else
    disp('Warnning ... the sign of x0 is not determined! Please check!');
end

% calculate x1 and y1
rotation_theta = (-wd-az0)*pi/180;
cos_theta = cos(rotation_theta);
sin_theta = sin(rotation_theta);
% x1 = cos_theta*x0 - sin_theta*y0;
% y1 = sin_theta*x0 + cos_theta*y0;
x1 = cos_theta*x0 + sin_theta*y0;
y1 = -sin_theta*x0 + cos_theta*y0;











function d = get_distance(user_lat,user_lon,lat,lon)
% sub function to calculate distance
% out put d is in [m]
R=earthRadius('m');%radius of the earth in meters
lat1=degtorad(user_lat);
lat2=degtorad(lat);
delta_lat=degtorad(lat-user_lat);
delta_lon=degtorad(lon-user_lon);
%a=(sin(delta_lat/2))*(sin(delta_lat/2))+(cos(lat1))*(cos(lat2))*(sin(delta_lon/2))*(sin(delta_lon/2));
a=(sin(delta_lat./2)).*(sin(delta_lat./2))+(cos(lat1)).*(cos(lat2)).*(sin(delta_lon./2)).*(sin(delta_lon./2));
c=2.*asin(sqrt(a));
d=R.*c;
