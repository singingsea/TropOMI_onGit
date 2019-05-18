function read_TROPOMI(file_nm)
% this is a simple function to read TROPOMI use MATLAB, currently, the
% functions used to read TROPOMI are phython codes
    
stlat =43.7810; % Downsview
stlon =-79.4680;
    
info =ncinfo(file_nm);

date_str = info.Groups.Name;

% read location info
lat = ncread(file_nm,[date_str '/PRODUCT/latitude']);
lon = ncread(file_nm,[date_str '/PRODUCT/longitude']);


lat_bounds = ncread(file_nm,[date_str '/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds']);
lon_bounds = ncread(file_nm,[date_str '/PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds']);

figure;hold all;
plot(lon_bounds,lat_bounds,'.');
plot(lon,lat,'.');


x0=111.3*(lon-stlon)*cos(stlat/180*3.1415926);
y0=111.3*(lat-stlat);

figure;hold all;
plot(x0,y0,'.');

    