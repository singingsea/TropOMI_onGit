function TropOMI = get_TropOMI_pixel_transport_location()
% this function need read in the extracted TropOMI data by
% "read_TropOMI_multidays.py", it needs ERA mearged TropOMI data
% this function add up/down wind locations based on wind
addpath(genpath('C:\Users\ZhaoX\Documents\MATLAB\matlab'));
%data_path = 'C:\Projects\TropOMI\data\NO2_output\OFFL\Downsview_ERA\';
%data_path = 'C:\Projects\TropOMI\data\NO2_output\OFFL\Egbert_ERA\';
%output_path = 'C:\Projects\TropOMI\data\NO2_output\OFFL\Egbert_ERA\transport\';

avg_wind = true;% if true, then average wind at frist three pressure levels 
pressure_level = '1000hPa';% only used, if avg_win = false;
delta_t = 1; % trave time in [hour]

% location of interest
site = 'Downsview';
%site = 'Egbert';
%site = 'FortMcKay';
%site = 'StGeorge';
%site = 'FortMcKay_suncrude';

add_2sites_windrotation = false;% if yes, we will also calculate the wind-roation pixel, based on azimuth of two sites
site2 = 'StGeorge';
%site2 = 'Downsview';
%site2 = 'FortMcKay_suncrude';
%site2 = 'FortMcKay';


data_path = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\'];

output_path = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport\'];
%output_path = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport_0dot5\'];
%output_path = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport_1hr_avgwind\'];
mkdir(output_path);


all_files = dir(data_path);% this is all files in the folder
j = 1;
for i=1:length(all_files)
    TF = length(all_files(i).name) == length(['TropOMI_at_' site 'yyyymmdd.csv']);% get only the files have all TropOMI pixels
    if TF
        all_pixel_files(j,:) = all_files(i).name;% get only the files have all TropOMI pixels
        j = j +1;
    end
end

% default lon/lat information for Pandora site
if strcmp(site,'Downsview')
    user_lat=43.7810; % Downsview
    user_lon=-79.4680;
elseif strcmp(site,'Egbert')
    user_lat=44.2300; 
    user_lon=-79.7800;
elseif strcmp(site,'FortMcKay')
    user_lat=57.1836;
    user_lon=-111.6400;
elseif strcmp(site,'StGeorge')
    user_lat=43.6605;
    user_lon=-79.39860;
elseif strcmp(site,'FortMcKay_suncrude')
    user_lat=57.0333;
    user_lon=-111.6167;  
end

if strcmp(site2,'Downsview')
    user_lat2=43.7810; % Downsview
    user_lon2=-79.4680;
elseif strcmp(site2,'Egbert')
    user_lat2=44.2300; 
    user_lon2=-79.7800;
elseif strcmp(site2,'FortMcKay')
    user_lat2=57.1836;
    user_lon2=-111.6400;
elseif strcmp(site2,'StGeorge')
    user_lat2=43.6605;
    user_lon2=-79.39860;
elseif strcmp(site2,'FortMcKay_suncrude')
    user_lat2=57.0333;
    user_lon2=-111.6167;    
end

for i = 1:length(all_pixel_files)
    filename = all_pixel_files(i,:);
    % use the wind-transport method --> for each pixel, we calculate its
    % new location based on given wind speed/direction and transportation
    % time
    data = calculate_pixel_up_down_wind_location(data_path,filename,user_lat,user_lon,pressure_level,delta_t,avg_wind);
    
    % use pixel_roation method ---> calcualte wind-rotated new location for
    % all pixels
    data = pixel_rotation(data,user_lat,user_lon);

    % use pixel_rotation method, but this time we rotate the wind to
    % specified up-down wind, that is defined by azimuth of two sites
    if add_2sites_windrotation
        data = pixel_rotation_2sites(data,user_lat,user_lon,user_lat2,user_lon2);
    end
    
    writetable(data,[output_path filename(1:end-4) '_transport.csv']);
    if i == 1
        data_all = data;
    else
        data_all = [data_all;data];
        data = [];
    end
end

% save data
TropOMI = data_all;
save([output_path 'TropOMI_transport'],'TropOMI');

%%
function data = pixel_rotation(data,user_lat,user_lon)
for i =1:height(data)
    %[lat_r,lon_r] = wind_rotation_v2(data.u_wind,data.v_wind,user_lat,user_lon,data.lat,data.lon);
    [lat_r,lon_r,x0,y0,x1,y1] = wind_rotation(data.u_wind(i,:),data.v_wind(i,:),user_lat,user_lon,data.lat(i,:),data.lon(i,:));    
    data.lat_r(i,:) = lat_r;
    data.lon_r(i,:) = lon_r;
    data.x0(i,:) = x0;
    data.y0(i,:) = y0;
    data.x1(i,:) = x1;
    data.y1(i,:) = y1;
end
data.lat_site1 = repmat(user_lat,height(data),1); % add rotation centre location
data.lon_site1 = repmat(user_lon,height(data),1); % add rotation centre location
%%
function data = pixel_rotation_2sites(data,user_lat,user_lon,user_lat2,user_lon2)
for i =1:height(data)
    %[lat_r,lon_r] = wind_rotation_v2(data.u_wind,data.v_wind,user_lat,user_lon,data.lat,data.lon);
    %[lat_r,lon_r,x0,y0,x1,y1,sites_distance] = wind_rotation_2sites(data.u_wind(i,:),data.v_wind(i,:),user_lat,user_lon,user_lat2,user_lon2,data.lat(i,:),data.lon(i,:));    

    [lat_r,lon_r,x0,y0,x1,y1,x0_st2,y0_st2,x1_st2,y1_st2] = wind_rotation_2sites_v3(data.u_wind(i,:),data.v_wind(i,:),user_lat,user_lon,user_lat2,user_lon2,data.lat(i,:),data.lon(i,:));    
    data.lat_r_2site(i,:) = lat_r;
    data.lon_r_2site(i,:) = lon_r;
    data.x0_2site(i,:) = x0;
    data.y0_2site(i,:) = y0;
    data.x1_2site(i,:) = x1;
    data.y1_2site(i,:) = y1;
    %data.sites_distance(i,:) = sites_distance;
    data.x0_st2_2site(i,:) = x0_st2;
    data.y0_st2_2site(i,:) = y0_st2;
    data.x1_st2_2site(i,:) = x1_st2;
    data.y1_st2_2site(i,:) = y1_st2;
end
data.lat_site2 = repmat(user_lat2,height(data),1); % add the 2nd site location
data.lon_site2 = repmat(user_lon2,height(data),1); % add the 2nd site location
%%
function data = calculate_pixel_up_down_wind_location(data_path,filename,user_lat,user_lon,pressure_level,delta_t,avg_wind)

data = importfile_v4([data_path filename]);

%% test to average the wind
if avg_wind
%     p = 650:50:1000;% pressure levels
%     w=exp(p.*1e-3)./sum(exp(p.*1e-3));% weight function used to average the wind
%     u = w(1).*data.u_650hPa + w(2).*data.u_700hPa + w(3).*data.u_750hPa + w(4).*data.u_800hPa + w(5).*data.u_850hPa + + w(6).*data.u_900hPa + w(7).*data.u_950hPa + w(8).*data.u_1000hPa;
%     v = w(1).*data.v_650hPa + w(2).*data.v_700hPa + w(3).*data.v_750hPa + w(4).*data.v_800hPa + w(5).*data.u_850hPa + + w(6).*data.v_900hPa + w(7).*data.v_950hPa + w(8).*data.v_1000hPa;
     u = (data.u_900hPa + data.u_950hPa + data.u_1000hPa)./3;
     v = (data.v_900hPa + data.v_950hPa + data.v_1000hPa)./3;
else
    % the following line is original wind based on pressure level
    eval(['u = data.u_' pressure_level ';']);% u wind in m/s! 
    eval(['v = data.v_' pressure_level ';']);% v wind in m/s! 
end
%%
%wd = atan2d(u,v);% wind direction [degree]
wd = atan2d(u,v) + 180;% wind direction [degree] <-- to follow the meterological defination; there is a 180 degree offset!
ws = hypot(u,v)/1000*60*60;% wind speed [km/hr]
wd = rem(360+wd, 360);% wind direction in range of [0 360]
wd_opp = rem(wd+180,360);% the opposite direction of wind! will be used to calculate the upwind location

transport_d = ws.*delta_t;% horizontal transport distance [km]

distUnits = 'km';% this unit is used for earth radius
% Convert input distance to earth degrees (Lat, Lon are typicaly given in degrees)
arclen = rad2deg(transport_d/earthRadius(distUnits)); % get the arc length of the horizontal transport
pixel_lat = data.lat;
pixel_lon = data.lon;
% [trans_lat_upw,trans_lon_upw] = reckon(pixel_lat, pixel_lon,arclen,wd_opp);% get the new lat and lon --> up wind
% [trans_lat_dnw,trans_lon_dnw] = reckon(pixel_lat, pixel_lon,arclen,wd);% get the new lat and lon--> down wind
[trans_lat_upw,trans_lon_upw] = reckon(pixel_lat, pixel_lon,arclen,wd);% get the new lat and lon --> up wind; note now we change the wd to follow the meterological defination, so the az is opsite of wind direction!!!!
[trans_lat_dnw,trans_lon_dnw] = reckon(pixel_lat, pixel_lon,arclen,wd_opp);% get the new lat and lon--> down wind

data.trans_lat_upw = trans_lat_upw;% the new lat of the pixel --> up wind direction
data.trans_lon_upw = trans_lon_upw;% the new lon of the pixel --> up wind direction
data.trans_lat_dnw = trans_lat_dnw;% the new lat of the pixel --> down wind direction
data.trans_lon_dnw = trans_lon_dnw;% the new lon of the pixel --> down wind direction

[arclen_upw,az_upw] = distance(user_lat,user_lon,trans_lat_upw,trans_lon_upw);% get the arc length and angle between site and "transported pixel" --> up wind direction
[arclen_dnw,az_dnw] = distance(user_lat,user_lon,trans_lat_dnw,trans_lon_dnw);% get the arc length and angle between site and "transported pixel" --> down wind direction
[arclen_2,az_2] = distance(user_lat,user_lon,data.lat,data.lon);% get the arc length and angle between site and "original pixel"
data.distence2stn_trans_upw = deg2rad(arclen_upw).*earthRadius(distUnits);% distance between station and "transported pixel" --> up wind direction
data.distence2stn_trans_dnw = deg2rad(arclen_dnw).*earthRadius(distUnits);% distance between station and "transported pixel" --> down wind direction
data.distence2stn_original = deg2rad(arclen_2).*earthRadius(distUnits);% distance between station and "original pixel"

% save the calculated wind information
data.u_wind = u;% [m/s]
data.v_wind = v;% [m/s]
data.windspeed = ws; % [km/hr]
data.winddirection = wd;% [degree]
data.transport_distance = transport_d; % [km]



%d = get_distance(user_lat,user_lon,trans_lat,trans_lon);


function d = get_distance(user_lat,user_lon,lat,lon)
% sub function to calculate distance
R=6371000;%radius of the earth in meters
lat1=degtorad(user_lat);
lat2=degtorad(lat);
delta_lat=degtorad(lat-user_lat);
delta_lon=degtorad(lon-user_lon);
%a=(sin(delta_lat/2))*(sin(delta_lat/2))+(cos(lat1))*(cos(lat2))*(sin(delta_lon/2))*(sin(delta_lon/2));
a=(sin(delta_lat./2)).*(sin(delta_lat./2))+(cos(lat1)).*(cos(lat2)).*(sin(delta_lon./2)).*(sin(delta_lon./2));
c=2.*asin(sqrt(a));
d=R.*c;


%%

function TropOMIatDownsview20180422 = importfile_v2(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIATDOWNSVIEW20180422 = IMPORTFILE(FILENAME) Reads data from text
%   file FILENAME for the default selection.
%
%   TROPOMIATDOWNSVIEW20180422 = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   TropOMIatDownsview20180422 = importfile('TropOMI_at_Downsview20180422.csv', 2, 531);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/02/26 11:22:43

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: datetimes (%{yyyy-MM-dd HH:mm:ss.SSS}D)
%   column43: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%{yyyy-MM-dd HH:mm:ss.SSS}D%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
TropOMIatDownsview20180422 = table(dataArray{1:end-1}, 'VariableNames', {'lon','lat','ground_pixel','scanline','qa','time','delta_time','no2_trop','no2_trop_err','amf_trop','amf','snow_cover','albedo_modis_nosnow','albedo_modis_snow','d','u_1000hPa','v_1000hPa','u_950hPa','v_950hPa','u_900hPa','v_900hPa','u_850hPa','v_850hPa','u_800hPa','v_800hPa','u_750hPa','v_750hPa','u_700hPa','v_700hPa','u_650hPa','v_650hPa','no2','no2_err','cloud_frac','no2_strat','no2_strat_err','amf_strat','snow_ice_flag','albedo_no2','albedo_cloud','time_seconds','utc_time','Datenum'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIatDownsview20180422.utc_time=datenum(TropOMIatDownsview20180422.utc_time);

%%
function TropOMIatDownsview20180421 = importfile_v3(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIATDOWNSVIEW20180421 = IMPORTFILE(FILENAME) Reads data from text
%   file FILENAME for the default selection.
%
%   TROPOMIATDOWNSVIEW20180421 = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   TropOMIatDownsview20180421 = importfile('TropOMI_at_Downsview20180421.csv', 2, 541);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/03/09 09:24:15

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: double (%f)
%   column43: double (%f)
%	column44: datetimes (%{yyyy-MM-dd HH:mm:ss.SSS}D)
%   column45: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%{yyyy-MM-dd HH:mm:ss.SSS}D%f%[^\n\r]';


%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
TropOMIatDownsview20180421 = table(dataArray{1:end-1}, 'VariableNames', {'lon','lat','ground_pixel','scanline','qa','time','delta_time','no2_trop','no2_trop_err','amf_trop','amf','snow_cover','albedo_modis_nosnow','albedo_modis_snow','ECCC_NO2','ECCC_AMF','d','u_1000hPa','v_1000hPa','u_950hPa','v_950hPa','u_900hPa','v_900hPa','u_850hPa','v_850hPa','u_800hPa','v_800hPa','u_750hPa','v_750hPa','u_700hPa','v_700hPa','u_650hPa','v_650hPa','no2','no2_err','cloud_frac','no2_strat','no2_strat_err','amf_strat','snow_ice_flag','albedo_no2','albedo_cloud','time_seconds','utc_time','Datenum'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIatDownsview20180421.utc_time=datenum(TropOMIatDownsview20180421.utc_time);


%%
function TropOMIatDownsview20180308 = importfile_v4(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIATDOWNSVIEW20180308 = IMPORTFILE(FILENAME) Reads data from text
%   file FILENAME for the default selection.
%
%   TROPOMIATDOWNSVIEW20180308 = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   TropOMIatDownsview20180308 = importfile('TropOMI_at_Downsview20180308.csv', 2, 15);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/03/15 10:44:35

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: double (%f)
%   column43: double (%f)
%	column44: double (%f)
%   column45: datetimes (%{yyyy-MM-dd HH:mm:ss.SSSSSS+00:00}D)
%	column46: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%{yyyy-MM-dd HH:mm:ss.SSSSSS+00:00}D%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
TropOMIatDownsview20180308 = table(dataArray{1:end-1}, 'VariableNames', {'lon','lat','ground_pixel','scanline','qa','time','delta_time','total_time','no2_trop','no2_trop_err','amf_trop','amf','snow_cover','albedo_modis_nosnow','albedo_modis_snow','ECCC_NO2','ECCC_AMF','d','u_1000hPa','v_1000hPa','u_950hPa','v_950hPa','u_900hPa','v_900hPa','u_850hPa','v_850hPa','u_800hPa','v_800hPa','u_750hPa','v_750hPa','u_700hPa','v_700hPa','u_650hPa','v_650hPa','no2','no2_err','cloud_frac','no2_strat','no2_strat_err','amf_strat','snow_ice_flag','albedo_no2','albedo_cloud','time_seconds','utc_time','Datenum'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIatDownsview20180308.utc_time=datenum(TropOMIatDownsview20180308.utc_time);



