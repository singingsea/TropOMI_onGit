function pixel_averaging()
%load('C:\Projects\TropOMI\data\NO2_output\OFFL\FortMcKay_suncrude_ERA\transport\TropOMI_transport.mat');
%load('C:\Projects\TropOMI\data\NO2_output\OFFL\FortMcKay_ERA\transport\TropOMI_transport.mat');
addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab\');
%site = 'Downsview';
site = 'FortMcKay';
method = 'original location';%'rotation','twosites rotation','original location'
save_fig = 1;
data_path = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport\'];
load([data_path 'TropOMI_transport.mat']);
plot_path = [data_path 'pixel_averaging\'];mkdir(plot_path);
data = TropOMI;
DU = 2.6870e+16;
if strcmp(method,'rotation')
    data.lon_pixel = data.lon_r;
    data.lat_pixel = data.lat_r;
elseif strcmp(method,'twosites rotation')
    data.lon_pixel = data.lon_r_2site;
    data.lat_pixel = data.lat_r_2site;
elseif strcmp(method,'original location')
    data.lon_pixel = data.lon;
    data.lat_pixel = data.lat;
end
% simple QC
TF_qa = data.qa < 1;
data(TF_qa,:) = [];

% wind speed TF
TF_ws = (data.windspeed < 30);
data(TF_ws,:) = [];

% lon_step = 0.05;
% lat_step = 0.025;
% lon_step = 0.05/3;
% lat_step = 0.025/3;
lon_step = 0.08;
lat_step = 0.08;
lons = min(data.lon_pixel):lon_step:max(data.lon_pixel);
lats = min(data.lat_pixel):lat_step:max(data.lat_pixel);
[lon_grids,lat_grids] = meshgrid(lons,lats);

mean_of_pixel = NaN(length(lats),length(lons));
mean_of_pixel_x1 = NaN(length(lats),length(lons));
mean_of_pixel_y1 = NaN(length(lats),length(lons));
for i=1:length(lons)
    for j =1:length(lats)
        pixel_lon = lons(i);
        pixel_lat = lats(j);
        lon_grid_Ledge = pixel_lon - lon_step/2;
        lon_grid_Redge = pixel_lon + lon_step/2;
        lat_grid_Uedge = pixel_lat + lat_step/2;
        lat_grid_Dedge = pixel_lat - lat_step/2;
        
        TF_lon = (data.lon_pixel >= lon_grid_Ledge) & (data.lon_pixel < lon_grid_Redge);
        TF_lat = (data.lat_pixel <= lat_grid_Uedge) & (data.lat_pixel > lat_grid_Dedge);
        TF = TF_lon & TF_lat;
        data_at_1pixel = data(TF,:);
        
        mean_of_pixel(j,i) = mean(data_at_1pixel.no2);
        mean_of_pixel_x1(j,i) = mean(data_at_1pixel.x1);
        mean_of_pixel_y1(j,i) = mean(data_at_1pixel.y1);
    end
end

%% pixel averaging
fig_name = ['pxiel_averaging_' site '_' method];
figure;hold all;title(method);
plot_sites(site);
if strcmp(method,'twosites rotation')
    scatter(data.lon_site2,data.lat_site2,'ks','MarkerFaceColor',[1 1 1]);
end
surf(lon_grids,lat_grids,zeros(size(lon_grids)),mean_of_pixel./DU,'AlphaData',0.5,'FaceAlpha',0.5,'edgecolor','none');
view(2);
colormap(jet);
caxis([0.1 0.3]);
xlabel('Lon.');ylabel('Lat.');
h=colorbar;
ylabel(h,'NO_2 VCD [DU]');
% ylim([56.2 58.2]);
% xlim([-113.5 -109.7]);

addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab\zoharby-plot_google_map-6a7991e');
plot_google_map('apiKey', 'AIzaSyBcdJDmLtps3qRBBtoufM3Tc3IsWirS-k0','maptype','satellite','language','fr','showLabels',1,'MapScale',1);
print_setting(0,save_fig,[plot_path fig_name]);
%% meshed to high res
fig_name = ['pxiel_averaging_' site '_' method '_meshed_highres'];
figure;hold all;title(method);
plot_sites(site);
[lon_grids_highres,lat_grids_highres] = meshgrid(min(data.lon_pixel):0.001:max(data.lon_pixel),min(data.lat_pixel):0.001:max(data.lat_pixel));
qC = griddata(lons,lats,mean_of_pixel./DU,lon_grids_highres,lat_grids_highres);
%mesh(lons,lats,qC);
surf(lon_grids_highres,lat_grids_highres,zeros(size(lon_grids_highres)),qC,'AlphaData',0.3,'FaceAlpha',0.3,'edgecolor','none');
view(2);

colormap(jet);
caxis([0.1 0.3]);
xlabel('Lon.');ylabel('Lat.');
h=colorbar;
ylabel(h,'NO_2 VCD [DU]');

plot_google_map('apiKey', 'AIzaSyBcdJDmLtps3qRBBtoufM3Tc3IsWirS-k0','maptype','satellite','language','fr','showLabels',1,'MapScale',1);
print_setting(0,save_fig,[plot_path fig_name]);
%% distance view
fig_name = ['pxiel_averaging_' site '_' method '_distanceview'];
figure;hold all;title(method);
quiver(data.x0(1),data.y0(1),data.u_wind(1),data.v_wind(1));
quiver(data.x1(1),data.y1(1),data.u_wind(1),data.v_wind(1));
if strcmp(method,'twosites rotation')
    scatter(data.x0_2site,data.y0_site2,'ks','MarkerFaceColor',[1 1 1]);
end
surf(mean_of_pixel_x1,mean_of_pixel_y1,zeros(size(mean_of_pixel_x1)),mean_of_pixel./DU,'AlphaData',0.5,'FaceAlpha',0.5,'edgecolor','none');
view(2);
colormap(jet);
caxis([0.1 0.3]);
xlabel('Distance [km]');ylabel('Distance [km]');
h=colorbar;
ylabel(h,'NO_2 VCD [DU]');
print_setting(0,save_fig,[plot_path fig_name]);
% ylim([56.2 58.2]);
% xlim([-113.5 -109.7]);

%%
function plot_sites(site)
if strcmp(site,'Downsview') || strcmp(site,'StGeorge') || strcmp(site,'Egbert')
    user_lat=43.7810; % Downsview
    user_lon=-79.4680;
    plot(user_lon,user_lat,'ko','MarkerFaceColor',[1 1 1]);  
    user_lat=44.2300;  % St George
    user_lon=-79.7800;
    plot(user_lon,user_lat,'ko','MarkerFaceColor',[1 1 1]);
    user_lat=43.6605;% Egbert
    user_lon=-79.39860;
    plot(user_lon,user_lat,'ko','MarkerFaceColor',[1 1 1]);
    
    user_lat=43+47/60; % UTSC
    user_lon=-(79 + 11/60);
    plot(user_lon,user_lat,'ko','MarkerFaceColor',[1 1 1]);
elseif strcmp(site,'FortMcKay')
    user_lat=57.1836;
    user_lon=-111.6400;
    plot(user_lon,user_lat,'ko','MarkerFaceColor',[1 1 1]);
end

        