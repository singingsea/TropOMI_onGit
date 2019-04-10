function check_rotation_statistic()
% this function can check the statistic of the rotation method
save_fig = 1;

site = 'Downsview';
Pandora_no = '103';
% site = 'FortMcKay';
% Pandora_no = '122';
% site = 'StGeorge';
% Pandora_no = '145';
% site = 'Egbert';
% Pandora_no = '108';
method = 'rotation';%'rotation','twosites_rotation','wind_driven'


save_fig = 1;
plot_path = ['C:\Projects\TropOMI\plots\Wind_transport\' Pandora_no '\' method '\statistic\'];
mkdir(plot_path);
  
TropOMI_data = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport\TropOMI_transport.mat'];
addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab');

load(TropOMI_data);
DU = 2.6870e+16;

%% QC
% TF = TropOMI.qa == 1;
% TropOMI = TropOMI(TF,:);
%%
fig_name = 'original_density_all';
figure;dscatter(TropOMI.x0,TropOMI.y0);
xlabel('x distance [km]');
ylabel('y distance [km]');
xlim([-120 120]);
ylim([-120 120]);
title('original position of TropOMI pixel');
colorbar;
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'rotated_density_all';
figure;dscatter(TropOMI.x1,TropOMI.y1);
xlabel('x distance [km]');
ylabel('y distance [km]');
xlim([-120 120]);
ylim([-120 120]);
title('rotated position of TropOMI pixel');
colorbar;
print_setting(0,save_fig,[plot_path fig_name]);

%TF = (abs(TropOMI.x0) <= 30) & (abs(TropOMI.y0) <= 30);
TF = (TropOMI.x0.^2 + TropOMI.y0.^2).^0.5 <=50;
data = TropOMI(TF,:);

%% original statistic
fig_name = 'original_density_selected';
figure;
hold all;
dscatter(data.x0,data.y0);
xlabel('x distance [km]');
ylabel('y distance [km]');
title('original position of TropOMI pixel');
colorbar;

dis = 7;
TF0 = (abs(data.x0)<=dis) & (abs(data.y0)<=dis);
n1 = sum(TF0);
plot([-dis,-dis],[-dis,dis],'b--');
plot([dis,dis],[-dis,dis],'b--');
plot([-dis,dis],[dis,dis],'b--');
plot([-dis,dis],[-dis,-dis],'b--');
text(dis,dis,num2str(n1),'Color','blue');
print_setting(0,save_fig,[plot_path fig_name]);

%% rotated statistic
fig_name = 'rotated_density_selected';
figure;
hold all;
dscatter(data.x1,data.y1);
xlabel('x distance [km]');
ylabel('y distance [km]');
title('rotated position of TropOMI pixel');
colorbar;

dis = 7;
TF1 = abs(data.x1)<=dis;
n1 = sum(TF1);
plot([-dis,-dis],[-50,50],'b--');
plot([dis,dis],[-50,50],'b--');
text(dis,dis,num2str(n1),'Color','blue');

TF1_y = (abs(data.y1)<=dis) & TF1;
n1 = sum(TF1_y);
plot([-dis,dis],[dis,dis],'b--');
plot([-dis,dis],[-dis,-dis],'b--');
text(dis,-dis,num2str(n1),'Color','blue');

dis = 14;
TF1_y = (abs(data.y1)<=dis) & TF1;
n1 = sum(TF1_y);
plot([-7,7],[dis,dis],'b--');
plot([-7,7],[-dis,-dis],'b--');
text(7,-dis,num2str(n1),'Color','blue');

dis = 21;
TF1_y = (abs(data.y1)<=dis) & TF1;
n1 = sum(TF1_y);
plot([-7,7],[dis,dis],'b--');
plot([-7,7],[-dis,-dis],'b--');
text(7,-dis,num2str(n1),'Color','blue');

dis = 28;
TF1_y = (abs(data.y1)<=dis) & TF1;
n1 = sum(TF1_y);
plot([-7,7],[dis,dis],'b--');
plot([-7,7],[-dis,-dis],'b--');
text(7,-dis,num2str(n1),'Color','blue');

dis = 35;
TF1_y = (abs(data.y1)<=dis) & TF1;
n1 = sum(TF1_y);
plot([-7,7],[dis,dis],'b--');
plot([-7,7],[-dis,-dis],'b--');
text(7,-dis,num2str(n1),'Color','blue');



dis = 14;
TF2 = abs(data.x1)<=dis;
n1 = sum(TF2);
plot([-dis,-dis],[-50,50],'r--');
plot([dis,dis],[-50,50],'r--');
text(dis,dis,num2str(n1),'Color','red');

dis = 21;
TF3 = abs(data.x1)<=dis;
n1 = sum(TF3);
plot([-dis,-dis],[-50,50],'k--');
plot([dis,dis],[-50,50],'k--');
text(dis,dis,num2str(n1),'Color','black');
print_setting(0,save_fig,[plot_path fig_name]);

%% location of three dis
fig_name = 'original_locations_in_three_bins';
figure;hold all;
plot(data.x0(TF3,:),data.y0(TF3,:),'.');
plot(data.x0(TF2,:),data.y0(TF2,:),'.');
plot(data.x0(TF1,:),data.y0(TF1,:),'.');
xlabel('x distance [km]');
ylabel('y distance [km]');
title('original position of TropOMI pixel (within bins)');
legend({'group 3','group 2','group 1'});
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'original_locations_group1_quiver';
figure;hold all;
quiver(data.x0(TF1,:),data.y0(TF1,:),data.u_wind(TF1,:),data.v_wind(TF1,:));
xlabel('x distance [km]');
ylabel('y distance [km]');
title('original position of TropOMI pixel (within bins)');
print_setting(0,save_fig,[plot_path fig_name]);


% figure;hold all;
% scatter(data.x0(TF1,:),data.y0(TF1,:),10,(data.no2_trop(TF1,:)+data.no2_strat(TF1,:)));
% xlabel('x distance [km]');
% ylabel('y distance [km]');
% title('original position of TropOMI pixel (within bins)');

pixel_averging_xy(data,'rotation',site);
h = gca;
fig_name = 'rotated_pixel_averaged';
title('pixel averaged (rotated)');
print_setting(0,save_fig,[plot_path fig_name]);

pixel_averging_xy(data,'original location',site);
h = gca;
fig_name = 'original_pixel_averaged';
title('pixel averaged (original)');
print_setting(0,save_fig,[plot_path fig_name]);

function pixel_averging_xy(data,method,site)
DU = 2.6870e+16;
%% pixel averageing
if strcmp(method,'rotation')
    data.x_pixel = data.x1;
    data.y_pixel = data.y1;
elseif strcmp(method,'twosites rotation')
    data.x_pixel = data.x1_2site;
    data.y_pixel = data.y1_2site;
elseif strcmp(method,'original location')
    data.x_pixel = data.x0;
    data.y_pixel = data.y0;
end


% x_step = 0.08;
% y_step = 0.08;
x_step = 3.5;
y_step = 3.5;
xs = min(data.x_pixel):x_step:max(data.x_pixel);
ys = min(data.y_pixel):y_step:max(data.y_pixel);
[lon_grids,lat_grids] = meshgrid(xs,ys);

mean_of_pixel = NaN(length(ys),length(xs));
mean_of_pixel_x1 = NaN(length(ys),length(xs));
mean_of_pixel_y1 = NaN(length(ys),length(xs));
for i=1:length(xs)
    for j =1:length(ys)
        pixel_lon = xs(i);
        pixel_lat = ys(j);
        lon_grid_Ledge = pixel_lon - x_step/2;
        lon_grid_Redge = pixel_lon + x_step/2;
        lat_grid_Uedge = pixel_lat + y_step/2;
        lat_grid_Dedge = pixel_lat - y_step/2;
        
        TF_lon = (data.x_pixel >= lon_grid_Ledge) & (data.x_pixel < lon_grid_Redge);
        TF_lat = (data.y_pixel <= lat_grid_Uedge) & (data.y_pixel > lat_grid_Dedge);
        TF = TF_lon & TF_lat;
        data_at_1pixel = data(TF,:);
        
        mean_of_pixel(j,i) = mean(data_at_1pixel.no2);
        mean_of_pixel_x1(j,i) = mean(data_at_1pixel.x1);
        mean_of_pixel_y1(j,i) = mean(data_at_1pixel.y1);
    end
end


fig_name = ['pxiel_averaging_' site '_' method '_distanceview'];
figure;hold all;title(method);
% quiver(data.x0(1),data.y0(1),data.u_wind(1),data.v_wind(1));
% quiver(data.x1(1),data.y1(1),data.u_wind(1),data.v_wind(1));
% if strcmp(method,'twosites rotation')
%     scatter(data.x0_2site,data.y0_site2,'ks','MarkerFaceColor',[1 1 1]);
% end
surf(mean_of_pixel_x1,mean_of_pixel_y1,zeros(size(mean_of_pixel_x1)),mean_of_pixel./DU,'AlphaData',1,'FaceAlpha',1,'edgecolor','none');
view(2);
colormap(jet);
caxis([0.15 0.35]);
xlabel('Distance [km]');ylabel('Distance [km]');
h=colorbar;
ylabel(h,'NO_2 VCD [DU]');
% print_setting(0,save_fig,[plot_path fig_name]);
% ylim([56.2 58.2]);
% xlim([-113.5 -109.7]);

figure;hold all;
[lon_grids_highres,lat_grids_highres] = meshgrid(min(data.x_pixel):2:max(data.x_pixel),min(data.y_pixel):2:max(data.y_pixel));
qC = griddata(xs,ys,mean_of_pixel./DU,lon_grids_highres,lat_grids_highres);
%mesh(lons,lats,qC);
surf(lon_grids_highres,lat_grids_highres,zeros(size(lon_grids_highres)),qC,'AlphaData',1,'FaceAlpha',1,'edgecolor','none');
view(2);
colormap(jet);
caxis([0.15 0.35]);
xlabel('Distance [km]');ylabel('Distance [km]');
h=colorbar;
ylabel(h,'NO_2 VCD [DU]');
