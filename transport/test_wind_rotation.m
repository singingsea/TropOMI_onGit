function test_wind_rotation()

%lons = -75:-0.5:-78;
%lats = 40:1:45;
%lats = 45:0.5:48;
lons = -78.5:-0.1:-79.5;
lats = 43.5:0.1:44.0;
%lons = 72.5;
%lats = 0.5;

    
[lons_grid,lats_grid]=meshgrid(lons,lats);

%stlon = -79.1;
%stlat = 43.8;
%stlon = -79;
%stlat = 43.6;


    stlat=43.7810; % Downsview
    stlon=-79.4680;
%     stlat2=43.6605;
%     stlon2=-79.39860;
% stlat2=43.700000000000000;
%     stlon2=-79.200000000000000;
    
    stlat2=43.700000000000000;
    stlon2=-79.200000000000000;

u = -25*4;
v = 10*4;

N = size(lons_grid);
for i_lons = 1:N(1)
    for j_lats = 1:N(2)
        try
        [lat_r(i_lons,j_lats),lon_r(i_lons,j_lats),x0(i_lons,j_lats),y0(i_lons,j_lats),x1(i_lons,j_lats),y1(i_lons,j_lats)] = wind_rotation(u,v,stlat,stlon,lats_grid(i_lons,j_lats),lons_grid(i_lons,j_lats));
        
        [lat_r2(i_lons,j_lats),lon_r2(i_lons,j_lats),x0_r2(i_lons,j_lats),y0_r2(i_lons,j_lats),x1_r2(i_lons,j_lats),y1_r2(i_lons,j_lats)] = wind_rotation_v2(u,v,stlat,stlon,lats_grid(i_lons,j_lats),lons_grid(i_lons,j_lats));
        
        %[lat_r2(i_lons,j_lats),lon_r2(i_lons,j_lats),x0_r2(i_lons,j_lats),y0_r2(i_lons,j_lats),x1_r2(i_lons,j_lats),y1_r2(i_lons,j_lats)] = wind_rotation_2sites(u,v,stlat,stlon,stlat2,stlon2,lats_grid(i_lons,j_lats),lons_grid(i_lons,j_lats));
        %[lat_r2(i_lons,j_lats),lon_r2(i_lons,j_lats),x0_r2(i_lons,j_lats),y0_r2(i_lons,j_lats),x1_r2(i_lons,j_lats),y1_r2(i_lons,j_lats)] = wind_rotation_2sites_v2(u,v,stlat,stlon,stlat2,stlon2,lats_grid(i_lons,j_lats),lons_grid(i_lons,j_lats));
        
        %[lat_r2(i_lons,j_lats),lon_r2(i_lons,j_lats),x0,y0,x1,y1] = wind_rotation_2sites(u,v,stlat,stlon,stlat2,stlon2,lats_grid(i_lons,j_lats),lons_grid(i_lons,j_lats))
        catch
            test;
        end
    end
end

lat = reshape(lats_grid,[numel(lats_grid),1]);
lon = reshape(lons_grid,[numel(lons_grid),1]);
lat_r = reshape(lat_r,[numel(lat_r),1]);
lon_r = reshape(lon_r,[numel(lon_r),1]);
lat_r2 = reshape(lat_r2,[numel(lat_r2),1]);
lon_r2 = reshape(lon_r2,[numel(lon_r2),1]);

% check lon and lat
figure;hold all;
scatter(stlon,stlat,'ks','filled');

scatter(stlon2,stlat2,'ms','filled');

scatter(lon,lat,'ro');
scatter(lon_r,lat_r,'go');
scatter(lon_r2,lat_r2,'bo');

scatter(lon(1),lat(1),'r.');
scatter(lon_r(1),lat_r(1),'g.');
scatter(lon_r2(1),lat_r2(1),'b.');

scatter(lon(numel(lon)),lat(numel(lat)),'r*');
scatter(lon_r(numel(lon_r)),lat_r(numel(lat_r)),'g*');
scatter(lon_r2(numel(lon_r2)),lat_r2(numel(lat_r2)),'b*');


%%
y0 = reshape(y0,[numel(y0),1]);
x0 = reshape(x0,[numel(x0),1]);
y1 = reshape(y1,[numel(y1),1]);
x1 = reshape(x1,[numel(x1),1]);
y1_r2 = reshape(y1_r2,[numel(y1_r2),1]);
x1_r2 = reshape(x1_r2,[numel(x1_r2),1]);

% check lon and lat
figure;hold all;
grid on;
%scatter([0,0],'ks','filled');

%scatter(stlon2,stlat2,'ms','filled');

scatter(x0,y0,'ro');
scatter(x1,y1,'go');
scatter(x1_r2,y1_r2,'bo');

scatter(x0(1),y0(1),'r.');
scatter(x1(1),y1(1),'g.');
scatter(x1_r2(1),y1_r2(1),'b.');
quiver(x0(1),y0(1),u,v);
quiver(x1_r2(1),y1_r2(1),u,v);

ii = 45;
scatter(x0(ii),y0(ii),'rx');
scatter(x1(ii),y1(ii),'gx');
scatter(x1_r2(ii),y1_r2(ii),'bx');

scatter(x0(numel(x0)),y0(numel(y0)),'r*');
scatter(x1(numel(x1)),y1(numel(y1)),'g*');
scatter(x1_r2(numel(x1_r2)),y1_r2(numel(y1_r2)),'b*');

        