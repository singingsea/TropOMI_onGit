function bin_by_wd()
type = '1site'; % '2sites'
%type = '2sites'; % '2sites'
addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab\');
method = 'rotation';%'rotation','twosites_rotation','wind_driven'
Pandora_no = 103;
%method = 'original location';%'rotation','twosites rotation','original location'
save_fig = 0;
if strcmp(method, 'rotation')
    data_path = ['C:\Projects\TropOMI\plots\Wind_transport\' num2str(Pandora_no) '\rotation\' ];
    load([data_path type 'TropOMI_Pandora_transport.mat']);
    plot_path = [data_path '\pixel_averaging\'];mkdir(plot_path);
elseif strcmp(method, 'twosites_rotation')
    data_path = ['C:\Projects\TropOMI\plots\Wind_transport\' num2str(Pandora_no) '\twosites_rotation\'];
    load([data_path type 'TropOMI_Pandora_transport.mat']);
    plot_path = [data_path '\pixel_averaging\2sites\'];mkdir(plot_path);
end




%load('C:\Projects\TropOMI\plots\Wind_transport\103\rotation\1siteTropOMI_Pandora_transport.mat');
%load('C:\Projects\TropOMI\plots\Wind_transport\145\rotation\1siteTropOMI_Pandora_transport.mat');
%load('C:\Projects\TropOMI\plots\Wind_transport\122\rotation\1siteTropOMI_Pandora_transport.mat');
%load('C:\Projects\TropOMI\plots\Wind_transport\108\rotation\1siteTropOMI_Pandora_transport.mat');
%load('C:\Projects\TropOMI\plots\Wind_transport\104\rotation\1siteTropOMI_Pandora_transport.mat');
%load('C:\Projects\TropOMI\plots\Wind_transport\122\twosites_rotation\2sitesTropOMI_Pandora_transport.mat');
data = C_all;

%% filter
%TF = (data.x1 <= 10) & (abs(data.y1) <= 10) & (abs(data.arrival_time) <= 1);
%TF = (data.x1 <= 7) & (abs(data.y1) <= 7*5) & (abs(data.arrival_time) <= 5);
%TF = (data.x1 <= 7) & (data.y1 <= 7*5) & (data.y1 > 0 ) & (abs(data.arrival_time) <= 5);
TF = (data.x1 <= 7) & (data.y1 >= - 7*5) & (data.y1 < 0 ) & (abs(data.arrival_time) <= 5);
% % %TF = (data.x1 <= 10) & (abs(data.y1) <=800) & (abs(data.arrival_time) <= 50);
data(~TF,:) = [];

%data.NO2_VCD = data.NO2_VCD.*0.7;
%data.NO2_VCD = data.NO2_VCD.*0.6;
%data.NO2_VCD = data.NO2_VCD_site2;

% data.NO2_VCD = data.NO2_VCD_site2;
% data.NO2_VCD_err = data.NO2_VCD_err_site2;

%% bin by wind direction
DU = 2.6870e+16;
N_bins = 12; % number of bins
%N_bins = 36; % number of bins
data_bins = NaN(height(data),N_bins);% this bin is to store TropOMI NO2
data_bins_amf = NaN(height(data),N_bins);% this bin is to store TropOMI NO2 amf
data_bins_no2_err = NaN(height(data),N_bins);% this bin is to store TropOMI NO2 err
data_bins_albedo_no2 = NaN(height(data),N_bins);% this bin is to store TropOMI no2 albedo
data_bins_ground_pixel = NaN(height(data),N_bins);% this bin is to store TropOMI ground_pixel
data_bins_no2_trop = NaN(height(data),N_bins);% this bin is to store TropOMI no2_trop
data_bins_amf_trop = NaN(height(data),N_bins);% this bin is to store TropOMI no2_trop amf
data_bins_no2_strat = NaN(height(data),N_bins);% this bin is to store TropOMI no2_start
data_bins_pandora = NaN(height(data),N_bins);% this bin is to store Pandora NO2 at site 1
data_bins_pandora2 = NaN(height(data),N_bins);% this bin is to store Pandora NO2 at site 2
for i = 1:N_bins
    wd_min = (i-1)*(360/N_bins);
    wd_max = (i)*(360/N_bins);
    TF = (data.winddirection>= wd_min) & (data.winddirection < wd_max);
    data_bins(TF,i) = data.no2(TF)./DU;
    data_bins_amf(TF,i) = data.amf(TF);
    data_bins_no2_err(TF,i) = data.no2_err(TF)./DU;
    data_bins_albedo_no2(TF,i) = data.albedo_no2(TF);
    data_bins_ground_pixel(TF,i) = data.ground_pixel(TF);
    data_bins_no2_trop(TF,i) = data.no2_trop(TF)./DU;
    data_bins_amf_trop(TF,i) = data.amf_trop(TF);
    data_bins_no2_strat(TF,i) = data.no2_strat(TF)./DU;
    N(i) = sum(TF);
    data_bins_pandora(TF,i) = data.NO2_VCD(TF);
    if strcmp(method, 'twosites_rotation')
        data_bins_pandora2(TF,i) = data.NO2_VCD_site2(TF);
    end
end

figure;fig_name = 'Mean_NO2_bin_by_wd_TropOMI';
boxplot(data_bins);
xticks([1 3 5 7 9 11]);
xticklabels({'15','75','135','195','255','315'});
ylabel('TropOMI NO_2 mean VCD [DU]');
xlabel('Wind direction [degree]');
ylim([0 1]);
print_setting(0,save_fig,[plot_path fig_name]);

figure;fig_name = 'Mean_NO2_bin_by_wd_Pandora';
boxplot(data_bins_pandora);    
xticks([1 3 5 7 9 11]);
xticklabels({'15','75','135','195','255','315'});
ylabel('Pandora NO_2 mean VCD [DU]');
xlabel('Wind direction [degree]');
ylim([0 1]);
print_setting(0,save_fig,[plot_path fig_name]);

if strcmp(method, 'twosites_rotation')
    figure;fig_name = 'Mean_NO2_bin_by_wd_Pandora_2ndsite';
    boxplot(data_bins_pandora2);    
    xticks([1 3 5 7 9 11]);
    xticklabels({'15','75','135','195','255','315'});
    ylabel('Pandora NO_2 mean VCD [DU]');
    xlabel('Wind direction [degree]');
    ylim([0 1]);
    print_setting(0,save_fig,[plot_path fig_name]);
end

%% major wind bin plot
figure;hold all;fig_name = 'Mean_NO2_bin_by_wd';
y=mean(data_bins,'omitnan');
y_err = std(data_bins,'omitnan')./(N).^0.5;
errorbar(1:N_bins,y,y_err);

y=mean(data_bins_no2_trop,'omitnan');
y_err = std(data_bins_no2_trop,'omitnan')./(N).^0.5;
errorbar(1:N_bins,y,y_err);

y=mean(data_bins_no2_strat,'omitnan');
y_err = std(data_bins_no2_strat,'omitnan')./(N).^0.5;
errorbar(1:N_bins,y,y_err);

y_pandora=mean(data_bins_pandora,'omitnan');
y_pandora_err = std(data_bins_pandora,'omitnan')./(N).^0.5;
errorbar(1:N_bins,y_pandora,y_pandora_err);

if strcmp(method, 'twosites_rotation')
    y_pandora2=mean(data_bins_pandora2,'omitnan');
    y_pandora2_err = std(data_bins_pandora2,'omitnan')./(N).^0.5;
    errorbar(1:N_bins,y_pandora2,y_pandora2_err);
end
% legend({'TropOMI',['Pandora ' num2str(Pandora_no)]});
legend({'TropOMI','TropOMI trop.','TropOMI start.',['Pandora ' num2str(Pandora_no)]});
xticks([1 3 5 7 9 11]);
xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
ylabel('NO_2 mean VCD [DU]');
xlabel('Wind direction');
ylim([0 1]);
print_setting(0,save_fig,[plot_path fig_name]);

%% scatter plot
fig_name = 'Mean_NO2_bin_by_wd_scatter';
y=mean(data_bins,'omitnan');
tf = isnan(y);
y = y(~tf);
y_pandora = y_pandora(~tf);
line_fits(y_pandora',y');
xlabel(['Pandora'  num2str(Pandora_no) ' mean VCD [DU]']);
ylabel('TropOMI mean VCD [DU]');
new_N_bins = sum(~tf);
scatter(y_pandora',y',10,(1:new_N_bins)','filled');
print_setting(0,save_fig,[plot_path fig_name]);

%% 2nd wind bin plot for amf
figure;hold all;fig_name = 'Mean_NO2_amf_albedo_bin_by_wd';
y=mean(data_bins_amf,'omitnan');
y_err = std(data_bins_amf,'omitnan')./(N).^0.5;
%plot(1:N_bins,y);
yyaxis left
errorbar(1:N_bins,y,y_err);
ylabel('NO_2 AMF');
ylim([1 3.1]);
yyaxis right
y_albedo_no2=mean(data_bins_albedo_no2,'omitnan');
y_albedo_no2_err = std(data_bins_albedo_no2,'omitnan')./(N).^0.5;
errorbar(1:N_bins,y_albedo_no2,y_albedo_no2_err);
legend({'TropOMI',['Pandora ' num2str(Pandora_no)]});
xticks([1 3 5 7 9 11]);
xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
ylabel('NO_2 albedo');
xlabel('Wind direction');
ylim([0 0.1]);
print_setting(0,save_fig,[plot_path fig_name]);


%% 3rd wind bin plot for pixel number and no2 err
figure;hold all;fig_name = 'pixel_number_and_no2err_bin_by_wd';
y=mean(data_bins_ground_pixel,'omitnan');
y_err = std(data_bins_ground_pixel,'omitnan')./(N).^0.5;
%plot(1:N_bins,y);
yyaxis left
errorbar(1:N_bins,y,y_err);
ylabel('mean TropOMI ground pixel number');
%ylim([1 3.1]);
yyaxis right
y2=mean(data_bins_no2_err,'omitnan');
y2_err = std(data_bins_no2_err,'omitnan')./(N).^0.5;
errorbar(1:N_bins,y2,y2_err);
legend({'TropOMI',['Pandora ' num2str(Pandora_no)]});
xticks([1 3 5 7 9 11]);
xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
ylabel('NO_2 err [DU]');
xlabel('Wind direction');
%ylim([0 0.1]);
print_setting(0,save_fig,[plot_path fig_name]);

%% 4th wind bin plot for trop. amf and trop. column
figure;hold all;fig_name = 'trop_amf_and_trop_column_by_wd';
y=mean(data_bins_amf_trop,'omitnan');
y_err = std(data_bins_amf_trop,'omitnan')./(N).^0.5;
%plot(1:N_bins,y);
yyaxis left
errorbar(1:N_bins,y,y_err);
ylabel('TropOMI NO_2 trop. AMF');
ylim([0 2]);
yyaxis right
y2=mean(data_bins_no2_trop,'omitnan');
y2_err = std(data_bins_no2_trop,'omitnan')./(N).^0.5;
errorbar(1:N_bins,y2,y2_err);
legend({'TropOMI',['Pandora ' num2str(Pandora_no)]});
xticks([1 3 5 7 9 11]);
xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
ylabel('NO_2 trop. column [DU]');
xlabel('Wind direction');
ylim([0 0.4]);
print_setting(0,save_fig,[plot_path fig_name]);

%% histogram
fig_name = 'histogram_of_ground_pixel';
figure;histogram(data.ground_pixel,'BinWidth',15);ylabel('f');xlabel('ground pixel');
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_NO2_err';
figure;histogram(data.no2_err,'BinWidth',0.01);ylabel('f');xlabel('TropOMI NO_2 err [DU]');
xlim([0 0.15]);
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_NO2_AMF';
figure;histogram(data.amf,'BinWidth',0.2);ylabel('f');xlabel('TropOMI NO_2 AMF');
xlim([1 4.5]);
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_NO2_Albedo';
figure;histogram(data.albedo_no2,'BinWidth',0.01);ylabel('f');xlabel('TropOMI NO_2 albedo');
xlim([0 0.1]);
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_qa';
figure;histogram(data.qa,'BinWidth',0.01);ylabel('f');xlabel('TropOMI qa');
print_setting(0,save_fig,[plot_path fig_name]);



