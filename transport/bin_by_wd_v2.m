function bin_by_wd_v2()
use_sum_no2 = true;
type = '1site'; % '2sites'
%type = '2sites'; % '2sites'
addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab\');
method = 'rotation';%'rotation','twosites_rotation','wind_driven'
% site = 'Downsview';
% Pandora_no = '104';
site = 'FortMcKay';
Pandora_no = '122';
% site = 'StGeorge';
% Pandora_no = '145';
% site = 'Egbert';
% Pandora_no = '108';

use_eccc_no2 = true;
%Pandora2_no = 122;
Pandora2_no = Pandora_no;
    
save_fig = 1;
if strcmp(method, 'rotation')
    if use_eccc_no2
        data_path = ['C:\Projects\TropOMI\plots\Wind_transport\v9_smooth2OMI\' num2str(Pandora_no) '\rotation_eccc\' ];
        data2_path = ['C:\Projects\TropOMI\plots\Wind_transport\v9_smooth2OMI\' num2str(Pandora2_no) '\rotation_eccc\' ];
    else
        data_path = ['C:\Projects\TropOMI\plots\Wind_transport\v9_smooth2OMI\' num2str(Pandora_no) '\rotation_knmi\' ];
        data2_path = ['C:\Projects\TropOMI\plots\Wind_transport\v9_smooth2OMI\' num2str(Pandora2_no) '\rotation_knmi\' ];
    end
    load([data_path type 'TropOMI_Pandora_transport.mat']);
    data = C_all;
    load([data2_path type 'TropOMI_Pandora_transport.mat']);
    data2 = C_all;
    %plot_path = [data_path '\pixel_averaging\'];mkdir(plot_path);
    %plot_path = [data_path '\pixel_averaging_newbins\'];mkdir(plot_path);
    %plot_path = [data_path '\pixel_averaging_comparetoOMI\'];mkdir(plot_path);
    plot_path = [data_path '\pixel_averaging_correctedWD\'];mkdir(plot_path);
elseif strcmp(method, 'twosites_rotation')
    data_path = ['C:\Projects\TropOMI\plots\Wind_transport\' num2str(Pandora_no) '\twosites_rotation\'];
    load([data_path type 'TropOMI_Pandora_transport.mat']);
    data = C_all;
    plot_path = [data_path '\pixel_averaging\2sites\'];mkdir(plot_path);
end

if use_sum_no2
    data.no2 = data.no2_strat+data.no2_trop;
    data2.no2 = data2.no2_strat+data2.no2_trop;
end

%% filter

% TF = (abs(data.x1) <= 5) & (abs(data.y1) <= 30) & (abs(data.arrival_time) <= 1);
TF = (abs(data.x1) <= 20) & (abs(data.y1) <= 30) & (abs(data.arrival_time) <= 1);
data(~TF,:) = [];

% TF = (abs(data2.x1) <= 5) & (abs(data2.y1) <= 30) & (abs(data2.arrival_time) <= 1);
TF = (abs(data2.x1) <= 20) & (abs(data2.y1) <= 30) & (abs(data2.arrival_time) <= 1);

data2(~TF,:) = [];

%data.NO2_VCD = data.NO2_VCD.*0.7;
%data.NO2_VCD = data.NO2_VCD.*0.6;
%data.NO2_VCD = data.NO2_VCD_site2;

% data.NO2_VCD = data.NO2_VCD_site2;
% data.NO2_VCD_err = data.NO2_VCD_err_site2;
data_bins = make_bins(data,method,plot_path,save_fig,Pandora_no);
combined_histogram(data,data2,Pandora_no,Pandora2_no,save_fig,plot_path);

%% generate date list
datelist = unique(fix(datenum(data.utc_time)));

%% save data
save([plot_path 'data_bins'],'data_bins');
save([plot_path 'datelist'],'datelist');

function data_bins = make_bins(data,method,plot_path,save_fig,Pandora_no)
%% bin by wind direction
DU = 2.6870e+16;
N_bins = 12; % number of bins
%N_bins = 36; % number of bins
data_bins = table;
data_bins.no2 = NaN(height(data),N_bins+1);% this bin is to store TropOMI NO2
data_bins.amf = NaN(height(data),N_bins+1);% this bin is to store TropOMI NO2 amf
data_bins.no2_err = NaN(height(data),N_bins+1);% this bin is to store TropOMI NO2 err
data_bins.albedo_no2 = NaN(height(data),N_bins+1);% this bin is to store TropOMI no2 albedo
data_bins.ground_pixel = NaN(height(data),N_bins+1);% this bin is to store TropOMI ground_pixel
data_bins.no2_trop = NaN(height(data),N_bins+1);% this bin is to store TropOMI no2_trop
data_bins.amf_trop = NaN(height(data),N_bins+1);% this bin is to store TropOMI no2_trop amf
data_bins.no2_strat = NaN(height(data),N_bins+1);% this bin is to store TropOMI no2_start
data_bins.pandora = NaN(height(data),N_bins+1);% this bin is to store Pandora NO2 at site 1
data_bins.pandora2 = NaN(height(data),N_bins+1);% this bin is to store Pandora NO2 at site 2

% data_bins.no2 = NaN(height(data),N_bins);% this bin is to store TropOMI NO2
% data_bins.amf = NaN(height(data),N_bins);% this bin is to store TropOMI NO2 amf
% data_bins.no2_err = NaN(height(data),N_bins);% this bin is to store TropOMI NO2 err
% data_bins.albedo_no2 = NaN(height(data),N_bins);% this bin is to store TropOMI no2 albedo
% data_bins.ground_pixel = NaN(height(data),N_bins);% this bin is to store TropOMI ground_pixel
% data_bins.no2_trop = NaN(height(data),N_bins);% this bin is to store TropOMI no2_trop
% data_bins.amf_trop = NaN(height(data),N_bins);% this bin is to store TropOMI no2_trop amf
% data_bins.no2_strat = NaN(height(data),N_bins);% this bin is to store TropOMI no2_start
% data_bins.pandora = NaN(height(data),N_bins);% this bin is to store Pandora NO2 at site 1
% data_bins.pandora2 = NaN(height(data),N_bins);% this bin is to store Pandora NO2 at site 2


% data_bins = NaN(height(data),N_bins);% this bin is to store TropOMI NO2
% data_bins_amf = NaN(height(data),N_bins);% this bin is to store TropOMI NO2 amf
% data_bins_no2_err = NaN(height(data),N_bins);% this bin is to store TropOMI NO2 err
% data_bins_albedo_no2 = NaN(height(data),N_bins);% this bin is to store TropOMI no2 albedo
% data_bins_ground_pixel = NaN(height(data),N_bins);% this bin is to store TropOMI ground_pixel
% data_bins_no2_trop = NaN(height(data),N_bins);% this bin is to store TropOMI no2_trop
% data_bins_amf_trop = NaN(height(data),N_bins);% this bin is to store TropOMI no2_trop amf
% data_bins_no2_strat = NaN(height(data),N_bins);% this bin is to store TropOMI no2_start
% data_bins_pandora = NaN(height(data),N_bins);% this bin is to store Pandora NO2 at site 1
% data_bins_pandora2 = NaN(height(data),N_bins);% this bin is to store Pandora NO2 at site 2
for i = 1:(N_bins+1)
    % wind bins start from 0 degree
%     wd_min = (i-1)*(360/N_bins);
%     wd_max = (i)*(360/N_bins);
%     TF = (data.winddirection >= wd_min) & (data.winddirection < wd_max);
    
    % wind bins use 0 degree as centre    
    if (i == 1) | (i == (N_bins + 1))
        TF = (data.winddirection <= 15) | (data.winddirection > 345);
    else
        wd_min = (i-1)*(360/N_bins)-15;
        wd_max = (i)*(360/N_bins)-15;
        TF = (data.winddirection >= wd_min) & (data.winddirection < wd_max);
    end

    data_bins.no2(TF,i) = data.no2(TF)./DU;
%    data_bins.amf(TF,i) = data.amf(TF);
%    data_bins.no2_err(TF,i) = data.no2_err(TF)./DU;
%    data_bins.albedo_no2(TF,i) = data.albedo_no2(TF);
%    data_bins.ground_pixel(TF,i) = data.ground_pixel(TF);
    data_bins.no2_trop(TF,i) = data.no2_trop(TF)./DU;
%    data_bins.amf_trop(TF,i) = data.amf_trop(TF);
    data_bins.no2_strat(TF,i) = data.no2_strat(TF)./DU;
    N(i) = sum(TF);
    data_bins.pandora(TF,i) = data.NO2_VCD(TF);
    if strcmp(method, 'twosites_rotation')
        data_bins.pandora2(TF,i) = data.NO2_VCD_site2(TF);
    end
end

figure;fig_name = 'Mean_NO2_bin_by_wd_TropOMI';
boxplot(data_bins.no2);% this is TropOMI NO2
xticks([1 3 5 7 9 11 13]);
%xticklabels({'15','75','135','195','255','315'});
xticklabels({'0','60','120','180','240','300','360'});
ylabel('TropOMI NO_2 mean VCD [DU]');
xlabel('Wind direction [degree]');
ylim([0 1]);
print_setting(0,save_fig,[plot_path fig_name]);

figure;fig_name = 'Mean_NO2_bin_by_wd_Pandora';
boxplot(data_bins.pandora);    
xticks([1 3 5 7 9 11 13]);
%xticklabels({'15','75','135','195','255','315'});
xticklabels({'0','60','120','180','240','300','360'});
ylabel('Pandora NO_2 mean VCD [DU]');
xlabel('Wind direction [degree]');
ylim([0 1]);
print_setting(0,save_fig,[plot_path fig_name]);

if strcmp(method, 'twosites_rotation')
    figure;fig_name = 'Mean_NO2_bin_by_wd_Pandora_2ndsite';
    boxplot(data_bins.pandora2);    
    xticks([1 3 5 7 9 11 13]);
    %xticklabels({'15','75','135','195','255','315'});
    xticklabels({'0','60','120','180','240','300','360'});
    ylabel('Pandora NO_2 mean VCD [DU]');
    xlabel('Wind direction [degree]');
    ylim([0 1]);
    print_setting(0,save_fig,[plot_path fig_name]);
end

%% major wind bin plot
figure;hold all;fig_name = 'Mean_NO2_bin_by_wd';
y=mean(data_bins.no2,'omitnan');
y_err = std(data_bins.no2,'omitnan')./(N).^0.5;
%errorbar(1:N_bins,y,y_err);
errorbar(1:(N_bins+1),y,y_err);

y=mean(data_bins.no2_trop,'omitnan');
y_err = std(data_bins.no2_trop,'omitnan')./(N).^0.5;
%errorbar(1:N_bins,y,y_err);
errorbar(1:(N_bins+1),y,y_err);

y=mean(data_bins.no2_strat,'omitnan');
y_err = std(data_bins.no2_strat,'omitnan')./(N).^0.5;
%errorbar(1:N_bins,y,y_err);
errorbar(1:(N_bins+1),y,y_err);

y_pandora=mean(data_bins.pandora,'omitnan');
y_pandora_err = std(data_bins.pandora,'omitnan')./(N).^0.5;
%errorbar(1:N_bins,y_pandora,y_pandora_err);
errorbar(1:(N_bins+1),y_pandora,y_pandora_err);

if strcmp(method, 'twosites_rotation')
    y_pandora2=mean(data_bins.pandora2,'omitnan');
    y_pandora2_err = std(data_bins.pandora2,'omitnan')./(N).^0.5;
    %errorbar(1:N_bins,y_pandora2,y_pandora2_err);
    errorbar(1:(N_bins+1),y_pandora2,y_pandora2_err);
end
% legend({'TropOMI',['Pandora ' num2str(Pandora_no)]});
legend({'TropOMI','TropOMI trop.','TropOMI start.',['Pandora ' num2str(Pandora_no)]});
xticks([1 3 5 7 9 11 13]);
%xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
xticklabels({'0^o','60^o','120^o','180^o','240^o','300^o','360^o'});
ylabel('NO_2 mean VCD [DU]');
xlabel('Wind direction');
ylim([0 1.2]);
print_setting(0,save_fig,[plot_path fig_name]);

%% major wind bin plot --> abs diff
figure;hold all;fig_name = 'Mean_NO2_bin_by_wd_abs_diff';
y=mean(data_bins.no2,'omitnan');% tropomi no2
y_err = std(data_bins.no2,'omitnan')./(N).^0.5;

y_pandora=mean(data_bins.pandora,'omitnan');% pandora no2
y_pandora_err = std(data_bins.pandora,'omitnan')./(N).^0.5;

y_diff = y_pandora - y;% absolute difference
y_diff_err = (y_pandora_err.^2 + y_err.^2).^0.5;% combined error

%errorbar(1:N_bins,y_diff,y_diff_err);
errorbar(1:(N_bins+1),y_diff,y_diff_err);

legend({['Pandora ' num2str(Pandora_no) ' - TropOMI']});
xticks([1 3 5 7 9 11 13]);
%xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
xticklabels({'0^o','60^o','120^o','180^o','240^o','300^o','360^o'});
ylabel('Pandora - TropOMI NO_2 VCD [DU]');
xlabel('Wind direction');
ylim([-0.35 0.35]);
print_setting(0,save_fig,[plot_path fig_name]);

%% major wind bin plot --> rel diff
figure;hold all;fig_name = 'Mean_NO2_bin_by_wd_rel_diff';
y=mean(data_bins.no2,'omitnan');% tropomi no2
y_err = std(data_bins.no2,'omitnan')./(N).^0.5;

y_pandora=mean(data_bins.pandora,'omitnan');% pandora no2
y_pandora_err = std(data_bins.pandora,'omitnan')./(N).^0.5;

y_rel_diff = (y_pandora - y)./((y_pandora + y)/2).*100;% absolute difference

%plot(1:N_bins,y_rel_diff,'.-');
plot(1:(N_bins+1),y_rel_diff,'.-');

legend({['Pandora ' num2str(Pandora_no) ' - TropOMI']});
xticks([1 3 5 7 9 11]);
%xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
xticklabels({'0^o','60^o','120^o','180^o','240^o','300^o','360^o'});
ylabel('\delta_r_e_l Pandora - TropOMI NO_2 VCD [%]');
xlabel('Wind direction');
ylim([-100 100]);
print_setting(0,save_fig,[plot_path fig_name]);

%% scatter plot
fig_name = 'Mean_NO2_bin_by_wd_scatter';
%y=mean(data_bins.no2,'omitnan');
y=mean(data_bins.no2(:,1:12),'omitnan');
tf = isnan(y);
y = y(~tf);
y_pandora = y_pandora(:,1:12);
y_pandora = y_pandora(~tf);
line_fits(y_pandora',y');
xlabel(['Pandora'  num2str(Pandora_no) ' mean VCD [DU]']);
ylabel('TropOMI mean VCD [DU]');
new_N_bins = sum(~tf);
scatter(y_pandora',y',10,(1:new_N_bins)','filled');
print_setting(0,save_fig,[plot_path fig_name]);

%% 2nd wind bin plot for amf
figure;hold all;fig_name = 'Mean_NO2_amf_albedo_bin_by_wd';
y=mean(data_bins.amf,'omitnan');
y_err = std(data_bins.amf,'omitnan')./(N).^0.5;

y2=mean(data_bins.amf_trop,'omitnan');
y2_err = std(data_bins.amf_trop,'omitnan')./(N).^0.5;
%plot(1:N_bins,y);
yyaxis left
%errorbar(1:N_bins,y,y_err);
%errorbar(1:N_bins,y2,y2_err);
errorbar(1:(N_bins+1),y,y_err);
errorbar(1:(N_bins+1),y2,y2_err);
ylabel('NO_2 AMF');
ylim([0.5 3.1]);

yyaxis right
y_albedo_no2=mean(data_bins.albedo_no2,'omitnan');
y_albedo_no2_err = std(data_bins.albedo_no2,'omitnan')./(N).^0.5;
%errorbar(1:N_bins,y_albedo_no2,y_albedo_no2_err);
errorbar(1:(N_bins+1),y_albedo_no2,y_albedo_no2_err);
%legend({'TropOMI',['Pandora ' num2str(Pandora_no)]});
legend({'TropOMI AMF','TropOMI trop. AMF','TropOMI albedo NO_2'});
xticks([1 3 5 7 9 11 13]);
%xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
xticklabels({'0^o','60^o','120^o','180^o','240^o','300^o','360^o'});
ylabel('NO_2 albedo');
xlabel('Wind direction');
ylim([0 0.1]);
print_setting(0,save_fig,[plot_path fig_name]);


%% 3rd wind bin plot for pixel number and no2 err
figure;hold all;fig_name = 'pixel_number_and_no2err_bin_by_wd';
y=mean(data_bins.ground_pixel,'omitnan');
y_err = std(data_bins.ground_pixel,'omitnan')./(N).^0.5;
%plot(1:N_bins,y);
yyaxis left
%errorbar(1:N_bins,y,y_err);
errorbar(1:(N_bins+1),y,y_err);
ylabel('mean TropOMI ground pixel number');
%ylim([1 3.1]);
yyaxis right
y2=mean(data_bins.no2_err,'omitnan');
y2_err = std(data_bins.no2_err,'omitnan')./(N).^0.5;
%errorbar(1:N_bins,y2,y2_err);
errorbar(1:(N_bins+1),y2,y2_err);
legend({'TropOMI mean ground pixel','TropOMI NO_2 err'});
xticks([1 3 5 7 9 11 13]);
%xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
xticklabels({'0^o','60^o','120^o','180^o','240^o','300^o','360^o'});
ylabel('NO_2 err [DU]');
xlabel('Wind direction');
%ylim([0 0.1]);
print_setting(0,save_fig,[plot_path fig_name]);

%% 4th wind bin plot for trop. amf and trop. column
figure;hold all;fig_name = 'trop_amf_and_trop_column_by_wd';
y=mean(data_bins.amf_trop,'omitnan');
y_err = std(data_bins.amf_trop,'omitnan')./(N).^0.5;
%plot(1:N_bins,y);
yyaxis left
%errorbar(1:N_bins,y,y_err);
errorbar(1:(N_bins+1),y,y_err);
ylabel('TropOMI NO_2 trop. AMF');
ylim([0 2]);
yyaxis right
y2=mean(data_bins.no2_trop,'omitnan');
y2_err = std(data_bins.no2_trop,'omitnan')./(N).^0.5;
%errorbar(1:N_bins,y2,y2_err);
errorbar(1:(N_bins+1),y2,y2_err);
legend({'TropOMI trop. AMF','TropOMI trop. NO_2'});
xticks([1 3 5 7 9 11 13]);
%xticklabels({'15^o','75^o','135^o','195^o','255^o','315^o'});
xticklabels({'0^o','60^o','120^o','180^o','240^o','300^o','360^o'});
ylabel('NO_2 trop. column [DU]');
xlabel('Wind direction');
ylim([0 0.4]);
print_setting(0,save_fig,[plot_path fig_name]);

%% histogram
function combined_histogram(data,data2,Pandora_no,Pandora2_no,save_fig,plot_path)
DU = 2.6870e+16;
fig_name = 'histogram_of_ground_pixel';
h1 = figure;hold all;
histogram(data.ground_pixel,'BinWidth',15);ylabel('f');xlabel('ground pixel');
histogram(data2.ground_pixel,'BinWidth',15);ylabel('f');xlabel('ground pixel');
legend({num2str(Pandora_no),num2str(Pandora2_no)});
xlim([0 450]);
ylim([0 250]);
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_NO2_err';
h2 = figure;hold all;
histogram(data.no2_err./DU,'BinWidth',0.01);ylabel('f');xlabel('TropOMI NO_2 err [DU]');
histogram(data2.no2_err./DU,'BinWidth',0.01);ylabel('f');xlabel('TropOMI NO_2 err [DU]');
legend({num2str(Pandora_no),num2str(Pandora2_no)});
xlim([0 0.2]);
ylim([0 250]);
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_NO2_AMF';
h3 = figure;hold all;
histogram(data.amf,'BinWidth',0.2);ylabel('f');xlabel('TropOMI NO_2 AMF');
histogram(data2.amf,'BinWidth',0.2);ylabel('f');xlabel('TropOMI NO_2 AMF');
legend({num2str(Pandora_no),num2str(Pandora2_no)});
xlim([0 4.5]);
ylim([0 250]);
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_NO2_AMF_trop';
h4 = figure;hold all;
histogram(data.amf_trop,'BinWidth',0.2);ylabel('f');xlabel('TropOMI NO_2 trop. AMF');
histogram(data2.amf_trop,'BinWidth',0.2);ylabel('f');xlabel('TropOMI NO_2 trop. AMF');
legend({num2str(Pandora_no),num2str(Pandora2_no)});
xlim([0 4.5]);
ylim([0 250]);
print_setting(0,save_fig,[plot_path fig_name]);

fig_name = 'histogram_of_TropOMI_NO2_Albedo';
h5 = figure;hold all;
histogram(data.albedo_no2,'BinWidth',0.01);ylabel('f');xlabel('TropOMI NO_2 albedo');
histogram(data2.albedo_no2,'BinWidth',0.01);ylabel('f');xlabel('TropOMI NO_2 albedo');
legend({num2str(Pandora_no),num2str(Pandora2_no)});
xlim([0 0.1]);
ylim([0 250]);
print_setting(0,save_fig,[plot_path fig_name]);



% fig_name = 'histogram_of_TropOMI_qa';
% h5 = figure;hold all;
% histogram(data.qa,'BinWidth',0.01);ylabel('f');xlabel('TropOMI qa');
% histogram(data2.qa,'BinWidth',0.01);ylabel('f');xlabel('TropOMI qa');
% legend({num2str(Pandora_no),num2str(Pandora2_no)});
% print_setting(0,save_fig,[plot_path fig_name]);



