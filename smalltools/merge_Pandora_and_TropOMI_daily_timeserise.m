function merge_Pandora_and_TropOMI_daily_timeserise()

save_fig = 1;
time_bin_Pandora = 30; % min --> pandora measurements will be averaged
distance = 10; % km --> range of TropOMI pixels will be averaged

output_path = ['C:\Projects\TropOMI\plots\classical\'];

mkdir(output_path);
% site = 'Downsview';
% Pandora_no = '103';
% site = 'FortMcKay';
% Pandora_no = '122';
% site = 'StGeorge';
% Pandora_no = '145';
% site = 'Egbert';
% Pandora_no = '108';

site1 = 'Downsview';
Pandora_no1 = '103';
[site1_h1,site1_h2,site1_h3] = plot_one_site(time_bin_Pandora,distance,site1,Pandora_no1);

site2 = 'StGeorge';
Pandora_no2 = '145';
[site2_h1,site2_h2,site2_h3] = plot_one_site(time_bin_Pandora,distance,site2,Pandora_no2);

site3 = 'Egbert';
Pandora_no3 = '108';
[site3_h1,site3_h2,site3_h3] = plot_one_site(time_bin_Pandora,distance,site3,Pandora_no3);

site4 = 'FortMcKay';
Pandora_no4 = '122';
[site4_h1,site4_h2,site4_h3] = plot_one_site(time_bin_Pandora,distance,site4,Pandora_no4);


fig_name = 'combined_timeserise';
t_start = datenum('2018-03-01');
t_end = datenum('2019-01-01');
h0 = figure;
h0_1 = subplot(4,1,1);
combin2subplot(site1_h1,h0_1);xlim([t_start t_end]);ylim([0 1.4]);
datetick('x','mm-dd','keeplimits');title(site1);
h0_2 = subplot(4,1,2);
combin2subplot(site2_h1,h0_2);xlim([t_start t_end]);ylim([0 1.4]);
datetick('x','mm-dd','keeplimits');title(site2);
h0_3 = subplot(4,1,3);
combin2subplot(site3_h1,h0_3);xlim([t_start t_end]);ylim([0 1.4]);
datetick('x','mm-dd','keeplimits');title(site3);
h0_4 = subplot(4,1,4);
combin2subplot(site4_h1,h0_4);xlim([t_start t_end]);ylim([0 1.4]);
datetick('x','mm-dd','keeplimits');title(site4);

print_setting(1/2, 0, [output_path fig_name]);
datetick('x','mm-dd','keeplimits');
print_setting(1/2, save_fig, [output_path fig_name]);


fig_name = 'combined_KNMI_vs_Pandora';
h0 = figure;title('KNMI vs. Pandora');
h0_1 = subplot(2,2,1);
combin2subplot(site1_h2,h0_1);xlim([0 1.4]);ylim([0 1.4]);title(site1);
xlabel('TROPOMI (KNMI) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
h0_2 = subplot(2,2,2);
combin2subplot(site2_h2,h0_2);xlim([0 1.4]);ylim([0 1.4]);title(site2);
xlabel('TROPOMI (KNMI) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
h0_3 = subplot(2,2,3);
combin2subplot(site3_h2,h0_3);xlim([0 1.4]);ylim([0 1.4]);title(site3);
xlabel('TROPOMI (KNMI) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
h0_4 = subplot(2,2,4);
combin2subplot(site4_h2,h0_4);xlim([0 1.4]);ylim([0 1.4]);title(site4);
xlabel('TROPOMI (KNMI) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
print_setting(1, save_fig, [output_path fig_name]);

fig_name = 'combined_ECCC_vs_Pandora';
h0 = figure;title('ECCC vs. Pandora');
h0_1 = subplot(2,2,1);
combin2subplot(site1_h3,h0_1);xlim([0 1.4]);ylim([0 1.4]);title(site1);
xlabel('TROPOMI (ECCC) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
h0_2 = subplot(2,2,2);
combin2subplot(site2_h3,h0_2);xlim([0 1.4]);ylim([0 1.4]);title(site2);
xlabel('TROPOMI (ECCC) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
h0_3 = subplot(2,2,3);
combin2subplot(site3_h3,h0_3);xlim([0 1.4]);ylim([0 1.4]);title(site3);
xlabel('TROPOMI (ECCC) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
h0_4 = subplot(2,2,4);
combin2subplot(site4_h3,h0_4);xlim([0 1.4]);ylim([0 1.4]);title(site4);
xlabel('TROPOMI (ECCC) NO_2 [DU]');ylabel('Pandora NO_2 [DU]');
print_setting(1, save_fig, [output_path fig_name]);

function combin2subplot(original_ax,target_ax)
axes_to_be_copied = findobj(original_ax,'type','axes', '-not', 'tag', 'legend');
% Identify the Legend
%legend_to_be_copied = findobj(original_ax,'type','axes', 'tag', 'legend');
legend_to_be_copied = findobj(original_ax,'type','legend', 'tag', 'legend');
% Identify the children of this axes 
chilred_to_be_copied = get(axes_to_be_copied,'children'); 
% Identify orientation of the axes 
[az,el] = view; 
% Copy the children of the axes 
copyobj(chilred_to_be_copied,target_ax);
% If there is a legend
if isfloat(legend_to_be_copied)
    copyobj(legend_to_be_copied, target_ax);
end
% Set the limits and orientation of the subplot as the original figure 
set(target_ax,'Xlim',get(axes_to_be_copied,'XLim')) 
set(target_ax,'Ylim',get(axes_to_be_copied,'YLim')) 
set(target_ax,'Zlim',get(axes_to_be_copied,'ZLim')) 
view(target_ax,[az,el]) 

function [h1,h2,h3] = plot_one_site(time_bin_Pandora,distance,site,Pandora_no,save_fig)

if nargin <6
    save_fig =1;
end

remove_snow_cover = false;


output_path = ['C:\Projects\TropOMI\plots\classical\' Pandora_no '\'];

mkdir(output_path);

Pandora_data = ['C:\Projects\Pandora\output\' Pandora_no '\Pandora' Pandora_no 's1_' site '_L2Tot_rnvs0p1-5.mat'];

  
TropOMI_data = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport\TropOMI_transport.mat'];

addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab');

load(Pandora_data);
load(TropOMI_data);
DU = 2.6870e+16;
    

%% QC for TropOMI
TF_qa = TropOMI.qa >=1;
TropOMI(~TF_qa,:) = [];

if remove_snow_cover
    TF_snow = TropOMI.snow_cover ==4;
    TropOMI(~TF_snow,:) = [];
end

%% test generate time stamp here by use "total_time"
% pivot_time = datetime('2010-01-01');
% TropOMI.utc_time_new = datetime(datestr(TropOMI.total_time./60./60./24 + pivot_time));
% 
% figure;
% plot(TropOMI.utc_time,TropOMI.utc_time_new,'.');
% xlabel('utc using day + delta time');
% ylabel('utc using total time');

%% select TropOMI

TF_distance = TropOMI.d <= distance.*1000;% note the d in TropOMI is in metre!
TropOMI(~TF_distance,:) = [];

TropOMI.ECCC_total = TropOMI.ECCC_NO2 + TropOMI.no2_strat;
TropOMI.KNMI_total = TropOMI.no2_trop + TropOMI.no2_strat;

TropOMI = table2timetable(TropOMI,'RowTimes','utc_time');

%%
% get TropOMI daily mean and std
TropOMI_mean = retime(TropOMI,'hourly','mean');
TropOMI_std= retime(TropOMI,'hourly',@std);

TF_nan = isnan(TropOMI_mean.lon);
TropOMI_mean(TF_nan,:) = [];
TropOMI_std(TF_nan,:) = [];

%% select Pandora measurements based on hourly averaged TropOMI time stamp
t_offset = minutes(time_bin_Pandora);
TropOMI_mean_time = TropOMI_mean.total_time./60./60./24 + datetime('2010-01-01');
TropOMI_mean.pairing_start_time = TropOMI_mean_time - t_offset;
TropOMI_mean.pairing_end_time = TropOMI_mean_time + t_offset;


Pandora = table2timetable(Pandora,'RowTimes','UTC');
for i = 1:height(TropOMI_mean)
    TF_paired= Pandora.UTC >= TropOMI_mean.pairing_start_time(i,:) & Pandora.UTC <= TropOMI_mean.pairing_end_time(i,:);
    Pandora_paired_1hour = Pandora(TF_paired,:);
    if i == 1
        Pandora_paired_all = Pandora_paired_1hour;
    else
        Pandora_paired_all = [Pandora_paired_all;Pandora_paired_1hour];
    end
    Pandora_paired_1hour = [];
end
Pandora_mean = retime(Pandora_paired_all,'hourly','mean');
Pandora_std = retime(Pandora_paired_all,'hourly',@std);

TF_nan = isnan(Pandora_mean.NO2_VCD);
Pandora_mean(TF_nan,:) = [];
Pandora_std(TF_nan,:) = [];

%% merge Pandora and TropOMI
TropOMI_mean.join_timestamp = datetime(TropOMI_mean.utc_time.Year,TropOMI_mean.utc_time.Month,TropOMI_mean.utc_time.Day,TropOMI_mean.utc_time.Hour,0,0);
TropOMI_std.join_timestamp = datetime(TropOMI_mean.utc_time.Year,TropOMI_mean.utc_time.Month,TropOMI_mean.utc_time.Day,TropOMI_mean.utc_time.Hour,0,0);

Pandora_mean.join_timestamp = datetime(Pandora_mean.UTC.Year,Pandora_mean.UTC.Month,Pandora_mean.UTC.Day,Pandora_mean.UTC.Hour,0,0);
Pandora_std.join_timestamp = datetime(Pandora_mean.UTC.Year,Pandora_mean.UTC.Month,Pandora_mean.UTC.Day,Pandora_mean.UTC.Hour,0,0);

mean_merged = innerjoin(TropOMI_mean,Pandora_mean,'LeftKeys','join_timestamp','RightKeys','join_timestamp');
std_merged = innerjoin(TropOMI_std,Pandora_std,'LeftKeys','join_timestamp','RightKeys','join_timestamp');





%%
fig_name = ['timeseries_' num2str(time_bin_Pandora) '_min_' num2str(distance) '_km'];
h1 = figure;hold all;
x = mean_merged.utc_time;
y = mean_merged.KNMI_total./DU;
y_err =  std_merged.KNMI_total./DU;
% plot(x,y,'.-');
errorbar(datenum(x),y,y_err);

x = mean_merged.utc_time;
y = mean_merged.ECCC_total./DU;
y_err = std_merged.ECCC_total./DU;
% plot(x,y,'.-');
errorbar(datenum(x),y,y_err);

x = mean_merged.utc_time;
y = mean_merged.NO2_VCD;
y_err = std_merged.NO2_VCD;
% plot(x,y,'.-');
errorbar(datenum(x),y,y_err);
ylabel('NO_2 VCD [DU]');
datetick('x','yyyy-mm-dd','keeplimits');
legend({'KNMI','ECCC',['Pandora' num2str(Pandora_no)]});
title(['+/-' num2str(time_bin_Pandora) ' min ; ' num2str(distance) ' km']);
print_setting(1/2, save_fig, [output_path fig_name]);

%%
fig_name = ['scatter_KNMI_' num2str(time_bin_Pandora) '_min_' num2str(distance) '_km'];
x = mean_merged.NO2_VCD;% Pandora
x_err = std_merged.NO2_VCD;
y = mean_merged.KNMI_total./DU;% KNMI TropOMI
y_err =  std_merged.KNMI_total./DU;

line_fits(x,y);
h2 = gca;
xlabel('Pandora NO_2 VCD [DU]');
ylabel('KNMI TROPOMI NO_2 VCD [DU]');


title(['+/-' num2str(time_bin_Pandora) ' min ; ' num2str(distance) ' km']);
print_setting(1/4, save_fig, [output_path fig_name]);

%%
fig_name = ['scatter_ECCC_' num2str(time_bin_Pandora) '_min_' num2str(distance) '_km'];
x = mean_merged.NO2_VCD;% Pandora
x_err = std_merged.NO2_VCD;
y = mean_merged.ECCC_total./DU;% ECCC TropOMI
y_err =  std_merged.ECCC_total./DU;

line_fits(x,y);
h3 = gca;
xlabel('Pandora NO_2 VCD [DU]');
ylabel('ECCC TROPOMI NO_2 VCD [DU]');

% tf = y_err ==0;
% y_err(tf,:) = 0.0001;
%[a, b, sigma_a, sigma_b, b_save] = york_fit(x',y',x_err',y_err');
%y_fits = a+b.*x;
%plot(x,y_fits,'-');
title(['+/-' num2str(time_bin_Pandora) ' min ; ' num2str(distance) ' km']);
print_setting(1/4, save_fig, [output_path fig_name]);
%% save the data


