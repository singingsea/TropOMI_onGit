function diff_by_radius()
% site = 'Downsview';
% Pandora_no = '103';
% site = 'FortMcKay';
% Pandora_no = '122';
% site = 'StGeorge';
% Pandora_no = '145';
site = 'Egbert';
Pandora_no = '108';
save_fig = 1;
load(['C:\Projects\TropOMI\plots\Wind_transport\v9\' num2str(Pandora_no) '\rotation_knmi\1siteTropOMI_Pandora_transport.mat']);
plot_path = ['C:\Projects\TropOMI\plots\Wind_transport\v9\' num2str(Pandora_no) '\rotation_analysis\'];
%plot_path = ['C:\Projects\TropOMI\plots\Wind_transport\v7\' num2str(Pandora_no) '\rotation_analysis_7hr\'];
mkdir(plot_path);
data = C_all;

use_ECCC_no2 = true;
diff_table_ECCC = get_diff(data,use_ECCC_no2);

use_ECCC_no2 = false;
diff_table_KNMI = get_diff(data,use_ECCC_no2);

figure;cmaps = jet(2);

fig_name = 'mean_difference';
subplot(2,1,1);
title(site);
hold all;
x = diff_table_KNMI.bin_right;
y = diff_table_KNMI.mean;
y_err = diff_table_KNMI.mean_err;
%y_err = diff_table_KNMI.std;
%errorbar(x,y,y_err);
boundedline(x,y, y_err, 'transparency', 0.5,'rx-','alpha');
x = diff_table_ECCC.bin_right;
y = diff_table_ECCC.mean;
y_err = diff_table_ECCC.mean_err;
%y_err = diff_table_ECCC.std;
%errorbar(x,y,y_err,'o-');
boundedline(x,y, y_err, 'transparency', 0.5,'bo-','alpha');
%legend({'KNMI','err of mean','ECCC','err of mean'});
legend({'KNMI - Pandora (err of mean)','KNMI - Pandora (mean)','ECCC - Pandora (err of mean)','ECCC - Pandora (mean)'});
ylabel('NO_2 VCD (TROPOMI - Pandora) [DU]');
xticks(diff_table_KNMI.bin_right);
xticklabels({'0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100'});
xlabel('Radius [km]');
ylim([-0.3 0.2]);


subplot(2,1,2);
hold all;
x = diff_table_KNMI.bin_right;
y = diff_table_KNMI.N;
plot(x,y,'rx-');
x = diff_table_ECCC.bin_right;
y = diff_table_ECCC.N;
plot(x,y,'bo-');
legend({'KNMI','ECCC'});
ylabel('No. of measurements');
xticks(diff_table_KNMI.bin_right);
xticklabels({'0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100'});
xlabel('Radius [km]');
ylim([0 500]);

print_setting(1,save_fig,[plot_path fig_name]);

function diff_table = get_diff(data,use_ECCC_no2)
DU = 2.6870e+16;

if use_ECCC_no2
    data.no2 = (data.ECCC_NO2 + data.no2_strat)./DU;
else
    data.no2 = (data.no2_trop + data.no2_strat)./DU;
end
%%
TF = (abs(data.x1) <= 5) & (abs(data.arrival_time) <=1);
data(~TF,:) = [];

%%
diff_table = table;
%distance_bins = 0:10:100;
distance_bins = 0:10:50;
for i = 1:(numel(distance_bins)-1)
    bin_left = distance_bins(i);
    bin_right = distance_bins(i+1);
    TF = (data.distence2stn_original >= bin_left) & (data.distence2stn_original < bin_right);
    data_1bin = data(TF,:);
    diff = data_1bin.no2 - data_1bin.NO2_VCD;
    diff_table.mean(i,1) = mean(diff);
    diff_table.std(i,1) = std(diff);
    diff_table.mean_err(i,1) = std(diff)./(sum(TF))^0.5;
    diff_table.N(i,1) = sum(TF);
    diff_table.bin_left(i,1) = bin_left;
    diff_table.bin_right(i,1) = bin_right;
    TF = [];data_1bin=[];
end



    
    

