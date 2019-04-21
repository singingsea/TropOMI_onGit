function diff_by_radius()
% site = 'Downsview';
% Pandora_no = '104';
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

%fig_name = 'mean_difference';
fig_name = 'status_difference';
subplot(2,2,1);
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
ylim([-0.2 0.2]);



subplot(2,2,2);
title(site);
hold all;
x = diff_table_KNMI.bin_right;
y = diff_table_KNMI.R;
yL = diff_table_KNMI.RL;
yU = diff_table_KNMI.RU;
errorbar(x,y,yL,yU,'rx-');
% plot(x,y,'rx-');
x = diff_table_ECCC.bin_right;
y = diff_table_ECCC.R;
yL = diff_table_ECCC.RL;
yU = diff_table_ECCC.RU;
errorbar(x,y,yL,yU,'bo-');
% plot(x,y,'bo-');
legend({'KNMI','ECCC'});
ylabel('R');
xticks(diff_table_KNMI.bin_right);
xticklabels({'0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100'});
xlabel('Radius [km]');
ylim([-0.5 1]);


subplot(2,2,3);
title(site);
hold all;
x = diff_table_KNMI.bin_right;
y = diff_table_KNMI.slop;
y_err = diff_table_KNMI.slop_err;
boundedline(x,y, y_err, 'transparency', 0.5,'rx-','alpha');
x = diff_table_ECCC.bin_right;
y = diff_table_ECCC.slop;
y_err = diff_table_ECCC.slop_err;
boundedline(x,y, y_err, 'transparency', 0.5,'bo-','alpha');
legend({'err of slop (KNMI)','slop (KNMI)','err or slop (ECCC)','slop (ECCC)'});
ylabel('slop');
xticks(diff_table_KNMI.bin_right);
xticklabels({'0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100'});
xlabel('Radius [km]');
ylim([0.5 1.5]);

subplot(2,2,4);
title(site);
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

print_setting(1,save_fig,[plot_path fig_name],'png');

function diff_table = get_diff(data,use_ECCC_no2)
DU = 2.6870e+16;

if use_ECCC_no2
    data.no2 = (data.ECCC_NO2 + data.no2_strat)./DU;% satellite no2, total column by sum start and trop
else
    data.no2 = (data.no2_trop + data.no2_strat)./DU;
end
%%
TF = (abs(data.x1) <= 5) & (abs(data.arrival_time) <=2);
data(~TF,:) = [];

%%
diff_table = table;
%distance_bins = 0:10:100;
distance_bins = 0:10:50;
for i = 1:(numel(distance_bins)-1)
    bin_left = distance_bins(i);
    bin_right = distance_bins(i+1);
    TF = (data.distence2stn_original >= bin_left) & (data.distence2stn_original < bin_right);
    data_1bin = data(TF,:);% data within one distance bin
    
    y = data_1bin.no2; % satellite total NO2
    x = data_1bin.NO2_VCD; % Pandora total NO2
    diff = y - x;% data of NO2 difference, within one distance bin

    [intercept,slop,slop_nlm,mdl_lm,mdl_nlm] = line_fits(x,y);
    slop = mdl_nlm.Coefficients.Estimate;
    slop_err = mdl_nlm.Coefficients.SE;
    %R = mdl_lm.Rsquared.Ordinary^0.5;
    [R,P,RL,RU] = corrcoef(x,y);
    R = R(1,2); RL = R - RL(1,2); RU = RU(1,2) - R;
    
    % create output table, which store the values, such as mean, std and
    % etc
    diff_table.mean(i,1) = mean(diff);
    diff_table.std(i,1) = std(diff);
    diff_table.mean_err(i,1) = std(diff)./(sum(TF))^0.5;
    diff_table.N(i,1) = sum(TF);
    diff_table.bin_left(i,1) = bin_left;
    diff_table.bin_right(i,1) = bin_right;
    
    diff_table.R(i,1) = R;R = [];
    diff_table.RL(i,1) = RL;RL = [];
    diff_table.RU(i,1) = RU;RU = [];
    diff_table.slop(i,1) = slop;slop = [];
    diff_table.slop_err(i,1) = slop_err;slop_err = [];
    
    TF = [];data_1bin=[];
end



    
    

