function daily_weekly_mean()
load('C:\Projects\TropOMI\plots\Wind_transport\v7\103\rotation_knmi\1siteTropOMI_Pandora_transport.mat');
DU = 2.6870e+16;
data = C_all;
data.no2_sum = data.no2_trop + data.no2_strat;
data.no2_sum_err = (data.no2_trop_err.^2 + data.no2_strat_err.^2).^0.5;

data.no2_eccc_sum = data.ECCC_NO2 + data.no2_strat;


data = table2timetable(data,'RowTimes','UTC');

%% QC
TF = abs(data.x1) <=5 & abs(data.y1) <=10 & abs(data.arrival_time) <=1;
data(~TF,:) = [];


%%
data_daily = retime(data,'daily','mean');
data_daily_std = retime(data,'daily',@std);
TF_NAN = isnan(data_daily.no2_sum);
data_daily(TF_NAN,:) = [];
data_daily_std(TF_NAN,:) = [];
%% plots
figure;hold all;
x = datenum(data_daily.UTC);
y = data_daily.no2_sum./DU;
y_err = data_daily_std.no2_sum./DU;
errorbar(x,y,y_err,'-.');
%plot(x,y,'.-');
datetick('x','yyyy-mm-dd','keeplimits');

x = datenum(data_daily.UTC);
y = data_daily.no2_eccc_sum./DU;
y_err = data_daily_std.no2_eccc_sum./DU;
errorbar(x,y,y_err,'.-');
%plot(x,y,'.-');
datetick('x','yyyy-mm-dd','keeplimits');

x = datenum(data_daily.UTC);
y = data_daily.NO2_VCD;
y_err = data_daily_std.NO2_VCD;
errorbar(x,y,y_err,'.-');
%plot(x,y,'.-');
datetick('x','yyyy-mm-dd','keeplimits');

legend({'KNMI','ECCC','Pandora'});


