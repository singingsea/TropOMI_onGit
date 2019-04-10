function [status,bins_all] = bin_by_x(method,data,x_name,y_name,type)
step = 0.1; % step in x
i_bins = 0:step:1.5;% steps in x

DU = 2.6870e+16;
switch type
    case 'original' % 'original', 'up wind', 'down wind'       
        TF = data.type == 0;
    case 'up wind'       
        TF = data.type == 1;
    case 'down wind'      
        TF = data.type == -1;
    case 'all'
        TF = ~isnan(data.type);
end
data = data(TF,:);
switch x_name
    case 'Pandora'        
        x = data.NO2_VCD;% Pandora --> 1st site
        x_label = 'Pandora VCD [DU]';
    case 'Pandora2'
        x = data.NO2_VCD_site2;% Pandora --> 2nd site
        x_label = 'Pandora VCD (site2) [DU]';
    case 'TropOMI'
        x = data.no2./DU;% TropOMI
        x_label = 'TropOMI VCD [DU]';
end

switch y_name
    case 'Pandora'
        y = data.NO2_VCD;% Pandora --> 1st site [DU]
        y_label = 'Pandora VCD [DU]';
    case 'Pandora2'
        y = data.NO2_VCD_site2;% Pandora --> 2nd site [DU]
        y_label = 'Pandora VCD (site2) [DU]';
    case 'TropOMI'
        y = data.no2./DU;% TropOMI [DU]
        y_label = 'TropOMI VCD [DU]';
end


bins_all = NaN(height(data),length(i_bins)-1);
for i = 1:length(i_bins)-1
    TF = (x >= i_bins(i)) & (x < i_bins(i+1));
    data_in1bin = y(TF,:); 
    N_in1bin(i,1) = sum(TF);
    if sum(TF)
        bins_all(1:length(data_in1bin),i) = data_in1bin;
        
    end
    data_in1bin = [];
end

status = table;
status.x_centre = (i_bins(1:end-1)+step./2)';
status.y_mean = mean(bins_all,1,'omitnan')';
status.N =  N_in1bin;
status.y_std = std(bins_all,1,'omitnan')';
status.y_std_err = status.y_std./(status.N).^0.5;
status.y_median = median(bins_all,1,'omitnan')';

fig_name = ['scatter_' method '_' y_name '_vs_' x_name '_' type(~isspace(type))];
x = status.x_centre;
y = status.y_mean;
x_err = status.y_std_err;
line_fits(x,y);   
errorbar(x,y,x_err);
%scatter(x,y,10,c,'filled');
ylabel(y_label);
xlabel(x_label);
title(type);

legend({[y_name ' vs. ' x_name],'linear fit','linear fit (no offset)','1-on-1'},'err of mean');
xlim([0 1.4]);
ylim([0 1.4]);

print_setting(1/4,save_fig,[output_path fig_name]);