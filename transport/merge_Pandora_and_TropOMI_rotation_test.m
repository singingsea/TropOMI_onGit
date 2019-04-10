function merge_Pandora_and_TropOMI_rotation_test()
% this function is based on the wind-driven method! 
% need define the transportation time [hr], which should be consistent with
% the transportation time used in the "get_TropOMI_pixel_transport_location()"
% site = 'Downsview';
% Pandora_no = '104';
site = 'FortMcKay';
Pandora_no = '122';
% site = 'StGeorge';
% Pandora_no = '145';
% site = 'Egbert';
% Pandora_no = '108';
method = 'wind_driven';%'rotation','twosites_rotation','wind_driven'
bin_by = 'travel_t';% ws, travel_t, travel_d, wd
remove_snow_cover = false;

site2 = 'StGeorge';% only used if method is 'twosites_rotation'
Pandora_no2 = '145';
% site2 = 'Downsview';
% Pandora_no2 = '103';
% site2 = 'FortMcKay';
% Pandora_no2 = '122';

delta_t = 1;% transportation time in [hour] --> make sure this is the same number use in "get_TropOMI_pixel_transport_location()"
use_sum_no2 = true;% if true, then TropOMI total column is the Nsum_v (Ntrop + Nstrat); if false, TropOMI total column is the N_v --> for details, see the TropOMI NO2 data description file
use_ECCC_no2 = false;
save_fig = 1;
if use_ECCC_no2
    output_path = ['C:\Projects\TropOMI\plots\Wind_transport\' Pandora_no '\' method '_eccc\'];
else
    output_path = ['C:\Projects\TropOMI\plots\Wind_transport\' Pandora_no '\' method '_knmi\'];
end
%output_path = ['C:\Projects\TropOMI\plots\Wind_transport_0dot5\' Pandora_no '\'];
%output_path = ['C:\Projects\TropOMI\plots\Wind_transport_1hr_avgwind\' Pandora_no '\'];
mkdir(output_path);
%Pandora_data = 'C:\Projects\Pandora\output\103\Pandora103s1_Downsview_L2Tot_rnvs0p1-5.mat';
%Pandora_data = 'C:\Projects\Pandora\output\104\Pandora104s1_Downsview_L2Tot_rnvs0p1-5.mat';
Pandora_data = ['C:\Projects\Pandora\output\' Pandora_no '\Pandora' Pandora_no 's1_' site '_L2Tot_rnvs0p1-5.mat'];
if strcmp(method,'twosites_rotation')
    Pandora_data2 = ['C:\Projects\Pandora\output\' Pandora_no2 '\Pandora' Pandora_no2 's1_' site2 '_L2Tot_rnvs0p1-5.mat'];
end


TropOMI_data = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport\TropOMI_transport.mat'];
%TropOMI_data = 'C:\Projects\TropOMI\data\NO2_output\OFFL\Downsview_ERA\transport\TropOMI_transport.mat';
%TropOMI_data = ['C:\Projects\TropOMI\data\NO2_output\OFFL\FortMcKay_suncrude_ERA\transport\TropOMI_transport.mat'];

%TropOMI_data = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport_0dot5\TropOMI_transport.mat'];
%TropOMI_data = ['C:\Projects\TropOMI\data\NO2_output\OFFL\' site '_ERA\transport_1hr_avgwind\TropOMI_transport.mat'];


addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab');
if strcmp(method,'twosites_rotation')
    load(Pandora_data2);
    Pandora2 = Pandora;
end
load(Pandora_data);
load(TropOMI_data);
DU = 2.6870e+16;
    
%% if use ECCC trop no2, the ECCC no2 is only trop!!!
if use_ECCC_no2 && use_sum_no2
    TropOMI.no2_trop = TropOMI.ECCC_NO2;
end
%% QC for TropOMI
if remove_snow_cover
    TF_snow = TropOMI.snow_cover ==4;
    TropOMI(TF_snow,:) = [];
end

%%
switch method
   
    case 'twosites_rotation'
        % rotation method for 2 sites
        C_all = pair_Pandora_TropOMI_rotation_2sites(Pandora,Pandora2,TropOMI);

        data = C_all;
        
        %TF_wd = (data.winddirection >= 165) & (data.winddirection <= 225);
%         TF_wd = (data.winddirection >= 127) & (data.winddirection <= 187);
%         data(~TF_wd,:) = [];
        %x_name = 'Pandora'; % 'Pandora','Pandora2','TropOMI'
        %y_name = 'TropOMI'; % 'Pandora','Pandora2','TropOMI'
        %c_name = 't1'; % 't1', 't2', 'd1', 'd2'
        %type = 'original';% 'original', 'down wind', 'up wind', 'all'
        method = '2sites';
        auto_scatter(method,data,'Pandora','TropOMI','wd','all',save_fig,output_path,use_sum_no2);           
        
        %auto_scatter(method,data,'Pandora','TropOMI','wd','original',save_fig,output_path);           
        %auto_scatter(method,data,'Pandora2','TropOMI','wd','original',save_fig,output_path);  

        %auto_scatter(method,data,'Pandora','Pandora2','wd','original',save_fig,output_path);
        auto_scatter(method,data,'Pandora','Pandora2','wd','up wind',save_fig,output_path,use_sum_no2);
        auto_scatter(method,data,'Pandora','Pandora2','wd','down wind',save_fig,output_path,use_sum_no2);

        auto_scatter(method,data,'Pandora','TropOMI','wd','up wind',save_fig,output_path,use_sum_no2);
        auto_scatter(method,data,'Pandora','TropOMI','wd','down wind',save_fig,output_path,use_sum_no2);
        auto_scatter(method,data,'Pandora2','TropOMI','wd','up wind',save_fig,output_path,use_sum_no2);
        auto_scatter(method,data,'Pandora2','TropOMI','wd','down wind',save_fig,output_path,use_sum_no2);

    
    case 'rotation'
        % rotation method for 1 site
%         distance = 7;
%         y_dis = 15;
         distance = 5;
         y_dis = 5;    
%          distance = 20;% smooth to OMI resolution
%          y_dis = 10;    

%          arrival_time_limit = 7;% hour
        arrival_time_limit = 1;% hour
        data = pair_Pandora_TropOMI_rotation(Pandora,TropOMI,distance,y_dis,arrival_time_limit);
    
        method = '1site';

        group_status = table;
        [group_status.intercept(1),group_status.slop(1),group_status.slop_nlm(1),group_status.R(1),group_status.N(1)]  ...
            = auto_scatter(method,data,'Pandora','TropOMI','d1','original',save_fig,output_path,use_sum_no2,['-' num2str(distance) 'km' num2str(y_dis) 'km']);
        descirption = [num2str(distance) 'kmX' num2str(y_dis) 'km'];
        group_status.description(1) = {descirption};
        bin_by_x(method,data,'Pandora','TropOMI','original',save_fig,output_path,use_sum_no2,descirption);       
        
        
        data = pair_Pandora_TropOMI_rotation(Pandora,TropOMI,distance,y_dis*2,arrival_time_limit);
        
        [group_status.intercept(2),group_status.slop(2),group_status.slop_nlm(2),group_status.R(2),group_status.N(2)]  ...
            = auto_scatter(method,data,'Pandora','TropOMI','d1','original',save_fig,output_path,use_sum_no2,['-' num2str(distance) 'km' num2str(y_dis*2) 'km']);     
        descirption = [num2str(distance) 'kmX' num2str(y_dis*2) 'km'];
        group_status.description(2) = {descirption};
        bin_by_x(method,data,'Pandora','TropOMI','original',save_fig,output_path,use_sum_no2,descirption);
        

        data = pair_Pandora_TropOMI_rotation(Pandora,TropOMI,distance,y_dis*3,arrival_time_limit);
        
        [group_status.intercept(3),group_status.slop(3),group_status.slop_nlm(3),group_status.R(3),group_status.N(3)]  ...
            = auto_scatter(method,data,'Pandora','TropOMI','d1','original',save_fig,output_path,use_sum_no2,['-' num2str(distance) 'km' num2str(y_dis*3) 'km']);   
        descirption = [num2str(distance) 'kmX' num2str(y_dis*3) 'km'];
        group_status.description(3) = {descirption};        
        %group_status.description(3) = {'7km*21km'};
        bin_by_x(method,data,'Pandora','TropOMI','original',save_fig,output_path,use_sum_no2,descirption);
        

        data = pair_Pandora_TropOMI_rotation(Pandora,TropOMI,distance,y_dis*4,arrival_time_limit);
        
        [group_status.intercept(4),group_status.slop(4),group_status.slop_nlm(4),group_status.R(4),group_status.N(4)]  ...
            = auto_scatter(method,data,'Pandora','TropOMI','d1','original',save_fig,output_path,use_sum_no2,['-' num2str(distance) 'km' num2str(y_dis*4) 'km']);   
        descirption = [num2str(distance) 'kmX' num2str(y_dis*4) 'km'];
        group_status.description(4) = {descirption};        
        %group_status.description(4) = {'7km*28km'};
        bin_by_x(method,data,'Pandora','TropOMI','original',save_fig,output_path,use_sum_no2,descirption);


        data = pair_Pandora_TropOMI_rotation(Pandora,TropOMI,distance,y_dis*5,arrival_time_limit);
        
        [group_status.intercept(5),group_status.slop(5),group_status.slop_nlm(5),group_status.R(5),group_status.N(5)]  ...
            = auto_scatter(method,data,'Pandora','TropOMI','d1','original',save_fig,output_path,use_sum_no2,['-' num2str(distance) 'km' num2str(y_dis*5) 'km']);   
        descirption = [num2str(distance) 'kmX' num2str(y_dis*5) 'km'];
        group_status.description(5) = {descirption};        
        %group_status.description(5) = {'7km*35km'};
        bin_by_x(method,data,'Pandora','TropOMI','original',save_fig,output_path,use_sum_no2,descirption);


        data = pair_Pandora_TropOMI_rotation(Pandora,TropOMI,distance,y_dis*6,arrival_time_limit);
        
        [group_status.intercept(6),group_status.slop(6),group_status.slop_nlm(6),group_status.R(6),group_status.N(6)]  ...
            = auto_scatter(method,data,'Pandora','TropOMI','d1','original',save_fig,output_path,use_sum_no2,['-' num2str(distance) 'km' num2str(y_dis*6) 'km']);   
        descirption = [num2str(distance) 'kmX' num2str(y_dis*6) 'km'];
        group_status.description(6) = {descirption};        
        %group_status.description(6) = {'7km*42km'};
        bin_by_x(method,data,'Pandora','TropOMI','original',save_fig,output_path,use_sum_no2,descirption);

 
        fig_name = ['group_status_' num2str(distance) 'km' num2str(y_dis) 'km'];
        figure;hold all;
        ax1 = subplot(3,1,1);        
        plot(group_status.slop_nlm,'.-');
        ylim([0.3 1.2]);
        xticks(ax1,1:6);
        ylabel('slop');
        ax2 = subplot(3,1,2);
        plot(group_status.R,'.-');
        ylim([0.3 1]);
        xticks(ax2,1:6);
        ylabel('R');
        ax3 = subplot(3,1,3);
        plot(group_status.N,'.-'); 
        ylim([0 1200]);
        ylabel('N');
        xticks(ax3,1:6);
        xticklabels(ax3,group_status.description);
        rotateXLabels(ax3,45);       
        print_setting(1,save_fig,[output_path fig_name]);
        
        C_all = data;% this is the table will be auto saved
        
    case 'wind_driven'    
        % wind-dirven method
        up_down_original = 'original';% 'original', 'up wind', 'down wind'
        fig_name = ['scatter_' up_down_original];
        C_closest_orginal = pair_Pandora_TropOMI(Pandora,TropOMI,delta_t,up_down_original,use_sum_no2);
        print_setting(1/4,save_fig,[output_path fig_name]);

        up_down_original = 'up wind';% 'original', 'up wind', 'down wind'
        fig_name = ['scatter_' up_down_original];
        C_closest_up_wind = pair_Pandora_TropOMI(Pandora,TropOMI,delta_t,up_down_original,use_sum_no2);
        print_setting(1/4,save_fig,[output_path fig_name]);

        up_down_original = 'down wind';% 'original', 'up wind', 'down wind'
        fig_name = ['scatter_' up_down_original];
        C_closest_down_wind = pair_Pandora_TropOMI(Pandora,TropOMI,delta_t,up_down_original,use_sum_no2);
        print_setting(1/4,save_fig,[output_path fig_name]);

        %
        C_all = [C_closest_orginal;C_closest_up_wind;C_closest_down_wind];
        fig_name = ['scatter_all'];
        if use_sum_no2
            y = (C_all.no2_trop + C_all.no2_strat)./DU;% TropOMI sum no2
        else
            y = C_all.no2./DU;% TropOMI total no2
        end
        x = C_all.NO2_VCD;% Pandora
        line_fits(x,y);
        %scatter(x,y,10,C_all.distance_4pair,'filled');
        %scatter(x,y,10,C_all.type,'filled');
        gscatter(x,y,C_all.type);
        ylabel('TropOMI NO_2 VCD [DU]');
        xlabel('Pandora NO_2 VCD [DU]');
        legend({'TropOMI vs. Pandora','linear fit','linear fit (no offset)','1-on-1','up wind pixel','original pixel','down wind pixel'});
        xlim([0 1.5]);
        ylim([0 1.5]);
        print_setting(1/4,save_fig,[output_path fig_name]);
end


%% save the data
save([output_path '\' method 'TropOMI_Pandora_transport.mat'],'C_all');

%%
function C_rotated_match_2sites = pair_Pandora_TropOMI_rotation_2sites(Pandora,Pandora2,TropOMI)
% simple QC
TF_qa = TropOMI.qa == 1;
TropOMI(~TF_qa,:) = [];

% filter by distance
distance = 15; % km
TF1 = abs(TropOMI.x1_2site) <= distance ;% pixel is within 10 km
%TF2 = abs(TropOMI.x1_st2_2site) <= distance ;% rotated station 2 is within 10 km
%TF = TF1 & TF2;
TF = TF1;
TropOMI = TropOMI(TF,:);

% calculate arrival time at site1
TropOMI.arrival_time = TropOMI.y1_2site./TropOMI.windspeed; %[hr] --> y1 is in [km], windspeed is in [km/hr]
% note --> negative y1 is down wind, so, negative arrival_time is also for
% downwind pixel

% calculate arrival time at site2
%TropOMI.arrival_time2 = (TropOMI.y1_2site - TropOMI.sites_distance)./TropOMI.windspeed; %[hr] --> y1 is in [km], windspeed is in [km/hr]
TropOMI.arrival_time2 = (TropOMI.y1_2site - TropOMI.y1_st2_2site)./TropOMI.windspeed; %[hr] --> y1 is in [km], windspeed is in [km/hr]

% QC based on transportation time, note that the wind field is hourly, so,
% it is only "accurate" if the transportation time is within 1hr
% TF_arrival = (abs(TropOMI.arrival_time) > 5) | (abs(TropOMI.arrival_time2) > 5);
% TropOMI(TF_arrival,:) = [];

% calculate the time stamp for each rotated pixel, that will be used to
% pair with Pandora
TropOMI.UTC_4pair_rotation = TropOMI.utc_time + hours(TropOMI.arrival_time);

% calculate the time stamp for each rotated pixel, that will be used to
% pair with Pandora
TropOMI.UTC_4pair_rotation_at_site2 = TropOMI.utc_time + hours(TropOMI.arrival_time2);

% for all TropOMI data, pair each satellite data with a closest (in time) Pandora
% measurement
time_tolerance = 15;% min    
TF = (abs(TropOMI.y1) <= 25) &(hours(abs(TropOMI.arrival_time)) <= minutes(40));
TropOMI(~TF,:)=[];
%idx_upwind = (hours(TropOMI.arrival_time) > minutes(time_tolerance)) & (TropOMI.y1_2site > 10);
%idx_downwind = (hours(TropOMI.arrival_time) < -minutes(time_tolerance)) & (TropOMI.y1_2site < -10);
%idx_original = (hours(abs(TropOMI.arrival_time)) <= minutes(time_tolerance)) & (abs(TropOMI.y1_2site) <= 10);

%idx_upwind = (hours(TropOMI.arrival_time) > 0);
%idx_downwind = (hours(TropOMI.arrival_time2) < 0);
%idx_original = (hours(TropOMI.arrival_time) <0) & (hours(TropOMI.arrival_time2) >0);

% idx_upwind = (TropOMI.arrival_time > 1);
% idx_downwind = (TropOMI.arrival_time2 < -1);
% idx_original = (TropOMI.arrival_time <= 1) & (TropOMI.arrival_time2 >= -1) & (abs(TropOMI.arrival_time - TropOMI.arrival_time2))<=1; 
idx_upwind = (TropOMI.y1_st2_2site > 0);
idx_downwind = (TropOMI.y1_st2_2site < 0);
%idx_original = abs(TropOMI.arrival_time2) <1 & abs(TropOMI.arrival_time) < 1 & abs(TropOMI.y1_2site) < 20 & abs(TropOMI.y1_2site - TropOMI.y1_st2_2site) < 20;
idx_original = [];


idx_original = (hours(abs(TropOMI.arrival_time)) <= minutes(time_tolerance)) & (abs(TropOMI.y1) <= 10);

TropOMI.type = NaN(height(TropOMI),1);
TropOMI.type(idx_upwind,:) = 1;
TropOMI.type(idx_downwind,:) = -1;
TropOMI.type(idx_original,:) = 0;

j = 1;
for i = 1:height(TropOMI)
    [min_value,min_index] = min(abs(Pandora.UTC-TropOMI.UTC_4pair_rotation(i)));
    if min_value <= minutes(time_tolerance)
        C_rotated_match(j,:) = [TropOMI(i,:),Pandora(min_index,:)];
        j = j +1;
    end
end

names = Pandora2.Properties.VariableNames;
for i = 1:length(names)
    Pandora2.Properties.VariableNames{names{1,i}} = [names{1,i} '_site2'];
end
    
  
j = 1;
for i = 1:height(C_rotated_match)
    [min_value,min_index] = min(abs(Pandora2.UTC_site2 - C_rotated_match.UTC_4pair_rotation_at_site2(i)));
    if min_value <= minutes(time_tolerance)
        C_rotated_match_2sites(j,:) = [C_rotated_match(i,:),Pandora2(min_index,:)];
        j = j +1;
    end
end

%%
function C_output = pair_Pandora_TropOMI_rotation(Pandora,TropOMI,distance,y_dis,arrival_time_limit)
% simple QC
TF_qa = TropOMI.qa == 1;
TropOMI(~TF_qa,:) = [];

% label Pandora
Pandora.meas_serieal_no = (1:height(Pandora))';

% filter by distance
%distance = 10; % km --> v1 and v2 used this
%distance = 7; % km --> v3 and v4
TF = abs(TropOMI.x1) <= distance ;
TropOMI = TropOMI(TF,:);

% calculate arrival time
TropOMI.arrival_time = TropOMI.y1./TropOMI.windspeed; %[hr] --> y1 is in [km], windspeed is in [km/hr]
% note --> negative y1 is down wind pixel, so, negative arrival_time is also for
% downwind pixel

% 2nd type arrival_time 
%TropOMI.arrival_time = (TropOMI.y1- (distance.^2-TropOMI.x1.^2).^0.5)./TropOMI.windspeed; 

% QC based on transportation time, note that the wind field is hourly, so,
% it is only "accurate" if the transportation time is within 1hr
%TF_arrival = abs(TropOMI.arrival_time) > 1;
%TropOMI(TF_arrival,:) = [];

% calculate the time stamp for each rotated pixel, that will be used to
% pair with Pandora
TropOMI.UTC_4pair_rotation = TropOMI.utc_time + hours(TropOMI.arrival_time);


% for all TropOMI data, pair each satellite data with a closest (in time) Pandora
% measurement
time_tolerance = 10;% [min]
%time_tolerance = 60;% [min]
%y_dis = 10;% distance in y direction  --> v1 and v2 used this
%y_dis = 7*5;% distance in y direction  --> v3
%y_dis = 100;% distance in y direction  --> v4
%idx_upwind = (hours(TropOMI.arrival_time) > minutes(time_tolerance)) & (TropOMI.y1 > y_dis);
%idx_downwind = (hours(TropOMI.arrival_time) < -minutes(time_tolerance)) & (TropOMI.y1 < -y_dis);
% idx_upwind = (hours(TropOMI.arrival_time) > minutes(time_tolerance)) & (TropOMI.y1 > y_dis);
% idx_downwind = (hours(TropOMI.arrival_time) < -minutes(time_tolerance)) & (TropOMI.y1 < -y_dis);
% idx_original = (hours(abs(TropOMI.arrival_time)) <= minutes(time_tolerance)) & (abs(TropOMI.y1) <= y_dis);
idx_upwind =  (TropOMI.y1 > y_dis) & (hours(abs(TropOMI.arrival_time)) <= hours(7));
idx_downwind = (TropOMI.y1 < -y_dis) & (hours(abs(TropOMI.arrival_time)) <= hours(7));
idx_original = (abs(TropOMI.y1) <= y_dis) & (hours(abs(TropOMI.arrival_time)) < hours(arrival_time_limit)) ;


TropOMI.type = NaN(height(TropOMI),1); % --> should include this!!!
TropOMI.type(idx_upwind,:) = 1;
TropOMI.type(idx_downwind,:) = -1;
TropOMI.type(idx_original,:) = 0;

j = 1;
for i = 1:height(TropOMI)
    [min_value,min_index] = min(abs(Pandora.UTC-TropOMI.UTC_4pair_rotation(i)));
    if min_value <= minutes(time_tolerance)
        C_rotated_match(j,:) = [TropOMI(i,:),Pandora(min_index,:)];
        j = j +1;
    end
end

% % remove duplicate matching ... 
% [~,ia,ic] = unique(C_rotated_match.meas_serieal_no);% Unique values
% [Ncount, ~, idxcount] = histcounts(ic,numel(ia));% count unique values
% idxkeep = Ncount(idxcount)>1;% Where is greater than one occurence
% C_nonunique = C_rotated_match(idxkeep,:);% Extract from C
% C_realunique = C_rotated_match(~idxkeep,:);% Extract from C
% 
% means_serial_need_loop = unique(C_nonunique.meas_serieal_no);
% for i = 1:numel(means_serial_need_loop)
%     TF = C_nonunique.meas_serieal_no == means_serial_need_loop(i,1);
%     sub_C_nonunique = C_nonunique(TF,:);
%     [~,index_min] = min(abs(sub_C_nonunique.arrival_time));% if oversampled, use the closest in time pixel
%     %[~,index_min] = min(abs(sub_C_nonunique.x1));% if oversampled, use the closest in spatial pixel
%     
%     C_selected(i,:) = sub_C_nonunique(index_min,:);
% end
% 
% C_output = [C_realunique;C_selected];
% C_output = sortrows(C_output,'meas_serieal_no');


C_output = C_rotated_match;
%%
function C_closest = pair_Pandora_TropOMI(Pandora,TropOMI,delta_t,up_down_original,use_sum_no2)
% simple QC
TF_qa = TropOMI.qa == 1;
TropOMI(~TF_qa,:) = [];

% filter by distance
distance = 10; % km
if strcmp(up_down_original,'original')
    TropOMI.UTC_4pair = TropOMI.utc_time;
    %TropOMI = table2timetable(TropOMI,'RowTimes','UTC');
    TF = TropOMI.distence2stn_original <= distance;
    TropOMI.distance_4pair = TropOMI.distence2stn_original;
    TropOMI.type = repmat(0,height(TropOMI),1);
elseif strcmp(up_down_original,'up wind')
    TropOMI.UTC_4pair = TropOMI.utc_time - hours(delta_t);
    %TropOMI = table2timetable(TropOMI,'RowTimes','UTC');
    TF = TropOMI.distence2stn_trans_upw <= distance;
    TropOMI.distance_4pair = TropOMI.distence2stn_trans_upw;
    TropOMI.type = repmat(-1,height(TropOMI),1);
elseif strcmp(up_down_original,'down wind')
    TropOMI.UTC_4pair = TropOMI.utc_time + hours(delta_t);
    %TropOMI = table2timetable(TropOMI,'RowTimes','UTC');
    TF = TropOMI.distence2stn_trans_dnw <= distance;
    TropOMI.distance_4pair = TropOMI.distence2stn_trans_dnw;
    TropOMI.type = repmat(1,height(TropOMI),1);
end
TropOMI(~TF,:) = [];

% for all TropOMI data, pair each satellite data with a closest (in time) Pandora
% measurement
time_tolerance = 10;% min    
j = 1;
for i = 1:height(TropOMI)
    [min_value,min_index] = min(abs(Pandora.UTC-TropOMI.UTC_4pair(i)));
    if min_value <= minutes(time_tolerance)
        C(j,:) = [TropOMI(i,:),Pandora(min_index,:)];
        j = j +1;
    end
end

% only keep the TropOMI data that is cloest to the site
Unique_Pandora_UTC = unique(C.UTC);
for i = 1:length(Unique_Pandora_UTC)
    TF_unique_Pandora = C.UTC == Unique_Pandora_UTC(i);
    C_sub = C(TF_unique_Pandora,:);
    [min_distance,min_index] = min(C_sub.distance_4pair);
    C_closest(i,:) = C_sub(min_index,:);
end


%%
DU = 2.6870e+16;
if use_sum_no2
    y = (C_closest.no2_trop + C_closest.no2_strat)./DU;% TropOMI sum no2
else
    y = C_closest.no2./DU;% TropOMI total no2
end
x = C_closest.NO2_VCD;% Pandora
% remove NaN
TF_NaN = isnan(y) | isnan(x);
y = y(~TF_NaN,:);
x = x(~TF_NaN,:);
color_code = C_closest.distance_4pair(~TF_NaN,:);

line_fits(x,y);
scatter(x,y,10,color_code,'filled');
ylabel('TropOMI NO_2 VCD [DU]');
xlabel('Pandora NO_2 VCD [DU]');

cc = colorbar;
cc.Label.String = 'Distance [km]';
title(up_down_original);
xlim([0 1.5]);
ylim([0 1.5]);


%%
function [intercept,slop,slop_nlm,R,N] = auto_scatter(method,data,x_name,y_name,c_name,type,save_fig,output_path,use_sum_no2,fig_nm_append)
if nargin == 9
    fig_nm_append ='';
    fig_title_append = fig_nm_append;
else
    fig_title_append = fig_nm_append;
    fig_nm_append = char(fig_nm_append);
    tf = isspace(fig_nm_append);
    fig_nm_append = fig_nm_append(1,~tf);
    
    fig_nm_append = strrep(fig_nm_append,'<=','ltoe');
    fig_nm_append = strrep(fig_nm_append,'>=','mtoe');
    fig_nm_append = strrep(fig_nm_append,'>','mt');
    fig_nm_append = strrep(fig_nm_append,'<','lt');

end
            
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

switch x_name
    case 'Pandora'
        x = data.NO2_VCD(TF,:);% Pandora --> 1st site
        x_label = 'Pandora VCD [DU]';
    case 'Pandora2'
        x = data.NO2_VCD_site2(TF,:);% Pandora --> 2nd site
        x_label = 'Pandora VCD (site2) [DU]';
    case 'TropOMI'
        if use_sum_no2
            x = (data.no2_trop(TF,:) + data.no2_strat(TF,:))./DU;% TropOMI sum no2
        else
            x = data.no2(TF,:)./DU;% TropOMI total no2
        end
        x_label = 'TropOMI VCD [DU]';
end

switch y_name
    case 'Pandora'
        y = data.NO2_VCD(TF,:);% Pandora --> 1st site [DU]
        y_label = 'Pandora VCD [DU]';
    case 'Pandora2'
        y = data.NO2_VCD_site2(TF,:);% Pandora --> 2nd site [DU]
        y_label = 'Pandora VCD (site2) [DU]';
    case 'TropOMI'
        if use_sum_no2
            y = (data.no2_trop(TF,:) + data.no2_strat(TF,:))./DU;% TropOMI sum no2 [DU]
        else
            y = data.no2(TF,:)./DU;% TropOMI total [DU]
        end
        y_label = 'TropOMI VCD [DU]';
end

switch c_name
    case 't1'
        c = data.arrival_time(TF,:);% [hr]
        c_label = 'transport time [hr]';
    case 't2'
        c = data.arrival_time2(TF,:);% [hr]
        c_label = 'transport time (site2) [hr]';
    case 'd1'
        c = data.y1(TF,:);% [km]
        c_label = 'transport distance [km]';
    case 'd2'
        c = data.y1_2site(TF,:);% [km]
        c_label = 'transport distance (site2) [km]';
    case 'wd'
        c = data.winddirection(TF,:);% [km]
        c_label = 'wind direction [degree]';
end

fig_name = ['scatter_' method '_' y_name '_vs_' x_name '_' type(~isspace(type)) '_codeby_' c_name '_' fig_nm_append];
% remove NaN
TF_NaN = isnan(x) | isnan(y);
x = x(~TF_NaN,1);
y = y(~TF_NaN,1);
c = c(~TF_NaN,1);

[intercept,slop,slop_nlm,mdl_lm,mdl_nlm] = line_fits(x,y);       
R = mdl_lm.Rsquared.Ordinary.^0.5;
N = mdl_lm.NumObservations  ;

scatter(x,y,10,c,'filled');
ylabel(y_label);
xlabel(x_label);
title([type fig_title_append]);
h = colorbar;
ylabel(h,c_label);
legend({[y_name ' vs. ' x_name],'linear fit','linear fit (no offset)','1-on-1',type});
xlim([0 1.4]);
ylim([0 1.4]);

print_setting(1/4,save_fig,[output_path fig_name]);

%%
function [status,bins_all] = bin_by_x(method,data,x_name,y_name,type,save_fig,output_path,use_sum_no2, name_ext)
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
        if use_sum_no2
            x = (data.no2_trop + data.no2_strat)./DU;% TropOMI sum no2
        else
            x = data.no2./DU;% TropOMI total no2
        end
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
        if use_sum_no2
            y = (data.no2_trop + data.no2_strat)./DU;% TropOMI sum no2 [DU]
        else
            y = data.no2./DU;% TropOMI total no2 [DU]
        end
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

figure;hold all;
fig_name = ['binned_scatter_' method '_' y_name '_vs_' x_name '_' type(~isspace(type)) '_' name_ext];
x = status.x_centre;
y = status.y_mean;
x_err = status.y_std_err;
%line_fits(x,y);   
errorbar(x,y,x_err,'.');
%scatter(x,y,10,c,'filled');
ylabel(y_label);
xlabel(x_label);
title(type);
legend({[y_name ' vs. ' x_name]});
%legend({[y_name ' vs. ' x_name],'linear fit','linear fit (no offset)','1-on-1'},'err of mean');
xlim([0 1.4]);
ylim([0 1.4]);

print_setting(1/4,save_fig,[output_path fig_name]);