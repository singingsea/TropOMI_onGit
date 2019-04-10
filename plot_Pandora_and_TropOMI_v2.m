function data = plot_Pandora_and_TropOMI_v2()

%site = 'Downsview';Pandora_no = '104';
%site = 'Egbert';Pandora_no = '108';
%site = 'StGeorge';Pandora_no = '145';
site = 'FortMcKay';Pandora_no = '122';
TropOMI_type = 'OFFL';
avg_type = 'closest';%'avg24', 'avg7'

ERA_merged = true;% use the ERA merged data or not

DU = 2.6870e+16;
addpath('C:\Users\ZhaoX\Documents\MATLAB\matlab');
save_fig = 1;
if ERA_merged == true
    plot_path = ['C:\Projects\TropOMI\plots\' TropOMI_type '\' site '\' avg_type '\' Pandora_no '_ERA\'];
    filename = ['C:\Projects\TropOMI\data\NO2_merged_with_Pandora_ERA\TropOMI_' TropOMI_type '_' site '_' avg_type '_Pandora' Pandora_no site '_simple.csv'];
    data = importfile_simple_ERA(filename);
    Pandora_filename = ['C:\Projects\TropOMI\data\NO2_merged_with_Pandora\Pandora' Pandora_no site '.csv'];
    Pandora = importfile_Pandora_raw(Pandora_filename);
else
    plot_path = ['C:\Projects\TropOMI\plots\' TropOMI_type '\' site '\' avg_type '\' Pandora_no '\'];
    filename = ['C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_' TropOMI_type '_' site '_' avg_type '_Pandora' Pandora_no site '_simple.csv'];
    data = importfile_merged_simple(filename);
    Pandora_filename = ['C:\Projects\TropOMI\data\NO2_merged_with_Pandora\Pandora' Pandora_no site '.csv'];
    Pandora = importfile_Pandora_raw(Pandora_filename);
end
mkdir(plot_path);
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_OFFL_Downsview_avg7km_Pandora103Downsview.csv';
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_OFFL_Downsview_closest_Pandora103Downsview.csv';
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_OFFL_Downsview_avg24km_Pandora103Downsview.csv'
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_NRTI_Downsview_avg7km_Pandora103Downsview.csv'
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_NRTI_Downsview_closest_Pandora103Downsview.csv'
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_NRTI_Downsview_avg24km_Pandora103Downsview.csv'
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_OFFL_FortMcKay_avg7km_Pandora122FortMcKay.csv';
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_OFFL_FortMcKay_closest_Pandora122FortMcKay.csv';
%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_NRTI_FortMcKay_closest_Pandora122FortMcKay.csv';

%filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\TropOMI_OFFL_Egbert_closest_Pandora108Egbert_simple.csv'

% try
%     data = importfile(filename);
% catch
%     data = importfile_closest(filename);
% end

% try
%     data = importfile_P122_avg(filename);
% catch
%     data = importfile_P122_closest(filename);
% end
    

%Pandora_filename = 'C:\Projects\TropOMI\data\NO2_merged_with_Pandora\Pandora103Downsview.csv';
%Pandora = importfile_Pandora_raw(Pandora_filename);

%data.TropOMI_NO2 = data.no2;
%data.Pandora_NO2 = data.Column8NitrogendioxidetotalverticalcolumnamountDobsonUnits9e99r;
%data.Pandora_NO2_err = data.Column9UncertaintyofnitrogendioxidetotalverticalcolumnamountDob;
%data.TropOMI_datetime = data.Datetime_coarse;
%data.Pandora_datetime = data.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86;

TF = isnan(data.Pandora_NO2);
data(TF,:) = [];

%% give a trick TF
%   TF_trick1 = abs(data.TropOMI_NO2./DU-data.Pandora_NO2) > 0.8;
%   TF_trick2 = data.d./1000 > 20;
%   TF_trick = TF_trick1 | TF_trick2;
%   data(TF_trick,:)=[];
%data(TF_trick1,:)=[];
%data(TF_trick2,:)=[];
u = data.u_1000hPa;
v = data.v_1000hPa;
wind_direction_a = atan2d(u,v);
data.wind_speed_1000hPa = hypot(u,v);% wind speed [m/s]
data.wind_direction_1000hPa = rem(360+wind_direction_a, 360);% wind direction in range of [0 360]

u = data.u_900hPa;
v = data.v_900hPa;
wind_direction_a = atan2d(u,v);
data.wind_speed_900hPa = hypot(u,v);
data.wind_direction_900hPa = rem(360+wind_direction_a, 360);


u = data.u_800hPa;
v = data.v_800hPa;
wind_direction_a = atan2d(u,v);
data.wind_speed_800hPa = hypot(u,v);
data.wind_direction_800hPa = rem(360+wind_direction_a, 360);


u = data.u_700hPa;
v = data.v_700hPa;
wind_direction_a = atan2d(u,v);
data.wind_speed_700hPa = hypot(u,v);
data.wind_direction_700hPa = rem(360+wind_direction_a, 360);

figure;
fig_name = 'Pandora_TropOMI_winfunction_1000hPa';
wind_direction = data.wind_direction_1000hPa.*(pi/180);
rho = data.wind_speed_1000hPa;
data.rel_diff = (data.TropOMI_NO2./DU - data.Pandora_NO2)./((data.Pandora_NO2 + data.TropOMI_NO2./DU)./2).*100;
C = data.rel_diff;
sz = data.Pandora_NO2.*300;
plot_table = table;
plot_table.wind_direction = wind_direction;
plot_table.rho = rho;
plot_table.sz = sz;
plot_table.C = C;
plot_table = sortrows(plot_table,'C','descend');
polarscatter([plot_table.wind_direction;280.*(pi/180)],[plot_table.rho;20],[plot_table.sz;1*300],[plot_table.C;0],'filled','MarkerFaceAlpha',.6);
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
colormap(jet);
colorbar;rlim([0 20]);caxis([-100 100]);
title(['TropOMI NO_2 and Pandora \Delta_r_e_l [%] Wind function (1000 hPa)']);
text(300.*(pi/180),34,[site]);
text(284.*(pi/180),30.5,['Scale:']);
text(280.*(pi/180),30,['NO_2 VCD = 1 DU ']);
text(276.*(pi/180),30,['\Delta_r_e_l = 0 %']);
print_setting(1/2,save_fig,[plot_path fig_name]);

figure;
fig_name = 'Pandora_TropOMI_winfunction_900hPa';
wind_direction = data.wind_direction_900hPa.*(pi/180);
rho = data.wind_speed_900hPa;
data.rel_diff = (data.TropOMI_NO2./DU - data.Pandora_NO2)./((data.Pandora_NO2 + data.TropOMI_NO2./DU)./2).*100;
C = data.rel_diff;
sz = data.Pandora_NO2.*300;
plot_table = table;
plot_table.wind_direction = wind_direction;
plot_table.rho = rho;
plot_table.sz = sz;
plot_table.C = C;
plot_table = sortrows(plot_table,'C','descend');
polarscatter([plot_table.wind_direction;280.*(pi/180)],[plot_table.rho;20],[plot_table.sz;1*300],[plot_table.C;0],'filled','MarkerFaceAlpha',.6);
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
colormap(jet);
colorbar;rlim([0 20]);caxis([-100 100]);
title(['TropOMI NO_2 and Pandora \Delta_r_e_l [%] Wind function (900 hPa)']);
text(300.*(pi/180),34,[site]);
text(284.*(pi/180),30.5,['Scale:']);
text(280.*(pi/180),30,['NO_2 VCD = 1 DU ']);
text(276.*(pi/180),30,['\Delta_r_e_l = 0 %']);
print_setting(1/2,save_fig,[plot_path fig_name]);


figure;
fig_name = 'Pandora_TropOMI_winfunction_800hPa';
wind_direction = data.wind_direction_800hPa.*(pi/180);
rho = data.wind_speed_800hPa;
data.rel_diff = (data.TropOMI_NO2./DU - data.Pandora_NO2)./((data.Pandora_NO2 + data.TropOMI_NO2./DU)./2).*100;
C = data.rel_diff;
sz = data.Pandora_NO2.*300;
plot_table = table;
plot_table.wind_direction = wind_direction;
plot_table.rho = rho;
plot_table.sz = sz;
plot_table.C = C;
plot_table = sortrows(plot_table,'C','descend');
polarscatter([plot_table.wind_direction;280.*(pi/180)],[plot_table.rho;20],[plot_table.sz;1*300],[plot_table.C;0],'filled','MarkerFaceAlpha',.6);
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
colormap(jet);
colorbar;rlim([0 20]);caxis([-100 100]);
title(['TropOMI NO_2 and Pandora \Delta_r_e_l [%] Wind function (800 hPa)']);
text(300.*(pi/180),34,[site]);
text(284.*(pi/180),30.5,['Scale:']);
text(280.*(pi/180),30,['NO_2 VCD = 1 DU ']);
text(276.*(pi/180),30,['\Delta_r_e_l = 0 %']);
print_setting(1/2,save_fig,[plot_path fig_name]);

figure;
fig_name = 'Pandora_TropOMI_winfunction_700hPa';
wind_direction = data.wind_direction_700hPa.*(pi/180);
rho = data.wind_speed_700hPa;
data.rel_diff = (data.TropOMI_NO2./DU - data.Pandora_NO2)./((data.Pandora_NO2 + data.TropOMI_NO2./DU)./2).*100;
C = data.rel_diff;
sz = data.Pandora_NO2.*300;
plot_table = table;
plot_table.wind_direction = wind_direction;
plot_table.rho = rho;
plot_table.sz = sz;
plot_table.C = C;
plot_table = sortrows(plot_table,'C','descend');
polarscatter([plot_table.wind_direction;280.*(pi/180)],[plot_table.rho;20],[plot_table.sz;1*300],[plot_table.C;0],'filled','MarkerFaceAlpha',.6);
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
colormap(jet);
colorbar;rlim([0 20]);caxis([-100 100]);
title(['TropOMI NO_2 and Pandora \Delta_r_e_l [%] Wind function (700 hPa)']);
text(300.*(pi/180),34,[site]);
text(284.*(pi/180),30.5,['Scale:']);
text(280.*(pi/180),30,['NO_2 VCD = 1 DU ']);
text(276.*(pi/180),30,['\Delta_r_e_l = 0 %']);
print_setting(1/2,save_fig,[plot_path fig_name]);

%% check the quality of merge
fig_name = 'TropOMI_vs_pandora_merge_quality';
data.TropOMI_datenum = datenum(data.TropOMI_datetime);
data.Pandora_datenum = datenum(data.Pandora_datetime);
line_fits(data.Pandora_datenum,data.TropOMI_datenum);
datetick('x','mmm-dd','keeplimits');
datetick('y','mmm-dd','keeplimits');
ylabel(['TropOMI measurement time']);
xlabel(['Pandora measurement time']);
print_setting(1/4,save_fig,[plot_path fig_name]);

%% check the quality of the data
fig_name = 'TropOMI_vs_Pandora_total_column_codebydis';
line_fits(data.Pandora_NO2,data.TropOMI_NO2./DU);
%H=errorbarxy(data.Pandora_NO2,data.TropOMI_NO2./DU,data.,ye,{'ko-', 'b', 'r'});
scatter(data.Pandora_NO2,data.TropOMI_NO2./DU,15,data.d./1000,'filled');
ylabel(['TropOMI NO_2 [DU]']);
xlabel(['Pandora NO_2 [DU]']);
title(['Color coded by distance [km]']);
colorbar;
xlim([0 1.1]);
ylim([0 1.1]);
print_setting(1/4,save_fig,[plot_path fig_name]);

%% check the quality of the data
fig_name = 'TropOMI_vs_Pandora_total_column_codebytime';
line_fits(data.Pandora_NO2,data.TropOMI_NO2./DU);
delta_t = (data.TropOMI_datenum - data.Pandora_datenum).*24.*60;
scatter(data.Pandora_NO2,data.TropOMI_NO2./DU,15,delta_t,'filled');
ylabel(['TropOMI NO_2 [DU]']);
xlabel(['Pandora NO_2 [DU]']);
title(['Color coded by delta time (TropOMI - Pandora)[min]']);
colorbar;
xlim([0 1.1]);
ylim([0 1.1]);
print_setting(1/4,save_fig,[plot_path fig_name]);

%% time series
fig_name = 'Pandora_vs_TropOMI_timeserise';
figure;hold all;
plot(Pandora.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86,Pandora.Column8NitrogendioxidetotalverticalcolumnamountDobsonUnits9e99r,'k.')
plot(data.TropOMI_datetime,data.TropOMI_NO2./DU,'rs');
plot(data.Pandora_datetime,data.Pandora_NO2,'bs');
ylabel(['NO_2 [DU]']);
legend({'Pandora','TropOMI','Pandora (conincident)'});
%legend({'TropOMI','Pandora (conincident)'});
print_setting(1/2,save_fig,[plot_path fig_name]);
ylim([0 1.8]);


%%
function TropOMIOFFLDownsviewavg7kmPandora103Downsview = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIOFFLDOWNSVIEWAVG7KMPANDORA103DOWNSVIEW = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   TROPOMIOFFLDOWNSVIEWAVG7KMPANDORA103DOWNSVIEW = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   TropOMIOFFLDownsviewavg7kmPandora103Downsview = importfile('TropOMI_OFFL_Downsview_avg7km_Pandora103Downsview.csv', 2, 19);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/06/20 10:50:01

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

dateFormats = {'yyyy-MM-dd HH:mm:ss', 'yyyyMMdd''T''HHmmss''Z'''};
dateFormatIndex = 1;
blankDates = cell(1,size(raw,2));
anyBlankDates = false(size(raw,1),1);
invalidDates = cell(1,size(raw,2));
anyInvalidDates = false(size(raw,1),1);
for col=[2,8]% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
    try
        dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[2,8]}, 'InputFormat', dateFormats{col==[2,8]}); %#ok<AGROW>
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, 'UniformOutput', false);
            dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[2,8]}, 'InputFormat', dateFormats{col==[2,8]}); %#ok<AGROW>
        catch
            dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<AGROW>
        end
    end
    
    dateFormatIndex = dateFormatIndex + 1;
    blankDates{col} = dataArray{col} == '';
    anyBlankDates = blankDates{col} | anyBlankDates;
    invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
    anyInvalidDates = invalidDates{col} | anyInvalidDates;
end
dates = dates(:,[2,8]);
blankDates = blankDates(:,[2,8]);
invalidDates = invalidDates(:,[2,8]);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43]);
rawStringColumns = string(raw(:, [7,44,45,46,47]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
TropOMIOFFLDownsviewavg7kmPandora103Downsview = table;
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Datenum = cell2mat(rawNumericColumns(:, 1));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Datetime_coarse = dates{:, 1};
TropOMIOFFLDownsviewavg7kmPandora103Downsview.no2_trop = cell2mat(rawNumericColumns(:, 2));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.no2_trop_err = cell2mat(rawNumericColumns(:, 3));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.d = cell2mat(rawNumericColumns(:, 4));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.no2 = cell2mat(rawNumericColumns(:, 5));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.UTC = rawStringColumns(:, 1);
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86 = dates{:, 2};
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column2Fractionaldayssince1Jan2000UTmidnightforcenterofmeasurem = cell2mat(rawNumericColumns(:, 6));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column3Effectivedurationofmeasurementinseconds = cell2mat(rawNumericColumns(:, 7));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column4Solarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 8));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column5Solarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 9));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column6Lunarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 10));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column7Lunarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 11));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column8NitrogendioxidetotalverticalcolumnamountDobsonUnits9e99r = cell2mat(rawNumericColumns(:, 12));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column9UncertaintyofnitrogendioxidetotalverticalcolumnamountDob = cell2mat(rawNumericColumns(:, 13));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column10Directnitrogendioxideairmassfactor = cell2mat(rawNumericColumns(:, 14));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column11Diffusecorrectionappliedbeforefittingateffectivefitting = cell2mat(rawNumericColumns(:, 15));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column12L2dataqualityflagfornitrogendioxide0assuredhighquality1 = cell2mat(rawNumericColumns(:, 16));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column13Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 17));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column14Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 18));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column15Fittingresultindex012noerrororwarning34warning4error = cell2mat(rawNumericColumns(:, 19));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column16Normalizedrmsofspectralfittingresidualsweightedwithmeas = cell2mat(rawNumericColumns(:, 20));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column17Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 21));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column18Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 22));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column19Climatologicalstationpressurembar = cell2mat(rawNumericColumns(:, 23));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column20Dataprocessingtypeindex = cell2mat(rawNumericColumns(:, 24));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column21Calibrationfileversion = cell2mat(rawNumericColumns(:, 25));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column22Calibrationfilevaliditystartingdate = cell2mat(rawNumericColumns(:, 26));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column23Level2Fitdataqualityflag0assuredhighquality1assuredmedi = cell2mat(rawNumericColumns(:, 27));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column24Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 28));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column25Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 29));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column26Level1dataqualityflag0assuredhighquality1assuredmediumq = cell2mat(rawNumericColumns(:, 30));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column27Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 31));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column28Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 32));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column29WavelengtheffectivetemperatureC999noeffectivetemperatur = cell2mat(rawNumericColumns(:, 33));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column30Estimatedaverageresidualstraylightlevelonlyvalidforstra = cell2mat(rawNumericColumns(:, 34));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column31Retrievedwavelengthshiftfromlevel1datanm9nowavelengthch = cell2mat(rawNumericColumns(:, 35));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column32Retrievedwavelengthshiftfromspectralfittingnm9nowavelen = cell2mat(rawNumericColumns(:, 36));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column33Integrationtimems = cell2mat(rawNumericColumns(:, 37));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column34Numberofdarkcountcycles = cell2mat(rawNumericColumns(:, 38));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column35Effectivepositionoffilterwheel10filterwheelnotused19are = cell2mat(rawNumericColumns(:, 39));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column36Effectivepositionoffilterwheel20filterwheelnotused19are = cell2mat(rawNumericColumns(:, 40));
TropOMIOFFLDownsviewavg7kmPandora103Downsview.instrument = rawStringColumns(:, 2);
TropOMIOFFLDownsviewavg7kmPandora103Downsview.location = rawStringColumns(:, 3);
TropOMIOFFLDownsviewavg7kmPandora103Downsview.time = rawStringColumns(:, 4);
TropOMIOFFLDownsviewavg7kmPandora103Downsview.LTC = rawStringColumns(:, 5);

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIOFFLDownsviewavg7kmPandora103Downsview.Datetime_coarse=datenum(TropOMIOFFLDownsviewavg7kmPandora103Downsview.Datetime_coarse);TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86=datenum(TropOMIOFFLDownsviewavg7kmPandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86);

%%
function TropOMIOFFLDownsviewclosestPandora103Downsview = importfile_closest(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIOFFLDOWNSVIEWCLOSESTPANDORA103DOWNSVIEW = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   TROPOMIOFFLDOWNSVIEWCLOSESTPANDORA103DOWNSVIEW = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   TropOMIOFFLDownsviewclosestPandora103Downsview = importfile('TropOMI_OFFL_Downsview_closest_Pandora103Downsview.csv', 2, 35);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/06/20 11:09:03

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,5,6,7,8,10,11,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

dateFormats = {'yyyy-MM-dd HH:mm:ss', 'yyyyMMdd''T''HHmmss''Z'''};
dateFormatIndex = 1;
blankDates = cell(1,size(raw,2));
anyBlankDates = false(size(raw,1),1);
invalidDates = cell(1,size(raw,2));
anyInvalidDates = false(size(raw,1),1);
for col=[12,14]% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
    try
        dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[12,14]}, 'InputFormat', dateFormats{col==[12,14]}); %#ok<AGROW>
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, 'UniformOutput', false);
            dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[12,14]}, 'InputFormat', dateFormats{col==[12,14]}); %#ok<AGROW>
        catch
            dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<AGROW>
        end
    end
    
    dateFormatIndex = dateFormatIndex + 1;
    blankDates{col} = dataArray{col} == '';
    anyBlankDates = blankDates{col} | anyBlankDates;
    invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
    anyInvalidDates = invalidDates{col} | anyInvalidDates;
end
dates = dates(:,[12,14]);
blankDates = blankDates(:,[12,14]);
invalidDates = invalidDates(:,[12,14]);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,2,3,5,6,7,8,10,11,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49]);
rawStringColumns = string(raw(:, [4,9,13,50,51,52,53]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
TropOMIOFFLDownsviewclosestPandora103Downsview = table;
TropOMIOFFLDownsviewclosestPandora103Downsview.lon = cell2mat(rawNumericColumns(:, 1));
TropOMIOFFLDownsviewclosestPandora103Downsview.lat = cell2mat(rawNumericColumns(:, 2));
TropOMIOFFLDownsviewclosestPandora103Downsview.qa = cell2mat(rawNumericColumns(:, 3));
TropOMIOFFLDownsviewclosestPandora103Downsview.utc_time = rawStringColumns(:, 1);
TropOMIOFFLDownsviewclosestPandora103Downsview.no2_trop = cell2mat(rawNumericColumns(:, 4));
TropOMIOFFLDownsviewclosestPandora103Downsview.no2_trop_err = cell2mat(rawNumericColumns(:, 5));
TropOMIOFFLDownsviewclosestPandora103Downsview.d = cell2mat(rawNumericColumns(:, 6));
TropOMIOFFLDownsviewclosestPandora103Downsview.no2 = cell2mat(rawNumericColumns(:, 7));
TropOMIOFFLDownsviewclosestPandora103Downsview.time_x = rawStringColumns(:, 2);
TropOMIOFFLDownsviewclosestPandora103Downsview.Datenum = cell2mat(rawNumericColumns(:, 8));
TropOMIOFFLDownsviewclosestPandora103Downsview.hour = cell2mat(rawNumericColumns(:, 9));
TropOMIOFFLDownsviewclosestPandora103Downsview.Datetime_coarse = dates{:, 1};
TropOMIOFFLDownsviewclosestPandora103Downsview.UTC = rawStringColumns(:, 3);
TropOMIOFFLDownsviewclosestPandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86 = dates{:, 2};
TropOMIOFFLDownsviewclosestPandora103Downsview.Column2Fractionaldayssince1Jan2000UTmidnightforcenterofmeasurem = cell2mat(rawNumericColumns(:, 10));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column3Effectivedurationofmeasurementinseconds = cell2mat(rawNumericColumns(:, 11));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column4Solarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 12));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column5Solarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 13));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column6Lunarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 14));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column7Lunarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 15));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column8NitrogendioxidetotalverticalcolumnamountDobsonUnits9e99r = cell2mat(rawNumericColumns(:, 16));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column9UncertaintyofnitrogendioxidetotalverticalcolumnamountDob = cell2mat(rawNumericColumns(:, 17));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column10Directnitrogendioxideairmassfactor = cell2mat(rawNumericColumns(:, 18));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column11Diffusecorrectionappliedbeforefittingateffectivefitting = cell2mat(rawNumericColumns(:, 19));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column12L2dataqualityflagfornitrogendioxide0assuredhighquality1 = cell2mat(rawNumericColumns(:, 20));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column13Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 21));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column14Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 22));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column15Fittingresultindex012noerrororwarning34warning4error = cell2mat(rawNumericColumns(:, 23));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column16Normalizedrmsofspectralfittingresidualsweightedwithmeas = cell2mat(rawNumericColumns(:, 24));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column17Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 25));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column18Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 26));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column19Climatologicalstationpressurembar = cell2mat(rawNumericColumns(:, 27));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column20Dataprocessingtypeindex = cell2mat(rawNumericColumns(:, 28));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column21Calibrationfileversion = cell2mat(rawNumericColumns(:, 29));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column22Calibrationfilevaliditystartingdate = cell2mat(rawNumericColumns(:, 30));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column23Level2Fitdataqualityflag0assuredhighquality1assuredmedi = cell2mat(rawNumericColumns(:, 31));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column24Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 32));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column25Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 33));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column26Level1dataqualityflag0assuredhighquality1assuredmediumq = cell2mat(rawNumericColumns(:, 34));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column27Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 35));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column28Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 36));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column29WavelengtheffectivetemperatureC999noeffectivetemperatur = cell2mat(rawNumericColumns(:, 37));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column30Estimatedaverageresidualstraylightlevelonlyvalidforstra = cell2mat(rawNumericColumns(:, 38));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column31Retrievedwavelengthshiftfromlevel1datanm9nowavelengthch = cell2mat(rawNumericColumns(:, 39));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column32Retrievedwavelengthshiftfromspectralfittingnm9nowavelen = cell2mat(rawNumericColumns(:, 40));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column33Integrationtimems = cell2mat(rawNumericColumns(:, 41));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column34Numberofdarkcountcycles = cell2mat(rawNumericColumns(:, 42));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column35Effectivepositionoffilterwheel10filterwheelnotused19are = cell2mat(rawNumericColumns(:, 43));
TropOMIOFFLDownsviewclosestPandora103Downsview.Column36Effectivepositionoffilterwheel20filterwheelnotused19are = cell2mat(rawNumericColumns(:, 44));
TropOMIOFFLDownsviewclosestPandora103Downsview.instrument = rawStringColumns(:, 4);
TropOMIOFFLDownsviewclosestPandora103Downsview.location = rawStringColumns(:, 5);
TropOMIOFFLDownsviewclosestPandora103Downsview.time_y = rawStringColumns(:, 6);
TropOMIOFFLDownsviewclosestPandora103Downsview.LTC = rawStringColumns(:, 7);

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIOFFLDownsviewclosestPandora103Downsview.Datetime_coarse=datenum(TropOMIOFFLDownsviewclosestPandora103Downsview.Datetime_coarse);TropOMIOFFLDownsviewclosestPandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86=datenum(TropOMIOFFLDownsviewclosestPandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86);

function Pandora103Downsview = importfile_Pandora_raw(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   PANDORA103DOWNSVIEW = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   PANDORA103DOWNSVIEW = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   Pandora103Downsview = importfile('Pandora103Downsview.csv', 2, 9444);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/06/20 11:36:07

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using the
% specified date format.
try
    dates{1} = datetime(dataArray{1}, 'Format', 'yyyyMMdd''T''HHmmss''Z''', 'InputFormat', 'yyyyMMdd''T''HHmmss''Z''');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{1} = cellfun(@(x) x(2:end-1), dataArray{1}, 'UniformOutput', false);
        dates{1} = datetime(dataArray{1}, 'Format', 'yyyyMMdd''T''HHmmss''Z''', 'InputFormat', 'yyyyMMdd''T''HHmmss''Z''');
    catch
        dates{1} = repmat(datetime([NaN NaN NaN]), size(dataArray{1}));
    end
end

dates = dates(:,1);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37]);
rawStringColumns = string(raw(:, [38,39,40,41]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
Pandora103Downsview = table;
Pandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86 = dates{:, 1};
Pandora103Downsview.Column2Fractionaldayssince1Jan2000UTmidnightforcenterofmeasurem = cell2mat(rawNumericColumns(:, 1));
Pandora103Downsview.Column3Effectivedurationofmeasurementinseconds = cell2mat(rawNumericColumns(:, 2));
Pandora103Downsview.Column4Solarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 3));
Pandora103Downsview.Column5Solarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 4));
Pandora103Downsview.Column6Lunarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 5));
Pandora103Downsview.Column7Lunarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 6));
Pandora103Downsview.Column8NitrogendioxidetotalverticalcolumnamountDobsonUnits9e99r = cell2mat(rawNumericColumns(:, 7));
Pandora103Downsview.Column9UncertaintyofnitrogendioxidetotalverticalcolumnamountDob = cell2mat(rawNumericColumns(:, 8));
Pandora103Downsview.Column10Directnitrogendioxideairmassfactor = cell2mat(rawNumericColumns(:, 9));
Pandora103Downsview.Column11Diffusecorrectionappliedbeforefittingateffectivefitting = cell2mat(rawNumericColumns(:, 10));
Pandora103Downsview.Column12L2dataqualityflagfornitrogendioxide0assuredhighquality1 = cell2mat(rawNumericColumns(:, 11));
Pandora103Downsview.Column13Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 12));
Pandora103Downsview.Column14Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 13));
Pandora103Downsview.Column15Fittingresultindex012noerrororwarning34warning4error = cell2mat(rawNumericColumns(:, 14));
Pandora103Downsview.Column16Normalizedrmsofspectralfittingresidualsweightedwithmeas = cell2mat(rawNumericColumns(:, 15));
Pandora103Downsview.Column17Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 16));
Pandora103Downsview.Column18Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 17));
Pandora103Downsview.Column19Climatologicalstationpressurembar = cell2mat(rawNumericColumns(:, 18));
Pandora103Downsview.Column20Dataprocessingtypeindex = cell2mat(rawNumericColumns(:, 19));
Pandora103Downsview.Column21Calibrationfileversion = cell2mat(rawNumericColumns(:, 20));
Pandora103Downsview.Column22Calibrationfilevaliditystartingdate = cell2mat(rawNumericColumns(:, 21));
Pandora103Downsview.Column23Level2Fitdataqualityflag0assuredhighquality1assuredmedi = cell2mat(rawNumericColumns(:, 22));
Pandora103Downsview.Column24Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 23));
Pandora103Downsview.Column25Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 24));
Pandora103Downsview.Column26Level1dataqualityflag0assuredhighquality1assuredmediumq = cell2mat(rawNumericColumns(:, 25));
Pandora103Downsview.Column27Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 26));
Pandora103Downsview.Column28Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 27));
Pandora103Downsview.Column29WavelengtheffectivetemperatureC999noeffectivetemperatur = cell2mat(rawNumericColumns(:, 28));
Pandora103Downsview.Column30Estimatedaverageresidualstraylightlevelonlyvalidforstra = cell2mat(rawNumericColumns(:, 29));
Pandora103Downsview.Column31Retrievedwavelengthshiftfromlevel1datanm9nowavelengthch = cell2mat(rawNumericColumns(:, 30));
Pandora103Downsview.Column32Retrievedwavelengthshiftfromspectralfittingnm9nowavelen = cell2mat(rawNumericColumns(:, 31));
Pandora103Downsview.Column33Integrationtimems = cell2mat(rawNumericColumns(:, 32));
Pandora103Downsview.Column34Numberofdarkcountcycles = cell2mat(rawNumericColumns(:, 33));
Pandora103Downsview.Column35Effectivepositionoffilterwheel10filterwheelnotused19are = cell2mat(rawNumericColumns(:, 34));
Pandora103Downsview.Column36Effectivepositionoffilterwheel20filterwheelnotused19are = cell2mat(rawNumericColumns(:, 35));
Pandora103Downsview.instrument = cell2mat(rawNumericColumns(:, 36));
Pandora103Downsview.location = rawStringColumns(:, 1);
Pandora103Downsview.time = rawStringColumns(:, 2);
Pandora103Downsview.UTC = rawStringColumns(:, 3);
Pandora103Downsview.LTC = rawStringColumns(:, 4);

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% Pandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86=datenum(Pandora103Downsview.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86);

%%
function TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay = importfile_P122_avg(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIOFFLFORTMCKAYAVG7KMPANDORA122FORTMCKAY = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   TROPOMIOFFLFORTMCKAYAVG7KMPANDORA122FORTMCKAY = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay = importfile('TropOMI_OFFL_FortMcKay_avg7km_Pandora122FortMcKay.csv', 2, 33);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/06/20 13:38:27

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

dateFormats = {'yyyy-MM-dd HH:mm:ss', 'yyyyMMdd''T''HHmmss''Z'''};
dateFormatIndex = 1;
blankDates = cell(1,size(raw,2));
anyBlankDates = false(size(raw,1),1);
invalidDates = cell(1,size(raw,2));
anyInvalidDates = false(size(raw,1),1);
for col=[2,8]% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
    try
        dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[2,8]}, 'InputFormat', dateFormats{col==[2,8]}); %#ok<AGROW>
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, 'UniformOutput', false);
            dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[2,8]}, 'InputFormat', dateFormats{col==[2,8]}); %#ok<AGROW>
        catch
            dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<AGROW>
        end
    end
    
    dateFormatIndex = dateFormatIndex + 1;
    blankDates{col} = dataArray{col} == '';
    anyBlankDates = blankDates{col} | anyBlankDates;
    invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
    anyInvalidDates = invalidDates{col} | anyInvalidDates;
end
dates = dates(:,[2,8]);
blankDates = blankDates(:,[2,8]);
invalidDates = invalidDates(:,[2,8]);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42]);
rawStringColumns = string(raw(:, [7,43,44,45]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay = table;
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Datenum = cell2mat(rawNumericColumns(:, 1));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Datetime_coarse = dates{:, 1};
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.no2_trop = cell2mat(rawNumericColumns(:, 2));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.no2_trop_err = cell2mat(rawNumericColumns(:, 3));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.d = cell2mat(rawNumericColumns(:, 4));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.no2 = cell2mat(rawNumericColumns(:, 5));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.UTC = rawStringColumns(:, 1);
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86 = dates{:, 2};
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column2Fractionaldayssince1Jan2000UTmidnightforcenterofmeasurem = cell2mat(rawNumericColumns(:, 6));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column3Effectivedurationofmeasurementinseconds = cell2mat(rawNumericColumns(:, 7));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column4Solarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 8));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column5Solarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 9));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column6Lunarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 10));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column7Lunarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 11));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column8NitrogendioxidetotalverticalcolumnamountDobsonUnits9e99r = cell2mat(rawNumericColumns(:, 12));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column9UncertaintyofnitrogendioxidetotalverticalcolumnamountDob = cell2mat(rawNumericColumns(:, 13));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column10Directnitrogendioxideairmassfactor = cell2mat(rawNumericColumns(:, 14));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column11L2dataqualityflagfornitrogendioxide0highquality1mediumq = cell2mat(rawNumericColumns(:, 15));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column12Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 16));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column13Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 17));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column14Fittingresultindex012noerrororwarning34warning4error = cell2mat(rawNumericColumns(:, 18));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column15Normalizedrmsofspectralfittingresidualsweightedwithmeas = cell2mat(rawNumericColumns(:, 19));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column16Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 20));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column17Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 21));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column18Climatologicalstationpressurembar = cell2mat(rawNumericColumns(:, 22));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column19Dataprocessingtypeindex = cell2mat(rawNumericColumns(:, 23));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column20Calibrationfileversionused = cell2mat(rawNumericColumns(:, 24));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column21Level2Fitdataqualityflag0highquality1mediumquality2lowq = cell2mat(rawNumericColumns(:, 25));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column22Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 26));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column23Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 27));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column24Level1dataqualityflag0highquality1mediumquality2lowqual = cell2mat(rawNumericColumns(:, 28));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column25Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 29));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column26Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 30));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column27WavelengtheffectivetemperatureC999noeffectivetemperatur = cell2mat(rawNumericColumns(:, 31));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column28Estimatedaverageresidualstraylightlevelonlyvalidforstra = cell2mat(rawNumericColumns(:, 32));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column29Retrievedwavelengthshiftfromlevel1datanm9nowavelengthch = cell2mat(rawNumericColumns(:, 33));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column30Retrievedwavelengthshiftfromspectralfittingnm9nowavelen = cell2mat(rawNumericColumns(:, 34));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column31Integrationtimems = cell2mat(rawNumericColumns(:, 35));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column32Numberofdarkcountcycles = cell2mat(rawNumericColumns(:, 36));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column33Effectivepositionoffilterwheel10filterwheelnotused19are = cell2mat(rawNumericColumns(:, 37));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column34Effectivepositionoffilterwheel20filterwheelnotused19are = cell2mat(rawNumericColumns(:, 38));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.instrument = cell2mat(rawNumericColumns(:, 39));
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.location = rawStringColumns(:, 2);
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.time = rawStringColumns(:, 3);
TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.LTC = rawStringColumns(:, 4);

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Datetime_coarse=datenum(TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Datetime_coarse);TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86=datenum(TropOMIOFFLFortMcKayavg7kmPandora122FortMcKay.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86);

function TropOMIOFFLFortMcKayclosestPandora122FortMcKay = importfile_P122_closest(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIOFFLFORTMCKAYCLOSESTPANDORA122FORTMCKAY = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   TROPOMIOFFLFORTMCKAYCLOSESTPANDORA122FORTMCKAY = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   TropOMIOFFLFortMcKayclosestPandora122FortMcKay = importfile('TropOMI_OFFL_FortMcKay_closest_Pandora122FortMcKay.csv', 2, 41);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/06/20 14:00:47

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,5,6,7,8,10,11,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

dateFormats = {'yyyy-MM-dd HH:mm:ss', 'yyyyMMdd''T''HHmmss''Z'''};
dateFormatIndex = 1;
blankDates = cell(1,size(raw,2));
anyBlankDates = false(size(raw,1),1);
invalidDates = cell(1,size(raw,2));
anyInvalidDates = false(size(raw,1),1);
for col=[12,14]% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
    try
        dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[12,14]}, 'InputFormat', dateFormats{col==[12,14]}); %#ok<AGROW>
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, 'UniformOutput', false);
            dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[12,14]}, 'InputFormat', dateFormats{col==[12,14]}); %#ok<AGROW>
        catch
            dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<AGROW>
        end
    end
    
    dateFormatIndex = dateFormatIndex + 1;
    blankDates{col} = dataArray{col} == '';
    anyBlankDates = blankDates{col} | anyBlankDates;
    invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
    anyInvalidDates = invalidDates{col} | anyInvalidDates;
end
dates = dates(:,[12,14]);
blankDates = blankDates(:,[12,14]);
invalidDates = invalidDates(:,[12,14]);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,2,3,5,6,7,8,10,11,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48]);
rawStringColumns = string(raw(:, [4,9,13,49,50,51]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
TropOMIOFFLFortMcKayclosestPandora122FortMcKay = table;
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.lon = cell2mat(rawNumericColumns(:, 1));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.lat = cell2mat(rawNumericColumns(:, 2));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.qa = cell2mat(rawNumericColumns(:, 3));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.utc_time = rawStringColumns(:, 1);
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.no2_trop = cell2mat(rawNumericColumns(:, 4));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.no2_trop_err = cell2mat(rawNumericColumns(:, 5));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.d = cell2mat(rawNumericColumns(:, 6));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.no2 = cell2mat(rawNumericColumns(:, 7));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.time_x = rawStringColumns(:, 2);
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Datenum = cell2mat(rawNumericColumns(:, 8));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.hour = cell2mat(rawNumericColumns(:, 9));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Datetime_coarse = dates{:, 1};
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.UTC = rawStringColumns(:, 3);
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86 = dates{:, 2};
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column2Fractionaldayssince1Jan2000UTmidnightforcenterofmeasurem = cell2mat(rawNumericColumns(:, 10));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column3Effectivedurationofmeasurementinseconds = cell2mat(rawNumericColumns(:, 11));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column4Solarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 12));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column5Solarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 13));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column6Lunarzenithangleforcenterofmeasurementindegree = cell2mat(rawNumericColumns(:, 14));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column7Lunarazimuthforcenterofmeasurementindegree0northincrease = cell2mat(rawNumericColumns(:, 15));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column8NitrogendioxidetotalverticalcolumnamountDobsonUnits9e99r = cell2mat(rawNumericColumns(:, 16));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column9UncertaintyofnitrogendioxidetotalverticalcolumnamountDob = cell2mat(rawNumericColumns(:, 17));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column10Directnitrogendioxideairmassfactor = cell2mat(rawNumericColumns(:, 18));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column11L2dataqualityflagfornitrogendioxide0highquality1mediumq = cell2mat(rawNumericColumns(:, 19));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column12Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 20));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column13Sumover2iusingthoseiforwhichthecorrespondingL2dataquali = cell2mat(rawNumericColumns(:, 21));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column14Fittingresultindex012noerrororwarning34warning4error = cell2mat(rawNumericColumns(:, 22));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column15Normalizedrmsofspectralfittingresidualsweightedwithmeas = cell2mat(rawNumericColumns(:, 23));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column16Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 24));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column17Expectednormalizedrmsofweightedspectralfittingresiduals = cell2mat(rawNumericColumns(:, 25));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column18Climatologicalstationpressurembar = cell2mat(rawNumericColumns(:, 26));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column19Dataprocessingtypeindex = cell2mat(rawNumericColumns(:, 27));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column20Calibrationfileversionused = cell2mat(rawNumericColumns(:, 28));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column21Level2Fitdataqualityflag0highquality1mediumquality2lowq = cell2mat(rawNumericColumns(:, 29));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column22Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 30));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column23Sumover2iusingthoseiforwhichthecorrespondingL2Fitdataqu = cell2mat(rawNumericColumns(:, 31));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column24Level1dataqualityflag0highquality1mediumquality2lowqual = cell2mat(rawNumericColumns(:, 32));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column25Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 33));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column26Sumover2iusingthoseiforwhichthecorrespondingL1dataquali = cell2mat(rawNumericColumns(:, 34));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column27WavelengtheffectivetemperatureC999noeffectivetemperatur = cell2mat(rawNumericColumns(:, 35));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column28Estimatedaverageresidualstraylightlevelonlyvalidforstra = cell2mat(rawNumericColumns(:, 36));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column29Retrievedwavelengthshiftfromlevel1datanm9nowavelengthch = cell2mat(rawNumericColumns(:, 37));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column30Retrievedwavelengthshiftfromspectralfittingnm9nowavelen = cell2mat(rawNumericColumns(:, 38));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column31Integrationtimems = cell2mat(rawNumericColumns(:, 39));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column32Numberofdarkcountcycles = cell2mat(rawNumericColumns(:, 40));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column33Effectivepositionoffilterwheel10filterwheelnotused19are = cell2mat(rawNumericColumns(:, 41));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column34Effectivepositionoffilterwheel20filterwheelnotused19are = cell2mat(rawNumericColumns(:, 42));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.instrument = cell2mat(rawNumericColumns(:, 43));
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.location = rawStringColumns(:, 4);
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.time_y = rawStringColumns(:, 5);
TropOMIOFFLFortMcKayclosestPandora122FortMcKay.LTC = rawStringColumns(:, 6);

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Datetime_coarse=datenum(TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Datetime_coarse);TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86=datenum(TropOMIOFFLFortMcKayclosestPandora122FortMcKay.Column1UTdateandtimeforcenterofmeasurementyyyymmddThhmmssZISO86);

function TropOMINRTIDownsviewclosestPandora103Downsviewsimple = importfile_merged_simple(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMINRTIDOWNSVIEWCLOSESTPANDORA103DOWNSVIEWSIMPLE =
%   IMPORTFILE(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%
%   TROPOMINRTIDOWNSVIEWCLOSESTPANDORA103DOWNSVIEWSIMPLE =
%   IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   TropOMINRTIDownsviewclosestPandora103Downsviewsimple = importfile('TropOMI_NRTI_Downsview_closest_Pandora103Downsview_simple.csv', 2, 35);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/06/22 10:36:04

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,4,5]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

dateFormats = {'yyyy-MM-dd HH:mm:ss', 'yyyyMMdd''T''HHmmss''Z'''};
dateFormatIndex = 1;
blankDates = cell(1,size(raw,2));
anyBlankDates = false(size(raw,1),1);
invalidDates = cell(1,size(raw,2));
anyInvalidDates = false(size(raw,1),1);
for col=[1,3]% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
    try
        dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[1,3]}, 'InputFormat', dateFormats{col==[1,3]}); %#ok<AGROW>
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, 'UniformOutput', false);
            dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[1,3]}, 'InputFormat', dateFormats{col==[1,3]}); %#ok<AGROW>
        catch
            dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<AGROW>
        end
    end
    
    dateFormatIndex = dateFormatIndex + 1;
    blankDates{col} = dataArray{col} == '';
    anyBlankDates = blankDates{col} | anyBlankDates;
    invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
    anyInvalidDates = invalidDates{col} | anyInvalidDates;
end
dates = dates(:,[1,3]);
blankDates = blankDates(:,[1,3]);
invalidDates = invalidDates(:,[1,3]);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [2,4,5]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
TropOMINRTIDownsviewclosestPandora103Downsviewsimple = table;
TropOMINRTIDownsviewclosestPandora103Downsviewsimple.TropOMI_datetime = dates{:, 1};
TropOMINRTIDownsviewclosestPandora103Downsviewsimple.TropOMI_NO2 = cell2mat(rawNumericColumns(:, 1));
TropOMINRTIDownsviewclosestPandora103Downsviewsimple.Pandora_datetime = dates{:, 2};
TropOMINRTIDownsviewclosestPandora103Downsviewsimple.Pandora_NO2 = cell2mat(rawNumericColumns(:, 2));
TropOMINRTIDownsviewclosestPandora103Downsviewsimple.d = cell2mat(rawNumericColumns(:, 3));

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMINRTIDownsviewclosestPandora103Downsviewsimple.TropOMI_datetime=datenum(TropOMINRTIDownsviewclosestPandora103Downsviewsimple.TropOMI_datetime);TropOMINRTIDownsviewclosestPandora103Downsviewsimple.Pandora_datetime=datenum(TropOMINRTIDownsviewclosestPandora103Downsviewsimple.Pandora_datetime);


%%
function TropOMIOFFLDownsviewclosestPandora104Downsviewsimple = importfile_simple_ERA(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TROPOMIOFFLDOWNSVIEWCLOSESTPANDORA104DOWNSVIEWSIMPLE =
%   IMPORTFILE(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%
%   TROPOMIOFFLDOWNSVIEWCLOSESTPANDORA104DOWNSVIEWSIMPLE =
%   IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   TropOMIOFFLDownsviewclosestPandora104Downsviewsimple = importfile('TropOMI_OFFL_Downsview_closest_Pandora104Downsview_simple.csv', 2, 154);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/12/04 14:53:21

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,4,5,6,7,8,9,10,11,12,13]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

dateFormats = {'yyyy-MM-dd HH:mm:ss', 'yyyyMMdd''T''HHmmss''Z'''};
dateFormatIndex = 1;
blankDates = cell(1,size(raw,2));
anyBlankDates = false(size(raw,1),1);
invalidDates = cell(1,size(raw,2));
anyInvalidDates = false(size(raw,1),1);
for col=[1,3]% Convert the contents of columns with dates to MATLAB datetimes using the specified date format.
    try
        dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[1,3]}, 'InputFormat', dateFormats{col==[1,3]}); %#ok<AGROW>
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{col} = cellfun(@(x) x(2:end-1), dataArray{col}, 'UniformOutput', false);
            dates{col} = datetime(dataArray{col}, 'Format', dateFormats{col==[1,3]}, 'InputFormat', dateFormats{col==[1,3]}); %#ok<AGROW>
        catch
            dates{col} = repmat(datetime([NaN NaN NaN]), size(dataArray{col})); %#ok<AGROW>
        end
    end
    
    dateFormatIndex = dateFormatIndex + 1;
    blankDates{col} = dataArray{col} == '';
    anyBlankDates = blankDates{col} | anyBlankDates;
    invalidDates{col} = isnan(dates{col}.Hour) - blankDates{col};
    anyInvalidDates = invalidDates{col} | anyInvalidDates;
end
dates = dates(:,[1,3]);
blankDates = blankDates(:,[1,3]);
invalidDates = invalidDates(:,[1,3]);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [2,4,5,6,7,8,9,10,11,12,13]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple = table;
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.TropOMI_datetime = dates{:, 1};
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.TropOMI_NO2 = cell2mat(rawNumericColumns(:, 1));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.Pandora_datetime = dates{:, 2};
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.Pandora_NO2 = cell2mat(rawNumericColumns(:, 2));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.d = cell2mat(rawNumericColumns(:, 3));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.u_1000hPa = cell2mat(rawNumericColumns(:, 4));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.v_1000hPa = cell2mat(rawNumericColumns(:, 5));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.u_900hPa = cell2mat(rawNumericColumns(:, 6));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.v_900hPa = cell2mat(rawNumericColumns(:, 7));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.u_800hPa = cell2mat(rawNumericColumns(:, 8));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.v_800hPa = cell2mat(rawNumericColumns(:, 9));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.u_700hPa = cell2mat(rawNumericColumns(:, 10));
TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.v_700hPa = cell2mat(rawNumericColumns(:, 11));

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.TropOMI_datetime=datenum(TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.TropOMI_datetime);TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.Pandora_datetime=datenum(TropOMIOFFLDownsviewclosestPandora104Downsviewsimple.Pandora_datetime);

