function TropOMI = reformat_Chris_TropOMI_data()
% simple function that can reformat Chris' TropOMI data to a similar format
% can be processed by transport code

% load Chris' data
%load('C:\Projects\TropOMI\data\NO2\OFFL\TropOMIknmidata_no2_GlobalSO2_011_600km_ECCC.mat');
load('C:\Projects\OMI\from_Chris\OMI_NO2_NA_cf0.5_2018_v3.1.mat');

TropOMI = table;

TropOMI.lon = data.lon';
TropOMI.lat = data.lat';
%TropOMI.qa = data.qual';
TropOMI.qa = repmat(1,height(TropOMI),1);

TropOMI.no2_trop = data.vcd0';
TropOMI.ECCC_NO2 = data.vcd';
TropOMI.no2_strat = data.vcds';
TropOMI.vza = data.vza';
TropOMI.cldfrac = data.cldfrac';

%TropOMI.snow_cover = data.snow';

TropOMI.u_1000hPa = data.u(1,:)';
TropOMI.v_1000hPa = data.v(1,:)';

TropOMI.u_950hPa = data.u(2,:)';	
TropOMI.v_950hPa = data.v(2,:)';	

TropOMI.u_900hPa = data.u(3,:)';	
TropOMI.v_900hPa = data.v(3,:)';	

TropOMI.u_850hPa = data.u(4,:)';	
TropOMI.v_850hPa = data.v(4,:)';	

TropOMI.u_800hPa = data.u(5,:)';	
TropOMI.v_800hPa = data.v(5,:)';	

TropOMI.u_750hPa = data.u(6,:)';	
TropOMI.v_750hPa = data.v(6,:)';	

% TropOMI.u_700hPa = NaN(numel(data.lon),1);		
% TropOMI.v_700hPa = NaN(numel(data.lon),1);		
% 
% TropOMI.u_650hPa = NaN(numel(data.lon),1);		
% TropOMI.v_650hPa = NaN(numel(data.lon),1);	

TropOMI.mdj = data.mjd';
[year,month,day,hour,minute,secs,ticks] = mjd2utc(TropOMI.mdj);
TropOMI.utc_time = datetime(year,month,day,hour,minute,secs,ticks);


tf = (TropOMI.lat >= 40) & (TropOMI.lat <= 50) & (TropOMI.lon >= -90) & (TropOMI.lon <= -70);
OMI = TropOMI(tf,:);
