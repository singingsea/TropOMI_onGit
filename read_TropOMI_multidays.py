# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 12:05:02 2018

@author: ZhaoX
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May  1 14:56:42 2018

@author: xiaoy
"""
import netCDF4
import numpy as np
import pandas as pd
import datetime

ERA_merged = True # if the TropOMI data has already get ERA wind data merged
#site = 'Downsview'
#site = 'StGeorge'
#site = 'Egbert'
site = 'FortMcKay'
#site = 'FortMcKay_suncrude'
TropOMI_datatype = 'OFFL'
#TropOMI_datatype = 'NRTI'
dis1 = 24e3
#dis1 = 7e3

if site == 'Downsview':
    user_lat=43.7810 # Downsview
    user_lon=-79.4680
    raw_data_area = 'Toronto'# this is the general area of TropOMI truncated data, note that we only have two areas now, one is Toronto (cover for Downsview, StGeogre, and Egbert) and the other one is OilSands (only cover FortMcKay)
elif site == 'StGeorge':
    user_lat=43.6605
    user_lon=-79.3986
    raw_data_area = 'Toronto'
elif site == 'Egbert':
    user_lat=44.2300 
    user_lon=-79.7800
    raw_data_area = 'Toronto'
elif site == 'FortMcKay':
    user_lat=57.1836
    user_lon=-111.6400
    raw_data_area = 'OilSands'
elif site == 'FortMcKay_suncrude':
    user_lat=57.0333
    user_lon=-111.6167
    raw_data_area = 'OilSands'
    
if ERA_merged == True:
    filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\' +  TropOMI_datatype + '\\' + raw_data_area + '_ERA\\'
    #filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\OFFL\\'
    output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\' + TropOMI_datatype + '\\' + site + '_ERA\\' # outp file path
else:
    filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\' +  TropOMI_datatype + '\\' + raw_data_area + '\\'
    output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\' + TropOMI_datatype + '\\' + site + '\\' # outp file path
    

#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\OFFL\\test1\\'    
#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\OFFL\\test\\'


#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\OFFL\\Egbert\\' # outp file path
#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\OFFL\\Downsview\\' # outp file path
#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\NRTI\\Downsview\\' # outp file path
#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\OFFL\\FortMcKay\\' # outp file path
#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\NRTI\\FortMcKay\\' # outp file path
#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\OFFL\\StGeorge\\' # outp file path
#output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\NRTI\\StGeorge\\' # outp file path

#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\OFFL\\Toronto\\'    # input file path
#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\NRTI\\Toronto\\'    # input file path
#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\OFFL\\OilSands\\'    # input file path
#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\NRTI\\OilSands\\'    # input file path

#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\NRTI\\'    # input file path
#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\OFFL\\Toronto\\'
#filepath = 'C:\\Projects\\TropOMI\\data\\NO2\\OFFL\\'


#site = 'Downsview'
#user_lat=43.7810 # Downsview
#user_lon=-79.4680

#site = 'StGeorge'
#user_lat=43.6605
#user_lon=-79.3986

#site = 'Egbert'
#user_lat=44.2300 
#user_lon=-79.7800

#site = 'FortMcKay'
#user_lat=57.1836
#user_lon=-111.6400

shelve_filename = output_path + 'TropOMI_' +TropOMI_datatype +'_'+ site + '_avg' + str(int(dis1/1000)) +'km.out'



#%%
def get_distance(user_lat,user_lon,lat,lon): 
    R=6371000#radius of the earth in meters
    lat1=np.radians(user_lat)
    lat2=np.radians(lat)
    delta_lat=np.radians(lat-user_lat)
    delta_lon=np.radians(lon-user_lon)
    a=(np.sin(delta_lat/2))*(np.sin(delta_lat/2))+(np.cos(lat1))*(np.cos(lat2))*(np.sin(delta_lon/2))*(np.sin(delta_lon/2))
    c=2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
    d=R*c
    return d

#%%
def get_TropOMI_single_day_v2(filename,output_path,site,user_lat,user_lon,dis1):
    abc = netCDF4.Dataset(filename, 'r')
    df_output = pd.DataFrame()
    days= abc.groups.keys() #these are the days of TROPomi data
    if len(days) > 0:
       
        i=1
        for day in days:
                df = pd.DataFrame()
                # read in variables in "PRODUCT"
                df['lon']=abc.groups[day].groups["PRODUCT"].variables["longitude"][:].ravel()
                df['lat']=abc.groups[day].groups["PRODUCT"].variables["latitude"][:].ravel()
                
                df['ground_pixel']=abc.groups[day].groups["PRODUCT"].variables["ground_pixel"][:].ravel() #'across-track dimension index' 'This coordinate variable defines the indices across track, from west to east; index starts at 0'
                df['scanline']=abc.groups[day].groups["PRODUCT"].variables["scanline"][:].ravel() # 'along-track dimension index' 'This coordinate variable defines the indices along track; index starts at 0'
                df['qa']=abc.groups[day].groups["PRODUCT"].variables["qa_value"][:].ravel() #'A continuous quality descriptor, varying between 0 (no data) and 1 (full quality data). Recommend to ignore data with qa_value < 0.5'
#                df['utc_time']=abc.groups[day].groups["PRODUCT"].variables["utc_time"][:].ravel()
                df['time']=abc.groups[day].groups["PRODUCT"].variables["time"][:].ravel()
                df['delta_time']=abc.groups[day].groups["PRODUCT"].variables["delta_time"][:].ravel()
                df['total_time']=abc.groups[day].groups["PRODUCT"].variables["total_time"][:].ravel() # this is generated by Debora! it is the time in seconds since 2010-01-01
                df['no2_trop']=abc.groups[day].groups["PRODUCT"].variables["nitrogendioxide_tropospheric_column"][:].ravel()*6.02214e+19
                df['no2_trop_err']=abc.groups[day].groups["PRODUCT"].variables["nitrogendioxide_tropospheric_column_precision"][:].ravel()*6.02214e+19
                df['amf_trop']=abc.groups[day].groups["PRODUCT"].variables["air_mass_factor_troposphere"][:].ravel()
                
                df['amf']=abc.groups[day].groups["PRODUCT"].variables["air_mass_factor_total"][:].ravel()
                df['snow_cover']=abc.groups[day].groups["PRODUCT"].variables["snow_cover"][:].ravel() # '1-sea, 2-land, 3-sea ice, 4-snow, 0-below equator'
                df['albedo_modis_nosnow']=abc.groups[day].groups["PRODUCT"].variables["modis_nosnow_albedo"][:].ravel() # 'MODIS Terra and Aqua BRDF/Albedo'
                df['albedo_modis_snow']=abc.groups[day].groups["PRODUCT"].variables["modis_snow_albedo"][:].ravel() # 'MODIS Terra and Aqua BRDF/Albedo'
                df['ECCC_NO2']=abc.groups[day].groups["PRODUCT"].variables["ECCC_NO2"][:].ravel()*6.02214e+19
                df['ECCC_AMF']=abc.groups[day].groups["PRODUCT"].variables["ECCC_AMF"][:].ravel()
                
                df['d'] = get_distance(user_lat,user_lon,df['lat'],df['lon'])
                
                if ERA_merged == True:
                    try:
                        # read in all ERA pressure levels
                        ERA_p_level=abc.groups[day].groups["PRODUCT"].variables["level"][:].ravel()
                        p_level_i = 0
                        for p_level in ERA_p_level:
                            # assign u wind from xxx hPa pressure to df
                            exec('df["u_' + str(int(p_level)) + 'hPa"] = abc.groups[day].groups["PRODUCT"].variables["u"][' + str(p_level_i) + ',:].ravel()' )
                            # assign v wind from xxx hPa pressure to df
                            exec('df["v_' + str(int(p_level)) + 'hPa"] = abc.groups[day].groups["PRODUCT"].variables["v"][' + str(p_level_i) + ',:].ravel()' )
                            p_level_i +=1
                    except:
                        print('Failed to read data from a file: ' + day)
                                

                
                # read in variables in "PRODUCT"/"SUPPORT_DATA"/"DETAILED_RESULTS"
                df['no2']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["DETAILED_RESULTS"].variables["nitrogendioxide_total_column"][:].ravel()*6.02214e+19
                df['no2_err']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["DETAILED_RESULTS"].variables["nitrogendioxide_total_column_precision"][:].ravel()*6.02214e+19
                df['cloud_frac']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["DETAILED_RESULTS"].variables["cloud_fraction_crb_nitrogendioxide_window"][:].ravel()# 'effective_cloud_area_fraction_assuming_fixed_cloud_albedo' 'cloud fraction at 439 nm for NO2 retrieval'
                df['no2_strat']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["DETAILED_RESULTS"].variables["nitrogendioxide_stratospheric_column"][:].ravel()*6.0221e+19 # 'stratospheric vertical column of nitrogen dioxide, derived from the TM5-MP vertical profiles'
                df['no2_strat_err']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["DETAILED_RESULTS"].variables["nitrogendioxide_stratospheric_column_precision"][:].ravel()*6.0221e+19
                df['amf_strat']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["DETAILED_RESULTS"].variables["air_mass_factor_stratosphere"][:].ravel()
                
                df['snow_ice_flag']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["INPUT_DATA"].variables["snow_ice_flag"][:].ravel() # 'snow-free_land sea_ice_1_percent sea_ice_2_percent sea_ice_3_percent sea_ice_4_percent sea_ice_5_percent sea_ice_6_percent sea_ice_7_percent sea_ice_8_percent sea_ice_9_percent sea_ice_10_percent sea_ice_11_percent sea_ice_12_percent sea_ice_13_percent sea_ice_14_percent sea_ice_15_percent sea_ice_16_percent sea_ice_17_percent sea_ice_18_percent sea_ice_19_percent sea_ice_20_percent sea_ice_21_percent sea_ice_22_percent sea_ice_23_percent sea_ice_24_percent sea_ice_25_percent sea_ice_26_percent sea_ice_27_percent sea_ice_28_percent sea_ice_29_percent sea_ice_30_percent sea_ice_31_percent sea_ice_32_percent sea_ice_33_percent sea_ice_34_percent sea_ice_35_percent sea_ice_36_percent sea_ice_37_percent sea_ice_38_percent sea_ice_39_percent sea_ice_40_percent sea_ice_41_percent sea_ice_42_percent sea_ice_43_percent sea_ice_44_percent sea_ice_45_percent sea_ice_46_percent sea_ice_47_percent sea_ice_48_percent sea_ice_49_percent sea_ice_50_percent sea_ice_51_percent sea_ice_52_percent sea_ice_53_percent sea_ice_54_percent sea_ice_55_percent sea_ice_56_percent sea_ice_57_percent sea_ice_58_percent sea_ice_59_percent sea_ice_60_percent sea_ice_61_percent sea_ice_62_percent sea_ice_63_percent sea_ice_64_percent sea_ice_65_percent sea_ice_66_percent sea_ice_67_percent sea_ice_68_percent sea_ice_69_percent sea_ice_70_percent sea_ice_71_percent sea_ice_72_percent sea_ice_73_percent sea_ice_74_percent sea_ice_75_percent sea_ice_76_percent sea_ice_77_percent sea_ice_78_percent sea_ice_79_percent sea_ice_80_percent sea_ice_81_percent sea_ice_82_percent sea_ice_83_percent sea_ice_84_percent sea_ice_85_percent sea_ice_86_percent sea_ice_87_percent sea_ice_88_percent sea_ice_89_percent sea_ice_90_percent sea_ice_91_percent sea_ice_92_percent sea_ice_93_percent sea_ice_94_percent sea_ice_95_percent sea_ice_96_percent sea_ice_97_percent sea_ice_98_percent sea_ice_99_percent sea_ice_100_percent permanent_ice snow mixed_pixels_at_coastlines suspect_ice_value corners ocean'
                df['albedo_no2']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["INPUT_DATA"].variables["surface_albedo_nitrogendioxide_window"][:].ravel() # 'Surface albedo in the NO2 fit window'
                df['albedo_cloud']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["INPUT_DATA"].variables["surface_albedo"][:].ravel() # 'Surface albedo in the cloud product'
                try:
                    # Debora's code has some issue, so the "time" is not right!!! stop use it
                    #df['time_seconds'] = df['time'] + df['delta_time']/1000
                    #df['utc_time']=pd.to_datetime(df['time_seconds'],unit='s',origin=pd.Timestamp('2010-01-01'))                    
                    #df['Datenum']=pd.to_numeric(df['utc_time'].values,errors='coerce')
                    
                    #df['time'] = df['utc_time'].apply(lambda x: datetime.datetime.strptime(x,'%Y-%m-%dT%H-%M-%S.%fZ'))
                    
                    df['time_seconds'] = df['time'] + df['delta_time']/1000            
                    #df['utc_time'] =  pd.to_datetime(df['delta_time'], unit='ms',origin=pd.Timestamp(day))
                    df['utc_time'] =  pd.to_datetime(df['delta_time'], unit='ms',origin=pd.Timestamp(day),utc=True)
                    df['Datenum']=pd.to_numeric(df['utc_time'].values,errors='coerce')
                except:
                    print(['failed in time convertion for ' + day])
                    #df['time'] =day
                    #df['Datenum'] = 0
                if i==1:
                    df_output = df
                    i+=1
                else:
                    df_output = pd.concat([df_output,df])
    else:
        day = 9999    
            
    if len(df_output) > 0:
        # quality control
        df_output = df_output[df_output.qa>0.5].copy()
        df_output = df_output[df['cloud_frac']<0.3].copy()
        #df_output = df_output[df_output.snow_ice_flag==0]
        if len(df_output) > 0:
            # save to csv
            df_output.to_csv(output_path + 'TropOMI_at_' + site + day + '.csv',index=False)
            
    return df_output, day
#%%
def average_within_dis1(df,site,day):
    # make average within dis1
    df_dis1 = df.copy()
    df_dis1['datetime'] = pd.to_datetime(df_dis1.time)
    df_dis1['hour'] = df_dis1.datetime.dt.hour
    df_dis1 = df_dis1[df_dis1.d<=dis1].copy()
    if len(df_dis1) > 0:
        for hour in pd.unique(df_dis1.hour):
            df_daily_dis1 = pd.DataFrame()
            df_daily_dis1['Datenum'] = [df_dis1[df_dis1.hour == hour]['Datenum'].mean()]
            df_daily_dis1['Datetime'] = pd.to_datetime(df_daily_dis1['Datenum'])
            try:
                df_daily_dis1['Datetime_coarse'] = df_daily_dis1['Datetime'].values.astype('datetime64[s]')
            except:
                df_daily_dis1['Datetime_coarse'] = df_daily_dis1['Datetime']    
            df_daily_dis1['no2_trop'] = df_dis1[df_dis1.hour == hour]['no2_trop'].mean()
            df_daily_dis1['no2_trop_err'] = df_dis1[df_dis1.hour == hour]['no2_trop_err'].mean()
            df_daily_dis1['d'] = df_dis1[df_dis1.hour == hour]['d'].mean()
            df_daily_dis1['no2'] = df_dis1[df_dis1.hour == hour]['no2'].mean()
            df_daily_dis1.reset_index(inplace=True)
            # save to csv
            df_daily_dis1.to_csv(output_path + 'TropOMI_at_' + site + day + '_hour'+str(int(hour)) + '_'+ str(int(dis1/1000)) + 'km.csv',index=False)
            # return data
            return df_daily_dis1

#%%
def get_closest_pixel(df,site,day):            
    # make closest pixel
    df_closest = df.copy()
    #df_closest['datetime'] = pd.to_datetime(df_closest.time)
    df_closest['datetime'] = df_closest['utc_time']
    df_closest['hour'] = df_closest.datetime.dt.hour
    
    if len(df_closest) > 0:
        for hour in pd.unique(df_closest.hour):
            df_hourly_closest = df_closest[df_closest.hour == hour].copy()
            min_dis = min(df_hourly_closest.d)
            
            df_hourly_closest_output = df_hourly_closest[df_hourly_closest.d == min_dis]
            #df_hourly_closest_output['Datetime'] = pd.to_datetime(df_hourly_closest_output['Datenum'])
            df_hourly_closest_output['Datetime_coarse'] = df_hourly_closest_output['datetime'].values.astype('datetime64[s]')
  

            df_hourly_closest_output.reset_index(inplace=True)
            df_hourly_closest_output.to_csv(output_path + 'TropOMI_at_' + site + day + '_hour'+str(int(hour)) +'_closest.csv',index=False)
    
    return df_hourly_closest_output               

          
#%%
def load_files(mypath):
    from os import listdir
    from os.path import isfile, join
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    onlyncfiles = []
    for f in onlyfiles:
        if f.find('.nc') != -1:
            onlyncfiles.append(f)
    return onlyncfiles

#%%

onlyncfiles = load_files(filepath)
first_valid_df_dis1 = True
first_valid_df_closest = True
for file_nm in onlyncfiles:
    # read in data as dataframe
    print('\n \n ')
    print('Read in data from file : \n' + filepath + file_nm)
    #df_daily,day = get_TropOMI_single_day(filepath + file_nm,output_path,site,user_lat,user_lon,dis1)
    df_daily,day = get_TropOMI_single_day_v2(filepath + file_nm,output_path,site,user_lat,user_lon,dis1) # this one included more variables from NC file,and also use cloud frac < 0.3 as filter
    df_daily_dis1 = pd.DataFrame()
    df_hourly_closest = pd.DataFrame()
    if len(df_daily) > 0:
        df_daily_dis1 = average_within_dis1(df_daily,site,day)
        df_hourly_closest = get_closest_pixel(df_daily,site,day)        
        
    if (df_daily_dis1 is None) == False:    
        if  (first_valid_df_dis1 == True):
            df_concat_dis1 = df_daily_dis1
            first_valid_df_dis1 = False
        else:
            df_concat_dis1 = pd.concat([df_concat_dis1,df_daily_dis1])
    

    if (df_hourly_closest is None) == False:    
        if  (first_valid_df_closest == True):
            df_concat_closest = df_hourly_closest
            first_valid_df_closest = False
        else:
            df_concat_closest = pd.concat([df_concat_closest,df_hourly_closest])
            
df_concat_dis1.drop('index',axis=1,inplace = True)            
df_concat_dis1 = df_concat_dis1[df_concat_dis1.Datenum != 0]# remove the ones we can't interp the timestamp! 
df_concat_dis1.sort_values(by='Datenum',inplace=True)
df_concat_dis1.to_csv(output_path + 'TropOMI_at_' + site + '_summary_dis' + str(dis1/1000) +'km.csv',index=False)

df_concat_closest.drop('index',axis=1,inplace = True) 
df_concat_closest = df_concat_closest[df_concat_closest.Datenum != 0]# remove the ones we can't interp the timestamp! 
df_concat_closest.sort_values(by='Datenum',inplace=True)
df_concat_closest.to_csv(output_path + 'TropOMI_at_' + site + '_summary_closest.csv',index=False)


# save data
print('\nSave dataframes:')
import shelve    
my_shelf = shelve.open(shelve_filename,'n') # 'n' for new
for key in dir():     
    if key.find('df_concat') != -1:
        try:
            my_shelf[key] = globals()[key]
            print(key + ' saved as dataframe! ')
        except TypeError:
            print('ERROR shelving: {0}'.format(key))
my_shelf.close()

# reload data
my_shelf = shelve.open(shelve_filename)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

