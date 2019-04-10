# -*- coding: utf-8 -*-
"""
Created on Tue May  1 14:56:42 2018

@author: xiaoy
"""
import netCDF4
import numpy as np
import pandas as pd

filename = 'C:\\Projects\\TropOMI\\data\\NO2\\tropomi_OilSands_NO2_OFFL_20180312.nc'
#filename = 'C:\\Projects\\TropOMI\\data\\NO2\\trop_overpass__NO2_____s20180308_e20180410_c20180430203837_Toronto.nc'
#filename = 'C:\\Projects\\TropOMI\\data\\NO2\\trop_overpass__NO2_____s20180308_e20180410_c20180430203856_Egbert.nc'
#filename = 'C:\\Projects\\TropOMI\\data\\NO2\\trop_overpass__NO2_____s20180308_e20180410_c20180430203843_McKay.nc'
output_path = 'C:\\Projects\\TropOMI\\data\\NO2_output\\'
site = 'Downsview'
user_lat=43.7810 # Downsview
user_lon=-79.4680

#site = 'Egbert'
#user_lat=44.2300 
#user_lon=-79.7800

#site = 'FortMcKay'
#user_lat=57.1836
#user_lon=-111.6400


dis1 = 24e3 #OMI pixel size 24 km
dis2 = 7e3 #TropOMI pixel size 7 km



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
abc = netCDF4.Dataset(filename, 'r')

days= abc.groups.keys() #these are the days of TROPomi data

data=["longitude","latitude","qa_value","nitrogendioxide_tropospheric_column"] 


i=1
for day in days:
        df = pd.DataFrame()
        df['lon']=abc.groups[day].groups["PRODUCT"].variables["longitude"][:].ravel()
        df['lat']=abc.groups[day].groups["PRODUCT"].variables["latitude"][:].ravel()
        df['qa']=abc.groups[day].groups["PRODUCT"].variables["qa_value"][:].ravel()
        df['no2_trop']=abc.groups[day].groups["PRODUCT"].variables["nitrogendioxide_tropospheric_column"][:].ravel()*6.02214e+19
        df['no2_trop_err']=abc.groups[day].groups["PRODUCT"].variables['nitrogendioxide_tropospheric_column_precision'][:].ravel()*6.02214e+19
        df['d'] = get_distance(user_lat,user_lon,df['lat'],df['lon'])
        df['time'] =day
        
        df['no2']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["DETAILED_RESULTS"].variables["nitrogendioxide_summed_total_column"][:].ravel()*6.02214e+19
        #df['snow_ice_flag']=abc.groups[day].groups["PRODUCT"].groups["SUPPORT_DATA"].groups["INPUT_DATA"].variables["snow_ice_flag"][:].ravel()
        
        if i==1:
            df_output = df
            i+=1
        else:
            df_output = pd.concat([df_output,df])
# quality control
df_output = df_output[df_output.qa>0.5]
#df_output = df_output[df_output.snow_ice_flag==0]
# save to csv
df_output.to_csv(output_path + 'TropOMI_at_' + site + '.csv',index=False)


df_dis1 = df_output.copy()
df_dis1 = df_dis1[df_dis1.d<=dis1]
df_daily_dis1 = pd.DataFrame()
df_daily_dis1['no2_trop'] = df_dis1.groupby(['time'])['no2_trop'].mean().copy()
df_daily_dis1['no2_trop_err'] = df_dis1.groupby(['time'])['no2_trop_err'].mean().copy()
df_daily_dis1['d'] = df_dis1.groupby(['time'])['d'].mean().copy()
df_daily_dis1['no2'] = df_dis1.groupby(['time'])['no2'].mean().copy()
df_daily_dis1.reset_index(inplace=True)
# save to csv
df_daily_dis1.to_csv(output_path + 'TropOMI_at_' + site + '_daily_dis1.csv',index=False)

df_dis2 = df_output.copy()
df_dis2 = df_dis2[df_dis2.d<=dis2]
df_daily_dis2 = pd.DataFrame()
df_daily_dis2['no2_trop'] = df_dis2.groupby(['time'])['no2_trop'].mean().copy()
df_daily_dis2['no2_trop_err'] = df_dis2.groupby(['time'])['no2_trop_err'].mean().copy()
df_daily_dis2['d'] = df_dis2.groupby(['time'])['d'].mean().copy()
df_daily_dis2['no2'] = df_dis2.groupby(['time'])['no2'].mean().copy()
df_daily_dis2.reset_index(inplace=True)
# save to csv
df_daily_dis2.to_csv(output_path + 'TropOMI_at_' + site + '_daily_dis2.csv',index=False)

            

        