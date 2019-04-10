# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 11:59:59 2018

@author: ZhaoX
"""
import netCDF4
import numpy as np
import pandas as pd

file = 'C:\\Projects\\TropOMI\\data\\NO2\\tropomi_OilSands_NO2_OFFL_20180311.nc'

abc = netCDF4.Dataset(file, 'r')
trop=abc.groups['PRODUCT'].variables.keys()
trop2=abc.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables.keys()

datas=pd.DataFrame(abc.groups['PRODUCT'].variables[trop[0]][:].ravel(),columns=[trop[0]])

for key in trop[1:]:
    print(key)
    datas[key]=abc.groups['PRODUCT'].variables[key][:].ravel()
for key in trop2:
    print(key)
    datas[key]=abc.groups['PRODUCT']. groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables[key][:].ravel()

#The time format is a bit tricky… I use this to deal with the time:
datas['time'] = datas['utc_time'].apply(lambda x: 
                                datetime.datetime.strptime(x,'%Y-%m-%dT%H-%M-%S.%fZ'))

datas['Datenum']=pd.to_numeric(datas['time'].values,errors='coerce')
