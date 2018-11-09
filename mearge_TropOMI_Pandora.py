# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:23:02 2018

@author: ZhaoX
"""

import pandas as pd

Pandora_no = '145'
#site = 'Downsview'
#site = 'Egbert'
site = 'StGeorge'
TropOMI_type = 'OFFL'
average_dis = '24' # TropOMI averaged distance in km

#Pandora_shelve_filename = '\\\\wdow05dtmibroh.ncr.int.ec.gc.ca\\GDrive\\Pandora\\103\\Blick\\L2\\Blick_L2.out'
#Pandora_shelve_filename = '\\\\wdow05dtmibroh.ncr.int.ec.gc.ca\\GDrive\\Pandora\\122\\Blick\\L2\\Blick_L2.out'
Pandora_shelve_filename = '\\\\wdow05dtmibroh.ncr.int.ec.gc.ca\\GDrive\\Pandora\\' + Pandora_no + '\\Blick\\L2\\Blick_L2.out'
#Pandora_shelve_filename = '\\\\wdow05dtmibroh.ncr.int.ec.gc.ca\\GDrive\\Pandora\\109\\Blick\\L2\\Blick_L2.out'
TropOMI_shelve_filepath = 'C:\\Projects\\TropOMI\\data\\NO2_output\\' + TropOMI_type + '\\' + site+ '\\'
#TropOMI_shelve_filepath = 'C:\\Projects\\TropOMI\\data\\NO2_output\\NRTI\\Downsview\\'
#TropOMI_shelve_filepath = 'C:\\Projects\\TropOMI\\data\\NO2_output\\OFFL\\FortMcKay\\'
#TropOMI_shelve_filepath = 'C:\\Projects\\TropOMI\\data\\NO2_output\\NRTI\\FortMcKay\\'
#TropOMI_shelve_filepath = 'C:\\Projects\\TropOMI\\data\\NO2_output\\NRTI\\StGeorge\\'

#TropOMI_shelve_filename = 'TropOMI_OFFL_Downsview_avg7km.out'
#TropOMI_shelve_filename = 'TropOMI_OFFL_Downsview_avg24km.out'
#TropOMI_shelve_filename = 'TropOMI_NRTI_Downsview_avg7km.out'
TropOMI_shelve_filename = 'TropOMI_' + TropOMI_type + '_' + site + '_avg' + average_dis + 'km.out'
#TropOMI_shelve_filename = 'TropOMI_OFFL_FortMcKay_avg7km.out'
#TropOMI_shelve_filename = 'TropOMI_NRTI_FortMcKay_avg7km.out'
#TropOMI_shelve_filename = 'TropOMI_NRTI_StGeorge_avg7km.out'
CSV_output_filepath = 'C:\\Projects\\TropOMI\\data\\NO2_merged_with_Pandora\\'

def open_shelf(shelve_filename):
    import shelve
    my_shelf = shelve.open(shelve_filename)
    
    for key in my_shelf:
        globals()[key]=my_shelf[key]
    my_shelf.close()  
    


# load Pandora dataframe
open_shelf(Pandora_shelve_filename)
# load TropOMI dataframe
open_shelf(TropOMI_shelve_filepath + TropOMI_shelve_filename)
exec('df_Pandora = Pandora' + Pandora_no + 's1_' + site + '_L2Tot_rnvs0p1.copy()') # only read NO2 data
#df_Pandora = Pandora104s1_Downsview_L2Tot_rnvs0p1.copy()
#df_Pandora = Pandora103s1_Downsview_L2Tot_rnvs0p1.copy()
#df_Pandora = Pandora122s1_FortMcKay_L2Tot_rnvs0p1.copy()
#df_Pandora = Pandora104s1_Downsview_L2Tot_rnvs0p1.copy()
#df_Pandora = Pandora109s1_StGeorge_L2Tot_rnvs0p1.copy()
df_TropOMI_dis1 = df_concat_dis1.copy()
df_TropOMI_closest = df_concat_closest.copy()


#%% BlickP1.5.7
# Pandora data filteration for P108 
df_Pandora = df_Pandora[df_Pandora['Column 10: Direct nitrogen dioxide air mass factor']<5]
df_Pandora = df_Pandora[(df_Pandora['Column 12: L2 data quality flag for nitrogen dioxide: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'] == 0) | (df_Pandora['Column 12: L2 data quality flag for nitrogen dioxide: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'] == 10)]
df_Pandora = df_Pandora[df_Pandora['Column 33: Integration time [ms]']<10000]
df_Pandora = df_Pandora[df_Pandora['Column 9: Uncertainty of nitrogen dioxide total vertical column amount [Dobson Units] based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no uncertainty could be retrieved']<0.2]
df_Pandora = df_Pandora[df_Pandora['Column 16: Normalized rms of spectral fitting residuals weighted with measured uncertainty, -9=fitting not successful or no uncertainty given']<0.05]

#%% BlickP1.5.2
# Pandora data filteration for P103 
#df_Pandora = df_Pandora[df_Pandora['Column 10: Direct nitrogen dioxide air mass factor']<5]
#df_Pandora = df_Pandora[df_Pandora['Column 12: L2 data quality flag for nitrogen dioxide: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality']<1]
#df_Pandora = df_Pandora[df_Pandora['Column 33: Integration time [ms]']<10000]
#df_Pandora = df_Pandora[df_Pandora['Column 9: Uncertainty of nitrogen dioxide total vertical column amount [Dobson Units] based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no uncertainty could be retrieved']<0.2]
#df_Pandora = df_Pandora[df_Pandora['Column 16: Normalized rms of spectral fitting residuals weighted with measured uncertainty, -9=fitting not successful or no uncertainty given']<0.05]

#%% BlickP 1.4
# Pandora data filteration for P122
#df_Pandora = df_Pandora[df_Pandora['Column 10: Direct nitrogen dioxide air mass factor']<5]
#df_Pandora = df_Pandora[df_Pandora['Column 11: L2 data quality flag for nitrogen dioxide: 0=high quality, 1=medium quality, 2=low quality']<1]
#df_Pandora = df_Pandora[df_Pandora['Column 31: Integration time [ms]']<10000]
#df_Pandora = df_Pandora[df_Pandora['Column 9: Uncertainty of nitrogen dioxide total vertical column amount [Dobson Units] based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no uncertainty could be retrieved']<0.2]
#df_Pandora = df_Pandora[df_Pandora['Column 15: Normalized rms of spectral fitting residuals weighted with measured uncertainty, -9=fitting not successful or no uncertainty given']<0.05]


#%%  merge and save averaged TropOMI and Pandora
df_TropOMI_dis1.set_index('Datetime',inplace = True)
df_TropOMI_dis1.index.tz_localize('utc')
df_TropOMI_dis1['UTC'] = df_TropOMI_dis1.index.tz_localize('utc')
df_TropOMI_dis1.sort_values(by=['UTC'],inplace=True)
x = pd.merge_asof(df_TropOMI_dis1 , df_Pandora, left_on = 'UTC', right_on='UTC',tolerance=pd.Timedelta('10min'))

instrument_names = pd.unique(x.instrument)
instrument_locations = pd.unique(x.location)
try:
    CSV_output_filename = TropOMI_shelve_filename[0:-4]+'_' + instrument_names[0]+ instrument_locations[0] + '.csv'
    instrument_name = instrument_names[0]
    instrument_location = instrument_locations[0]
except:
    CSV_output_filename = TropOMI_shelve_filename[0:-4]+'_' + instrument_names[1]+ instrument_locations[1] + '.csv'
    instrument_name = instrument_names[1]
    instrument_location = instrument_locations[1]
    
x.to_csv(CSV_output_filepath+CSV_output_filename,index=False)# general output

x_simple = pd.DataFrame()
x_simple['TropOMI_datetime'] = x.Datetime_coarse.copy()
x_simple['TropOMI_NO2'] = x.no2.copy()
x_simple['Pandora_datetime'] = x['Column 1: UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'].copy()
x_simple['Pandora_NO2'] = x['Column 8: Nitrogen dioxide total vertical column amount [Dobson Units], -9e99=retrieval not successful'].copy()
x_simple['d'] = x.d.copy()
x_simple.to_csv(CSV_output_filepath+CSV_output_filename[0:-4] + '_simple.csv',index=False)# output for plotting tools (used in matlab)


#%%  merge and save averaged TropOMI and Pandora
df_TropOMI_closest.set_index('datetime',inplace = True)
df_TropOMI_closest.index.tz_localize('utc')
df_TropOMI_closest['UTC'] = df_TropOMI_closest.index.tz_localize('utc')
df_TropOMI_closest.sort_values(by=['UTC'],inplace=True)
x = pd.merge_asof(df_TropOMI_closest , df_Pandora, left_on = 'UTC', right_on='UTC',tolerance=pd.Timedelta('10min'))

CSV_output_filename = TropOMI_shelve_filename[0:TropOMI_shelve_filename.find('avg')] +'closest'+'_' +instrument_name + instrument_location + '.csv'
x.to_csv(CSV_output_filepath+CSV_output_filename,index=False)

x_simple = pd.DataFrame()
x_simple['TropOMI_datetime'] = x.Datetime_coarse.copy()
x_simple['TropOMI_NO2'] = x.no2.copy()
x_simple['Pandora_datetime'] = x['Column 1: UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'].copy()
x_simple['Pandora_NO2'] = x['Column 8: Nitrogen dioxide total vertical column amount [Dobson Units], -9e99=retrieval not successful'].copy()
x_simple['d'] = x.d.copy()
x_simple.to_csv(CSV_output_filepath+CSV_output_filename[0:-4] + '_simple.csv',index=False)# output for plotting tools (used in matlab)


#%% save filtered Pandora
CSV_output_filename = instrument_name + instrument_location + '.csv'
df_Pandora.to_csv(CSV_output_filepath+CSV_output_filename,index=False)
