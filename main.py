# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 17:27:52 2018

@author: Daixin Zhao

"""

import ftplib # ftp utils
import os
from tqdm import tqdm # add processing bar
import numpy as np
import pandas as pd
from netCDF4 import Dataset # netcdf data processing
import utm # Bidirectional UTM-WGS84 converter
import scipy as sp
import scipy.io as spio
import datetime


#%%
"""
=============================================================================

Download executable for Jason2 s_gdr.

This executable is used to download a single track of altimetry *.nc data
to target directory from NOAA FTP.

Example usage:
    Input Parameters:
        target_dir = '<User's directory>'
        cycle_start = starting cycle number
        cycle_end = ending cycle number
        track = a single track number

    Output:
        Download to the target directory + '/jason2/gdr/s_gdr/'
        The subdirectories are same as the directories in FTP
        
=============================================================================
"""

def download_onetrack(target_dir = None, cycle_start = 0, cycle_end = 327, track = None):
    if track < 10:
        t1 = '0'
    else:
        t1 = ''
    if track < 100:
        t2 = '0'
    else:
        t2 =''
    filematch = '*_'+str(t1)+str(t2)+str(track)+'_*.nc'# find files match sinlge number of track
    for i in tqdm(range(cycle_start, cycle_end), ascii=True, desc='Downloading...'):
        if i == 304:  # in noaa database have missing cycle 304, ignore
            i = i + 1
        if i < 10:
            i1 = '0'
        else:
            i1 = ''
        if i < 100:
            i2 = '0'
        else:
            i2 = ''
        name = 'cycle'+str(i1)+str(i2)+str(i) # number of cycle
        if not os.path.exists(target_dir+'/jason2/gdr/s_gdr/'+name):
            os.makedirs(target_dir+'/jason2/gdr/s_gdr/'+name) # make directory
        save_dir = target_dir+'/jason2/gdr/s_gdr/'+name
        ftp=ftplib.FTP('ftp.nodc.noaa.gov', 'anonymous', 'anonymous')
        ftp.cwd('pub/data.nodc/jason2/gdr/s_gdr/' + name) # connect to ftp
        for filename in ftp.nlst(filematch):
            target_filename = os.path.join(save_dir, os.path.basename(filename))
            with open(target_filename, 'wb') as fhandle:
                ftp.retrbinary('RETR %s' % filename, fhandle.write) # download files
        ftp.quit() # close ftp connection
    return

#%%

"""
=============================================================================

File read executable for Jason2 s_gdr.

This executable is used to read *.nc data, generate variables list and data
matrix for later processing.

Example usage:
    Input Parameters:
        fname = read-in file with netCDF4.Dataset('filepath + *.nc', 'r')
                                                    'r' -- readable only

    Output Parameters:
        header: DataFrame({Variable names, Value<1 or None>})
        data: Numpy array for further process
        
=============================================================================
"""

def read_jason2_PH(fname):
    # Set lists for later use in dataframe 'ncData'
    Variables, Attname, Attval, DataMatrix, MatrixShape = ([] for x in range(5))
    
    # Get corresponding Attribute Units(Attval), Data Matrix(DataMatrix)
    # When variables have Attribute Units(Attval), set Attribute Name(Attname) as 'units'
    # Trace Matrix Shape(MatrixShape) for further calculation
    
    for i in fname.variables:
        Variables.append(i) # Save all Variables
        try:
            attval = fname.variables[i].units # When Variables have units
        except AttributeError:
            attval = None # When Variables don't have units, set Attval as None
        Attval.append(attval) # Save Attribute Units(Attval)
    
        for j in [attval]: # Looking for values in Attribute Units(Attval)
            if j == None:
                attname = None # Attribute have no units, set Attname as None
            else:
                attname = 'units' # Attribute have units, set Attname as 'units'
            Attname.append(attname) # Save Attribute Name(Attname)
    
        get_data_matrix = fname.variables[i][:]
        DataMatrix.append(get_data_matrix) # Save Data Matrix(DataMatrix)
        MatrixShape.append(get_data_matrix.shape) # Save Matrix Shape(MatrixShape)
    
    cols = ['Variables','AttributeValue','AttributeName','DataMatrix','MatrixShape'] # columns order
    ncData = pd.DataFrame({'Variables': Variables,
                           'AttributeName': Attname,
                           'AttributeValue': Attval,
                           'DataMatrix': DataMatrix,
                           'MatrixShape': MatrixShape}, columns = cols)
    
    data_row_num = ncData.at[0, 'MatrixShape'][0] # Number of row in data matrix
    data = np.empty([data_row_num,179]) # create final data matrix
    header_usage = [] # creat space for recording variable usage
    
    s = 0 # set for write data in final matrix
    for i in range(1, 180): 
        header = ncData[:i]['Variables'] # header: Variables
        m = ncData.at[i-1, 'MatrixShape'] # get matrix rows m, n(columns) is not used
    
        if i == 1: # set c = length of time matrix
            c = m
    
        if m == c: # choose the data, which rows match time matrix
            header_usage.append(1) # mark used as 1
            data[:, s] = ncData.at[i-1, 'DataMatrix']
            s = s + 1
        else:
            header_usage.append(np.nan) # mark unused as np.nan
    
    header_usage = pd.Series(header_usage) # change header_usage to pandas.Series
    header = pd.concat([header, header_usage], axis=1) # concatenate with header

    return header, data

#%%
    
"""
=============================================================================

Concatenate all useable tracks with their data matirx.

This executable is used to concatenate useable tracks with their data matrix,
within the chosen window and individual track.

Example usage:
    Input Parameters:
        latmin = Minimum latitude of window
        latmax = Maximum latitude of window
        lonmin = Minimum longitude of window
        lonmax = Maximum longitude of window
        track = single track number
        JA2_dir = '<User's data directory>' # Please give full access
                                            # e.g.: 'C:/JASON2/JASON_2_PH/'

    Output Parameters:
        mat: List contains numpy array for later time series processing

Notes:
    <1> If encounter 'UnboundLocalError: local variable referenced before assignment'
    initialize global variables TmSri & TmSric before function 'corrected_time_series'
    announce 'global TmSri&TmSric' in function itself
    
    <2> Please put 'EGM2008_5.mat' file and main.py in the same directory    
    
=============================================================================
"""

# Please refer to Notes <1>

#TmSri = [] # Time Series
#TmSric = [] # Corrected Time Series

def altprocess(vlat, vlon, SR, latmin, latmax, lonmin, lonmax, track, JA2_dir = None):
    mat = []  
    if track < 10:
        t1 = '0'
    else:
        t1 = ''
    if track < 100:
        t2 = '0'
    else:
        t2 =''
    # interation over cycle 1 to 327, using tqdm to add processing bar
    for c in tqdm(range(1, 328), ascii=True, desc='Processing All Cycles'):
        if c < 10:
            c1 = '0'
        else:
            c1 = ''
        if c < 100:
            c2 = '0'
        else:
            c2 = ''
    # file name for match: JA2_IPH_2PTP + cycle number _ track number _
        id_track = 'JA2_IPH_2PTP'+str(c1)+str(c2)+str(c)+'_'+str(t1)+str(t2)+str(track)+'_' 
        name = 'cycle_'+str(c1)+str(c2)+str(c) # cycle number
        path = os.path.join(JA2_dir + name) # go to data path
        d = os.listdir(path) # list contents
        data_dir = JA2_dir + name + '/'
        data_cycle = []
        for i in range(0, len(d)):
            if id_track in d[i]: # match file name
                fname = Dataset(data_dir+d[i], 'r')
                header, data = read_jason2_PH(fname)
                m, _ = data.shape

                for j in range(0,m):
                    # find data within the window
                    if data[j,2]<latmax and data[j,2]>latmin and data[j,1]<lonmax and data[j,1]>lonmin:
                        data_cycle.append(data[j,:])
                y = np.array(data_cycle)
        mat.append(y)
        y = None # after append reset y as None
        
    mat = [x for x in mat if x is not None] # eliminate None elements
    mat = [x for x in mat if len(x) != 0] # eliminate empty elements

    E1, N1, _, _ = utm.from_latlon(vlat, vlon)
    
    header_used = [] # get header list
    for n in range(len(header)):
        if header[0][n] == 1:
            header_used.append(header['Variables'][n])

    if mat:
#        global TmSri # modify global variable TmSri
        TmSri = [] # Time Series
        for i in range(0, len(mat)):
            s, _ = mat[i].shape
            if s > 0:
                for j in range(0, s):
                    E2, N2, _, _ = utm.from_latlon(mat[i][j][2], mat[i][j][1])
                    Dist = np.sqrt(np.square(N1 - N2) + np.square(E1 - E2))
                    if Dist <= SR: # if distance shorter than search radious
                        if len(mat[i][j]) == 179:
                            TmSri.append(mat[i][j])
    TmSri = np.array(TmSri)
     
    if TmSri.size:
#        global TmSric # modify global variable TmSric
        TmSric = [] # Corrected Time Series
        rowone = [TmSri[0]]
        TmSric = pd.DataFrame(np.concatenate((rowone, TmSri), axis=0), columns=header_used)
        # drop duplicates rows if column 'time', 'lat' and 'lon' are same
        TmSric = TmSric.drop_duplicates(subset = ['time','lat','lon']) 
#    
#    # get geoid height, the scipy interpolation might take 3 - 5 minutes
#    # read EGM2008_5.mat, please refer to Notes <2>
#    EGM = spio.loadmat('EGM2008_5.mat') 
#    # change XX, YY, GH to 1-dimensional for further interpolation
#    XX = EGM['XX'].flatten()
#    YY = EGM['YY'].flatten()
#    GH = EGM['GH'].flatten()
#    points = np.column_stack((XX,YY))
#    
#    if TmSric.size:
#        latitude = TmSric['lat'].as_matrix()
#        longitude = TmSric['lon'].as_matrix()
#
#    print('Interpolating griddata from EGM2008_5...')
#    print('This might take a while, normally 3 - 5 minutes...')
#    height = sp.interpolate.griddata(points, GH, (longitude, latitude), method='linear')
#    TmSric['Geoid height'] = pd.Series(height, index = TmSric.index)
#    print('Interpolation terminated')
    
    return TmSric

#%%

"""
=============================================================================

Concatenate all useable tracks with their data matirx.

This executable is used to concatenate useable tracks with their data matrix,
within the chosen window and individual track.

Example usage:
    Input Parameters:
        latmin = Minimum latitude of window
        latmax = Maximum latitude of window
        lonmin = Minimum longitude of window
        lonmax = Maximum longitude of window
        track = single track number
        JA2_dir = '<User's data directory>' # Please give full access
                                            # e.g.: 'C:/JASON2/JASON_2_PH/'

    Output Parameters:
        mat: List contains numpy array for later time series processing

Notes:
    <1> If encounter 'UnboundLocalError: local variable referenced before assignment'
    initialize global variables TmSri & TmSric before function 'corrected_time_series'
    announce 'global TmSri&TmSric' in function itself
    
    <2> Please put 'EGM2008_5.mat' file and main.py in the same directory    
    
=============================================================================
"""

def time_series(TmSric): 
    if TmSric.size:
        m, _ = TmSric.shape
        TmSric = TmSric.reset_index(drop=True) # reset dataframe index 
        
        # initialization
        Alt_1hz,Range_ku,Range_c,Range_oce3_ku,Range_oce3_c=(np.zeros((m,1)) for x in range(5))
        Range_red3_c,Range_ice3_ku,Range_ice3_c,Range_red3_ku=(np.zeros((m,1)) for x in range(4))
        Cor = np.zeros((m,8)) # Matrix of correction parameters 
        
        for i in range(0, m):    
            # Altitude, 1Hz
            Alt_1hz[i,0] = TmSric['alt'][i]   
            
            ### 1Hz, Corrections ###    
            # Inverse barometric correction
            if TmSric['inv_bar_corr'][i] != 32767 and TmSric['inv_bar_corr'][i] != -32768 and TmSric['inv_bar_corr'][i] != 2.147483647e+09:
                InvBar_ku = TmSric['inv_bar_corr'][i]
            else:
                InvBar_ku = 0
    
            # Sea state bias
            if TmSric['sea_state_bias_ku'][i] != 32767 and TmSric['sea_state_bias_ku'][i] != -32768 and TmSric['sea_state_bias_ku'][i] != 2.147483647e+09:
                SeSbias_ku = TmSric['sea_state_bias_ku'][i]
            else:
                SeSbias_ku = 0
    #        if TmSric['sea_state_bias_ku'][i] != -32768 and TmSric['sea_state_bias_ku'][i] != 2.147483647e+09:
    #            SeSbias_ku = TmSric['sea_state_bias_ku'][i]
    ##        elif TmSric['sea_state_bias_ku'][i] >= 32767:
    ##            SeSbias_ku = 3.2767
    #        else:
    #            SeSbias_ku = 0
    
            # Ionospheric correction
            if TmSric['iono_corr_alt_ku'][i] != 32767 and TmSric['iono_corr_alt_ku'][i] != -32768 and TmSric['iono_corr_alt_ku'][i] != 2.147483647e+09:
                IonCor_ku = TmSric['iono_corr_alt_ku'][i]
            else:
                IonCor_ku = 3.2767
    #            IonCor_ku = 0
    
            # Ocean tide
            if TmSric['ocean_tide_sol1'][i] != 32767 and TmSric['ocean_tide_sol1'][i] != -32768 and TmSric['ocean_tide_sol1'][i] != 2.147483647e+09:
                OcTide_ku = TmSric['ocean_tide_sol1'][i]
            else:
                OcTide_ku = 0
    
            # Pole Tide
            if TmSric['pole_tide'][i] != 32767 and TmSric['pole_tide'][i] != -32768 and TmSric['pole_tide'][i] != 2.147483647e+09:
                PoTide_ku = TmSric['pole_tide'][i]
            else:
                PoTide_ku = 0
    
            # Earth tide
            if TmSric['solid_earth_tide'][i] != 32767 and TmSric['solid_earth_tide'][i] != -32768 and TmSric['solid_earth_tide'][i] != 2.147483647e+09:
                ETide_ku = TmSric['solid_earth_tide'][i]
            else:
                ETide_ku = 0
    
            # Wet tropospheric correction
            if TmSric['model_wet_tropo_corr'][i] != 32767 and TmSric['model_wet_tropo_corr'][i] != -32768 and TmSric['model_wet_tropo_corr'][i] != 2.147483647e+09:
                WTCor_ku = TmSric['model_wet_tropo_corr'][i]
    #        if TmSric['rad_wet_tropo_corr'][i] != 32767 and TmSric['rad_wet_tropo_corr'][i] != -32768 and TmSric['rad_wet_tropo_corr'][i] != 2.147483647e+09:
    #            WTCor_ku = TmSric['rad_wet_tropo_corr'][i]
            else:
                WTCor_ku = 0
    
            # Dry tropospheric correction
            if TmSric['model_dry_tropo_corr'][i] != 32767 and TmSric['model_dry_tropo_corr'][i] != -32768 and TmSric['model_dry_tropo_corr'][i] != 2.147483647e+09:
                DTCor_ku = TmSric['model_dry_tropo_corr'][i]
            else:
                DTCor_ku = 0
                
            # Performing corrections
#            Correction_ku = InvBar_ku + SeSbias_ku + IonCor_ku + OcTide_ku + PoTide_ku + ETide_ku + WTCor_ku + DTCor_ku + TmSric['Geoid height'][i]
            Correction_ku = InvBar_ku + SeSbias_ku + IonCor_ku + OcTide_ku + PoTide_ku + ETide_ku + WTCor_ku + DTCor_ku + TmSric['geoid_EGM2008'][i]
            # Correction matrix
            Cor[i,:] = [InvBar_ku, SeSbias_ku, IonCor_ku, OcTide_ku, PoTide_ku, ETide_ku, WTCor_ku, DTCor_ku]
    
            # range correction
            Range_ku[i, 0] = TmSric['range_ku'][i] + Correction_ku
            Range_c[i, 0] = TmSric['range_c'][i] + Correction_ku
            Range_oce3_ku[i, 0] = TmSric['range_oce3_ku'][i] + Correction_ku
            Range_oce3_c[i, 0] = TmSric['range_oce3_c'][i] + Correction_ku
            Range_red3_ku[i, 0] = TmSric['range_red3_ku'][i] + Correction_ku
            Range_red3_c[i, 0] = TmSric['range_red3_c'][i] + Correction_ku
            Range_ice3_ku[i, 0] = TmSric['range_ice3_ku'][i] + Correction_ku
            Range_ice3_c[i, 0] = TmSric['range_ice3_c'][i] + Correction_ku
            
        # Sea Surface Height
        SSH_ku = Alt_1hz - Range_ku
        SSH_c = Alt_1hz - Range_c
        SSH_oce3_ku = Alt_1hz - Range_oce3_ku
        SSH_oce3_c = Alt_1hz - Range_oce3_c
        SSH_red3_ku = Alt_1hz - Range_red3_ku
        SSH_red3_c = Alt_1hz - Range_red3_c
        SSH_ice3_ku = Alt_1hz - Range_ice3_ku
        SSH_ice3_c = Alt_1hz - Range_ice3_c
        
        
        # Creating Time Series Dataframe
        
        # initialization
        s = 0 # set for write data in time series dataframe
        ind = [0] # initialize index
        y, mo, d = (np.zeros((m,1)) for x in range(3))
        TSyear, TSmonth, TSday, TStotal, SSH_ku_mean, SSH_ku_median, SSH_ku_std=([] for x in range(7))
        SSH_c_mean, SSH_c_median, SSH_c_std, SSH_oce3_ku_mean, SSH_oce3_ku_median, SSH_oce3_ku_std=([] for x in range(6))
        SSH_oce3_c_mean, SSH_oce3_c_median, SSH_oce3_c_std, SSH_red3_ku_mean, SSH_red3_ku_median=([] for x in range(5))
        SSH_red3_ku_std, SSH_red3_c_mean, SSH_red3_c_median, SSH_red3_c_std, SSH_ice3_ku_mean=([] for x in range(5))
        SSH_ice3_ku_median, SSH_ice3_ku_std, SSH_ice3_c_mean, SSH_ice3_c_median, SSH_ice3_c_std=([] for x in range(5))
        
        # get first date
        a0 = datetime.timedelta(seconds=TmSric['time'][0]) + datetime.datetime(2000, 1, 1)
        y[0] = a0.timetuple()[0]
        mo[0] = a0.timetuple()[1]
        d[0] = a0.timetuple()[2]
    
        for t in range(1,m):
            
            a = datetime.timedelta(seconds=TmSric['time'][t]) + datetime.datetime(2000, 1, 1)
            y[t] = a.timetuple()[0]
            mo[t] = a.timetuple()[1]
            d[t] = a.timetuple()[2]
            # append data under these conditions
            if d[t,0] != d[t-1,0] or t+1 == m:
                s = s + 1 
                ind.append(t) # append index
                
                # append time lists
                TStotal.append(datetime.timedelta(seconds=TmSric['time'][t-1]) + datetime.datetime(2000, 1, 1))
                TSyear.append(y[t-1,0])
                TSmonth.append(mo[t-1,0])
                TSday.append(d[t-1,0])
    
                # SSH ku band
                SSH_ku_mean.append(np.nanmean(SSH_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_ku_median.append(np.nanmedian(SSH_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_ku_std.append(np.std(SSH_ku[ind[s-1]:ind[s]-1], axis = 0))
                
                # SSH c band 
                SSH_c_mean.append(np.nanmean(SSH_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_c_median.append(np.nanmedian(SSH_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_c_std.append(np.std(SSH_c[ind[s-1]:ind[s]-1], axis = 0))
                
                # SSH oce3 ku band
                SSH_oce3_ku_mean.append(np.nanmean(SSH_oce3_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_oce3_ku_median.append(np.nanmedian(SSH_oce3_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_oce3_ku_std.append(np.std(SSH_oce3_ku[ind[s-1]:ind[s]-1], axis = 0))
                
                #SSH oce3 c band
                SSH_oce3_c_mean.append(np.nanmean(SSH_oce3_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_oce3_c_median.append(np.nanmedian(SSH_oce3_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_oce3_c_std.append(np.std(SSH_oce3_c[ind[s-1]:ind[s]-1], axis = 0))
                
                #SSH red3 ku band
                SSH_red3_ku_mean.append(np.nanmean(SSH_red3_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_red3_ku_median.append(np.nanmedian(SSH_red3_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_red3_ku_std.append(np.std(SSH_red3_ku[ind[s-1]:ind[s]-1], axis = 0))
                
                #SSH red3 c band
                SSH_red3_c_mean.append(np.nanmean(SSH_red3_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_red3_c_median.append(np.nanmedian(SSH_red3_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_red3_c_std.append(np.std(SSH_red3_c[ind[s-1]:ind[s]-1], axis = 0))
                
                #SSH ice3 ku band
                SSH_ice3_ku_mean.append(np.nanmean(SSH_ice3_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_ice3_ku_median.append(np.nanmedian(SSH_ice3_ku[ind[s-1]:ind[s]-1], axis = 0))
                SSH_ice3_ku_std.append(np.std(SSH_ice3_ku[ind[s-1]:ind[s]-1], axis = 0))
                
                #SSH ice3 c band
                SSH_ice3_c_mean.append(np.nanmean(SSH_ice3_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_ice3_c_median.append(np.nanmedian(SSH_ice3_c[ind[s-1]:ind[s]-1], axis = 0))
                SSH_ice3_c_std.append(np.std(SSH_ice3_c[ind[s-1]:ind[s]-1], axis = 0))
                
                
        # reshape list as 1-dimensional array for dataframe        
        SSH_ku_mean = np.asarray(SSH_ku_mean).reshape((-1))
        SSH_ku_median = np.asarray(SSH_ku_median).reshape((-1))
        SSH_ku_std = np.asarray(SSH_ku_std).reshape((-1))
        
        SSH_c_mean = np.asarray(SSH_c_mean).reshape((-1))
        SSH_c_median = np.asarray(SSH_c_median).reshape((-1))
        SSH_c_std = np.asarray(SSH_c_std).reshape((-1))
        
        SSH_oce3_ku_mean = np.asarray(SSH_oce3_ku_mean).reshape((-1))
        SSH_oce3_ku_median = np.asarray(SSH_oce3_ku_median).reshape((-1))
        SSH_oce3_ku_std = np.asarray(SSH_oce3_ku_std).reshape((-1))
         
        SSH_oce3_c_mean = np.asarray(SSH_oce3_c_mean).reshape((-1))
        SSH_oce3_c_median = np.asarray(SSH_oce3_c_median).reshape((-1))
        SSH_oce3_c_std = np.asarray(SSH_oce3_c_std).reshape((-1))
         
        SSH_red3_ku_mean = np.asarray(SSH_red3_ku_mean).reshape((-1))
        SSH_red3_ku_median = np.asarray(SSH_red3_ku_median).reshape((-1))
        SSH_red3_ku_std = np.asarray(SSH_red3_ku_std).reshape((-1))
         
        SSH_red3_c_mean = np.asarray(SSH_red3_c_mean).reshape((-1))
        SSH_red3_c_median = np.asarray(SSH_red3_c_median).reshape((-1))
        SSH_red3_c_std = np.asarray(SSH_red3_c_std).reshape((-1))
         
        SSH_ice3_ku_mean = np.asarray(SSH_ice3_ku_mean).reshape((-1))
        SSH_ice3_ku_median = np.asarray(SSH_ice3_ku_median).reshape((-1))
        SSH_ice3_ku_std = np.asarray(SSH_ice3_ku_std).reshape((-1))
         
        SSH_ice3_c_mean = np.asarray(SSH_ice3_c_mean).reshape((-1))
        SSH_ice3_c_median = np.asarray(SSH_ice3_c_median).reshape((-1))
        SSH_ice3_c_std = np.asarray(SSH_ice3_c_std).reshape((-1))     
        
        # columns order
        cols = ['Date Detail','Year','Month','Day','SSH_ku_mean','SSH_ku_median','SSH_ku_std','SSH_c_mean',
                'SSH_c_median','SSH_c_std','SSH_oce3_ku_mean','SSH_oce3_ku_median','SSH_oce3_ku_std',
                'SSH_oce3_c_mean','SSH_oce3_c_median','SSH_oce3_c_std','SSH_red3_ku_mean','SSH_red3_ku_median',
                'SSH_red3_ku_std','SSH_red3_c_mean','SSH_red3_c_median','SSH_red3_c_std','SSH_ice3_ku_mean',
                'SSH_ice3_ku_median','SSH_ice3_ku_std','SSH_ice3_c_mean','SSH_ice3_c_median','SSH_ice3_c_std']
        
        # create time series dataframe
        TS = pd.DataFrame({'Date Detail':TStotal,'Year':TSyear,'Month':TSmonth,'Day':TSday,'SSH_ku_mean':SSH_ku_mean,
                           'SSH_ku_median':SSH_ku_median,'SSH_ku_std':SSH_ku_std,'SSH_c_mean':SSH_c_mean,
                           'SSH_c_median':SSH_c_median,'SSH_c_std':SSH_c_std,'SSH_oce3_ku_mean':SSH_oce3_ku_mean,
                           'SSH_oce3_ku_median':SSH_oce3_ku_median,'SSH_oce3_ku_std':SSH_oce3_ku_std,
                           'SSH_oce3_c_mean':SSH_oce3_c_mean,'SSH_oce3_c_median':SSH_oce3_c_median,
                           'SSH_oce3_c_std':SSH_oce3_c_std,'SSH_red3_ku_mean':SSH_red3_ku_mean,
                           'SSH_red3_ku_median':SSH_red3_ku_median,'SSH_red3_ku_std':SSH_red3_ku_std,
                           'SSH_red3_c_mean':SSH_red3_c_mean,'SSH_red3_c_median':SSH_red3_c_median,
                           'SSH_red3_c_std':SSH_red3_c_std,'SSH_ice3_ku_mean':SSH_ice3_ku_mean,
                           'SSH_ice3_ku_median':SSH_ice3_ku_median,'SSH_ice3_ku_std':SSH_ice3_ku_std,
                           'SSH_ice3_c_mean':SSH_ice3_c_mean,'SSH_ice3_c_median':SSH_ice3_c_median,
                           'SSH_ice3_c_std':SSH_ice3_c_std}, columns = cols) 
    
    return TS

#%%  
 
TmSric = altprocess(8.0002, 7.749, 1000, 7.9, 8.1, 7.6, 7.9, 20, JA2_dir = 'D:/JASON2/JASON_2_PH/')    
TS = time_series(TmSric)    
TS = TS.drop_duplicates(subset=['Year','Month','Day'], keep=False)

##def TStcor(TS): # drop duplicates and get mean values from TS, return TSc
#_, n = TS.shape
#
#TSc = TS.drop_duplicates(subset=['Year','Month','Day'], keep=False)
#d = TS[['Year', 'Month', 'Day']]
#TSc_res = np.zeros((TSc.shape))
#
#for i in range(0,len(TSc)):
#    if 
## return TSc

def Res_data(data):  

    m, n = data.shape
    residual = np.zeros((m, n))
    
    if data.size:
        f1 = data.index[data['Month'] == 1].tolist(); f2 = data.index[data['Month'] == 2].tolist() 
        f3 = data.index[data['Month'] == 3].tolist(); f4 = data.index[data['Month'] == 4].tolist()
        f5 = data.index[data['Month'] == 5].tolist(); f6 = data.index[data['Month'] == 6].tolist()
        f7 = data.index[data['Month'] == 7].tolist(); f8 = data.index[data['Month'] == 8].tolist()
        f9 = data.index[data['Month'] == 9].tolist(); f10 = data.index[data['Month'] == 10].tolist()
        f11 = data.index[data['Month'] == 11].tolist(); f12 = data.index[data['Month'] == 12].tolist()
        Jan1 = np.zeros((1, n-2)); Feb1 = np.zeros((1, n-2)); Mar1 = np.zeros((1, n-2))
        Apr1 = np.zeros((1, n-2)); May1 = np.zeros((1, n-2)); Jun1 = np.zeros((1, n-2))
        Jul1 = np.zeros((1, n-2)); Aug1 = np.zeros((1, n-2)); Sep1 = np.zeros((1, n-2))
        Oct1 = np.zeros((1, n-2)); Nov1 = np.zeros((1, n-2)); Dec1 = np.zeros((1, n-2))
        M_Data = np.zeros((12, n-2)); Mts = np.zeros((m,1))
        
        for i in range(2,n):
    
            Jan1[0, i-2] = np.nanmean(data.iloc[f1,i]); Feb1[0, i-2] = np.nanmean(data.iloc[f2,i]);
            Mar1[0, i-2] = np.nanmean(data.iloc[f3,i]); Apr1[0, i-2] = np.nanmean(data.iloc[f4,i]);
            May1[0, i-2] = np.nanmean(data.iloc[f5,i]); Jun1[0, i-2] = np.nanmean(data.iloc[f6,i]);
            Jul1[0, i-2] = np.nanmean(data.iloc[f7,i]); Aug1[0, i-2] = np.nanmean(data.iloc[f8,i]);
            Sep1[0, i-2] = np.nanmean(data.iloc[f9,i]); Oct1[0, i-2] = np.nanmean(data.iloc[f10,i]);
            Nov1[0, i-2] = np.nanmean(data.iloc[f11,i]); Dec1[0, i-2] = np.nanmean(data.iloc[f12,i]);
            
        M_Data[0,:] = Jan1; M_Data[1,:] = Feb1; M_Data[2,:] = Mar1; M_Data[3,:] = Apr1;
        M_Data[4,:] = May1; M_Data[5,:] = Jun1; M_Data[6,:] = Jul1; M_Data[7,:] = Aug1;
        M_Data[8,:] = Sep1; M_Data[9,:] = Oct1; M_Data[10,:] = Nov1; M_Data[11,:] = Dec1;
        
        residual[:,0] = data['Year']
        residual[:,1] = data['Month']
        for i in range(0, m):
            if n <= 3: 
                residual[i,2:n] = data.iloc[i,2:n] - M_Data[int(data['Month'][i])-1, 0]
                Mts[i, 0] = M_Data[int(data['Month'][i])-1, 0]
            else:
                residual[i,2:n] = data.iloc[i,2:n] - M_Data[int(data['Month'][i])-1, 1]
                Mts[i, 0] = M_Data[int(data['Month'][i])-1, 1]

##else:
##    for i in range(0, 11):
##        M_Data[i,1] = i
    return residual, M_Data, Mts

#
##def Outcor(TS,col,kr,rep):
kr = 2.9
rep = 1
TScor = TS
level = TS['SSH_ice3_ku_median']  # level = TS.iloc[:,col]
_, M, _ = Res_data(TS)
non_empty = 1
while non_empty != 0:
    level_resi, _, _ = Res_data(TScor[['Year','Month','SSH_ice3_ku_median']].copy())
    o1 = np.where(level_resi[:,2] >= np.nanmedian(level_resi[:,2])+kr*np.nanstd(level_resi[:,2]))
    o2 = np.where(level_resi[:,2] <= np.nanmedian(level_resi[:,2])-kr*np.nanstd(level_resi[:,2]))
    o1 = np.asarray(o1)
    o2 = np.asarray(o2)
    o = np.concatenate((o1,o2), axis=1)
    try: # when no values found
        f_ind1 = np.argmax(level_resi[o,0]) # might be a problem
        level.iloc[o[0,f_ind1]] = np.nan
        TScor['SSH_ice3_ku_median'] = level
    except ValueError:
        non_empty = 0
    

if rep == 1:
    _, M, _ = Res_data(TScor)
#    f_ind2 = TScor['SSH_ice3_ku_median'].index[TScor['SSH_ice3_ku_median'].apply(np.isnan)]
    f_ind2 = pd.isnull(TScor).any(1).nonzero()[0]
    for i in range(0, len(f_ind2)):
        TScor['SSH_ice3_ku_median'][f_ind2[i]] = M[int(TScor['Month'][f_ind2[i]])-1,21]
                                                                            # column-2? = 
                                                                            # why column-2
#return TScor
                                                                            
                                                                            
#TS_t = Outcor(TS, 19, 2.9, 1);
#s = find(TS_t(:,20) > 0.7);
#s0 = find(TS_t(:,20) < 0.7);
#TS_t(s,20) = nanmean(TS_t(s0,20));
#TS_final = [TS_t(:,5) TS_t(:,19:20)];

TS_final = TScor[['Date Detail','SSH_ice3_ku_median','SSH_ice3_ku_std']].copy()

#%%
import matplotlib.pyplot as plt
height_res = TS_final['SSH_ice3_ku_median'].as_matrix()
time_res = TS_final['Date Detail'].as_matrix()
plt.plot(time_res,height_res)
plt.show()
