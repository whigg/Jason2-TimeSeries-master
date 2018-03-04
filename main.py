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

This executable is used to read *.nc data, generate used variables list and data
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

def read_jason2_PH_nc(fname):
    # Set lists for later use in forming dataframe 'ncData'
    Variables, Attname, Attval, DataMatrix, MatrixShape = ([] for x in range(5))

#####################################################################
# Loop through all variables in netCDF data, 
# get corresponding Attribute Units(Attval), Data Matrix(DataMatrix)
# When variables have Attribute Units(Attval), set Attribute Name(Attname) as 'units'
# Trace Matrix Shape(MatrixShape) for further calculation and bug checking
#####################################################################

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

    ncData = pd.DataFrame({'Variables': Variables,
                           'AttributeName': Attname,
                           'AttributeValue': Attval,
                           'DataMatrix': DataMatrix,
                           'MatrixShape': MatrixShape})

#  Possible bug: e.g. shape of 'amplitude_20hz_ku_mle3' in Matlab 20x2963, Python 2963x20
#  Solved by check m = shape_matrix[0] and n != 20

    data_row_num = ncData.at[0, 'MatrixShape'][0] # Number of row in final data matrix
    data = np.empty([data_row_num,138]) # create empty space for final data matrix
    header_usage = [] # creat empty space for recording variable usage
    s = 0 # set for write in final data matrix iteration

    for i in range(1, 180): 
        header = ncData[:i]['Variables'] # header: Variable names in this iteration
        try:
            shape_matrix = ncData.at[i-1, 'MatrixShape'] # get matrix shape m, n
            m, n = shape_matrix
        except ValueError: # 1D matrix set n as 1
            m = shape_matrix[0]
            n = 1

        if i == 1: # set c = length of time matrix
            c = m

        if m == c and n != 20: # choose the data, we don't use 20hz data here
            header_usage.append(1) # mark used as 1
            data[:, s] = ncData.at[i-1, 'DataMatrix']
            s = s + 1
        else:
            header_usage.append(None) # mark unused as None

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

    Output Parameters:
        mat: List contains numpy array for later time series processing
        
=============================================================================
"""

def JA2_PH_crt(latmin, latmax, lonmin, lonmax, track):
    mat = []
    if lonmin < 0:
        lonmin = lonmin + 360
    if lonmax < 0:
        lonmax = lonmax + 360

    if track < 10:
        t1 = '0'
    else:
        t1 = ''
    if track < 100:
        t2 = '0'
    else:
        t2 =''
    for c in tqdm(range(0, 327), ascii=True, desc='Processing All Cycles'):
        if c == 304:
            c = c + 1  # dataset have no cycle304
        if c < 10:
            c1 = '0'
        else:
            c1 = ''
        if c < 100:
            c2 = '0'
        else:
            c2 = ''
        id_track = 'JA2_GPS_2PdP'+str(c1)+str(c2)+str(c)+'_'+str(t1)+str(t2)+str(track)+'_'
        name = 'cycle'+str(c1)+str(c2)+str(c) # cycle number
        path = os.path.join('D:/jason2/gdr/s_gdr/' + name)
        d = os.listdir(path)
        data_dir = 'D:/jason2/gdr/s_gdr/' + name + '/'
        data_cycle = []
        for i in range(0, len(d)):
            if id_track in d[i]:
                fname = Dataset(data_dir+d[i], 'r')
                _, data = read_jason2_PH_nc(fname)
                m, _ = data.shape

                for j in range(0,m):
                    if data[j,1]<latmax and data[j,1]>latmin and data[j,2]<lonmax and data[j,2]>lonmin:
                        # Lat varies from around 80 to 248, Lon varies from around -66 to 66
                        # Tourian might got it wrong with column number of 'mat' matrix
                        # The second column should be Lat, third should be Lon
                        # in python data[:,1] = lat, data[:,2] = lon
                        # in matlab data[:,2] = lat, data[:,3] = lon
                        data_cycle.append(data[j,:])
                y = np.array(data_cycle)
        mat.append(y)
        y = None
    mat = [x for x in mat if x is not None]
    mat = [x for x in mat if len(x) != 0]

    return mat

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

    Output Parameters:
        mat: List contains numpy array for later time series processing
        
=============================================================================
"""

def correcting_TS(mat, vlat, vlon, SR):
    E1, N1, _, _ = utm.from_latlon(vlat, vlon)
    if mat:
        TmSri = [] # Time Series
        for i in range(0, len(mat)):
            s, _ = mat[i].shape
            if s > 0:
                for j in range(0, s):
                    E2, N2, _, _ = utm.from_latlon(mat[i][j][1], mat[i][j][2]-360)
                    Dist = np.sqrt(np.square(N1 - N2) + np.square(E1 - E2))
                    if Dist <= SR:
                        if len(mat[i][j]) == 138:
                            TmSri.append(mat[i][j])
        TmSri = np.array(TmSri)
       
    if TmSri.size:
        rowone = [TmSri[0]]
        TmSric = np.concatenate((rowone, TmSri), axis=0)
        TmSric = np.unique(TmSric, axis=0)
        # the comparation should within the first 3 columns   
         
#def geoid_height():
    EGM = sp.io.loadmat('EGM2008_5.mat')
    XX = EGM['XX']
    YY = EGM['YY']
    GH = EGM['GH']
    #XX = np.genfromtxt('EGM2008_5_XX.csv', delimiter=',')
    #YY = np.genfromtxt('EGM2008_5_YY.csv', delimiter=',')
    #GH = np.genfromtxt('EGM2008_5_GH.csv', delimiter=',')
    #TmSri = correcting_TS()
    if TmSric.size:
        lat = TmSric[:,1]
        lon = TmSric[:,2]
        if np.where(lon > 180):
            lon = lon - 360
#    height = scipy.interpolate.griddata((EGM['XX'],EGM['YY']), EGM['GH'], (lon, lat), method='cubic')
#    TmSri = np.concatenate()
    return TmSric

#%%
    
def geophysical_cor(TmSric): 
    if TmSric.size:
        m, _ = TmSric.shape
        Alt_1hz,Range_ku,Range_c,Range_oce3_ku,Range_oce3_c=(np.zeros((m,1)) for x in range(5))
        Range_red3_c,Range_ice3_ku,Range_ice3_c,Range_red3_ku=(np.zeros((m,1)) for x in range(4))
        for i in range(0, m):
            # Altitude, 1Hz
            Alt_1hz[i,0] = TmSric[i,25]
            # 1Hz, Corrections
            # Inverse barometric correction
            if TmSric[i, 98] != 32767 and TmSric[i, 98] != -32768 and TmSric[i, 98] != 2.147483647000000e+05:
                InvBar_ku = TmSric[i, 98]
            else:
                InvBar_ku = 0
    
            # Sea state bias
            if TmSric[i, 50] != 32767 and TmSric[i, 50] != -32768 and TmSric[i, 50] != 2.147483647000000e+05:
                SeSbias_ku = TmSric[i, 50]
            else:
                SeSbias_ku = 0
    
            # Ionospheric correction
            if TmSric[i, 46] != 32767 and TmSric[1, 46] != -32768 and TmSric[i, 46] != 2.147483647000000e+05:
                IonCor_ku = TmSric[i, 46]
            else:
                IonCor_ku = 0
    
            # Ocean tide
            if TmSric[i, 100] != 32767 and TmSric[i, 100] != -32768 and TmSric[i, 100] != 2.147483647000000e+05:
                OcTide_ku = TmSric[i, 100]
            else:
                OcTide_ku = 0
    
            # Polar Tide
            if TmSric[i, 109] != 32767 and TmSric[i, 109] != -32768 and TmSric[i, 109] != 2.147483647000000e+05:
                PoTide_ku = TmSric[i, 109]
            else:
                PoTide_ku = 0
    
            # Earth tide
            if TmSric[i, 108] != 32767 and TmSric[i, 108] != -32768 and TmSric[i, 108] != 2.147483647000000e+05:
                ETide_ku = TmSric[i, 108]
            else:
                ETide_ku = 0
    
            # Wet tropospheric correction
            if TmSric[i, 41] != 32767 and TmSric[i, 41] != -32768 and TmSric[i, 41] != 2.147483647000000e+05:
                WTCor_ku = TmSric[i, 41]
            else:
                WTCor_ku = 0
    
            # Dry tropospheric correction
            if TmSric[i, 39] != 32767 and TmSric[i, 39] != -32768 and TmSric[i, 39] != 2.147483647000000e+05:
                DTCor_ku = TmSric[i, 39]
            else:
                DTCor_ku = 0
                
            # Performing corrections
    #        Correction_ku = InvBar_ku + SeSbias_ku + IonCor_ku + OcTide_ku + PoTide_ku + ETide_ku + WTCor_ku + DTCor_ku + TmSric(i,180)
            Correction_ku=InvBar_ku+SeSbias_ku+IonCor_ku+OcTide_ku+PoTide_ku+ETide_ku+WTCor_ku+DTCor_ku
    #        Cor[i,:] = [InvBar_ku, SeSbias_ku, IonCor_ku, OcTide_ku, PoTide_ku, ETide_ku, WTCor_ku, DTCor_ku]
    
            Range_ku[i, 0] = TmSric[i, 27] + Correction_ku
            Range_c[i, 0] = TmSric[i, 28] + Correction_ku
            Range_oce3_ku[i, 0] = TmSric[i, 29] + Correction_ku
            Range_oce3_c[i, 0] = TmSric[i, 30] + Correction_ku
            Range_red3_ku[i, 0] = TmSric[i, 31] + Correction_ku
            Range_red3_c[i, 0] = TmSric[i, 32] + Correction_ku
            Range_ice3_ku[i, 0] = TmSric[i, 33] + Correction_ku
            Range_ice3_c[i, 0] = TmSric[i, 34] + Correction_ku
            
        # Sea Surface Height
    
        SSH_ku = Alt_1hz - Range_ku
        SSH_c = Alt_1hz - Range_c
        SSH_oce3_ku = Alt_1hz - Range_oce3_ku
        SSH_oce3_c = Alt_1hz - Range_oce3_c
        SSH_red3_ku = Alt_1hz - Range_red3_ku
        SSH_red3_c = Alt_1hz - Range_red3_c
        SSH_ice3_ku = Alt_1hz - Range_ice3_ku  # Used
        SSH_ice3_c = Alt_1hz - Range_ice3_c
        
        # Creating Time Matrix
        
    #    TSyear,TSmonth,TSday,TSyears,TSsec,TSku_mean,TSku_median,TSku_std=([] for x in range(8))
    #    TSc_mean,TSc_median,TSc_std,TSoce3ku_mean,TSoce3ku_median,TSoce3ku_std=([] for x in range(6))
    #    TSoce3c_mean,TSoce3c_median,TSoce3c_std,TSred3ku_mean,TSred3ku_median,TSred3ku_std=([] for x in range(6))
    #    TSred3c_mean,TSred3c_median,TSred3c_std,TSice3ku_mean,TSice3ku_median,TSice3ku_std=([] for x in range(6))
    #    TSice3c_mean,TSice3c_median,TSice3c_std=([] for x in range(3))

        s = 0
        ind = [0]
        TSyear,TSmonth,TSday,TStotal,SSH_ice3_ku_median,SSH_ice3_ku_std=([] for x in range(6))
        y, mo, d = (np.zeros((m,1)) for x in range(3))
        a0 = datetime.timedelta(seconds=TmSric[0,0]) + datetime.datetime(2000,1,1)
        y[0] = a0.timetuple()[0]
        mo[0] = a0.timetuple()[1]
        d[0] = a0.timetuple()[2]
        
        for t in range(1,m):
            a = datetime.timedelta(seconds=TmSric[t,0]) + datetime.datetime(2000,1,1)
            y[t] = a.timetuple()[0]
            mo[t] = a.timetuple()[1]
            d[t] = a.timetuple()[2]
            if d[t,0] != d[t-1,0] or t+1 == m:
                ind.append(t)
                TSyear.append(y[t-1,0])
                TSmonth.append(mo[t-1,0])
                TSday.append(d[t-1,0])
                
                s = s + 1
                TStotal.append(datetime.timedelta(seconds=TmSric[t-1,0]) + datetime.datetime(2000,1,1))
                SSH_ice3_ku_median.append(np.nanmedian(SSH_ice3_ku[ind[s-2]], axis=0))
                SSH_ice3_ku_std.append(np.std(SSH_ice3_ku[ind[s-2]], axis=0))
        cols = ['Date Detail','Year','Month','Day','SSH_ice3_ku_median','SSH_ice3_ku_std']
        TS = pd.DataFrame({'Year': TSyear, 
                                 'Month': TSmonth, 
                                 'Day': TSday, 
                                 'Date Detail': TStotal, 
                                 'SSH_ice3_ku_median': SSH_ice3_ku_median,
                                 'SSH_ice3_ku_std': SSH_ice3_ku_std}, columns = cols)
        TS = TS.set_index('Date Detail')
    #            days = datetime.timedelta(seconds=TmSric[t-1,0]) + datetime.datetime(2000,1,1) - datetime.datetime(int(y[t-1,0]),1,1)
    #            if y[t-1,0]==1996 or y[t-1,0]==2000 or y[t-1,0]==2004 or y[t-1,0]==2008 or y[t-1,0]==2012 or y[t-1,0]==2016:
    #                le = 366
    #            else:
    #                le = 365   
    #            TSyears.append(days/le+TSyear[s-1]) 

#%%