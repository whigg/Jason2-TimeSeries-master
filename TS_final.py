# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 12:48:22 2018

@author: Daixin
"""
import main
#from netCDF4 import Dataset
import utm
import numpy as np
import pandas as pd
import datetime

#header, data = main.read_jason2_PH_nc(Dataset('JA2_IPH_2PTP001_020_20080712_190809_20080712_200421.nc','r'))
vlat = 8.0002
vlon = 7.749
SR = 60000
#
#vlat = 65.902
#vlon = -129.045
#SR = 1000000
E1, N1, _, _ = utm.from_latlon(vlat, vlon)
TmSri = [] # Time Series
TmSric = [] # Time Series corrected

def time_series():
#    mat, header = JA2_PH_crt(60,65,-125,-120,228)
#    mat, header = JA2_PH_crt(8.95,9.2,-1.23,-1,46)
    mat, header = main.JA2_PH_crt(7.9,8.1,7.6,7.9,20)
#    mat, header = JA2_PH_crt(65,66,-130,-127,228)
    header_used = []
    for n in range(len(header)):
        if header[0][n] == 1:
            header_used.append(header['Variables'][n])

    if mat:
        global TmSri # modify global var TmSri, 
        #avoid 'UnboundLocalError: local variable referenced before assignment'
        for i in range(0, len(mat)):
            s, _ = mat[i].shape
            if s > 0:
                for j in range(0, s):
                    E2, N2, _, _ = utm.from_latlon(mat[i][j][2], mat[i][j][1])
                    Dist = np.sqrt(np.square(N1 - N2) + np.square(E1 - E2))
                    if Dist <= SR:
                        if len(mat[i][j]) == 179:
                            TmSri.append(mat[i][j])
    TmSri = np.array(TmSri)
     
    if TmSri.size:
        global TmSric
        rowone = [TmSri[0]]
        TmSric = pd.DataFrame(np.concatenate((rowone, TmSri), axis=0), columns=header_used)
        TmSric = TmSric.drop_duplicates(subset = ['time','lat','lon'])
  
    return TmSric, header_used, mat
#    return TmSri, header_used, mat, header
TmSric, header_used, mat = time_series()
#%%
## def geoid_height(): # this might take 3 - 5 minutes.
#import scipy 
#import scipy.io as spio
#EGM = spio.loadmat('EGM2008_5.mat')
#XX = EGM['XX']
#YY = EGM['YY']
#GH = EGM['GH']
#XX = XX.flatten()
#YY = YY.flatten()
#GH = GH.flatten()
#points = np.column_stack((XX,YY))
#
#if TmSric.size:
#    lat = TmSric['lat'].as_matrix()
#    lon = TmSric['lon'].as_matrix()
##    if np.where(lon > 180):
##        lon = lon - 360
#
#    print('Interpolating griddata from EGM2008_5...')
#    print('This might take a while...')
#    height = scipy.interpolate.griddata(points, GH, (lon, lat), method='linear')
#    TmSric['Geoid height'] = pd.Series(height, index = TmSric.index)
#    print('Interpolation terminated')
#    return

#%%

if TmSric.size:
    m, _ = TmSric.shape
    TmSric = TmSric.reset_index(drop=True)

    Alt_1hz,Range_ku,Range_c,Range_oce3_ku,Range_oce3_c=(np.zeros((m,1)) for x in range(5))
    Range_red3_c,Range_ice3_ku,Range_ice3_c,Range_red3_ku,Correction_ku=(np.zeros((m,1)) for x in range(5))
    Cor = np.zeros((m,8)) # Correction parameters matrix
    for i in range(0, m):
        
        # Altitude, 1Hertz
        Alt_1hz[i,0] = TmSric['alt'][i]
        
        ### 1Hertz, Corrections ###
        
#         Inverse barometric correction
        if TmSric['inv_bar_corr'][i] != 32767 and TmSric['inv_bar_corr'][i] != -32768 and TmSric['inv_bar_corr'][i] != 2.147483647000000e+09:
            InvBar_ku = TmSric['inv_bar_corr'][i]
        else:
            InvBar_ku = 0

#         Sea state bias
        if TmSric['sea_state_bias_ku'][i] != 32767 and TmSric['sea_state_bias_ku'][i] != -32768 and TmSric['sea_state_bias_ku'][i] != 2.147483647000000e+09:
            SeSbias_ku = TmSric['sea_state_bias_ku'][i]
        else:
            SeSbias_ku = 0
#        if TmSric['sea_state_bias_ku'][i] != -32768 and TmSric['sea_state_bias_ku'][i] != 2.147483647000000e+09:
#            SeSbias_ku = TmSric['sea_state_bias_ku'][i]
##        elif TmSric['sea_state_bias_ku'][i] >= 32767:
##            SeSbias_ku = 3.2767
#        else:
#            SeSbias_ku = 0

        # Ionospheric correction
        if TmSric['iono_corr_alt_ku'][i] != 32767 and TmSric['iono_corr_alt_ku'][i] != -32768 and TmSric['iono_corr_alt_ku'][i] != 2.147483647000000e+09:
            IonCor_ku = TmSric['iono_corr_alt_ku'][i]
        else:
            IonCor_ku = 3.2767
#            IonCor_ku = 0

        # Ocean tide
        if TmSric['ocean_tide_sol1'][i] != 32767 and TmSric['ocean_tide_sol1'][i] != -32768 and TmSric['ocean_tide_sol1'][i] != 2.147483647000000e+09:
            OcTide_ku = TmSric['ocean_tide_sol1'][i]
        else:
            OcTide_ku = 0

        # Pole Tide
        if TmSric['pole_tide'][i] != 32767 and TmSric['pole_tide'][i] != -32768 and TmSric['pole_tide'][i] != 2.147483647000000e+09:
            PoTide_ku = TmSric['pole_tide'][i]
        else:
            PoTide_ku = 0

        # Earth tide
        if TmSric['solid_earth_tide'][i] != 32767 and TmSric['solid_earth_tide'][i] != -32768 and TmSric['solid_earth_tide'][i] != 2.147483647000000e+09:
            ETide_ku = TmSric['solid_earth_tide'][i]
        else:
            ETide_ku = 0

        # Wet tropospheric correction
        if TmSric['model_wet_tropo_corr'][i] != 32767 and TmSric['model_wet_tropo_corr'][i] != -32768 and TmSric['model_wet_tropo_corr'][i] != 2.147483647000000e+09:
            WTCor_ku = TmSric['model_wet_tropo_corr'][i]
#        if TmSric['rad_wet_tropo_corr'][i] != 32767 and TmSric['rad_wet_tropo_corr'][i] != -32768 and TmSric['rad_wet_tropo_corr'][i] != 2.147483647000000e+09:
#            WTCor_ku = TmSric['rad_wet_tropo_corr'][i]
        else:
            WTCor_ku = 0

        # Dry tropospheric correction
        if TmSric['model_dry_tropo_corr'][i] != 32767 and TmSric['model_dry_tropo_corr'][i] != -32768 and TmSric['model_dry_tropo_corr'][i] != 2.147483647000000e+09:
            DTCor_ku = TmSric['model_dry_tropo_corr'][i]
        else:
            DTCor_ku = 0
            
        # Performing corrections
#        Correction_ku = InvBar_ku + SeSbias_ku + IonCor_ku + OcTide_ku + PoTide_ku + ETide_ku + WTCor_ku + DTCor_ku + TmSric['Geoid height'][i]
        Correction_ku=InvBar_ku+SeSbias_ku+IonCor_ku+OcTide_ku+PoTide_ku+ETide_ku+WTCor_ku+DTCor_ku+TmSric['geoid_EGM2008'][i]
        Cor[i,:] = [InvBar_ku, SeSbias_ku, IonCor_ku, OcTide_ku, PoTide_ku, ETide_ku, WTCor_ku, DTCor_ku]


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
    
    # Creating Time Matrix
    
#    TSyear,TSmonth,TSday,TSyears,TSsec,TSku_mean,TSku_median,TSku_std=([] for x in range(8))
#    TSc_mean,TSc_median,TSc_std,TSoce3ku_mean,TSoce3ku_median,TSoce3ku_std=([] for x in range(6))
#    TSoce3c_mean,TSoce3c_median,TSoce3c_std,TSred3ku_mean,TSred3ku_median,TSred3ku_std=([] for x in range(6))
#    TSred3c_mean,TSred3c_median,TSred3c_std,TSice3ku_mean,TSice3ku_median,TSice3ku_std=([] for x in range(6))

#    TSice3c_mean,TSice3c_median,TSice3c_std=([] for x in range(3))

    s = 0
    ind = [0]
    TSyear,TSmonth,TSday,TStotal,SSH_ice3_ku_mean,SSH_ice3_ku_median,SSH_ice3_ku_std=([] for x in range(7))
    y, mo, d = (np.zeros((m,1)) for x in range(3))
    a0 = datetime.timedelta(seconds=TmSric['time'][0]) + datetime.datetime(2000,1,1)
    y[0] = a0.timetuple()[0]
    mo[0] = a0.timetuple()[1]
    d[0] = a0.timetuple()[2]
#    sec = []
    for t in range(1,m):
        a = datetime.timedelta(seconds=TmSric['time'][t]) + datetime.datetime(2000,1,1)
        y[t] = a.timetuple()[0]
        mo[t] = a.timetuple()[1]
        d[t] = a.timetuple()[2]
        if d[t,0] != d[t-1,0] or t+1 == m:
            s = s + 1
            ind.append(t)
            TSyear.append(y[t-1,0])
            TSmonth.append(mo[t-1,0])
            TSday.append(d[t-1,0])

#            sec = datetime.timedelta(seconds=TmSri[t-1,0]).total_seconds()
            
            TStotal.append(datetime.timedelta(seconds=TmSric['time'][t-1]) + datetime.datetime(2000,1,1))
            
            SSH_ice3_ku_mean.append(np.nanmean(SSH_ice3_ku[ind[s-1]:ind[s]-1], axis = 0))
            SSH_ice3_ku_median.append(np.nanmedian(SSH_ice3_ku[ind[s-1]:ind[s]-1], axis=0))
            SSH_ice3_ku_std.append(np.std(SSH_ice3_ku[ind[s-1]:ind[s]-1], axis=0))
            
            
    cols = ['Date Detail','Year','Month','Day','SSH_ice3_ku_mean','SSH_ice3_ku_median','SSH_ice3_ku_std']
    TS = pd.DataFrame({'Year': TSyear, 
                       'Month': TSmonth, 
                       'Day': TSday, 
                       'Date Detail': TStotal,
                       'SSH_ice3_ku_mean': SSH_ice3_ku_mean,
                       'SSH_ice3_ku_median': SSH_ice3_ku_median,
                       'SSH_ice3_ku_std': SSH_ice3_ku_std}, columns = cols)
    
#%%    
#    TS = TS.set_index('Date Detail')
#            days = datetime.timedelta(seconds=TmSric[t-1,0]) + datetime.datetime(2000,1,1) - datetime.datetime(int(y[t-1,0]),1,1)
#            if y[t-1,0]==1996 or y[t-1,0]==2000 or y[t-1,0]==2004 or y[t-1,0]==2008 or y[t-1,0]==2012 or y[t-1,0]==2016:
#                le = 366
#            else:
#                le = 365   
#            TSyears.append(days/le+TSyear[s-1])       
    
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
        TScor['SSH_ice3_ku_median'][f_ind2[i]] = M[int(TScor['Month'][f_ind2[i]])-1,3]
                                                                            # column-2? = 3
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
