# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 11:59:58 2018

@author: Daixin
"""
import numpy as np
import utm
from JA2_PH_crt import JA2_PH_crt
#import scipy 
import scipy.io as spio
import datetime
import pandas as pd
#from tqdm import tqdm

vlat = 65.4362
vlon = -127.5606
SR = 60000
E1, N1, _, _ = utm.from_latlon(vlat, vlon)
#_, header = JA2_PH_crt(60,65,-125,-120,228)
def correcting_TS():
    mat, header = JA2_PH_crt(60,65,-125,-120,228)
    header_df = []
    for n in range(len(header)):
        if header[0][n] == 1:
            header_df.append(header['Variables'][n])

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
        TmSric = pd.DataFrame(np.concatenate((rowone, TmSri), axis=0), columns=header_df)
        TmSric = TmSric.drop_duplicates(subset = ['time','lat','lon'])
    
    return TmSric
        
# def geoid_height():
EGM = spio.loadmat('EGM2008_5.mat')
#
XX = EGM['XX']
YY = EGM['YY']
GH = EGM['GH']
#XX = np.genfromtxt('EGM2008_5_XX.csv', delimiter=',')
#YY = np.genfromtxt('EGM2008_5_YY.csv', delimiter=',')
#GH = np.genfromtxt('EGM2008_5_GH.csv', delimiter=',')
TmSric = correcting_TS()
if TmSric.size:
    lat = TmSric['lat'].as_matrix()
    lon = TmSric['lon'].as_matrix()
    if np.where(lon > 180):
        lon = lon - 360
#    height = scipy.interpolate.griddata((EGM['XX'],EGM['YY']), EGM['GH'], (lon, lat), method='cubic')
#    TmSri = np.concatenate()

#def geophysical_cor(): 
if TmSric.size:
    m, _ = TmSric.shape
    TmSric = TmSric.reset_index(drop=True)
#    Alt_1hz = np.zeros((m,1))
    Alt_1hz,Range_ku,Range_c,Range_oce3_ku,Range_oce3_c=(np.zeros((m,1)) for x in range(5))
    Range_red3_c,Range_ice3_ku,Range_ice3_c,Range_red3_ku=(np.zeros((m,1)) for x in range(4))
    for i in range(0, m):
        # Altitude, 1Hertz
        Alt_1hz[i,0] = TmSric['qual_inst_corr_1hz_sig0_c'][i]
        # 1Hertz, Corrections
        # Inverse barometric correction
        if TmSric['agc_c'][i] != 32767 and TmSric['agc_c'][i] != -32768 and TmSric['agc_c'][i] != 2.147483647000000e+05:
            InvBar_ku = TmSric['agc_c'][i]
        else:
            InvBar_ku = 0

        # Sea state bias
        if TmSric['interp_flag_ocean_tide_sol2'][i] != 32767 and TmSric['interp_flag_ocean_tide_sol2'][i] != -32768 and TmSric['interp_flag_ocean_tide_sol2'][i] != 2.147483647000000e+05:
            SeSbias_ku = TmSric['interp_flag_ocean_tide_sol2'][i]
        else:
            Sesbias_ku = 0

        # Ionospheric correction
        if TmSric['interp_flag_tb'][i] != 32767 and TmSric['interp_flag_tb'][i] != -32768 and TmSric['interp_flag_tb'][i] != 2.147483647000000e+05:
            IonCor_ku = TmSric['interp_flag_tb'][i]
        else:
            InoCor_ku = 0

        # Ocean tide
        if TmSric['agc_rms_c'][i] != 32767 and TmSric['agc_rms_c'][i] != -32768 and TmSric['agc_rms_c'][i] != 2.147483647000000e+05:
            OcTide_ku = TmSric['agc_rms_c'][i]
        else:
            OcTide_ku = 0

        # Polar Tide
        if TmSric['tb_187'][i] != 32767 and TmSric['tb_187'][i] != -32768 and TmSric['tb_187'][i] != 2.147483647000000e+05:
            PoTide_ku = TmSric['tb_187'][i]
        else:
            PoTide_ku = 0

        # Earth tide
        if TmSric['off_nadir_angle_wf_ku'][i] != 32767 and TmSric['off_nadir_angle_wf_ku'][i] != -32768 and TmSric['off_nadir_angle_wf_ku'][i] != 2.147483647000000e+05:
            ETide_ku = TmSric['off_nadir_angle_wf_ku'][i]
        else:
            ETide_ku = 0

        # Wet tropospheric correction
        if TmSric['ecmwf_meteo_map_avail'][i] != 32767 and TmSric['ecmwf_meteo_map_avail'][i] != -32768 and TmSric['ecmwf_meteo_map_avail'][i] != 2.147483647000000e+05:
            WTCor_ku = TmSric['ecmwf_meteo_map_avail'][i]
        else:
            WTCor_ku = 0

        # Dry tropospheric correction
        if TmSric['orb_state_flag_diode'][i] != 32767 and TmSric['orb_state_flag_diode'][i] != -32768 and TmSric['orb_state_flag_diode'][i] != 2.147483647000000e+05:
            DTCor_ku = TmSric['orb_state_flag_diode'][i]
        else:
            DTCor_ku = 0
            
        # Performing corrections
#        Correction_ku = InvBar_ku + SeSbias_ku + IonCor_ku + OcTide_ku + PoTide_ku + ETide_ku + WTCor_ku + DTCor_ku + TmSric(i,180)

        Correction_ku=InvBar_ku+SeSbias_ku+IonCor_ku+OcTide_ku+PoTide_ku+ETide_ku+WTCor_ku+DTCor_ku

#        Cor[i,:] = [InvBar_ku, SeSbias_ku, IonCor_ku, OcTide_ku, PoTide_ku, ETide_ku, WTCor_ku, DTCor_ku]

        Range_ku[i, 0] = TmSric['qual_rad_1hz_tb238'][i] + Correction_ku
        Range_c[i, 0] = TmSric['qual_rad_1hz_tb340'][i] + Correction_ku
        Range_oce3_ku[i, 0] = TmSric['rad_averaging_flag'][i] + Correction_ku
        Range_oce3_c[i, 0] = TmSric['rad_land_frac_187'][i] + Correction_ku
        Range_red3_ku[i, 0] = TmSric['rad_land_frac_238'][i] + Correction_ku
        Range_red3_c[i, 0] = TmSric['rad_land_frac_340'][i] + Correction_ku
        Range_ice3_ku[i, 0] = TmSric['alt_state_flag_oper'][i] + Correction_ku
        Range_ice3_c[i, 0] = TmSric['alt_state_flag_c_band'][i] + Correction_ku
        
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
#    TSyear = []
#    TSmonth = []
#    TSday = []
    s = 0
    ind = [0]
    TSyear,TSmonth,TSday,TStotal,SSH_ice3_ku_median,SSH_ice3_ku_std=([] for x in range(6))
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
            ind.append(t)
            TSyear.append(y[t-1,0])
            TSmonth.append(mo[t-1,0])
            TSday.append(d[t-1,0])
            
            s = s + 1
#            sec = datetime.timedelta(seconds=TmSri[t-1,0]).total_seconds()
            TStotal.append(datetime.timedelta(seconds=TmSric['time'][t-1]) + datetime.datetime(2000,1,1))
            SSH_ice3_ku_median.append(np.nanmedian(SSH_ice3_ku[ind[s-2]], axis=0))
            SSH_ice3_ku_std.append(np.std(SSH_ice3_ku[ind[s-2]], axis=0))
    cols = ['Date Detail','Year','Month','Day','SSH_ice3_ku_median','SSH_ice3_ku_std']
    TS = pd.DataFrame({'Year': TSyear, 
                       'Month': TSmonth, 
                       'Day': TSday, 
                       'Date Detail': TStotal, 
                       'SSH_ice3_ku_median': SSH_ice3_ku_median,
                       'SSH_ice3_ku_std': SSH_ice3_ku_std}, columns = cols)
#    TS = TS.set_index('Date Detail')
#            days = datetime.timedelta(seconds=TmSric[t-1,0]) + datetime.datetime(2000,1,1) - datetime.datetime(int(y[t-1,0]),1,1)
#            if y[t-1,0]==1996 or y[t-1,0]==2000 or y[t-1,0]==2004 or y[t-1,0]==2008 or y[t-1,0]==2012 or y[t-1,0]==2016:
#                le = 366
#            else:
#                le = 365   
#            TSyears.append(days/le+TSyear[s-1])       

#def TStcor(TS):
#_, n = TS.shape
TS = TS.drop_duplicates(subset=['Year','Month','Day'], keep=False)
#test = np.nanmean(TS['SSH_ice3_ku_median'])
# ?? nanmean 'SSH_ice3_ku_median' and 'SSH_ice3_ku_std' ??

          

#def Res_data(data):  change TS to data after finish the function!!
m, n = TS.shape
residual = np.zeros((m, n))
if not TS.size:
    f1 = TS.index[TS['Month'] == 1].tolist(); f2 = TS.index[TS['Month'] == 2].tolist() 
    f3 = TS.index[TS['Month'] == 3].tolist(); f4 = TS.index[TS['Month'] == 4].tolist()
    f5 = TS.index[TS['Month'] == 5].tolist(); f6 = TS.index[TS['Month'] == 6].tolist()
    f7 = TS.index[TS['Month'] == 7].tolist(); f8 = TS.index[TS['Month'] == 8].tolist()
    f9 = TS.index[TS['Month'] == 9].tolist(); f10 = TS.index[TS['Month'] == 10].tolist()
    f11 = TS.index[TS['Month'] == 11].tolist(); f12 = TS.index[TS['Month'] == 12].tolist()
    Jan1 = np.zeros((1, n-2)); Feb1 = np.zeros((1, n-2)); Mar1 = np.zeros((1, n-2))
    Apr1 = np.zeros((1, n-2)); May1 = np.zeros((1, n-2)); Jun1 = np.zeros((1, n-2))
    Jul1 = np.zeros((1, n-2)); Aug1 = np.zeros((1, n-2)); Sep1 = np.zeros((1, n-2))
    Oct1 = np.zeros((1, n-2)); Nov1 = np.zeros((1, n-2)); Dec1 = np.zeros((1, n-2))
    M_Data = np.zeros((12, n-2))
    
    for i in range(0,n):
        Jan1[0, i] = 
    M_Data = [Jan1,Feb1,Mar1,Apr1,May1,Jun1,Jul1,Aug1,Sep1,Oct1,Nov1,Dec1]
    
    residual[:,0] = TS['Year']
    residual[:,1] = TS['Month']
    for i in range(0, m):
#        residual[i,2:n] = TS.iloc[i,2:n] - M_Data[TS['Month'][i]]

        Mts[i, 1] = M_Data[TS['Month'][i]]

else:
    for i in range(0, 11):
        M_Data[i,1] = i
#return residual, M_Data, Mts



#def Outcor(TS):
TScor = TS
level = TS['qual_inst_corr_1hz_range_ku_mle3']    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
 


       