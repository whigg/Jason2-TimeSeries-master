## This function process satellite altimetry data for a given coordinate to
## obtain water level time series
##
# function [TS, mat,TS_final,Cor] = Altprocess(vlat,vlon,SR,latmin, latmax, lonmin, lonmax, track)

##  input:
##       vlat: Virtual station's latitude [deg]
##       vlon: virtual station's longitude [deg]
##       SR: search radius for virtual station [m]
##       ofndr: applying off nadir correction 1, not applying 0
##       latmin: minimum latitude [deg]
##       latmax: maximum latitude  [deg]
##       lonmin: minimum longitude [deg]
##       lonmax: maximum longitude [deg]
##       track: track number []
##
##  output:
##       TS_final:
##               Column 1: year
##               Column 2: month
##               Column 3: surface water height from combined retracking
##               approach
##               Column 4: error of estimated water level

##        Time:
##               Time of measurements
##        mat: processed mat file within the defined box
#%%
from netCDF4 import Dataset  # handle netCDF files
import netCDF4
import numpy as np  # array processing
import os  # deal with directory
from tqdm import tqdm  # add processing bar
import utm
import ftplib

#%%

#################  Old version ########################
#def JA2_PH_crt(latmin, latmax, lonmin, lonmax, track):
#    mat = []
#    if track < 10:
#        t1 = '0'
#    else:
#        t1 = ''
#    if track < 100:
#        t2 = '0'
#    else:
#        t2 =''
#    k = 0
#    for c in tqdm(range(1, 327)):
#        B = []
#        s = 0
#        if c == 304:
#            c = c + 1  # NOAA ftp has no cycle304
#        if c < 10:
#            c1 = '0'
#        else:
#            c1 = ''
#        if c < 100:
#            c2 = '0'
#        else:
#            c2 = ''
#                # name from NOAA ftp  # cycle number                # track number
#        id_track = 'JA2_GPS_2PdP'+str(c1)+str(c2)+str(c)+'_'+str(t1)+str(t2)+str(track)+'_'
#            # name from GIS database  # cycle number                # track number
#       # id_track = 'JA2_IPH_2PTP'+str(i1)+str(i2)+str(i)+'_'+str(t1)+str(t2)+str(track)+'_'
#        name = 'cycle'+str(c1)+str(c2)+str(c) # cycle number
#        path = os.path.join('D:/jason2/gdr/s_gdr' + name)
#        d = os.listdir(path)
#
#        # """
#        # Encountered problem: cant acquire data from GIS server
#        #
#        # # path = os.path.join('//129.69.12.99/gis/Altimetry data/JASON2/JASON_2_PH/' + name)
#        # # d = os.listdir(path)
#
#        # Alternative: try to fetch from NOAA ftp
#        # But cannot read single *.nc directly from FTP
#        # ftp = ftplib.FTP('ftp.nodc.noaa.gov', 'anonymous', 'anonymous') # login as anonymous
#        # ftp.cwd('pub/data.nodc/jason2/gdr/s_gdr/' + name)    # change to target directory
#        # d = ftp.nlst() # list the contents
#        # """
#
#        for i in range(2, len(d)):

#            if id_track in d[i]:
#                fname = Dataset(d[i], 'r')
#                _, data = read_jason2_PH_nc(fname)
#                m, n = data.shape
#                for j in range(0,m):
#                    if data[j,3]<latmax and data[j,3]>latmin and data[j,2]<lonmax and data[j,2]>lonmin:
#                        B[s,:] = data[j,:]
#                        s = s + 1
#
#                if not B:
#                    mat[k, 1] = i
#                    mat[k, 2] = B
#                    k = k + 1
#
#    ftp.quit() # close connection with FTP
#    return mat


from tqdm import tqdm
import os
from netCDF4 import Dataset
from read_jason2_PH_nc import read_jason2_PH_nc
#import numpy as np

def JA2_PH_crt(latmin, latmax, lonmin, lonmax, track):
    mat = []
    B = []
    i_res = []
    B_res = []

    if track < 10:
        t1 = '0'
    else:
        t1 = ''
    if track < 100:
        t2 = '0'
    else:
        t2 =''
    for c in tqdm(range(0, 327)):
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
        for i in range(0, len(d)):
            if id_track in d[i]:
                fname = Dataset(data_dir+d[i], 'r')
                _, data = read_jason2_PH_nc(fname)
                m, _ = data.shape
                for j in range(0,m):
                    #  why ignore the first 2 ncData (bzw. Start from 3)
                    #  By python should change from 3 to 2, i set it as 0
                    #  Tourian wrote:     for j=3:length(d)
                    if data[j,1]<latmax and data[j,1]>latmin and data[j,2]<lonmax and data[j,2]>lonmin:
                        # Lat varies from around 80 to 248, Lon varies from around -66 to 66
                        # Tourian might got it wrong in column number of data matrix
                        # The second column should be Lat, third should be Lon
                        # in python data[:,1] = lat, data[:,2] = lon
                        B.append(data[j,:])
                if B:
                    i_res.append(c)
                    B_res.append(B)
#    B = np.array(B)
    mat = [i_res, B_res]
#    return mat, B
    return mat
#x = JA2_PH_crt(65, 66, -128, -127, 225)
#z = x[1]
#z1 = z[1]
#z11 = z[1][1]
#%%

##finding 20 Hz measurements inside the virtual station

E1, N1, _, _ = utm.from_latlon(vlat*pi/180, vlon*pi/180)
# utm.from_latlon(Latitude, longitude)
# return (Easting, Northing, # zone, zone letter)
m, n = mat.shape

def correcting_TS():
    if not len(mat):
        pp = 1
        # ic = 0
        TmSri = []
        for i in range(1, m):
             s1, s2 = np.shape(mat[i, 2]) # mat[i, 2].shape
             if s1 > 1:
                 for j in range(1, s1):  # possible range(0, s1)
                                            # mat{i,2}(j,3)          mat{1,2}(j,2)
                     E2, N2, _, _ = utm.from_latlon(mat[i, 2][j, 3]*pi/180, mat[i, 2][j, 2]*pi/180)
                     Dist = np.sqrt(np.square(N1 - N2) + np.square(E1 - E2))
                     if Dist <= SR:
                         #  length(mat{i,2}(j,:))
                         if len(mat[i, 2][j, :]) == 179:
                             TmSri[pp, :] = mat[i, 2][j, :]

                         pp = pp + 1
                     #     ic = 1
                     # if ic == 1:
                         # pp = pp + 1

                     # ic = 0


        if not len(TmSri):
            m, n = np.shape(TmSri) # TmSri.shape
            g = 2
            TmSric[1, :] = TmSri(1, :)
            for i in range(2, m):
                '''    ???????????
                    whats the output of f and how to set the condition'''
                f = [f for f if ]
                # [f for f in a if z > 2]

                if len(f) == 0:
                    TmSric[g,:] = TmSri[i,:]
                    g = g + 1
        else:
            TmSric = []
    else:
        TmSric = []

    return TmSric
##[N1, E1] = ell2utm(vlat*pi/180, vlon*pi/180)
##[m n] = size(mat)
##if isempty(mat)==0
##        pp=1;ic=0;
##
##        TmSri=[];
##        for i=1:m
##            [s1 s2]=size(mat{i,2});
##            if s1>1
##
##                for j=1:s1
##                    [N2,E2]=ell2utm(mat{i,2}(j,3)*pi/(180),mat{i,2}(j,2)*pi/(180));
##                    Dist=sqrt((N1-N2)^2+(E1-E2)^2);
##                    if Dist<=SR
##                        if length(mat{i,2}(j,:))==179
##                            TmSri(pp,:)=mat{i,2}(j,:);
##                        end
##
##                        ic=1;
##                    end
##                    if ic==1
##                        pp=pp+1;
##                    end
##                    ic=0;
##                end
##            end
##        end
##
##
##
##        %% correcting the time series
##        if isempty(TmSri)==0
##            [m n]=size(TmSri);
##            g=2;
##            TmSric(1,:)=TmSri(1,:);
##            for i=2:m
##                f=find(TmSric(:,1)==TmSri(i,1) & TmSric(:,2)==TmSri(i,2) & TmSric(:,3)==TmSri(i,3));
##                if isempty(f)==1
##                    TmSric(g,:)=TmSri(i,:);
##                    g=g+1;
##                end
##            end
##        else
##            TmSric=[];
##        end
##    else
##        TmSric=[];
##end
##




##Geoid height computation
from scipy.interpolate import griddata
def Geoid_height():
    if not len(TmSric):
        lat = TmSric(:, 3)
        lon = TmSric(:, 2)
        f = find(lon > 180)
        lon(f) = lon(f) - 360
        '''
                 ???????               why grid this one'''
        TmSric(:, 180) = griddata(#points, values, xi, method='linear', fill_value=nan, rescale=False)

    return
##if isempty(TmSric)==0
##lat=TmSric(:,3)
##lon=TmSric(:,2)
##f=find(lon>180)
##lon(f)=lon(f)-360
##TmSric(:,180) = griddata(XX(:),YY(:),GH(:),lon,lat)
##return

## geopysical correction
def geophysical_correction():
    if not len(TmSric):
        m, n = TmSric.shape

        for i in range (1, m):
            # altitude, 1Hertz
            Alt_1hz[i, 1] = TmSric[i, 26]
            # 20Hertz
            # Alt_20hz[i, 1:20] = TmSric[i, 154:173]/(10**4)

            # Ocean Range  %%%%% Ku %%%%%
            # 1 Hertz
            # Corrections
            # Inverse barometric correction
            if TmSric[i, 99] != 32767 and TmSric[i, 99] != -32768 and TmSric[i, 99] != 2.147483647000000e+05:
                InvBar_ku = TmSric[i, 99]
            else:
                InvBar_ku = 0

            # Sea state bias
            if TmSric[i, 51] != 32767 and TmSric[i, 51] != -32768 and TmSric[i, 51] != 2.147483647000000e+05:
                SeSbias_ku = TmSric[i, 51]
            else:
                Sesbias_ku = 0

            # Ionospheric correction
            if TmSric[i, 47] != 32767 and TmSric[1, 47] != -32768 and TmSric[i, 47] != 2.147483647000000e+05:
                IonCor_ku = TmSric[i, 47]
            else:
                InoCor_ku = 0

            # Ocean tide
            if TmSric[i, 101] != 32767 and TmSric[i, 101] != -32768 and TmSric[i, 101] != 2.147483647000000e+05:
                OcTide_ku = TmSric[i, 101]
            else:
                OcTide_ku = 0

            # Polar Tide
            if TmSric[i, 110] != 32767 and TmSric[i, 110] != -32768 and TmSric[i, 110] != 2.147483647000000e+05:
                PoTide_ku = TmSric[i, 110]
            else:
                PoTide_ku = 0

            # Earth tide
            if TmSric[i, 109] != 32767 and TmSric[i, 109] != -32768 and TmSric[i, 109] != 2.147483647000000e+05:
                ETide_ku = TmSric[i, 109]
            else:
                ETide_ku = 0

            # Wet tropospheric correction
            if TmSric[i, 42] != 32767 and TmSric[i, 42] != -32768 and TmSric[i, 42] != 2.147483647000000e+05:
                WTCor_ku = TmSric[i, 42]
            else:
                WTCor_ku = 0

            # Dry tropospheric correction
            if TmSric[i, 40] != 32767 and TmSric[i, 40] != -32768 and TmSric[i, 40] != 2.147483647000000e+05:
                DTCor_ku = TmSric[i, 42]
            else:
                DTCor_ku = 0

            # Performing corrections
            Correction_ku = InvBar_ku + SeSbias_ku + IonCor_ku + OcTide_ku +
                            PoTide_ku + ETide_ku + WTCor_ku + DTCor_ku + TmSric(i,180)
            Cor[i,:] = [InvBar_ku, SeSbias_ku, IonCor_ku, OcTide_ku, PoTide_ku, ETide_ku, WTCor_ku, DTCor_ku]

            Range_ku[i, 1] = TmSric[i, 28] + Correction_ku
            Range_c[i, 1] = TmSric[i, 29] + Correction_ku
            Range_oce3_ku[i, 1] = TmSric[i, 30] + Correction_ku
            Range_oce3_c[i, 1] = TmSric[i, 31] + Correction_ku
            Range_red3_ku[i, 1] = TmSric[i, 32] + Correction_ku
            Range_red3_c[i, 1] = TmSric[i, 33] + Correction_ku
            Range_ice3_ku[i, 1] = TmSric[i, 34] + Correction_ku
            Range_ice3_c[i, 1] = TmSric[i, 35] + Correction_ku

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
        s = 1
        ind = 1

        a =
        y =
        mo =
        d =
        for i in range(2, m):
            a =
            y[i, 1] =
            mo[i, 1] =
            d[i, 1] =
            if d[i, 1] != d[i-1, 1] or i == m:
                s = s + 1
                ind[s, 1] = i
                TS[s - 1, 1] = y[i - 1, 1]
                TS[s - 1, 2] = mo[i - 1, 1]
                TS[s - 1, 3] = d[i - 1, 1]

                days =
                if y[i-1, 1] == 1996 or y[i-1, 1] == 2000 or y[i-1, 1] == 2004 or y[i-1, 1] == 2008 or y[i-1, 1] == 2012 or y[i-1, 1] == 2016:
                    le = 366
                else:
                    le = 365

                TS[s-1, 4] = days / le + TS[s-1, 1]
                TS[s-1, 5] = TmSric[i-1, 1]/86400 + data


    #     '''
    #     % Creating Time matrix
    #     s=1;
    #     ind=1;
    #
    #     a=datestr(TmSric(1,1)/86400+datenum(2000,1,1),'yyyymmdd');
    #     y=str2num(a(1:4));
    #     mo=str2num(a(5:6));
    #     d=str2num(a(7:8));
    #     for i=2:m
    #         a=datestr(TmSric(i,1)/86400+datenum(2000,1,1),'yyyymmdd');
    #         y(i,1)=str2num(a(1:4));
    #         mo(i,1)=str2num(a(5:6));
    #         d(i,1)=str2num(a(7:8));
    #         if d(i,1)~=d(i-1,1) | i==m
    #             s=s+1;
    #             ind(s,1)=i;
    #             TS(s-1,1)=y(i-1,1);
    #             TS(s-1,2)=mo(i-1,1);
    #             TS(s-1,3)=d(i-1,1);
    #
    #             days=TmSric(i-1,1)/86400+datenum(2000,1,1)-datenum(y(i-1,1),1,1);
    #             if y(i-1,1)==2008 || y(i-1,1)==1996 ||y(i-1,1)==2000 ||y(i-1,1)==2004 ||y(i-1,1)==2012 || y(i-1,1)==2016
    #                 le=366;
    #             else
    #                 le=365;
    #             end
    #             TS(s-1,4)=days/le+TS(s-1,1);
    #             TS(s-1,5)=TmSric(i-1,1)/86400+datenum(2000,1,1);
    #
    #
    #
    #             TS(s-1,6)=nanmean(SSH_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,7)=nanmedian(SSH_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,8)=std(SSH_ku(ind(s-1):ind(s)-1));
    #
    #             TS(s-1,9)=nanmean(SSH_c(ind(s-1):ind(s)-1));
    #             TS(s-1,10)=nanmedian(SSH_c(ind(s-1):ind(s)-1));
    #             TS(s-1,11)=std(SSH_c(ind(s-1):ind(s)-1));
    #
    #             TS(s-1,12)=nanmean(SSH_oce3_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,13)=nanmedian(SSH_oce3_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,14)=std(SSH_oce3_ku(ind(s-1):ind(s)-1));
    #
    #             TS(s-1,12)=nanmean(SSH_oce3_c(ind(s-1):ind(s)-1));
    #             TS(s-1,13)=nanmedian(SSH_oce3_c(ind(s-1):ind(s)-1));
    #             TS(s-1,14)=std(SSH_oce3_c(ind(s-1):ind(s)-1));
    #
    #             TS(s-1,15)=nanmean(SSH_red3_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,16)=nanmedian(SSH_red3_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,17)=std(SSH_red3_ku(ind(s-1):ind(s)-1));
    #
    #             TS(s-1,15)=nanmean(SSH_red3_c(ind(s-1):ind(s)-1));
    #             TS(s-1,16)=nanmedian(SSH_red3_c(ind(s-1):ind(s)-1));
    #             TS(s-1,17)=std(SSH_red3_c(ind(s-1):ind(s)-1));
    #
    #             TS(s-1,18)=nanmean(SSH_ice3_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,19)=nanmedian(SSH_ice3_ku(ind(s-1):ind(s)-1));
    #             TS(s-1,20)=std(SSH_ice3_ku(ind(s-1):ind(s)-1));
    #
    #             TS(s-1,21)=nanmean(SSH_ice3_c(ind(s-1):ind(s)-1));
    #             TS(s-1,22)=nanmedian(SSH_ice3_c(ind(s-1):ind(s)-1));
    #             TS(s-1,23)=std(SSH_ice3_c(ind(s-1):ind(s)-1));
    #
    #         end
    #     end
    #     TS=TStcor(TS);
    #     TS_t=Outcor(TS,19,2.9,1);
    #     s=find(TS_t(:,20)>0.7);
    #     s0=find(TS_t(:,20)<0.7);
    #     TS_t(s,20)=nanmean(TS_t(s0,20));
    #     TS_final=[TS_t(:,5) TS_t(:,19:20)];
    # else
    #     TS_final=[];
    #     TS=[];
    # end
    #     '''







##
######################### MAIN FUNCTION HERE ######################
##'''
def Altprocess(vlat, vlon, SR, latmin, latmax, lonmin, lonmax, track):
##    loading EGM 2008 provide height from Geoid ????
   if lonmin < 0:
       lonmin = 360 + lonmin
   # else:
   #     pass
   if lonmax < 0:
       lonmax = 360 + lonmax
   # else:
   #     pass
   mat = JA2_PH_crt(latmin, latmax, lonmin, lonmax, track)
   # return mat, TS_final, Cor
##'''
####################################################################
