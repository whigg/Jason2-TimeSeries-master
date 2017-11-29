## This function process satellite altimetry data for a given coordinate to
## obtain water level time series
##
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


from netCDF4 import Dataset
import netCDF4
import numpy as np
import os
from tqdm import tqdm



def netcdf_reader():
    


    return ncDaata, ncid



def read_jasion2_PH_nc():
    
    
    return

DATA_DIR = '/home/geodasie/Desktop/Data/'
'''                                         225, 228'''
def JA2_PH_crt(latmin, latmax, lonmin, lonmax, track):
    mat = []
    if track < 10:
        t1 = '0'
    else:
        t1 = ''
    if track < 100:
        t2 = '0'
    else:
        t2 =''

    k = 1
    for i in tqdm(range(1, 320)):
        B = []
        s = 1

        if i < 10:
            i1 = '0'
        else:
            i1 = ''
        if i < 100:
            i2 = '0'
        else:
            i2 = ''

        id_track = 'JA2_cycle'+str(i1)+str(i2)+str(i)+'_track'+str(t1)+str(t2)+str(track)

        if i < 10:
            name1 = '0'
        else:
            name1 =''
        if i < 100:
            nn = '0'
        else:
            nn = ''
        name2 = str(nn)+str(i)
        name = 'cycle'+str(name1)+str(name2)
        path = os.path.join('/Desktop/Data/'+name)
        d = os.listdir(path)
        for j in range(3, len(d)):
            fndr = re.match(, id_track)
            if not len(fndr):
                fname =
                [header, data] = read_jason2_PH_nc(fname)

        
    return mat
##
######################### MAIN FUNCTION HERE ######################
##'''
##def Altprocess(vlat, vlon, SR, latmin, latmax, lonmin, lonmax, track):   
####    loading EGM 2008 provide height from Geoid    
##    if lonmin < 0:
##        lonmin = 360 + lonmin
##    else:
##        pass
##    if lonmax < 0:
##        lonmax = 360 + lonmax
##    else:
##        pass
####    mat = JA2_PH_crt(latmin, latmax, lonmin, lonmax, track)
##    return mat, TS_final, Cor
##'''
####################################################################
##finding 20 Hz measurements inside the virtual station
    
def ell2utm(vlat, vlon):
     return N1, E1
    
def finding20Hz():
    if not len(seq):
        pp = 1
        ic = 0
        TmSri = []
        for i in range(m):
             [s1, s2] = np.shape(mat[i, 2])
             if s1 > 1:
                 for j in range(s1):
                     N2, E2 = ell2utm(mat[i, 2][j, 3]*pi/180, mat[i, 2][j, 2]*pi/180)
                     Dist = np.sqrt(np.square(N1 - N2) + np.square(E1 - E2))
                     if Dist <= SR:
                         if len(mat[i, 2][j, :]) == 179:
                             TmSri(pp, :) = mat[i, 2][j, :]
                         
                         ic = 1
                     if ic == 1:
                         pp := pp + 1
                     ic = 0
    return

def correcting_TS():
    if not len(TmSri):
        m, n = np.shape(TmSri)
        g = 2
        TmSric(1, :) = TmSri(1, :)
        for i in range(2, m):
            f = 
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
        TmSric(:, 180) = griddata(#points, values, xi, method='linear', fill_value=nan, rescale=False)

    return
##if isempty(TmSric)==0
##lat=TmSric(:,3)
##lon=TmSric(:,2)
##f=find(lon>180)
##lon(f)=lon(f)-360
##TmSric(:,180) = griddata(XX(:),YY(:),GH(:),lon,lat)
##return

## geopysical correction:
