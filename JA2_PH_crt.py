# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 22:38:38 2018

@author: Daixin
"""

from tqdm import tqdm
import os
from netCDF4 import Dataset
from read_jason2_PH_nc import read_jason2_PH_nc
import numpy as np

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
                header, data = read_jason2_PH_nc(fname)
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

    return mat, header

#mat = JA2_PH_crt(60, 65, -125, -120, 228)