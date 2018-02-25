# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 22:38:38 2018

@author: Daixin
"""

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
#    for c in tqdm(range(0, 327)):
    for c in tqdm(range(300, 327)):
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
                    if data[j,1]<latmax and data[j,1]>latmin and data[j,2]<lonmax and data[j,2]>lonmin:
#                for i in range(0,m):
#                    if data[i,1]<latmax and data[i,1]>latmin and data[i,2]<lonmax and data[i,2]>lonmin:
                        # Tourian might got wrong with lat%lon column?!?
                        B.append(data[i,:])
                if B:
                    i_res.append(c)
                    B_res.append(B)
#    B = np.array(B)
    mat = [i_res, B_res]
#    return mat, B
    return mat

#mat = JA2_PH_crt(65,66,-128,-127,225)
#x = JA2_PH_crt(65, 66, -128, -127, 225)
#x, y = JA2_PH_crt(82, 83, -67, -66, 225)
#z = x[1]
#z1 = z[1]
#z11 = z[1][1]
