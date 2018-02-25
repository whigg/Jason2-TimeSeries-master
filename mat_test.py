# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 22:55:32 2018

@author: Daixin
"""
import os
from netCDF4 import Dataset
from read_jason2_PH_nc import read_jason2_PH_nc
import pandas as pd
import numpy as np

latmin = 82
latmax = 83
lonmin = -67
lonmax = -66
track = 225
B = []

mat = []
k = 0


i_res = []
B_res = []


dir = 'D:/jason2/gdr/s_gdr/cycle000/'
id_track = 'JA2_GPS_2PdP000_225_'
path = os.path.join('D:/jason2/gdr/s_gdr/cycle000')
d = os.listdir(path)
for i in range(0, len(d)):
#    if id_track in d[i]:
#        x = d[i]
#        fname = Dataset(dir+d[i], 'r')
#        print(fname.file_format)
    if id_track in d[i]:
        fname = Dataset(dir+d[i], 'r')
        _, data = read_jason2_PH_nc(fname)
        m, n = data.shape
        for j in range(0,m):
            if data[j,2]<latmax and data[j,2]>latmin and data[j,1]<lonmax and data[j,1]>lonmin:
                B.append(data[j,:])



        B = np.array(B)
        if B.size:
            mat
#            i_res.append(i)
#            B_res.append(B)

#            mat[k,0] = i
#            mat[k,1] = B
#            k = k + 1
mat = pd.DataFrame({'i': i_res, 'B': B_res})
