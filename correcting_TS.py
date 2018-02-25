# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 11:59:58 2018

@author: Daixin
"""
import numpy as np
from math import pi
import utm
from JA2_PH_crt import JA2_PH_crt
from tqdm import tqdm

vlat = 65.4362
vlon = -127.5606
SR = 6000
E1, N1, _, _ = utm.from_latlon(vlat*pi/180, vlon*pi/180)
#def correcting_TS():    
mat = JA2_PH_crt(60,65,-125,-120,228)
mat_data = mat[1]
TmSri = []
TmSric = []
if mat:
    for i in tqdm(range(0, len(mat_data))):
        s = len(mat_data[i])
        if s > 0:
            for j in range(0, s):
#                x = mat_data[i][j][1]
                E2, N2, _, _ = utm.from_latlon(mat_data[i][j][1]*pi/180, mat_data[i][j][2]*pi/180)
                Dist = np.sqrt(np.square(N1 - N2) + np.square(E1 - E2))
                if Dist <= SR:
                    if len(mat_data[i][j]) == 138:
                        TmSri.append(mat_data[i][j])
            
    TmSri = np.array(TmSri)
    if TmSri.size:
        m, _ = TmSri.shape
        g = 2
        TmSric.append(TmSri[0])
        for i in range(1, m):
            TmSric[]