# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:14:25 2018

@author: Daixin
"""

import scipy.io
import numpy as np

data = scipy.io.loadmat("EGM2008_5.mat")

for i in data:
	if '__' not in i and 'readme' not in i:
		np.savetxt((i+".csv"),data[i],delimiter=',')