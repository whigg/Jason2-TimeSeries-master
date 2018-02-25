# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:13:45 2018

@author: wosleepy
"""

import numpy as np
import pandas as pd
from netCDF4 import Dataset

# Read NetCDF3 data, 'r' -- readable only
"""
 careful with fname. TO DO #
"""
#fname = Dataset('./nc_data/225/cycle000/JA2_GPS_2PdP000_225_20080710_211338_20080710_220951.nc', 'r')

def read_jason2_PH_nc(fname):
    # Set lists for later use in forming dataframe 'ncData'
    Variables, Attname, Attval, DataMatrix, MatrixShape = [], [], [], [], []

    #####################################################################
    # Loop through all variables in netCDF data, get corresponding Attribute Units(Attval), Data Matrix(DataMatrix)
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
#  Solved by check m = shape_matrix[0] & n != 20

    data_row_num = ncData.at[0, 'MatrixShape'][0] # Number of row in final data matrix
    data = np.empty([data_row_num,138]) # create empty space for final data matrix
    header_usage = [] # creat empty space for recording variable usage
    s = 0 # set for write in final data matrix iteration

    for i in range(1, 180): # Tourian didn’t say why we only take the these 179 variables
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


#def Altprocess():
#    x = Dataset('./nc_data/225/cycle000/JA2_GPS_2PdP000_225_20080710_211338_20080710_220951.nc', 'r')
#    y1, y2 = read_jason2_PH_nc(x)
#    print(y1.shape, y2.shape)


#a = 1
#Altprocess()
