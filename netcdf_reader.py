from netCDF4 import Dataset
import netCDF4
import numpy as np
##import pandas as pd

data = Dataset('JA2_GPS_2PdP000_225_20080710_211338_20080710_220951.nc','r')
ncData = []
ncData.append(list(data.variables.keys())) # first row of ncData

for i in data.variables:
#     print(i)
    attval = []
    try:
        attval = data.variables[i].units
    except AttributeError:
        attval = None

    for j in [attval]:
        attname = []
        if j == None:
            attname = None
        else:
            attname = 'units'
        print([i, attname, attval])

print(len(data.variables.keys()))
for i in range(1, len(data.variables.keys())):
    sc = 1
    off = 0

    for j in range():

for i in data.variables:
    print(data.variables[i].shape)
# #     print(i)
#     shape = []
#     try:
#         sc = data.variables[i].scale_factor

#     except AttributeError:
#         sc = 1

#     try:
#         off = data.variables[i].add_offset
#     except AttributeError:
#         off = 0
#     print([sc, off])
#     a = data.variables[i].shape * sc + off
#     print(shape)

#     for j in [attval]:
#         attname = []
#         if j == None:
#             attname = None
#         else:
#             attname = 'units'

#         try:
#         attval = data.variables[i].scale_factor
#     except AttributeError:
#         attval = None
#     print(attval)

##print(data.variables.values())
##for i in data.variables:
##    print(i, data.variables[i].units, data.variables[i].shape)
##    d = np.array(data.variables['pole_tide'], dtype=type(data.variables))
##    print(d[:, 26, 36])

##print(data.variables.values())
##print(data.variables.type())


for i in data.variables:
    sc = 1  # scale_factor
    off = 0   # add_offset

    if attname == 'units':
        ncData[i, 2] = attname
        ncData[i, 3] = attval


    if attname == 'scale_factor':
        sc = attval


    if attname == 'add_offset':
        off = attval


##    mat = pd.DataFrame(np.random.randn(), index = [i], columns = [data.variables[i].units])
    ncData[i, 1] = i
    attname = None
    ncData[i, 4] = len(data.variables[i].shape)
##    mat = [i, data.variables[i].units, data.variables[i].shape]
    print(attval)



'''
<class 'netCDF4._netCDF4.Variable'>
int32 alt(time)
    _FillValue: 2147483647
    long_name: 1 Hz altitude of satellite
    standard_name: height_above_reference_ellipsoid
    units: m
    quality_flag: orb_state_flag_rest or orb_state_flag_diode
    add_offset: 1300000.0
    scale_factor: 0.0001
    coordinates: lon lat
    comment: Altitude of satellite above the reference ellipsoid. Associated quality flag is orb_state_flag_diode for the OGDR products, orb_state_flag_rest for the IGDR and GDR products
unlimited dimensions:
current shape = (2963,)
filling off
'''
