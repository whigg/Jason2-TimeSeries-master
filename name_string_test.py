# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 14:12:56 2018

@author: Daixin
"""
import ftplib
from tqdm import tqdm
from netCDF4 import Dataset
import os
#from read_jason2_PH_nc import read_jason2_PH_nc

#ftp = ftplib.FTP('ftp.nodc.noaa.gov', 'anonymous', 'anonymous')
##print(ftp.login()) # login as anonymous
#ftp.cwd('pub/data.nodc/jason2/gdr/s_gdr/cycle001')    # change to target directory
#d = ftp.nlst() # list the contents
#x = d[1]
#print(x)


def JA2_PH_crt(track):
    fname = []
    if track < 10:
        t1 = '0'
    else:
        t1 = ''
    if track < 100:
        t2 = '0'
    else:
        t2 =''
    for i in tqdm(range(1, 6)):

#        B = []
#        s = 0

        if i < 10:
            i1 = '0'
        else:
            i1 = ''
        if i < 100:
            i2 = '0'
        else:
            i2 = ''
                # name from NOAA ftp  # cycle number                # track number
        id_track = 'JA2_GPS_2PdP'+str(i1)+str(i2)+str(i)+'_'+str(t1)+str(t2)+str(track)+'_'
        # id_track = 'JA2_IPH_2PTP'+str(i1)+str(i2)+str(i)+'_'+str(t1)+str(t2)+str(track)+'_'
                # name from GIS database  # cycle number                # track number
        if i < 10:
            name1 = '0'
        else:
            name1 = ''
        if i < 100:
            nn = '0'
        else:
            nn = ''
        name2 = str(nn)+str(i)
        name = 'cycle'+str(name1)+str(name2) # cycle number

#        """
#        Encountered problem: cant acquire data from GIS server
#
#        # path = os.path.join('//129.69.12.99/gis/Altimetry data/JASON2/JASON_2_PH/' + name)
#        # d = os.listdir(path)
#        """
#        path = os.path.join('ftp://ftp.nodc.noaa.gov/pub/data.nodc/jason2/gdr/s_gdr/' + name)
#        d = os.listdir(path)
#        # Alternative: try to fetch directly from NOAA ftp
        ftp = ftplib.FTP('ftp.nodc.noaa.gov', 'anonymous', 'anonymous') # login as anonymous
        ftp.cwd('pub/data.nodc/jason2/gdr/s_gdr/' + name)    # change to target directory
        d = ftp.nlst() # list the contents
#        fname.append(d)
        for j in range(2, len(d)):
                    #  why ignore the first 2 ncData (bzw. Start from 3)
                    #  By python should change from 3 to 2
                    #  Tourian wrote:     for j=3:length(d)

            # fndr = re.search(id_track, d[j]) # possible bug area (re.match)
            #
            # if not len(fndr):
            #     fname = Dataset(d[j], 'r')
            if id_track in d[j]:
                fname = Dataset(d[j], 'r')
#                fname.append(d[j])
#    ftp.quit()
    return fname

out = JA2_PH_crt(24)
#%%
# input = fname
input = Dataset('./nc_data/225/cycle000/JA2_GPS_2PdP000_225_20080710_211338_20080710_220951.nc', 'r')
_, data = read_jason2_PH_nc(input)
m, n = data.shape





#%%
