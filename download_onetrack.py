# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 17:27:52 2018

@author: Daixin Zhao

=============================================================================

Download executable for Jason2 s_gdr.

This executable is used to download a single track of altimetry *.nc data
to target directory from NOAA FTP.

Example usage:
    Input Parameters:
        target_dir = '<User's directory>'
        cycle_start = starting cycle number
        cycle_end = ending cycle number
        track = a single track number

    Output:
        Download to the target directory + '/jason2/gdr/s_gdr/'
        The subdirectories are same as the directories in FTP

"""

import ftplib
import os
from tqdm import tqdm # add processing bar


def download_onetrack(target_dir = None, cycle_start = 0, cycle_end = 327, track = None):
    if track < 10:
        t1 = '0'
    else:
        t1 = ''
    if track < 100:
        t2 = '0'
    else:
        t2 =''
    filematch = '*_'+str(t1)+str(t2)+str(track)+'_*.nc'# find files match sinlge number of track
    for i in tqdm(range(cycle_start, cycle_end)):
        if i == 304:  # in noaa database have missing cycle 304, ignore
            i = i + 1
        if i < 10:
            i1 = '0'
        else:
            i1 = ''
        if i < 100:
            i2 = '0'
        else:
            i2 = ''
        name = 'cycle'+str(i1)+str(i2)+str(i) # number of cycle
        if not os.path.exists(target_dir+'/jason2/gdr/s_gdr/'+name):
            os.makedirs(target_dir+'/jason2/gdr/s_gdr/'+name) # make directory
        save_dir = target_dir+'/jason2/gdr/s_gdr/'+name
        ftp=ftplib.FTP('ftp.nodc.noaa.gov', 'anonymous', 'anonymous')
        ftp.cwd('pub/data.nodc/jason2/gdr/s_gdr/' + name) # connect to ftp
        for filename in ftp.nlst(filematch):
            target_filename = os.path.join(save_dir, os.path.basename(filename))
            with open(target_filename, 'wb') as fhandle:
                ftp.retrbinary('RETR %s' % filename, fhandle.write) # download files
        ftp.quit() # close ftp connection
    return
