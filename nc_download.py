from netCDF4 import Dataset
import netCDF4
from ftplib import FTP
import ftplib
import os
##from download_SGDR import download_ftp_tree

# connect to noaa ftp #

##ftp = FTP('ftp.nodc.noaa.gov')
##print(ftp.login())  # print login status
##ftp.cwd('pub/data.nodc/jason2/gdr/s_gdr/cycle009')
##print(ftp.retrlines('LIST')) # retrieve list
##download_ftp_tree(ftp, 'pub/data.nodc/jason2/gdr/s_gdr/', '/home/geodasie/Desktop/Data/')

# NetCDF download #

##version 1
##filenames = ftp.nlst()
##print(filenames)
##
##for filename in filenames:
##    local_filename = os.path.join('/home/geodasie/Desktop/Data/', filename)
##    file = open(local_filename, 'wb')
##    ftp.retrlines('RETR' + filename, file.write)
##
##    file.close()
##
##ftp.quit() # go gentle


##version 2
##os.chdir('/home/geodasie/Desktop/Data/')
##ls = ftp.nlst()
##count = len(ls)
##curr = 0
##print('found {} files'.format(count))
##
##for fn in ls:
##    curr  += 1
##    print('Processing file {} ... {} of {} ...'.format(fn, curr, count))
##    ftp.retrlines('RETR' + fn, open(fn, 'wb').write)
##
##ftp.quit()
##print('Download complete')

"""
##version 3   find specific file
##filematch = '*_228_*.nc'  # e.g. track number 228
##target_dir = '/home/geodasie/Desktop/Data'
##
##for filename in ftp.nlst(filematch):
##    target_file_name = os.path.join(target_dir, os.path.basename(filename))
##    with open(target_file_name, 'wb') as fhandle:
##        ftp.retrbinary('RETR %s' % filename, fhandle.write)
##
"""

##    maybe download and extract cifar10

##    dest_directory = FLAGS.data_dir
##    if not os.path.exists(data_directory):
##        os.makedirs(dest_directory)
##    filename = DATA_URL.split('/')[-1]
##    filepath = os.path.join(dest_directory, filename)
##    if not os.path.exists(filepath):
##        def _progress(count, block_size, total_size):
##            sys.stdout.write('\r>> Downloading %s %.1f%%' % (filename,
##                                                             float(count * block_size) / float(total_size) * 100.0))
##            sys.stdout.flush()
##        filepath, _ = urllib.request.urlretrieve(DATA_URL, filepath, _progress)
##        print()
##        statinfo = os.stat(filepath)
##        print('Successfully download', filename, statinfo.st_size, 'bytes.')
##    extracted_dir_path = os.path.join(dest_directory, 'S_GDR_cycle***')
##    if not os.path.exists(extracted_dir_path):
##        tarfile.open(filepath, 'r:gz').extractall(dest_directory)
##

'''
## track number: 225 228
## download netCDF from 225 228 from FTP
'''

## data properties check
'''                                               JA2_GPS_2PdP #cycle_#track_ *.nc   '''
##data = Dataset('/home/geodasie/Desktop/Data/225/JA2_GPS_2PdP000_225_20080710_211338_20080710_220951.nc')
##print(data.file_format)
##print(data.dimensions.keys()) # '# of measurments in the file, # of elementary measurements, # of waveform samples'
##print(data.dimensions['wvf_ind'])
##print(data.Conventions)
##print(data.variables['surface_type'])
