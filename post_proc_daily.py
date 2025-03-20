# ---------------------------- #
# Post processing binary files #
# ---------------------------- #

print('------------------')
print('STARTING POST-PROC')
print('------------------')


import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from glob import glob
import sys,os
import numpy as np

import datetime as dt

from netCDF4 import Dataset
import xarray as xr

print('here')

# Get Arguments
# -------------
idir = sys.argv[1]
odir = sys.argv[2]
dini = sys.argv[3]
freq = float(sys.argv[4])
nday = int(sys.argv[5])
fout = sys.argv[6]
fmask = sys.argv[7]
fcoord = sys.argv[8]

date_ini = dt.datetime.strptime(dini,'%Y%m%d')


# Load coordinates & mask
# -----------------------
ds = xr.open_dataset(fcoord)
lon = np.array(ds['longitude'])
lat = np.array(ds['latitude'])
depth = np.array(ds['depth'])
ds.close()

ds = Dataset(fmask,'r')
mask = np.array(ds['tmask'])
mask = np.array(mask,dtype=float)

ds.close()

sz,sy,sx = mask.shape

mask[np.where(mask==0)] = np.nan
mask[np.where(mask==1)] = 0 



# Load variables definitions
# --------------------------
xls = pd.ExcelFile(fout) 
sheet = xls.parse(0)


# Get initial N iteration
f = sorted(glob(idir+'/P1l.*.data'))[0]
f = os.path.basename(f)
n_ini = f.replace('P1l.','')
n_ini = n_ini.replace('.data','')

print('n_ini',n_ini)


print('Processing')
# Loop on variables
# -----------------
for i in range(0,sheet['Model name'].shape[0]):

 # Output variable
 if sheet['Output'][i] == 1:

    # Get meta data
    mname = sheet['Model name'][i]
    vname = sheet['NetCDF name'][i]
    ftag = sheet['File tag'][i]
    lname = sheet['Long name'][i]
    sname = sheet['Standard name'][i]
    units = sheet['Units'][i]
    info = sheet['Info'][i]

    print(vname)

    # Loop on date tag
    # ----------------
    for day in range(0,nday):

      full_data = []

      dnow = date_ini + dt.timedelta(days=day)
      jd = dnow.toordinal()
      dnow = dnow.strftime('%Y%m%d')

      # Loop on hours
      # -------------
      for hour in range(0,1):

         n = int(day*freq*24+float(n_ini)+hour*freq)
         n = str(n).zfill(10)

         # -> Special case chlorophyll
         if vname == 'chl':
           data = np.zeros(sz*sy*sx)
           for p in range(1,5):
              data += np.fromfile(idir+'/P'+str(p)+'l.'+str(n)+'.data',dtype='float32')

         # -> Special case phytoplankton as carbon
         elif vname == 'phyc':
           data = np.zeros(sz*sy*sx)
           for p in range(1,5):
              data += np.fromfile(idir+'/P'+str(p)+'c.'+str(n)+'.data',dtype='float32')

         # -> Special case Net primary production
         elif vname == 'nppv':
           ruPPYc = np.fromfile(idir+'/ruPPYc.'+str(n)+'.data',dtype='float32')
           resPPYc = np.fromfile(idir+'/resPPYc.'+str(n)+'.data',dtype='float32')
           exR2acP1 = np.fromfile(idir+'/exR2acP1.'+str(n)+'.data',dtype='float32')
           exR2acP2 = np.fromfile(idir+'/exR2acP2.'+str(n)+'.data',dtype='float32')
           exR2acP3 = np.fromfile(idir+'/exR2acP3.'+str(n)+'.data',dtype='float32')
           exR2acP4 = np.fromfile(idir+'/exR2acP4.'+str(n)+'.data',dtype='float32')

           data = ruPPYc - resPPYc - exR2acP1 - exR2acP2  - exR2acP3 - exR2acP4

         # -> Special case zooplankton
         elif vname == 'zooc':
           data = np.zeros(sz*sy*sx)
           for i in range(3,7):
             data += np.fromfile(idir+'/Z'+str(i)+'c.'+str(n)+'.data',dtype='float32')
             data = data/12.


         # -> Generic case
         else:
           # Load data
           data = np.fromfile(idir+'/'+mname+'.'+str(n)+'.data',dtype='float32')

         # Reshape 2D
         if vname != 'mld':
           data = np.reshape(data,(sz,sy,sx))
           data = data + mask # Faster than np.where

         # Reshape 3D
         else:
           data = np.reshape(data,(sy,sx))
           data = data + mask[0,:,:] 
           data = data.squeeze()

         full_data.append(data)

      full_data = np.array(full_data)

      # Write NetCDF file
      # -----------------
      fname = odir+'/'+dnow+'_h-OGS--'+ftag+'-MITgcmBFM-pilot8-b'+dini+'_fc-v01.nc'


      # Initialize if no exist
      if not os.path.exists(fname):

        # Create file
        dataset = Dataset(fname,'w',format='NETCDF4_CLASSIC')

        # Create dimension
        dtime = dataset.createDimension('time',1)
        ddepth = dataset.createDimension('depth',sz)
        dlat = dataset.createDimension('latitude',sy)
        dlon = dataset.createDimension('longitude',sx)


        # Create variables & set attributes
        vlon = dataset.createVariable('longitude',np.float32,('longitude'))
        vlon.units = "degrees east"
        vlon.standard_name = "longitude"
        vlon.long_name = "longitude"
        vlon.valid_min = np.amin(lon)
        vlon.valid_max = np.amax(lon)
        vlon.axis = "X"
        vlon[:] = lon

        vlat = dataset.createVariable('latitude',np.float32,('latitude'))
        vlat.units = "degrees north"
        vlatstandard_name = "longitude"
        vlat.long_name = "longitude"
        vlat.valid_min = np.amin(lat)
        vlat.valid_max = np.amin(lat)
        vlat.axis = "Y"
        vlat.long_name = "latitude"
        vlat[:] = lat

        vdepth = dataset.createVariable('depth',np.float32,('depth'))
        vdepth.standard_name = "depth"
        vdepth.positive = "down"
        vdepth.axis = "Z"
        vdepth.units = "meters"
        vdepth.valid_min  = np.amin(depth)
        vdepth.valid_max  = np.amax(depth)
        vdepth.long_name = "Depth"
        vdepth[:] = depth

        vtime = dataset.createVariable('time',np.float32,('time'))
        vtime.units = "hours since 1900-01-01 00:00:00"
        dori = dt.datetime(1900,1,1).toordinal()
        hori = (jd - dori )*24

        time = hori
        vtime[:] = time

      else:
 
        # Open file
        dataset = Dataset(fname,'a',format='NETCDF4_CLASSIC')

      # Write variable to NetCDF
      # ------------------------

      # 3D
      if vname != 'mld':
        vdata = dataset.createVariable(vname,np.float32,('time','depth','latitude','longitude'),zlib=True)
      # 2D
      else:
        vdata = dataset.createVariable(vname,np.float32,('time','latitude','longitude'),zlib=True)


      # Attributes
      vdata.units = units
      vdata.long_name = lname
      vdata.standard_name = sname

      if info != 'none':
         vdata.info = info
  
      vdata[:] = full_data
      dataset.close()



print('-------------')
print('POST-PROC END')
print('-------------')

