# ---------------------------- #
# Post processing binary files #
# ---------------------------- #


import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from glob import glob
import sys,os
import numpy as np

import datetime as dt

from netCDF4 import Dataset
import xarray as xr


# Load coordinates & mask
# -----------------------
ds = xr.open_dataset('coord.nc')
lon = np.array(ds['longitude'])
lat = np.array(ds['latitude'])
depth = np.array(ds['depth'])
ds.close()

ds = Dataset('mask.nc','r')
mask = np.array(ds['tmask'])
mask = np.array(mask,dtype=float)

ds.close()

sz,sy,sx = mask.shape

mask[np.where(mask==0)] = np.nan
mask[np.where(mask==1)] = 0 


# Load variables definitions
# --------------------------
vlist = xls = pd.ExcelFile(r"var_to_output.xls") 
sheet = xls.parse(0)

odir = 'NC'


# Get Arguments
# -------------
idir = sys.argv[1]
odir = sys.argv[2]
dini = sys.argv[3]
freq = float(sys.argv[4])
nday = float(sys.argv[5])


date_ini = dt.datetime.strptime(dini,'%Y%m%d')



# Get initial N iteration
f = sorted(glob('P1l.*.data'))[0]
n_ini = f.replace('P1l.','')
n_ini = n_ini.replace('.data','')


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


    # Loop on date tag
    for day in range(0,nday):
   
      full_data = []

      dnow = date_ini + dt.timedelta(days=day)
      jd = dnow.toordinal()
      dnow = dnow.strftime('%Y%m%d')

      # Loop on hours
      for hour in range(0,24):

         print(hour)

         n = int(float(n_ini)+hour*freq)
         n = str(n).zfill(10)

         # Special case chlorophyll
         if vname == 'chl':
           data = np.zeros(sz*sy*sx)
           for p in range(1,5):
              data += np.fromfile('P'+str(p)+'l.'+str(n)+'.data',dtype='float32')

         # Special case Net primary production
         elif vname == 'nppv':
           ruPPYc = np.fromfile('ruPPYc.'+str(n)+'.data',dtype='float32')
           resPYc = np.fromfile('resPPYc.'+str(n)+'.data',dtype='float32')
           exR2acP1 = np.fromfile('exR2acP1.'+str(n)+'.data',dtype='float32')
           exR2acP2 = np.fromfile('exR2acP2.'+str(n)+'.data',dtype='float32')
           exR2acP3 = np.fromfile('exR2acP3.'+str(n)+'.data',dtype='float32')
           exR2acP4 = np.fromfile('exR2acP4.'+str(n)+'.data',dtype='float32')

           data = ruPPYc - resPPYc - exR2acP1 - exR2acP2  - exR2acP3 - exR2acP4

         # Special case zooplankton
         elif vname == 'zooc':
           data = np.zeros(sz*sy*sx)
           for i in range(3,7):
             data += np.fromfile('Z'+str(i)+'c.'+str(n)+'.data',dtype='float32')


         # Generic case
         else:
           # Load data
           data = np.fromfile(mname+'.'+str(n)+'.data',dtype='float32')

         data = np.reshape(data,(sz,sy,sx))
         data = data + mask # Faster then np.where

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
        dtime = dataset.createDimension('time',24)
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
        vdepth.valid_min  = 0.75
        vdepth.valid_max  = 2927.323
        vdepth.long_name = "Depth"
        vdepth[:] = depth

        vtime = dataset.createVariable('time',np.float32,('time'))
        vtime.units = "hours since 1900-01-01 00:00:00"
        dori = dt.datetime(1900,1,1).toordinal()
        hori = (jd - dori )*24

        time = range(hori,hori+24)
        vtime[:] = time

      else:
 
        # Open file
        dataset = Dataset(fname,'a',format='NETCDF4_CLASSIC')

      # Write variable to NetCDF
      # ------------------------
      vdata = dataset.createVariable(vname,np.float32,('time','depth','latitude','longitude'),fill_value=-9999)
      vdata.units = units
      vdata.long_name = lname
      vdata.standard_name = sname

      if info != 'none':
         vdata.info = info
  
      vdata[:] = full_data
      dataset.close()

      exit()

