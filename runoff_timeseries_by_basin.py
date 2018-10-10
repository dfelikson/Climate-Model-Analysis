#!/usr/bin/env python

import os, sys, argparse, datetime
#from subprocess import call
import numpy as np
#np.seterr(invalid='ignore') # ignore invalid value warnings
from netCDF4 import Dataset
#from pyproj import Proj

#from scipy.interpolate import interp2d
from scipy import ndimage
from scipy.io import savemat

#import string

from RACMOgridAndStatsGlobals import *
from RasterClipperFunctions import *
from RACMOutilities import *

import time
#import progressbar

import pickle

#sys.path.append('/Users/dfelikso/Research/Software/ScriptsAndUtilities/pythonModules')
#import raster

#from uncertainties import unumpy as unp

# --- Setup ---
# Setup: RACMO
ncfilename = '/disk/staff/gcatania/polar/Arctic/data/RACMO/RACMO2.3/originalData/downscaled/runoff.1958-2017.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc'
clipfile = '/home/student/denis/ISMIP6/GrIS_Tidewater_Basins/GrIS_Tidewater_basins.mat.shp'
basinfilter = 'basin=*'
basinfilter = 'basin=527'

time_var = 'time';
x_var = 'lon';
y_var = 'lat';
runoff_var = 'runoffcorr'

# Setup: MAR

# Setup: general
oversampleGrid = False
oversampleZoom = 2.


# --- Processing ---
# A warning about how this script works
print(" ")
print("\033[93mWARNING! This script uses the exterior of each polygon! \033[0m") 
print(" ")

# Read NetCDF file{{{
print("reading netcdf file: " + ncfilename)
ncfile = Dataset(ncfilename, 'r')
t = ncfile.variables['time']
x = ncfile.variables[x_var][:]
y = ncfile.variables[y_var][:]
runoff = ncfile.variables['runoffcorr']

geoTransform = [x[0], x[1]-x[0], 0., y[-1], 0., y[0]-y[1]]
if oversampleGrid:
   geoTransform[1] = geoTransform[1] / oversampleZoom
   geoTransform[5] = geoTransform[5] / oversampleZoom
#}}}

# Create basin mask {{{
if not checkForField(clipfile, basinfilter.split('=')[0]):
   print('ERROR in ' + __file__ + ': Attribute "' + basinfilter.split('=')[0] + '" not found in ' + clipfile)
   sys.exit()

print("creating basin masks")
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(clipfile, 0)
layer = dataSource.GetLayer()

if '*' in basinfilter:
   basinValues = getUniqueFieldValues(clipfile, basinfilter.split('=')[0])
else:
   basinValues = [basinfilter.split('=')[1]]

if oversampleGrid:
   runoff_nrows = int(runoff.shape[1]*oversampleZoom)
   runoff_ncols = int(runoff.shape[2]*oversampleZoom)
else:
   runoff_nrows = int(runoff.shape[1])
   runoff_ncols = int(runoff.shape[2])

maskArray = np.zeros( (len(basinValues), runoff_nrows, runoff_ncols) )

for iValue, basinValue in enumerate(basinValues):
   AF = basinfilter.split('=')[0] + '=' + str(basinValue)
   layer.SetAttributeFilter(AF)
   for feature in layer:
      geometry = feature.GetGeometryRef()
      ring = geometry.GetGeometryRef(0)
      xClip = list()
      yClip = list()
      for ixy, xy in enumerate(ring.GetPoints()):
         xClip.append(xy[0])
         yClip.append(xy[1])
         
      maskArray_feature = clipImage(np.ones( (runoff_nrows, runoff_ncols) ), xClip, yClip, geoTransform)
      maskArray_feature[np.isnan(maskArray_feature)] = 0.
      maskArray[iValue,:,:] = maskArray[iValue,:,:] + maskArray_feature
#}}}

# Initialize output dict
#class runoffClass:
#   time = dict()
#   runoff = dict()
#rC = runoffClass()
runoffDict = dict()
runoffDict['time']   = {'time': t[:], 'units': t.units}
runoffDict['runoff'] = {'units': runoff.units}

# Mask at each timestep for all basins {{{
print("masking at each timestep for each basin")
tStart = time.time()
for iValue, basinValue in enumerate(basinValues):
   runoffDict['runoff'][basinValue] = np.empty( len(t) )
   #for itime in progressbar.progressbar( range(0, runoff.shape[0]) ):
   for itime in range(0, runoff.shape[0]):
      if oversampleGrid:
         runoff_timeSlice = ndimage.zoom(runoff[itime,:,:], oversampleZoom)
      else:
         runoff_timeSlice = runoff[itime,:,:]
      
      runoffSum = 0.
      runoffSum = np.sum( maskArray[iValue,:,:] * runoff_timeSlice)
      runoffDict['runoff'][basinValue][itime] = runoffSum

elapsed = time.time() - tStart
print('Elapsed: {:5.2f} sec'.format(elapsed))
#}}}

print("Saving to pickle and matfile")
pickle.dump(runoffDict, open('runoffDict.p', 'wb'))
savemat('runoffDict.mat', runoffDict)

