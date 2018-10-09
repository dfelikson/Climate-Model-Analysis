#!/usr/bin/env python

import os, sys, argparse, datetime
#from subprocess import call
import numpy as np
#np.seterr(invalid='ignore') # ignore invalid value warnings
from netCDF4 import Dataset
#from pyproj import Proj

#from scipy.interpolate import interp2d
from scipy import ndimage

#import string

from RACMOgridAndStatsGlobals import *
from RasterClipperFunctions import *
from RACMOutilities import *

import time
import progressbar

import pickle

#sys.path.append('/Users/dfelikso/Research/Software/ScriptsAndUtilities/pythonModules')
#import raster

#from uncertainties import unumpy as unp

ncfilename = '/Users/dfelikso/Research/Data/RACMO/RACMO2.3/originalData/runoff_WJB_int.2015.BN_1958_2013_1km.DD.nc'
clipfile = '/Users/dfelikso/Research/Projects/ISMIP6/GrIS_Tidewater_Basins/GrIS_Tidewater_basins.mat.shp'
attributefilter = 'basin=*'
attributefilter = 'basin=527'

oversampleGrid = False
oversampleZoom = 2.

print(" ")
print("\033[93mWARNING! This script uses the exterior of each polygon! \033[0m") 
print(" ")

# Read NetCDF file{{{
print("reading netcdf file: " + ncfilename)
ncfile = Dataset(ncfilename, 'r')
t = ncfile.variables['time']
x = ncfile.variables['x'][:]
y = ncfile.variables['y'][:]
runoff = ncfile.variables['runoffcorr']

geoTransform = [x[0], x[1]-x[0], 0., y[-1], 0., y[0]-y[1]]
if oversampleGrid:
   geoTransform[1] = geoTransform[1] / oversampleZoom
   geoTransform[5] = geoTransform[5] / oversampleZoom
#}}}

# Create basin mask {{{
if not checkForField(clipfile, attributefilter.split('=')[0]):
   print('ERROR in ' + __file__ + ': Attribute "' + attributefilter.split('=')[0] + '" not found in ' + clipfile)
   sys.exit()

print("creating basin masks")
driver = ogr.GetDriverByName("ESRI Shapefile")
dataSource = driver.Open(clipfile, 0)
layer = dataSource.GetLayer()

if '*' in attributefilter:
   attributeValues = getUniqueFieldValues(clipfile, attributefilter.split('=')[0])
else:
   attributeValues = [attributefilter.split('=')[1]]

if oversampleGrid:
   runoff_nrows = int(runoff.shape[1]*oversampleZoom)
   runoff_ncols = int(runoff.shape[2]*oversampleZoom)
else:
   runoff_nrows = int(runoff.shape[1])
   runoff_ncols = int(runoff.shape[2])

maskArray = np.zeros( (len(attributeValues), runoff_nrows, runoff_ncols) )

for iValue, attributeValue in enumerate(attributeValues):
   AF = attributefilter.split('=')[0] + '=' + str(attributeValue)
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
runoffDict = dict()
runoffDict['time'] = t[:].data
for attributeValue in attributeValues:
   runoffDict[attributeValue] = np.empty( len(t) )

# Mask at each timestep for all basins {{{
print("masking at each timestep for each basin")
#for itime in progressbar.progressbar( range(0, runoff.shape[0]) ):
for itime in range(0, runoff.shape[0]):
   #sys.stdout.write('{:03.0f} '.format(t[itime]))
   #sys.stdout.flush()

   if oversampleGrid:
      runoff_timeSlice = ndimage.zoom(runoff[itime,:,:], oversampleZoom)
   else:
      runoff_timeSlice = runoff[itime,:,:]
   
   runoffSum = 0.
   for iValue, attributeValue in enumerate(attributeValues):
      runoffSum = np.sum( maskArray[iValue,:,:] * runoff_timeSlice)
      runoffDict[attributeValue][itime] = runoffSum
      #sys.stdout.write('{:16.3f} '.format(runoffSum))
      #sys.stdout.flush()

   #sys.stdout.write('\n')
   #sys.stdout.flush()
#}}}

import pdb; pdb.set_trace()
pickle.dump(runoffDict, open('runoffDict', 'wb'))

