#!/bin/env python
# -----------
# Cheat sheet
# -----------
# To generate time series of downsampled runoff, spatially-integrated by catchment:
# python RACMOgridAndStats.py runoff_downscaled --nointerpolate --spatialresolution 1000 --clipfile '/disk/staff/gcatania/polar/Arctic/data/GIMP_DEM/drainagebasins/gimpbasinspoly.shp' --attributefilter 'basin=1.0' --temporalresolution daily --startdate 2014 01 01 --enddate 2014 12 31 --stats

# To generate a TIFF of the SMB anomaly from the downscaled product over the entire ice sheet:
# /home/student/denis/CentralWestGrISGlaciers/Analysis/climate/RACMO/RACMOgridAndStats.py smb_downscaled --nointerpolate --spatialresolution 1000 --anomaly --meanstartdate 1971 01 01 --meanenddate 1988 12 31 --startdate 1985 08 01 --enddate 2014 07 31 --temporalresolution monthly

# Interpolation and spatial units
#  1) Input in lat/lon
#     a) Output as lat/lon
#           Output units will not be area-weighted.
#     b) Output as projected
#           Project and multiply by the grid size of the projection (user-specified).
#  2) Input as projected
#     a) Output same as input projection
#           Multiply by the grid size of the projection (user-specified).
#     b) Output in a different projection
#           Not supported yet

# -------------
# Configuration
# -------------
# Debug flag
debug = False
    
# --------------
# Import modules
# --------------
import os, sys, argparse, datetime
from subprocess import call
import numpy as np
np.seterr(invalid='ignore') # ignore invalid value warnings
from netCDF4 import Dataset
from pyproj import Proj

from scipy.interpolate import interp2d

import string

from RACMOgridAndStatsGlobals import *
from RasterClipperFunctions import *
from RACMOutilities import *

sys.path.append('/home/student/denis/ScriptsAndUtilities/pythonModules')
import raster

from uncertainties import unumpy as unp

# For debugging
import csv

# ----------------------
# Command line arguments
# ----------------------
parser = argparse.ArgumentParser()
# RACMO variable and files
parser.add_argument('variable', type=str, help='RACMO variable desired',
               choices=['snowfall_daily', 'runoff_daily', 'runoff_downscaled', 'smb_daily', 'smb_monthly', 'smb_downscaled'])
parser.add_argument('--RACMOnetCDFdirectory', type=str, help='RACMO directory', default=None)

# Coordinates to project onto
parser.add_argument('--nointerpolate', action='store_true', help='do not interpolate; use "x" and "y" variables in ncfile as coordinates')
parser.add_argument('--coordinates', type=str, help='coordinates to project data', default='stere',
                     choices=['latlon','stere','utm'])
parser.add_argument('--spatialresolution', type=float, help='resolution of spatial coordinates in projection specified',
                     default=10000.0)
parser.add_argument('--interpmethod', type=str, help='interpolation method', default='griddata',
                     choices=['griddata','kriging'])
# Bounding box or clipfile
parser.add_argument('--boundingbox', type=four_floats, help='bounding box for catchment in selected coordinates [xmin xmax ymin ymax]')
parser.add_argument('--clipfile', type=str, help='file containing clip path', default='')
parser.add_argument('--attributefilter', type=str, help='attribute filter used to create a region from several polygons', default='')
parser.add_argument('--icemask', action='store_true', help='use the RACMO ice mask to select values only over ice')
parser.add_argument('--mask', type=str, help='use a user-specified mask to select values only over ice')
parser.add_argument('--maskval', type=float, help='value indicating ice in the mask', default=1.)

# Output
parser.add_argument('--outputdir', type=str, help='output directory', default='.')
parser.add_argument('--regionname', type=str, help='name of region for output file naming', default='regional')
parser.add_argument('--stats', action='store_true', help='output spatially integrated values in CSV file')
parser.add_argument('--tiffs', action='store_true', help='output maps of variable in GeoTIFF format')
parser.add_argument('--pngs', action='store_true', help='output maps of variable in PNG format')
parser.add_argument('--matfile', action='store_true', help='output MAT file with grids')
parser.add_argument('--csvfile', action='store_true', help='output CSV file of variable at original point locations')
parser.add_argument('--errortiff', action='store_true', help='output GeoTIFF with uncertainties (currently only works for SMB anomaly calculation)')
# Error on SMB is assumed to be 10% of the SMB
smb_error_percentage = 0.1

# Start/end times and time step
parser.add_argument('--temporalresolution', type=str, help='desired temporal resolution', default='weekly',
                     choices=['daily', 'weekly', 'monthly', 'yearly'])
parser.add_argument('--startdate', type=int, nargs=3, help='start year, month, day', default=[2000, 1, 1])
parser.add_argument('--enddate',   type=int, nargs=3, help='end year, month, day (inclusive)', default=[2000, 1, 31])
parser.add_argument('--outputlast', action='store_true', help='force output of last summation')

# Special calculations
parser.add_argument('--anomaly', action='store_true', help='calculate anomaly from a mean')
parser.add_argument('--meanstartdate', type=int, nargs=3, help='start year, month, day for calculation of mean field',
      default=[1961, 1, 1])
parser.add_argument('--meanenddate', type=int, nargs=3, help='end year, month, day for calculation of mean field',
      default=[1990, 12, 31])

args = parser.parse_args()

if args.RACMOnetCDFdirectory is not None:
   RACMOnetCDFdirectory = args.RACMOnetCDFdirectory
else:
   if 'downscaled' in args.variable:
      RACMOnetCDFdirectory = RACMOnetCDFdirectory + '/downscaled'

startyear = args.startdate[0]
startmonth = args.startdate[1]
startday = args.startdate[2]

endyear = args.enddate[0]
endmonth = args.enddate[1]
endday = args.enddate[2]

if args.anomaly:
   meanstartyear = args.meanstartdate[0]
   meanstartmonth = args.meanstartdate[1]
   meanstartday = args.meanstartdate[2]

   meanendyear = args.meanenddate[0]
   meanendmonth = args.meanenddate[1]
   meanendday = args.meanenddate[2]

# ----------
# Processing
# ----------
# Delete MAT file
if args.matfile:
   filename = args.outputdir + '/' + args.regionname + '_racmo' + '.mat'
   if os.path.isfile(filename):
      os.remove(filename)
   
# Read RACMO netCDF file(s)
print "reading RACMO netCDF file(s)"
if args.variable == 'snowfall_daily':
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_snowfall_daily_1991-2000.nc','r')
   variable = ncfile.variables['snowfall'][:,0,:,:]
   date_bnds = ncfile.variables['date_bnds'][:]
   start_date_str = str(date_bnds.min())
   year0 = int(start_date_str[0:4])
   month0 = int(start_date_str[4:6])
   day0 = int(start_date_str[6:8])
   epochdt = datetime.datetime(year0,month0,day0,0,0,0)

   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_snowfall_daily_2001-2010.nc','r')
   variable = np.vstack((variable, ncfile.variables['snowfall'][:,0,:,:]))
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_snowfall_daily_2011-2013.nc','r')
   variable = np.vstack((variable, ncfile.variables['snowfall'][:,0,:,:]))
   
   inputCoordinates = 'latlon'
   inputTemporalResolution = 'daily'
   variableShort = 'snowfall'

   units = 'kg m^-2 s^-1'
   
elif args.variable == 'runoff_daily':
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_runoff_daily_1991-2000.nc','r')
   variable = ncfile.variables['runoff'][:,0,:,:]
   date_bnds = ncfile.variables['date_bnds'][:]
   start_date_str = str(date_bnds.min())
   year0 = int(start_date_str[0:4])
   month0 = int(start_date_str[4:6])
   day0 = int(start_date_str[6:8])
   epochdt = datetime.datetime(year0,month0,day0,0,0,0)

   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_runoff_daily_2001-2010.nc','r')
   variable = np.vstack((variable, ncfile.variables['runoff'][:,0,:,:]))
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_runoff_daily_2011-2013.nc','r')
   variable = np.vstack((variable, ncfile.variables['runoff'][:,0,:,:]))

   lon = ncfile.variables['lon'][:]
   lat = ncfile.variables['lat'][:]

   inputCoordinates = 'latlon'
   inputTemporalResolution = 'daily'
   variableShort = 'runoff'

   units = 'kg m^-2 s^-1'
   
elif args.variable == 'runoff_downscaled':
   ncfilename = RACMOnetCDFdirectory + '/runoff_WJB_int.' + str(startyear) + '.BN_1958_2013_1km.DD.nc'
   ncfile = Dataset(ncfilename, 'r')
   variable = ncfile.variables['runoffcorr'][:,:,:]
   for year in range(startyear+1,endyear+1):
      ncfilename = RACMOnetCDFdirectory + '/runoff_WJB_int.' + str(year) + '.BN_1958_2013_1km.DD.nc'
      ncfile = Dataset(ncfilename, 'r')
      variable = np.vstack((variable, ncfile.variables['runoffcorr'][:,:,:]))

   year0 = int(startyear)
   month0 = int(01)
   day0 = int(01)
   epochdt = datetime.datetime(year0,month0,day0,0,0,0)

   # TBD: I don't understand what the following variables actually represent. Mins/maxes don't
   # make sense. Emailed Brice to figure these out. For now, it doesn't matter, since we'll be using the
   # x and y projected coordinates.
   lon = ncfile.variables['LON'][:]
   lat = ncfile.variables['LAT'][:]
   
   x = ncfile.variables['x'][:]
   y = ncfile.variables['y'][:]

   inputTemporalResolution = 'daily'
   variableShort = 'runoff'

   units = 'kg day^-1 m^-2'
   
elif args.variable == 'smb_daily':
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_smb_daily_1991-2000.nc','r')
   variable = ncfile.variables['smb'][:,0,:,:]
   date_bnds = ncfile.variables['date_bnds'][:]
   start_date_str = str(date_bnds.min())
   year0 = int(start_date_str[0:4])
   month0 = int(start_date_str[4:6])
   day0 = int(start_date_str[6:8])
   epochdt = datetime.datetime(year0,month0,day0,0,0,0)
   
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_smb_daily_2001-2010.nc','r')
   variable = np.vstack((variable, ncfile.variables['smb'][:,0,:,:]))
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_smb_daily_2011-2013.nc','r')
   variable = np.vstack((variable, ncfile.variables['smb'][:,0,:,:]))
   
   lon = ncfile.variables['lon'][:]
   lat = ncfile.variables['lat'][:]

   inputTemporalResolution = 'daily'
   variableShort = 'smb'
   
   units = 'kg m^-2 s^-1'
   
elif args.variable == 'smb_monthly':
   ncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_smb_monthly_1958-2013.nc','r')
   variable = ncfile.variables['smb'][:,:,:]
   #date_bnds = ncfile.variables['time'][:]
   #start_date_str = str(date_bnds.min())
   year0 = int(1958)
   month0 = int(1)
   day0 = int(1)
   epochdt = datetime.datetime(year0,month0,day0,0,0,0)
   
   lon = ncfile.variables['LON'][:]
   lat = ncfile.variables['LAT'][:]

   inputTemporalResolution = 'monthly'
   variableShort = 'smb'

   units = 'mmWE'
   
elif args.variable == 'smb_downscaled':
   if args.anomaly:
      readstartyear = np.min((startyear,meanstartyear))
      readendyear = np.max((endyear,meanendyear))
   else:
      readstartyear = startyear
      readendyear = endyear

   ncfilename = RACMOnetCDFdirectory + '/SMB_rec_WJB_int.' + str(readstartyear) + '.BN_1958_2013_1km.DD.nc'
   ncfile = Dataset(ncfilename, 'r')

   # TBD: I don't understand what the following variables actually represent. Mins/maxes don't
   # make sense. Emailed Brice to figure these out. For now, it doesn't matter, since we'll be using the
   # x and y projected coordinates.
   lon = ncfile.variables['LON'][:]
   lat = ncfile.variables['LAT'][:]
   
   x = ncfile.variables['x'][:]
   y = ncfile.variables['y'][:]

   variable = np.empty([(readendyear-readstartyear+1)*12, y.shape[0], x.shape[0]])
   for year in range(readstartyear,readendyear+1):
      print "reading year " + str(year)
      ncfilename = RACMOnetCDFdirectory + '/SMB_rec_WJB_int.' + str(year) + '.BN_1958_2013_1km.DD.nc'
      ncfile = Dataset(ncfilename, 'r')
      variableyear = ncfile.variables['SMB_rec'][:,:,:]

      startIdx = (year-readstartyear)*12
      for month in range(1,13):
         dstart = (datetime.date(year, month,   1) - datetime.date(year, 1,   1)).days
         if month == 12: dend   = (datetime.date(year+1, 1, 1) - datetime.date(year, 1,   1)).days - 1
         else:           dend   = (datetime.date(year, month+1, 1) - datetime.date(year, 1,   1)).days - 1
         
         variablemonthsum = np.sum(variableyear[dstart:dend,:,:],axis=0)
         variable[startIdx+month-1,:,:] = variablemonthsum

   year0 = int(readstartyear)
   month0 = int(01)
   day0 = int(01)
   epochdt = datetime.datetime(year0,month0,day0,0,0,0)

   #inputTemporalResolution = 'daily'
   inputTemporalResolution = 'monthly'
   variableShort = 'SMB'

   units = 'mmWE day^-1'

# Read the icemask
if args.icemask:
   if args.variable == 'smb_downscaled' or args.variable == 'runoff_downscaled':
      maskncfilename = RACMOnetCDFdirectory + '/Topo_icemask_1km_average_CW.nc'
      maskncfile = Dataset(maskncfilename, 'r')
      icemask = maskncfile.variables['icemask'][:]

      ds = gdal.Open('NETCDF:"'+ RACMOnetCDFdirectory + '/Topo_icemask_1km_average_CW.nc'+'":icemask')
      icemask_geoTransform = np.round(ds.GetGeoTransform())
      xicemask = icemask_geoTransform[0] + icemask_geoTransform[1] * np.arange(0,icemask.shape[1])
      yicemask = icemask_geoTransform[3] + icemask_geoTransform[5] * np.arange(0,icemask.shape[0])
   else:
      maskncfile = Dataset(RACMOnetCDFdirectory + '/ZGRN11_Masks.nc','r')
      icemask = maskncfile.variables['IceSheetMask'][0][0][:][:]

elif args.mask:
   mask = raster.readRasterBandAsArray(args.mask, 1)
   maskGeoTransform = raster.getCoordinates(args.mask, 1)

# a test
# import pdb; pdb.set_trace()
# variable = variable[:,2300:2400,200:300]
# variableError = smb_error_percentage * np.abs(variable)
# uvariable = unp.uarray(variable, variableError)

################################
# Setup the projection and units
################################
area = 1
if args.coordinates == 'utm':
   p = Proj(proj='utm',zone=10,ellps='WGS84')
   area = args.spatialresolution * args.spatialresolution

   units = units.replace(' m^-2', '')

elif args.coordinates == 'stere':
   p = Proj(proj='stere',lat_0=90,lat_ts=70,lon_0=-45,ellps='WGS84')
   area = args.spatialresolution * args.spatialresolution
   
   units = units.replace(' m^-2', '')

   
# Setup bounding box
iStep = args.spatialresolution

if args.boundingbox:
   xmin = args.boundingbox[0]
   xmax = args.boundingbox[1]
   ymin = args.boundingbox[2]
   ymax = args.boundingbox[3]
elif args.clipfile:
   # Shapefile
   (xClip, yClip) = basinUnionPolygon(args.clipfile, args.attributefilter)
   xmin = int(xClip.min() - iStep*2)
   xmax = int(xClip.max() + iStep*2)
   ymin = int(yClip.min() - iStep*2)
   ymax = int(yClip.max() + iStep*2)
else:
   xmin = np.min(x)
   xmax = np.max(x)
   ymin = np.min(y)
   ymax = np.max(y)

# Setup regular grid
print "setting up grid"
nx = (xmax - xmin) / iStep
ny = (ymax - ymin) / iStep

xi = np.linspace(xmin,xmax,num=nx+1)
yi = np.linspace(ymin,ymax,num=ny+1)

geoTransform = (xi.min()-iStep/2, iStep, 0, yi.max()+iStep/2, 0, -iStep)

xi, yi = np.meshgrid(xi, yi)

# Setup x's and y's at which RACMO data is posted
if args.coordinates == 'latlon':
   x = lon
   y = lat
elif args.nointerpolate:
   x, y = np.meshgrid(x, y)
   xi = x
   yi = y
   geoTransform = (xi.min()-iStep/2, iStep, 0, yi.max()+iStep/2, 0, -iStep)
else:
   x,y = p(lon,lat)

x = np.ravel(x)
y = np.ravel(y)

# Setup an array which will be used to clip the RACMO field
#  Here, we do nearneighbor interpolation from the icemask to each
#  data point. This is done because the icemask provided does not
#  have the same spatial extent as the RACMO data.
#if args.nointerpolate:
#   icemaski = icemask
#else:
if args.icemask:
   f = interp2d(xicemask, yicemask, icemask)
   maskArrayi = np.round(f(xi,yi))
else:
   maskArrayi = np.ones(xi.shape)

if args.clipfile:
   maskArrayi = np.flipud(clipImage(np.flipud(maskArrayi), xClip, yClip, geoTransform))

if args.mask:
   maskValsi = np.nan * np.ones( xi.shape )
   for irow in range(0,xi.shape[0]):
      for icol in range(0,xi.shape[1]):
         maskValsi[irow,icol] = raster.sampleRasterAtPoint(mask, maskGeoTransform, xi[irow,icol], yi[irow,icol], method='nearest')
   maskArrayi = np.where(maskValsi == args.maskval, maskArrayi, 0.)

# Setup the indexes
startdt = datetime.datetime(startyear,startmonth,startday,0,0,0)
enddt = datetime.datetime(endyear,endmonth,endday,0,0,0)

if inputTemporalResolution == 'daily':
   startIdx = (startdt - epochdt).days
   endIdx   = (enddt - epochdt).days

elif inputTemporalResolution == 'monthly':
   startIdx = monthdelta(epochdt, startdt)
   endIdx = monthdelta(epochdt, enddt) + 1
   
if debug: print 'startdate  = ' + str(args.startdate)
if debug: print 'startIdx = ' + str(startIdx)
if debug: print 'endIdx   = ' + str(endIdx)

# Anomaly
if args.anomaly:   
   meanstartdt = datetime.datetime(meanstartyear,meanstartmonth,meanstartday,0,0,0)
   meanenddt   = datetime.datetime(meanendyear,meanendmonth,meanendday,0,0,0)

   if inputTemporalResolution == 'daily':
      meanstartIdx = (meanstartdt - epochdt).days;     meanendIdx = (meanenddt - epochdt).days + 1
   elif inputTemporalResolution == 'monthly':
      meanstartIdx = monthdelta(epochdt, meanstartdt); meanendIdx = monthdelta(epochdt, meanenddt) + 1

   # Find mean field
   variableAnomalyMean = np.mean(variable[meanstartIdx:meanendIdx,:,:],axis=0)
   if args.clipfile or args.mask:
      variableAnomalyMean = np.where(maskArrayi > 0.5, variableAnomalyMean, np.nan)

   # SMB anomaly error
   # a text
   # if args.errortiff and args.variable == 'smb_downscaled': uvariableAnomalyMean = uvariable[meanstartIdx:meanendIdx,:,:].mean(axis=0)
   variableAnomalyMeanError = np.sqrt(np.sum(np.square((1./float(meanendIdx-meanstartIdx))*smb_error_percentage*variable[meanstartIdx:meanendIdx,:,:]),axis=0))

# Write header of stats file
if args.stats:
   filename = args.outputdir + '/' + args.regionname + '_' + args.variable + '_sum_' + args.temporalresolution + '.txt'
   f = open(filename, 'w')
   f.write('year, day of year, ' + variableShort + ' [' + units + ']\n')
   f.close()

# Sum the [variable] over the desired [temporalresolution]
variableSum = np.zeros(variable[1].shape)
if args.clipfile or args.mask:
   variableSum = np.where(maskArrayi > 0.5, np.zeros(maskArrayi.shape), np.nan)

print "summing " + args.variable + " at " + args.temporalresolution + " intervals"
print " from " + str(startIdx) + " to " + str(endIdx)

# if inputTemporalResolution == 'daily' and args.temporalresolution == 'daily':
#    resetEvery = 1  
# elif inputTemporalResolution == 'daily' and args.temporalresolution == 'weekly':
#    resetEvery = 7
# elif inputTemporalResolution == 'daily' and args.temporalresolution == 'yearly':
#    # TBD: This is hard-coded to 365 but needs to be 366 on leap years.
#    resetEvery = 365
# elif inputTemporalResolution == 'monthly' and args.temporalresolution == 'monthly':
#    resetEvery = 1
# else:
#    print "invalid temporal resolution specified: " + args.temporalresolution
#    sys.exit()

# Loop
dateCounter = startdt
year = dateCounter.year
dayOfYear = dateCounter.timetuple().tm_yday

# resetCounter = 1
variableSumiIntegratedVector = []
dnVector = []

if args.pngs or args.tiffs or args.stats or args.matfile or debug:
   for iIdx in range(startIdx, endIdx+1):
       year = dateCounter.year
       month = dateCounter.month
       day = dateCounter.day
       dayOfYear = dateCounter.timetuple().tm_yday  
       
       # Sum the variable
       if debug: print 'summing idx: ' + str(iIdx) + ', year: ' + str(year) + ', month: ' + str(month) + ', day: ' + str(day)
       variableSum = variableSum + variable[iIdx][:][:]
       
       # Check for reset
       resetFlag = False
       import pdb; pdb.set_trace()
       if args.temporalresolution == 'daily':
          if (dateCounter - datePrevious).days == 1: resetFlag = True
       #if args.temporalresolution == 'weekly':
       #   monday1 = (dateCounter - timedelta(days=datePrevious.weekday()))
       #   monday2 = (dateCounter - timedelta(days=datePrevious.weekday()))
       #   if (monday2 - monday1) % 7 == 0: resetFlag = True
       # if args.temporalresolution == 'monthly':
       #    if (dateCounter.year - datePrevious.year) * 12 + dateCounter.month - datePrevious.month % 1 == 0: resetFlag = True
       if args.temporalresolution == 'yearly':
          if relativedelta(end_date, start_date).years

       if resetCounter == resetEvery or (args.outputlast and iIdx == endIdx):    
           # Anomaly
           if args.anomaly:
              variableSum = variableSum - variableAnomalyMean
            
           # Interpolate to regularly spaced grid
           if args.nointerpolate:
              variableSumi = variableSum
           else:
              variableSum = np.ravel(variableSum)
              validIdx = ~np.isnan(variableSum)
              variableSumi = interpolate(args, (x[validIdx], y[validIdx]), variableSum[validIdx], (xi, yi))
           
           # HEAVY DEBUGGING #
           #f = open(args.regionname + 'original_gridpoints.csv','w')
           #np.savetxt(f, np.c_[x, y, np.ravel(variableSum)], fmt="%16.5f %16.5f %16.5f")
           #f.close()
           #outputGeoTIFF(args, year, dayOfYear,variableSumi,geoTransform)
           # HEAVY DEBUGGING #
            
           # Crop to clip path
           if args.clipfile or args.mask:
              #import pdb; pdb.set_trace()
              variableSumi = np.where(maskArrayi == 1, variableSumi, np.nan)
              #variableSumi = clipImage(variableSumi, xClip, yClip, geoTransform)
           #if debug: temp = args.regionname; args.regionname = 'DEBUG_clipped'; outputGeoTIFF(args, year, dayOfYear,variableSumi,geoTransform); args.regionname = temp
           
           # Multiply variable within each grid cell by the area of that cell
           # and divide by the number of time units we're summing over for an
           # average value of the variable over the time span
           variableSumi = area * variableSumi / resetCounter
           
           # Integrate over region
           #dateStr = '%4d%03d' % (year, dayOfYear)
           #dt = datetime.datetime.strptime(dateStr, '%Y%j')
           #dnVector.append(datetime2matlabdn(dt))
           variableSumiIntegrated = np.nansum(variableSumi)
           variableSumiIntegratedVector.append(variableSumiIntegrated)
   
           #if debug: np.savetxt('variableSum.csv', np.c_[x, y, variableSum], delimiter=',', fmt='%16.3f %16.3f %16.14f')
           
           #if debug: temp = args.regionname; args.regionname = 'DEBUG_region'; outputGeoTIFF(args, year, dayOfYear,variableSumi,geoTransform); args.regionname = temp
   
           # Write output
           if args.pngs:
              outputPNG(args, year, dayOfYear, variableSumi, geoTransform)
           if args.tiffs:
              outputGeoTIFF(args, year, dayOfYear, np.flipud(variableSumi), geoTransform)
           if args.stats:
              outputStats(args, year, dayOfYear, variableSumiIntegrated)
           if args.matfile:
              outputMAT(args, year, dayOfYear, variableSumi, geoTransform)
           if args.csvfile:
              variableSumi = np.ravel(variableSumi)
              outputCSV(args, year, dayOfYear, xi, yi, variableSumi)
              
           # Reset all variables
           # resetCounter = 0
           variableSum = np.zeros(variable[:][:].shape)
           variableSum = np.where(maskArrayi > 0.5, np.zeros(maskArrayi.shape), np.nan)
           #f = open(args.regionname + 'original_gridpoints.csv','w')
           #np.savetxt(f, np.c_[x, y, np.ravel(variableSum)], fmt="%16.5f %16.5f %16.5f")
           #f.close()
           variableSumiIntegratedVector = []
   
       # Increment the resetCounter and date
       datePrevious = dateCounter
       # resetCounter = resetCounter + 1
       if inputTemporalResolution == 'daily':
          dateCounter = dateCounter + datetime.timedelta(days=1)
       if inputTemporalResolution == 'monthly':
          dateCounter = add_months(dateCounter,1)
       

# Anomaly - output the total sum minus the average
#  NOTE: Anomaly is typically calculated for SMB monthly, which is a total mm W.E. and not a rate,
#        so there's no scaling by time below like there is in the loop above.
#{{{
if args.anomaly:
   variableAnomalyMean = np.ravel(variableAnomalyMean) # [mm W.E. / month]
   #from above: variableAnomalyMeanError
   
   # Sum
   variableSum = np.ravel(np.sum(variable[startIdx:endIdx,:,:],axis=0)) # [mm W.E.]
   # a test
   # if args.errortiff and args.variable == 'smb_downscaled': uvariableSum = uvariable[startIdx:endIdx,:,:].sum(axis=0)
   variableSumError = smb_error_percentage * np.sqrt(np.sum(np.square(variable[startIdx:endIdx,:,:]),axis=0))

   # Anomaly
   variableAnomaly = variableSum - variableAnomalyMean * (endIdx-startIdx) # [mm W.E.]
   # a test
   # if args.errortiff and args.variable == 'smb_downscaled': uvariableAnomaly = np.ravel(uvariableSum - uvariableAnomalyMean * (endIdx-startIdx))
   variableAnomalyError = np.ravel(np.sqrt( variableSumError**2 + (endIdx-startIdx)**2*variableAnomalyMeanError**2 ))
   
   # Interpolate to regularly spaced grid
   validIdx = ~np.isnan(variableAnomalyMean)
   variableSumi         = interpolate(args, (x[validIdx], y[validIdx]), variableSum[validIdx],     (xi, yi))
   variableAnomalyMeani = interpolate(args, (x[validIdx], y[validIdx]), variableAnomalyMean[validIdx], (xi, yi))
   variableAnomalyi     = interpolate(args, (x[validIdx], y[validIdx]), variableAnomaly[validIdx], (xi, yi))

   if args.errortiff and args.variable == 'smb_downscaled':
      variableAnomalyErrori= interpolate(args, (x[validIdx], y[validIdx]), variableAnomalyError[validIdx], (xi, yi))
   
   # DEBUG:
   #import pdb; pdb.set_trace()
   #lon = np.ravel(lon)
   #lat = np.ravel(lat)
   #filename = 'icemask.txt'
   #f = open(filename, 'w')
   #for i, v in enumerate(np.ravel(icemask)):
   #   printString = '%10.6f, %10.6f, %10.6f, %10.6f, %16.6f\n' % (lon[i], lat[i], x[i], y[i], v)
   #   f.write(printString)
   #f.close()
   #import pdb; pdb.set_trace()
   
   # Crop to clip path
   if args.clipfile or args.mask:
      #variableSumi = np.where(maskArrayi == 1, variableSumi, np.nan)
      variableSumi = clipImage(variableSumi, xClip, yClip, geoTransform)
      #variableAnomalyi = np.where(maskArrayi == 1, variableAnomalyi, np.nan)
      variableAnomalyi = clipImage(variableAnomalyi, xClip, yClip, geoTransform)
        
   # This needs to be replaced
   outputGeoTIFF2(args, 'sum',      variableSumi, geoTransform)
   outputGeoTIFF2(args, 'mean',     variableAnomalyMeani, geoTransform)
   outputGeoTIFF2(args, 'anomaly',  variableAnomalyi, geoTransform)

   if args.errortiff and args.variable == 'smb_downscaled':
      # a test
      # variableAnomalyError = unp.std_devs(uvariableAnomaly)
      outputGeoTIFF2(args, 'anomalyerror', variableAnomalyErrori, geoTransform)
#}}}

