#!/bin/env python

# -------------
# Configuration
# -------------
import os, sys, datetime, calendar
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

import argparse

# Interpolation
from scipy.interpolate import griddata
#from pykrige.ok import OrdinaryKriging
import numpy as np
#import pykrige.kriging_tools as kt

# GDAL
from osgeo import gdal, osr, gdalconst

def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + datetime.timedelta(days = 366)
   frac = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac

def outputPNG(args, year, dayOfYear, variableSumi, geoTransform):
   print('writing PNG for year ' + str(year) + ', DOY ' + str(dayOfYear))
   filename = args.outputdir + '/' + args.regionname + '_' + args.variable + '_sum_year%4d_doy%03d.png' % (year, dayOfYear)
   plt.imshow(np.flipud(variableSumi), vmin=0, vmax=10, cmap='Reds')
   plt.colorbar()
   plt.title(args.regionname + ' ' + str(year) + ' DOY ' + str(dayOfYear) + ': ' + args.variable)
   plt.savefig(filename)
   plt.close()
      
def outputGeoTIFF(args, year, dayOfYear, variableSumi, geoTransform):
   print('writing GTiff for year ' + str(year) + ', DOY ' + str(dayOfYear))
   filename = args.outputdir + '/' + args.regionname + '_' + args.variable + '_sum_year%4d_doy%03d.tif' % (year, dayOfYear)
   
   cols = variableSumi.shape[1]
   rows = variableSumi.shape[0]

   driver = gdal.GetDriverByName('GTiff')
   outRaster = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
   outRaster.SetGeoTransform(geoTransform)
   outBand = outRaster.GetRasterBand(1)

   arrayout = np.where(~np.isnan(variableSumi), variableSumi, -9999)

   outBand.WriteArray(arrayout)
   outBand.SetNoDataValue(-9999)
   outRasterSRS = osr.SpatialReference()
   outRasterSRS.ImportFromEPSG(3413)
   outRaster.SetProjection(outRasterSRS.ExportToWkt())
   outBand.FlushCache()

def outputGeoTIFF2(args, label, variableSumi, geoTransform):
   filename = args.outputdir + '/' + args.regionname + '_' + args.variable + '_' + label + '.tif'
   
   cols = variableSumi.shape[1]
   rows = variableSumi.shape[0]

   driver = gdal.GetDriverByName('GTiff')
   outRaster = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
   outRaster.SetGeoTransform(geoTransform)
   outBand = outRaster.GetRasterBand(1)

   arrayout = np.where(~np.isnan(variableSumi), variableSumi, -9999)

   outBand.WriteArray(arrayout)
   outBand.SetNoDataValue(-9999)
   outRasterSRS = osr.SpatialReference()
   outRasterSRS.ImportFromEPSG(3413)
   outRaster.SetProjection(outRasterSRS.ExportToWkt())
   outBand.FlushCache()
   
def outputStats(args, year, dayOfYear, variableSumiIntegrated):
   print('writing stats for year ' + str(year) + ', DOY ' + str(dayOfYear))
   filename = args.outputdir + '/' + args.regionname + '_' + args.variable + '_sum_' + args.temporalresolution + '.txt'
   f = open(filename, 'a')
   printString = '%4d, %03d, %12.3f\n' % (year, dayOfYear, variableSumiIntegrated)
   f.write(printString)
   f.close()
   
def outputMAT(args, year, dayOfYear, variableSumi, geoTransform):
   print("writing MAT file")
   filename = args.outputdir + '/' + args.regionname + '_racmo' + '.mat'
   
   # Load mat file
   if os.path.isfile(filename):
      matContents = sio.loadmat(filename)
      racmo = matContents['racmo']
      racmo = racmo[0,0]
   else:
      racmo = {}
   
   # Append to dictionary   
   if os.path.isfile(filename):
      variableStacked = np.dstack( (racmo[args.variable], variableSumi) )
      racmo[args.variable] = variableStacked
      yearStacked = np.vstack( (racmo['year'], year) )
      racmo['year'] = yearStacked
      dayOfYearStacked = np.vstack( (racmo['dayOfYear'], dayOfYear) )
      racmo['dayOfYear'] = dayOfYearStacked
   else:
      racmo[args.variable] = variableSumi
      racmo['year'] = year
      racmo['dayOfYear'] = dayOfYear
   
   racmo['geoTransform'] = geoTransform
   
   # Save
   sio.savemat(filename, {'racmo': racmo})

def outputCSV(args, year, dayOfYear, x, y, variableSum):
   print("writing CSV file")
   filename = args.outputdir + '/' + args.regionname + '_' + args.variable + '_' + str(year) + '_' + str(dayOfYear) + '.csv'
   f = open(filename, 'w')
   np.savetxt(f, np.c_[x, y, variableSum], delimiter=',', fmt='%16.5f %16.5f %16.5f')
   f.close()

# Time vector functions
def monthdelta(d1, d2):
   delta = 0
   while True:
      mdays = calendar.monthrange(d1.year, d1.month)[1]
      d1 += datetime.timedelta(days=mdays)
      if d1 <= d2:
         delta += 1
      else:
         break
   return delta

def add_months(sourcedate,months):
   month = sourcedate.month - 1 + months
   year = sourcedate.year + month / 12
   month = month % 12 + 1
   day = min(sourcedate.day,calendar.monthrange(year,month)[1])
   return datetime.date(year,month,day)

def interpolate(args, coords, data, coordsi):
   if args.interpmethod == 'griddata':
      #datai = griddata(coords, data, coordsi, method='linear')
      datai = griddata(coords, data, coordsi, method='nearest')
      datai = np.flipud(datai)
   if args.interpmethod == 'kriging':
      # Crop the data to bounding box, otherwise kriging takes very long
      idxCrop = (coords[0] >= args.boundingbox[0]) & (coords[0] <= args.boundingbox[1]) & (coords[1] >= args.boundingbox[2]) & (coords[1] <= args.boundingbox[3]) & (~np.isnan(data))
                       
      dataCropped = data[idxCrop]
      x = coords[0][idxCrop]
      y = coords[1][idxCrop]

      OK = OrdinaryKriging(x, y, dataCropped, variogram_model='exponential',
                           verbose=False, enable_plotting=False)
      import pdb; pdb.set_trace()
      datai, ss = OK.execute('grid', np.ravel(coordsi[0]), np.ravel(coordsi[1]))
   
   return datai

def four_floats(value):
    values = value.split()
    if len(values) != 4:
        raise argparse.ArgumentError
    values = map(float, values)
    return values
