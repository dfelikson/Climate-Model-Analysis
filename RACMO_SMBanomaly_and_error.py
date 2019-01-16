#!/usr/bin/env python

import os, sys, argparse, datetime
from subprocess import call
import numpy as np
#np.seterr(invalid='ignore') # ignore invalid value warnings
from netCDF4 import Dataset
from pyproj import Proj

from scipy.interpolate import interp2d

import string

import copy

from RACMOgridAndStatsGlobals import *
from RasterClipperFunctions import *
from RACMOutilities import *

sys.path.append('/home/student/denis/ScriptsAndUtilities/pythonModules')
import raster


# Setup
RACMO_directory = os.environ['RACMO_downscaled_DIR']
precip_netcdf = RACMO_directory + '/precip.1958-2017.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc'
runoff_netcdf = RACMO_directory + '/runoff.1958-2017.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc'
smb_netcdf    = RACMO_directory + '/smb_rec.1958-2017.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc'

output_directory = '/disk/staff/gcatania/polar/Arctic/data/RACMO/RACMO2.3/RACMO2.3p2anomaly'

# For present-day climate ... use the mean tiff:
meanstartyear = 2000; meanstartmonth = 9; meanstartday =  1;
meanendyear   = 2015; meanendmonth   = 8; meanendday   = 31;
startyear = 2000; startmonth = 9; startday =  1;
endyear   = 2015; endmonth   = 8; endday   = 31;

# For present-day climate anomaly from mean 1971-1988:
#meanstartyear = 1971; meanstartmonth = 9; meanstartday =  1;
#meanendyear   = 1988; meanendmonth   = 8; meanendday   = 31;
#startyear = 2000; startmonth = 9; startday =  1;
#endyear   = 2015; endmonth   = 8; endday   = 31;

# For stable climate anomaly from mean 1971-1988:
#meanstartyear = 1971; meanstartmonth = 9; meanstartday =  1;
#meanendyear   = 1988; meanendmonth   = 8; meanendday   = 31;
#startyear = 1971; startmonth = 9; startday =  1;
#endyear   = 1988; endmonth   = 8; endday   = 31;

startdate_string = datetime.date(startyear, startmonth, startday).strftime('%Y%b%d')
enddate_string = datetime.date(endyear, endmonth, endday).strftime('%Y%b%d')
meanstartdate_string = datetime.date(meanstartyear, meanstartmonth, meanstartday).strftime('%Y%b%d')
meanenddate_string = datetime.date(meanendyear, meanendmonth, meanendday).strftime('%Y%b%d')

precip_error = 0.1
runoff_error = 0.2


# Process
print("reading precip netcdf")
ncfilename = precip_netcdf
ncfile = Dataset(ncfilename, 'r')
precip = ncfile.variables['precipcorr']
x = ncfile.variables['lon'][:]
y = ncfile.variables['lat'][:]
xstep = x[1]-x[0]
ystep = y[1]-y[0]
geoTransform = (x.min()-xstep/2, xstep, 0, y.max()+ystep/2, 0, -ystep)

print("reading runoff netcdf")
ncfilename = runoff_netcdf
ncfile = Dataset(ncfilename, 'r')
runoff = ncfile.variables['runoffcorr']

print("reading smb netcdf")
ncfilename = smb_netcdf
ncfile = Dataset(ncfilename, 'r')
smb = ncfile.variables['SMB_rec']

epochdt = datetime.datetime(1958,1,15,0,0,0)
startdt = datetime.datetime(startyear,startmonth,startday,0,0,0)
enddt   = datetime.datetime(endyear,endmonth,endday,0,0,0)


# SMB anomaly
print('calculating smb anomaly')

# mean
meanstartdt = datetime.datetime(meanstartyear,meanstartmonth,meanstartday,0,0,0)
meanenddt   = datetime.datetime(meanendyear,meanendmonth,meanendday,0,0,0)
meanstartIdx = monthdelta(epochdt, meanstartdt); meanendIdx = monthdelta(epochdt, meanenddt)
smbMean    = copy.deepcopy(np.mean(smb[meanstartIdx:meanendIdx,:,:],axis=0)) # [mm W.E. / month]
precipMean = copy.deepcopy(np.mean(precip[meanstartIdx:meanendIdx,:,:],axis=0))

# # DEBUG: mean of annual means is no different than mean over all months
# years = np.arange(meanstartyear,meanendyear+1)
# smbMean_check = np.zeros( (len(years), smbMean.shape[0], smbMean.shape[1]) )
# for iyear,year in enumerate(np.arange(meanstartyear,meanendyear)):
#    meanstartdt = datetime.datetime(year,meanstartmonth,meanstartday,0,0,0)
#    meanenddt   = datetime.datetime(year+1,meanendmonth,meanendday,0,0,0)
#    meanstartIdx = monthdelta(epochdt, meanstartdt); meanendIdx = monthdelta(epochdt, meanenddt)
#    smbMean_check[iyear,:,:] = np.mean(smb[meanstartIdx:meanendIdx,:,:],axis=0)
# plt.imshow( np.flipud( (np.mean(smbMean_check,axis=0) - smbMean)/smbMean ) )

# sum
startIdx = monthdelta(epochdt, startdt); endIdx   = monthdelta(epochdt, enddt) + 1;
smbSum = copy.deepcopy(np.sum(smb[startIdx:endIdx,:,:],axis=0))  # [mm W.E.]

# anomaly
smbAnomaly = copy.deepcopy(smbSum - smbMean * (endIdx-startIdx)) # [mm W.E.]

smbSum         = np.flipud(copy.deepcopy(smbSum))
smbAnomalyMean = np.flipud(copy.deepcopy(smbMean))
smbAnomaly     = np.flipud(copy.deepcopy(smbAnomaly))

# write tifs
output_filename = output_directory + '/GrIS_smb_downscaled_sum_' + startdate_string + '-' + enddate_string + '_fromMean_' + meanstartdate_string + '-' + meanenddate_string + '_dh.tif'
print('writing: ' + output_filename)
raster.writeArrayAsRasterBand(output_filename, geoTransform, smbSum / 917., -9999.)

output_filename = output_directory + '/GrIS_smb_downscaled_mean_' + startdate_string + '-' + enddate_string + '_fromMean_' + meanstartdate_string + '-' + meanenddate_string + '_dh.tif'
print('writing: ' + output_filename)
raster.writeArrayAsRasterBand(output_filename, geoTransform, smbAnomalyMean / 917., -9999.)

output_filename = output_directory + '/GrIS_smb_downscaled_anomalySum_' + startdate_string + '-' + enddate_string + '_fromMean_' + meanstartdate_string + '-' + meanenddate_string + '_dh.tif'
print('writing: ' + output_filename)
raster.writeArrayAsRasterBand(output_filename, geoTransform, smbAnomaly / 917., -9999.)

# SMB error (calculated from percent errors on annual precipitation and annual runoff)
print('calculating smb error')
smbMeanError = np.zeros(smbMean.shape)
for year in np.arange(meanstartyear,meanendyear):
   dt1 = datetime.datetime(year,startmonth,startday,0,0,0)
   dt2 = datetime.datetime(year+1,endmonth,endday,0,0,0)
   idx1 = monthdelta(epochdt, dt1)
   idx2 = monthdelta(epochdt, dt2)
   
   precip_annual = np.sum(precip[idx1:idx2,:,:],axis=0)
   runoff_annual = np.sum(runoff[idx1:idx2,:,:],axis=0)
   precip_annual_error = precip_error * precip_annual
   runoff_annual_error = runoff_error * runoff_annual
   smbMeanError = smbMeanError + precip_annual_error**2 + runoff_annual_error**2

smbMeanError = (1./(meanendyear-meanstartyear+1)) * np.sqrt(smbMeanError)

smbCumulativeError = np.zeros(smbSum.shape)
for year in np.arange(startyear,endyear):
   dt1  = datetime.datetime(year,startmonth,startday,0,0,0)
   dt2  = datetime.datetime(year+1,endmonth,endday,0,0,0)
   idx1 = monthdelta(epochdt, dt1)
   idx2 = monthdelta(epochdt, dt2)
   
   precip_annual = np.sum(precip[idx1:idx2,:,:],axis=0)
   runoff_annual = np.sum(runoff[idx1:idx2,:,:],axis=0)
   precip_annual_error = precip_error * precip_annual
   runoff_annual_error = runoff_error * runoff_annual
   smbCumulativeError = smbCumulativeError + precip_annual_error**2 + runoff_annual_error**2

smbCumulativeError = np.sqrt(smbCumulativeError)

smbAnomalySumError = np.sqrt(smbCumulativeError**2 + (endyear-startyear+1)**2 * smbMeanError**2)
smbAnomalySumError = np.flipud(copy.deepcopy(smbAnomalySumError))

# write tifs
output_filename = output_directory + '/GrIS_smb_downscaled_anomalySumError_' + startdate_string + '-' + enddate_string + '_fromMean_' + meanstartdate_string + '-' + meanenddate_string + '_dh.tif'
print('writing: ' + output_filename)
raster.writeArrayAsRasterBand(output_filename, geoTransform, smbAnomalySumError / 917., -9999.)

# # account for range in density (750 +/- 150)
# output_filename = output_directory + '/GrIS_smb_downscaled_anomalySumError_' + startdate_string + '-' + enddate_string + '_fromMean_' + meanstartdate_string + '-' + meanenddate_string + '_densityError_dh.tif'
# print('writing: ' + output_filename)
# raster.writeArrayAsRasterBand(output_filename, geoTransform, smbAnomalySumError / 917 + smbAnomaly / 4500., -9999.)
# 
