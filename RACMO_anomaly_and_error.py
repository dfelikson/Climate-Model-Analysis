
   variableAnomalyMean = np.ravel(variableAnomalyMean) # [mm W.E. / month]
   
   # Sum
   variableSum = np.ravel(np.sum(variable[startIdx:endIdx,:,:],axis=0)) # [mm W.E.]
   variableSumError = smb_error_percentage * np.sqrt(np.sum(np.square(variable[startIdx:endIdx,:,:]),axis=0))

   # Anomaly
   variableAnomaly = variableSum - variableAnomalyMean * (endIdx-startIdx) # [mm W.E.]
   # a test
   # if args.errortiff and args.variable == 'smb_downscaled': uvariableAnomaly = np.ravel(uvariableSum - uvariableAnomalyMean * (endIdx-startIdx))
   variableAnomalyError = np.ravel(np.sqrt( variableSumError**2 + (endIdx-startIdx)**2*variableAnomalyMeanError**2 ))
   
   # Interpolate to regularly spaced grid
   if args.nointerpolate:
      variableSumi         = np.flipud(np.reshape(variableSum, xi.shape))
      variableAnomalyMeani = np.flipud(np.reshape(variableAnomalyMean, xi.shape))
      variableAnomalyi     = np.flipud(np.reshape(variableAnomaly, xi.shape))
   else:
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

