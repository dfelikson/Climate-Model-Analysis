import shapely
from shapely.ops import cascaded_union
from shapely.geometry import Polygon

import operator
from osgeo import gdal, gdalnumeric, ogr, osr
from PIL import Image, ImageDraw

import numpy as np

# For debugging
import scipy.io as sio
import sys
import matplotlib.pyplot as plt



def getUniqueFieldValues(shapefile, fieldName): #{{{
   driver = ogr.GetDriverByName("ESRI Shapefile")
   dataSource = driver.Open(shapefile, 0)
   layer = dataSource.GetLayer()
   
   values_list = list()
   for feature in layer:
      values_list.append(feature.GetField(fieldName))

   return sorted(set(values_list))
   #}}}

def checkForField(shapefile, fieldName): #{{{
   driver = ogr.GetDriverByName("ESRI Shapefile")
   dataSource = driver.Open(shapefile, 0)
   layer = dataSource.GetLayer()
   fieldNames = [field.name for field in layer.schema]
   if fieldName in fieldNames:
       return True
   else:
       return False
   #}}}

# This function creates a shapely polygon which is the union
# of several polygons within the given shapefile filtered
# based on attributeFilter.
def basinUnionPolygon(shapefile, attributeFilter):
   driver = ogr.GetDriverByName("ESRI Shapefile")
   dataSource = driver.Open(shapefile, 0)
   layer = dataSource.GetLayer()

   if attributeFilter:
      layer.SetAttributeFilter(attributeFilter)

   polys = []
   #for i in range(layer.GetFeatureCount()):
   for feature in layer:
      #feature = layer.GetFeature(i)
      wkt = feature.geometry().ExportToWkt()
      polys.append(shapely.from_wkt(wkt))

   polyUnion = cascaded_union(polys)

   # DEBUG
   #for poly in polys:
   #   plt.plot(np.array(poly.exterior.xy)[0,:],np.array(poly.exterior.xy)[1,:],'b')   
   #plt.plot(np.array(polyUnion.exterior.xy)[0,:],np.array(polyUnion.exterior.xy)[1,:],'r')
   #plt.show()
   
   if isinstance(polys[0], shapely.geometry.multipolygon.MultiPolygon):
      polyUnionArray = np.array(polyUnion.convex_hull.exterior.xy)
   else:
      polyUnionArray = np.array(polyUnion.exterior.xy)

   return (polyUnionArray[0,:], polyUnionArray[1,:])

# This function will convert the rasterized clipper shapefile to a
# mask for use within GDAL.    
def imageToArray(i):
    """
    Converts a Python Imaging Library array to a gdalnumeric image.
    """
    a=gdalnumeric.fromstring(i.tobytes(),'b')
    a.shape=i.im.size[1], i.im.size[0]
    return a

def arrayToImage(a):
    """
    Converts a gdalnumeric array to a Python Imaging Library Image.
    """
    i=Image.fromstring('L',(a.shape[1],a.shape[0]),
            (a.astype('b')).tostring())
    return i
     
def world2Pixel(geoMatrix, x, y):
  """
  Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
  the pixel location of a geospatial coordinate 
  """
  ulX = geoMatrix[0]
  ulY = geoMatrix[3]
  xDist = geoMatrix[1]
  yDist = geoMatrix[5]
  rtnX = geoMatrix[2]
  rtnY = geoMatrix[4]
  pixel = int((x - ulX) / xDist)
  line = int((ulY - y) / xDist)
  return (pixel, line) 

def clipImage(imgArray, xClip, yClip, imgGeoTransform):
  """
  Map points to pixels for drawing the field boundary on a blank
  8-bit, black and white, mask image.
  """
  pixels = []
  xClipPix = []
  yClipPix = []
  for i in range(len(xClip)):
    pixels.append(world2Pixel(imgGeoTransform, xClip[i], yClip[i]))
    xClipPix.append(pixels[i][0])
    yClipPix.append(pixels[i][1])

  pxHeight, pxWidth = imgArray.shape
  rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
  rasterize = ImageDraw.Draw(rasterPoly)
  rasterize.polygon(pixels, 0)
  mask = imageToArray(rasterPoly)

  """
  Clip the image using the mask
  """
  clip = gdalnumeric.choose(mask, (imgArray, np.nan))
  
  # DEBUG
  #import pdb; pdb.set_trace()
  #plt.imshow(imgArray)
  #plt.plot(xClipPix,yClipPix,'r')
  #plt.show()
  
  #import pdb; pdb.set_trace()
  #sio.savemat('test.mat',{'imgArray':imgArray, 'pixels':pixels})
  
  return clip
  
def histogram(a, bins=range(0,256)):
  """
  Histogram function for multi-dimensional array.
  a = array
  bins = range of numbers to match 
  """
  fa = a.flat
  n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)
  n = gdalnumeric.concatenate([n, [len(fa)]])
  hist = n[1:]-n[:-1] 
  return hist

def stretch(a):
  """
  Performs a histogram stretch on a gdalnumeric array image.
  """
  hist = histogram(a)
  im = arrayToImage(a)   
  lut = []
  for b in range(0, len(hist), 256):
    # step size
    step = reduce(operator.add, hist[b:b+256]) / 255
    # create equalization lookup table
    n = 0
    for i in range(256):
      lut.append(n / step)
      n = n + hist[i+b]
  im = im.point(lut)
  return imageToArray(im)
