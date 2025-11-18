#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Miscellaneous.

"""
__author__ = "François Métivier"
__copyright__ = "Copyright 2024"
__license__ = "GPL"
__version__ = "1.0"


#############################
#
# Librairies
#
#############################
import sys 
sys.path.append("..")

import sys


# Add path to pyFreeFem
machine='home'
if machine=='lab':
    sys.path.append('/home/francois/Nextcloud/src/pyFreeFem-master/')
else:
    sys.path.append('/home/metivier/Nextcloud/src/pyFreeFem-master/')

import pyFreeFem as pyff

import numpy as np
import time
import geopandas
import pandas

# from osgeo import gdal
# from osgeo import ogr
import cartopy.crs as ccrs
import cartopy.io.shapereader as shapereader

import glob
import shapely as shp

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

import triangle as tr

from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import spsolve

from copy import deepcopy

from geojson import LineString, Feature, FeatureCollection, dump
from pyproj import Transformer



def list_errors():
    """Small function to get the description of python errors.
    """
    import os
    import errno
    from pprint import pprint

    pprint({i: os.strerror(i) for i in sorted(errno.errorcode)})


def dump_geojson(cs, fname, epsg_ori='EPSG:3297', epsg_out='EPSG:4326'):
    """dumps a Matplotlib QuadContourSet of piezometric levels into a geojson file.

    ..Note:
    Geojson Dumps a series of lon, lat listrings whereas folium uses lat, lon format

    Parameters
    ----------
    cs : QuadContourSet
        objects returned by plt.contour() and containing the piezometric iso-levels
    fname : string
        output file name (without extension)
    epsg_ori : str, optional
        input EPSG projection of the piezometic level coordinates, by default 'EPSG:3297'
    epsg_out : str, optional
        output EPSG projection, by default 'EPSG:4326' hence lat lon
    """


    transformer = Transformer.from_crs(epsg_ori,epsg_out)
    
    features = []
    for i in range(len(cs.collections)):
        try:
            pp = cs.collections[i].get_paths()[0] # retrieve contours
            y, x = transformer.transform(pp.vertices[:,0], pp.vertices[:,1]) # get x, y from contours in lat lon
            dat = [[x[i],y[i]] for i in range(len(x))] # build linestring
                    
            LS = LineString(dat)

            features.append(Feature(geometry=LS, properties={"Piezometric level": cs.levels[i], 'Units': 'meters above sea level'}))
        except:
            pass

    feature_collection = FeatureCollection(features)

    with open('%s.geojson' % fname, 'w') as f:
        dump(feature_collection, f, indent=4)
