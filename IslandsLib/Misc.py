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

from osgeo import gdal
from osgeo import ogr
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




def list_errors():
    """Small function to get the description of python errors.
    """
    import os
    import errno
    from pprint import pprint

    pprint({i: os.strerror(i) for i in sorted(errno.errorcode)})

