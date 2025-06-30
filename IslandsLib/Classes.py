#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Classes used for FEM

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


class river:
    """container for rivers and by extension the boundaries of the domain, I can see the comments of some people..., :))
    """

    def __init__(self, rname, x, y, rtype="b", downstream=True):
        self.rname = rname
        self.x = x
        self.y = y
        self.z = x*0  # initialise à 0
        self.closed = False
        self.nodes = np.arange(len(x))
        self.downstream = downstream
        self.flux_nul = False  # si vrai impose une condition de flux nul
        self.segments = []
        self.rtype = rtype  # border type


class topo:
    """container for DEMs.
    """

    def __init__(self, dname, x, y, z):
        self.dname = dname
        self.x = x
        self.y = y
        self.z = z

