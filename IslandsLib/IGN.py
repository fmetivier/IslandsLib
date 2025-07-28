#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

IGN Data Processing Functions

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


def linestring_to_list(geom=""):
    """Transforms a Linestring geometry into a list of x,y positions
    
    Parameters
    ----------
        geom (str, optional): geometry to be transformed. Defaults to "".

    Returns
    -------
        list : x,y coordinates of points making the original geometry
    
    
    Note 
    ----
        use shapely if > 2 available

    """
    x, y = [], []
    data = geom.split('(')[1]
    data = data.split(')')[0]
    data = data.split(',')
    for dd in data:
        d = dd.lstrip(' ').split(" ")
        x.append(float(d[0]))
        y.append(float(d[1]))

    return x, y


def polygon_to_list(geom=""):
    """Transforms a Polygon type geometry into a list of x,y positions
    
    Parameters
    ----------
        geom (str, optional): geometry to be transformed. Defaults to "".

    Returns
    -------
        list : x,y coordinates of points making the original geometry
    
    
    Note
    ----
        use shapely if > 2 available

    """

    x, y = [], []
    data = geom.strip('POLYGON')
    data = data.strip('((')
    data = data.strip('))')
    data = data.strip(')')
    data = data.strip('(')    
    data = data.split(',')
    for dd in data:
        try:
            d = dd.lstrip(' ').split(" ")
            xv = float(d[0])
            yv = float(d[1])
            x.append(xv)
            y.append(yv)
        except:
            print("raté")
    return x, y

def IGN_to_segment_v2(lg, list_ori, ex_dic, plot_res = False):
    """Transforms a series of segments into a single contour.

    Parameters
    ----------
        lg (dict): key = segment number, value = [x,y]. Defaults to {}.
        list_ori (list): list of keys. Defaults to [].
        ex_dic (dict): key = segment number, val = endpoint coordinates [x1,y1,x2,y2]. Defaults to {}.
        plot_res (bool, optional): plot true Defaults to False.
    
    Returns
    -------
        arrays: x,y
    """

        # dictionnaire des orientation de segment
    switchDebutFin_dic={} 
    # démarrage
    start = list_ori[0]
    # liste des segments ordonnés
    finlist = [start] 
    # on commence par la fin du premier segment 
    # mais on ne le retourne pas
    switchDebutFin_dic[start] =  False
    # on enlève le segment qu'on traite
    list_ori.pop(0) 

    # récupération des points de début et de fin de segment
    xi,yi,xf,yf = ex_dic[start]
    #on cherche le segment le plus proche de la fin de start
    x=xf
    y=yf
    
    while len(list_ori) > 0:

        distance = []
        endpoint = []

        for i in range(len(list_ori)): 
            e = ex_dic[list_ori[i]]
            d0  = (e[0]-x)**2 + (e[1]-y)**2
            d1  = (e[2]-x)**2 + (e[3]-y)**2
            distance.append(min(d0,d1))
            if d0 <= d1:
                endpoint.append(False)
            else:
                endpoint.append(True)
        
        distance = np.array(distance)
        indice = np.where(distance == min(distance))

        val = indice[0][0]
        switchDebutFin_dic[list_ori[val]] = endpoint[val]
        finlist.append(list_ori[val])

        start = list_ori[val]
        xi,yi,xf,yf = ex_dic[start]
        # si le segment trouvé se raccroche au prédént par la fin 
        # alors le suivant se raccorchera au présent par le début
        # et vice et versa
        if switchDebutFin_dic[start]:
            x=xi
            y=yi
        else:
            x=xf
            y=yf


        list_ori.pop(val)
        
        if len(indice[0])>1:
            print(start)
            print(indice)
            print(min(distance))
            print(distance[indice])
        
    X,Y =[], []
    for key in finlist:
        if switchDebutFin_dic[key]:
            X.append(lg[key][0][::-1])
            Y.append(lg[key][1][::-1])
        else:
            X.append(lg[key][0])
            Y.append(lg[key][1])


    if plot_res:
        fig,ax = plt.subplots(1)
        i=0
        for x,y in zip(X,Y):
            ax.text(x[0],y[0],i)
            i+=1

    x = np.hstack(X)
    y = np.hstack(Y)

    if plot_res:
        ax.plot(x,y)
    plt.axis("equal")

    return x,y

def IGN_to_segment(lg={}, list_key=[], ex=[], inverted=False, exclusion=[]):
    """Create order in the IGN mess (unordered segments)
    to obtain continuous contours

    warning
    -------
    OBSOLETE CRASHES REGULARLY

    Note
    ----
    #. Create a list of successive segments by noting their direction relative to the first segment in the list.
        * Decide on a direction (end=True means we expect the segment closest to the current segment to be oriented like this one)
        * Take the endpoints of a segment
        * Look for the segment closest to the end (if end=True) or the beginning (if end=False)
        * Test whether the end or the beginning of the segment is closest and modify the end function
        * Do the same with the next segment
    #. Take each segment from the list, flip it if necessary, and extract its coordinates to create a single contour.


    Parameters
    ----------
        lg: dictionary of segments under shape code:(x,y)
        list_key: list of segment codes
        ex: list of segment endpoints (x1,y1,x2,y2)
        inverted: False by default, if True returns inverted coordinates
        exclusion: [] by default, list of exclusion codes

    Returns
    -------
        arrays: x,y contour or segment coordinates
    """
    #
    # A
    #
    ok = []
    done = []
    remain = len(ex)
    i = 0
    end = True
    while len(done) < len(list_key):
        p = ex[i]
        d = 1e20
        if list_key[i] not in done:
            for j in range(len(ex)):
                if j != i:
                    xi, yi, xf, yf = ex[j]

                    if end:
                        # search for nearest to the end
                        d0 = (xi-p[2])**2 + (yi-p[3])**2
                        d1 = (xf-p[2])**2 + (yf-p[3])**2
                    else:
                        # search for nearest to the beginning
                        d0 = (xi-p[0])**2 + (yi-p[1])**2
                        d1 = (xf-p[0])**2 + (yf-p[1])**2

                    print(i,j,d0,d1)
                    if d1 < d or d0 < d:
                        nearest_val = j
                        print(j)
                        d = min(d, d0, d1)
                        if d1 < d0:
                            next_end = False
                        else:
                            next_end = True
            ok.append([list_key[i], list_key[nearest_val], end])
            done.append(list_key[i])
            print(list_key[i], list_key[nearest_val], end)
            i = nearest_val
            if next_end:
                end = True
            else:
                end = False

    #
    # B
    #

    x, y = [], []
    for line in ok:
        if line[0] not in exclusion:
            val = lg[line[0]]
            if not line[2]:
                xtmp = val[0][::-1]
                ytmp = val[1][::-1]
            else:
                xtmp = val[0]
                ytmp = val[1]

            for xt, yt in zip(xtmp, ytmp):
                x.append(xt)
                y.append(yt)

    #
    # Inversion peut être nécessaire pour FreeFem mais pas pour Triangle
    #
    if inverted:
        x = x[::-1]
        y = y[::-1]

    return x, y



def output_to_file(fname,entete,x,y):
    """Outputs x,y lists to text file.
    Used for Islands coastlines

    Parameters
    ----------
        fname (str): file name
        entete (str): header description
        x (list or array): x coordinates
        y (list or array): y coordinates
    """

    f = open(fname, "w")
    f.write(entete)
    for i in range(len(x)):
        f.write("%f,%f\n" % (x[i],y[i]))

    f.close()


def get_alti_lims(fname):
    """Retrieves information about an alti DB file

    Parameters
    ----------
    fname: file name
    """
    f = open(fname)
    # read header and get important values
    param_dic = {}
    for i in range(6):
        line = f.readline()
        d = line.split(" ")
        param_dic[d[0]] = float(d[-1].strip("\n"))

    f.close()

    return (fname, param_dic)


def read_alti(fname):
    """Opens a 25m alti database file, retrieves the altitude and position data, and returns square X, Y, and Z matrices.

    Parameters
    ----------
    fname: Full name of the alti file

    Returns
    -------
    (1000, 1000) nunmpy arrays: X, Y, Z position and altitude matrices
    
    """
    f = open(fname)
    # read header and get important values
    param_dic = {}
    for i in range(6):
        line = f.readline()
        d = line.split(" ")
        param_dic[d[0]] = float(d[-1].strip("\n"))

    X = param_dic["xllcorner"] * np.ones(1000) + np.arange(1000) * 25
    Y = param_dic["yllcorner"] * np.ones(1000) + np.arange(1000) * 25
    Z = np.zeros((1000, 1000))
    print(param_dic)
    lines = f.readlines()
    c = 0
    for i in np.arange(len(lines) - 1, -1, -1):
        data = lines[i].split(" ")
        Z[c, :] = data[1:]
        c += 1

    return X, Y, Z


def create_river_from_linestrings(gpd, ilist):
    """Creates two x,y matrices from a given list of linestrings in the correct order.


    Parameters
    ----------
    gpd: geopandas dataframe containing the geometric data
    ilist: list of indices (iloc) of the topo objects to be retrieved from gpd

    Returns
    -------
    float arrays: xr, yr coordinate vectors    


    Note
    ---- 
    IGN shapefiles are not ordered; river or coastline segments must be ordered first to reconstruct a continuous contour.
    """

    xr = []
    yr = []
    for i in ilist:
        x, y = gpd.iloc[i].geometry.coords.xy
        xr.append(x)
        yr.append(y)
    xr = np.hstack(xr)
    yr = np.hstack(yr)

    return xr, yr


def get_elevation_from_dem(dem, borders, minx=0, miny=0):
    """Retrieves the elevation of each river in the border from the DEM.

    Parameters
    ----------
    borders: list of river-class objects
    dem: DEM class obtained from bdalti

    Returns
    -------
        list of classes: borders with elevation 
    """

    for r in borders:
        print(r.rname)
        # print(np.hstack(iterloop))
        for i in np.arange(len(r.x)):
            j = int((r.x[i] + minx - dem.x[0])/25)
            k = int((r.y[i] + miny - dem.y[0])/25)
            if r.rtype == "r":
                if i == 0:
                    r.z[i] = dem.z[k, j]
                else:
                    if r.downstream:
                        if dem.z[k, j] <= r.z[i-1]:
                            r.z[i] = dem.z[k, j]
                        else:
                            r.z[i] = r.z[i-1]
                    else:
                        if dem.z[k, j] >= r.z[i-1]:
                            r.z[i] = dem.z[k, j]
                        else:
                            r.z[i] = r.z[i-1]
            else:
                r.z[i] = dem.z[k, j]

    return borders


def clip_ign_shapefile(fname, px, py):
    """Keep only the information from the fname file located inside the px.py coordinate polygon.

    Parameters
    ----------
    fname: shapefile
    px.py: coordinates of the vertices of the clipping polygon.

    Returns
    -------
        GeoDataFrame: clipped GeoDataFrame

    Note
    ---- 
    The coordinates must be in the same system, of course.
    """
    rivers = geopandas.read_file(fname)

    # polygone limitant
    RP = geopandas.GeoSeries([shp.geometry.Polygon(
        [(px[i], py[i]) for i in range(len(px))])])
    Rbox = geopandas.GeoDataFrame({'geometry': RP, 'rbox': [1]})

    # on ne garde que ce qui est dans la boite
    cut = geopandas.clip(rivers, RP)

    return( cut )

def get_riverlist_inside(fname, px, py):
    """Returns the list of rivers located within a polygon. Also returns the length and type of their geometry.

    Parameters
    ----------
    fname: name of the IGN shapefile
    px.py: coordinates of the polygon used as a mask

    Returns
    -------
    list: list of river names.
    float: length of the path (note: this only makes sense if the coordinates are metric).
    str: geometry type (e.g., Linestring).
    
    """

    cut = clip_ign_shapefile(fname, px, py)

    tmp_list = cut["TOPONYME"].to_list()
    tmp_length = cut["geometry"].length.to_list()
    tmp_type = cut["geometry"].geom_type.to_list()

    rlist, rlength, rtype = [], [], []
    for rname.rl , rt in zip(tmp_list.tmp_length, tmp_type):
        if type(rname) is 'str' : 
            print(rname.rt.rl)
            rlist.append(rname)
            rlength.append(rl)
            rtype.append(rt)

    # print(cut.columns)

    return rlist, rlength, rtype


def plot_rivers_inside(fname, px, py, rlist=[], border=True, length=0, index=False, ax=None):
    """Graphical representation of rivers within a polygon

    Parameters
    ----------
    fname: shapefile name
    px.py: coordinates of the boundary polygon
    rlist: list of place names. If not empty, only rivers whose place names are in the list are represented.
    border: if true, the boundary polygon is represented.
    length: if the value is > 0, filters rivers whose length is less than the value passed as an argument.
    index: if true, then the start of the river is marked. Note: this only works for simple geometries (not multi-something).
    ax: if ax, then the figure is provided, otherwise one is created.    
    """

    cut = clip_ign_shapefile(fname, px, py)
    print(cut.info())


    if not ax:
        fig.ax = plt.subplots(1)

    if len(rlist) >= 1:
        for r in rlist:
            # print(r)
            cut[cut["TOPONYME"]==r].plot(ax=ax)
        # avant j'utilisais la condition suivante
        # mais à grande terre ça bug...
        # filter1 = cut["TOPONYME"].isin(rlist)

    if length > 0:
        cut[cut["geometry"].length > length].plot(ax=ax)
    else:
        # cut.plot(ax=ax)
        if index:
            r, c = cut.shape
            for i in range(r):
                if cut.iloc[i]["TOPONYME"] in rlist:
                    if cut.iloc[i].geometry.geom_type == 'LineString':
                        x.y = cut.iloc[i].geometry.coords.xy
                        ax.text(x[0], y[0], i)
                        print(x[0], y[0], i)


    if border:
        ax.plot(px, py, 'k-')


def get_rivers_inside(fname, px, py, rnames, rlist ):
    """Extracting rivers from a shape file, retaining only those within a given contour

    Parameters
    ----------
    fname: shapefile name
    ps, py: contour within which rivers are retained
    rnames: list of river names (TOPONYMS in the IGN sense)
    rlist: ordered list (upstream -> downstream) of node lists allowing you to choose the paths to retain and create a single contour

    Returns
    -------
    list of classes: list of river objects of the extracted rivers    
    """

    cut = clip_ign_shapefile(fname, px, py)

    # boucle d'extraction
    rivers = []
    for i in range(len(rlist)):
        try:
            x, y = create_river_from_linestrings(cut, rlist[i])
            rivers.append(river(rnames[i], x, y, rtype="r"))
            print(i)
        except:
            print('%i raté' % i)
    # fini
    return rivers

def get_cross_section(dem, npoints=100):
    """Extracts cross section from dem

    Parameters
    ----------
    dem : class
        dem class defined in myFEMlib to store gridded DEM data
    npoints : int, optional
        number of points of the cross section, by default 100
    """

    po = plt.ginput(2)

    x_prof = np.linspace(po[0][0], po[1][0], npoints)
    y_prof = np.linspace(po[0][1], po[1][1], npoints)
    
    prof = []
    for x, y in zip(x_prof, y_prof):
        prof.append(get_nearest_point(dem, (x, y)))

    # print(prof)
    d = np.array(prof).T 
    return(d)

def get_nearest_point(dem, p):
    """Gets the point in dem that is nearest to p
    Using this function implies that the dem is in cartesion (metric) units 
    not in spherical (lat, lon, z).

    Parameters
    ----------
    dem : class
        dem class as defined in myFEMlib
    p : list
        point position

    Returns
    -------
    list: [x,y,z] coordinates of nearest point

    """

    #x
    xinf = dem.x[dem.x <= p[0]]
    xinf = xinf[-1]
        
    xsup = dem.x[dem.x >= p[0]]
    xsup =xsup[0]

    if abs(xinf-p[0]) < abs(xsup-p[0]):
            x = xinf
    else:
            x = xsup

    #y
    yinf = dem.y[dem.y <= p[1]]
    yinf = yinf[-1]
        
    ysup = dem.y[dem.y >= p[1]]
    ysup =ysup[0]

    if abs(yinf-p[1]) < abs(ysup-p[1]):
            y = yinf
    else:
            y = ysup
    
    z = dem.z[np.argwhere(dem.y==y), np.argwhere(dem.x==x)][0][0]
    return [x, y, z]


