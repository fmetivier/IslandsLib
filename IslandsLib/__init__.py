#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IslandLib
"""

__author__ = "François Métivier"
__copyright__ = "Copyright 2025"
__license__ = "CC-By-4.0"
__version__ = "1.0"

__all__ = [
    "IslandLib",
    "Classes",
    "IGN",
    "FEM",
    "Misc"
]

#############################
#
# Librairies
#
#############################

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

from IslandsLib.Classes import *
from IslandsLib.IGN import *
from IslandsLib.FEM import *  
from IslandsLib.Misc import *  




def IslandLens(islands = 'Desirade', fname = "../data/Contours/Atlantic/Guadeloupe/Desirade.txt", ttype = "pq33a1000",  fi = 1e-4 , sub_sampling = 20, clockwise=True, lakes = [], plot = True):
    """Poisson Resolution for an Island (closed contour with null potential)


    Parameters
    ----------
        islands: str, optional
            Island name. Defaults to 'Desirade'.
        fname: str, optional 
            Coastline contour data file. Defaults to "../../data/Desirade.txt".
        ttype: str, optional 
            Triangle triangulation parameters. Defaults to "pq33a1000".
        fi: float, optional
            value of the Poisson coefficient (dimensionless). Defaults to 1e-4
        sub_sampling: int, optional
            sampling of the contour (1 = all points are preserved, 10 = one point every 10 points is conserved)
        lakes: boolean, None
            to be used as a dic when lake contours are to be included
        plot: boolean, True
            if True plots the figures

    Returns
    -------
        array: :math:`u = \phi^2` solution of poisson equation 
        Th object: Th domain of the solution 
        2D array: X positions on a grid
        2D array: Y positions on a grid
        2D array: Zm = masked water table height on a (X,Y) grid
        object: itp regular grid interpolator for Z
    """

    # 
    # read contour file
    #
    print(islands)
    f = open(fname,"r")

    lines = f.readlines()
    xi, yi = [], []
    for line in lines:
        if line[0] != "#": # skip header
            l = line.strip('\n').split(',')
            xi.append( float(l[0]) )
            yi.append( float(l[1]) )


    if not clockwise: # invert rotation
        xi = np.array(xi[::-1])
        yi = np.array(yi[::-1])

    # if subsampling
    xi = xi[::sub_sampling]
    yi = yi[::sub_sampling]

    if plot:
        fig,ax = plt.subplots(1)
        ax.plot(xi,yi)
        for i in range(len(yi)//100):
            ax.text(xi[100*i],yi[100*i],i*100)

    b = [river(islands, xi, yi)]

    lakes_list = []
    if lakes:
        for lake in lakes:
            print(lake[0])
            f = open(lake[1],"r")

            lines = f.readlines()
            f.close()
            x, y = [], []
            for line in lines:
                if line[0] != "#": # skip header
                    l = line.strip('\n').split(',')
                    x.append( float(l[0]) )
                    y.append( float(l[1]) )


            if not clockwise: # invert rotation
                x = np.array(x[::-1])
                y = np.array(y[::-1])

            # if subsampling
            x = x[::sub_sampling]
            y = y[::sub_sampling]

            lakes_list.append(river(lake[0],x,y))

    borders, vertices, segments, holes = prepare_boundaries(borders = b, lakes = lakes_list, rivers = [], for_FF=True)
    

    ##############################################################
    #
    # performs triangulation
    # 1) with triangle
    # 2) with Freefem
    #
    ##############################################################

    # ttype = 'p'
    xv, yv, mesh, borders = create_triangle_mesh(
            borders = borders, vertices = vertices, segments = segments,  holes = holes, plot = plot, ttype = ttype )


    print("Creating FreeFem mesh...")
    Th, Th_boundaries = create_FreeFem_mesh(
        xv, yv, mesh, adapt = False, borders = borders, plot = plot)

    # sur une ile le contour est à 0
    # on a lake it is given as a parameter.
    bcs = {}
    for T in Th_boundaries:
        if T.rname == islands:
            bcs[T.rname] = 0 * T.x
        if lakes:
            for lake in lakes:    
                if T.rname == lake[0]:
                    bcs[T.rname] = lake[2]*T.x

    for key, val in bcs.items():
        print(key, val)

    u, Th = solve_fem( Th = Th, bcs = bcs, fi = fi) 

    #####################################
    # Plot the potential and streamlines
    #####################################
    
    X, Y, Z, x, y, dx, dy, itp = create_grid_from_u(Th, u)
    fig, ax = plt.subplots(1)
    fig, ax = grid_plot_u(fig, ax, Th, X, Y, Z)
    p = add_mask(ax, xi, yi, outer=True, extent=[min(x),max(x),min(y),max(y)])
    Zm = mask_data(X,Y,Z,p)
    Zm = Z.copy()



        # plt.savefig("PhiPlot.svg", bbox_inches='tight')
    plt.savefig("%s.pdf" % (islands), bbox_inches='tight')

    plt.close()
    

    # end     
    return u, Th, X, Y, Zm, dx, dy, itp


def IslandBalance(islands = "Désirade", Z = None , dx = 10, dy = 10, R = 0.0025, K = 1, por = 0.25):
    """given the matric Z of water table elevation calculates, lens area, lens water volume and volumetric recharge. 
    Outputs results to csv file

    Parameters
    ----------
    islands : str, optional
        Island name, Defaults to Desirade
    Z : array of floats
        elevation of the water table, defaults to None 
    dx : float
        x-distance between to points in meters, defautls to 10 m
    dy : float
        y-distance between to points in meters, defautls to 10 m
    R : float
        recharge in meters / day. Defaults to 0.0025 m/d
    K : float
        Hydraulic conductivity in meters /day. Defaults to 1 m/d
    por : float
        porosity (dimensionless). Defaults to 0.025

    Returns
    -------
    float: Volume of lens V
    float: Area of lens S
    float: Yearly volumetric Recharge Rec 
        
    """

    h = np.mean(Z[Z>0])

    #Volume of water stored in the lens
    V = np.sum(np.nan_to_num(Z)) * por * dx * dy * 41
    
    # Area of the lens 
    S = len(Z[Z>0]) * dx * dy

    # yearly Volume of  recharge 
    Rech = S * R * 365.25
    
    print("Volume", V, "Surface", S, "Recharge", Rech)

    oname = "output_%s.csv" % (islands)
    f = open(oname,"w")
    f.write("Island,R (m/d),K (m/d),por, $h_{moy}$ (m), Volume (m$^3$),Surface (m$^2$),Recharge (m$^3$/yr)\n")
    f.write("%s,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e\n" % (islands,R,K,por,h,V,S,Rech))
    f.close()

    return V,S,Rech

if __name__ == '__main__':

    u, Th, X, Y, Zm, dx, dy, itp = IslandLens()
    IslandBalance(islands = "Désirade", Z = Zm)

    plt.show()