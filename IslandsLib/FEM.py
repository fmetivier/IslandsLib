#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Finite Element Functions

* Mesh generation using Triangle and pyFreeFem
* Functions to solve Laplace or poisson equations
* Graphic functions

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

from IslandsLib.Classes import *




def create_outer_vbn(borders, for_FF=False):
    """Create vertices, segments, and nodes from a list of borders
    for triangulation

    Parameters
    ----------
    borders: list 
        list of river objects that form the outer boundary of the area
    for_FF: boolean 
        if true, add an end node for FreeFem

    Returns
    -------
    list of floats: vertices = list of contour vertices
    list of ints: segments = list of contour segments defined by their end nodes
    list of ints: nodes = contour nodes


    Note
    ----
    * I Remove the last segment that doesn't correspond to anything and add a closing segment.
    * A flag should be added in case the contour is naturally closed.
    * be careful to close the borders or else triangle won't succeed in  segments creation and everything crashes.
    """
    vertices = []
    segments = []
    nodes = []

    start = 0

    tot_nodes = 0
    for b in borders:
        tot_nodes += len(b.x)

    for b in borders:
        print("adding border ", b.rname)
        bnodes = []
        for i in range(len(b.x)):
            vertices.append((b.x[i], b.y[i]))
            segments.append([start + i, start + i + 1])
            bnodes.append(start + i)
        if for_FF:
            if start+i+1 < tot_nodes:
                bnodes.append(start+i+1)
            else:
                bnodes.append(0)
        b.nodes = bnodes
        nodes.append(bnodes)
        start += len(b.x)

    segments.pop(-1)
    segments.append([len(segments), 0])

    return vertices, segments, nodes

def add_lakes_vbn(vertices, segments, lakes, len_xpt):
    """add lakes as  internal boundaries and create lake nodes

    Parameters
    ----------
    vertices : list
        list of vertices
    segments : list
        list of segments
    lakes : list of rivers
        island's lakes
    len_xpt: int
        original length of outer border (island contour)

    Returns
    -------
    list
        vertices
    list
        segments
    list
        lake nodes
    list
        lake center coordinates for holes calculation in triangle
    """
    
    lnodes = []
    holes  = []
    for lake in lakes:
        nodes = []
        for i in range(len(lake.x)):
            vertices.append((lake.x[i], lake.y[i]))
            # Freefem indices begin at 1 !
            segments.append([len_xpt + i, len_xpt + i + 1])
            nodes.append(len_xpt + i)

        segments.pop(-1)
        segments.append([len_xpt + len(lake.x) - 1, len_xpt])

        nodes.append(len_xpt)
        print("==============%s nodes===========" % (lake.rname))
        print(nodes)
        lnodes.append(nodes)

        lake.nodes = nodes
        holes.append([np.mean(lake.x), np.mean(lake.y)])
        len_xpt += len(lake.x)

    print("Done with lakes")

    return vertices, segments, lnodes, holes

def add_inner_vbn(vertices, segments, rivers):
    """Add inner boundaries (rivers) to the list of vertices
    and segments and  create river nodes

    Parameters
    ----------
    vertices: list
        outer contour vertices
    segments: list
        outer contour segments
    rivers: list of rivers
        inner rivers

    Returns
    -------
    lists of real numbers: vertices = list of contour vertices and inner rivers
    list of integers: segments = list of contour segments and inner rivers defined by their endpoint nodes
    list of integers: rnodes = contour nodes and inner rivers
    
    """
    rnodes = []
    for i in range(len(rivers)):
        print("adding river ", rivers[i].rname)
        if i == 0:  # first river get existing segments list length
            bl = len(segments)
        elif i > 0:  # for other rivers retrieve last segment to cut from next river
            segments.pop(-1)
        nodes = []
        for j in range(len(rivers[i].x)):
            vertices.append((rivers[i].x[j], rivers[i].y[j]))
            segments.append([bl + j, bl + j + 1])  # Freefem
            nodes.append(bl+j)

        rivers[i].nodes = nodes
        bl += len(nodes)  # keep the good number of nodes
        rnodes.append(nodes)
        segments.pop(-1)

    return vertices, segments, rnodes


def is_inside(node, check_list):
    """checks if a node is inside a list

    Parameters
    ----------
     node: int
        the number of a node of the triangulation
     check_list: list of int
        a list of node numbers

    Returns
    -------
    bool: True/False
    
    """

    if node in check_list:
        return True
    else:
        return False


def create_FreeFem_mesh(xv, yv, mesh, borders, adapt=False, plot=False):
    """Creates the FreeFem mesh from the mesh obtained with triangle

    Parameters
    ----------
    xv, yv: arrays 
        coordinates of the triangle nodes
    mesh: triangle mesh
    borders: dict 
        external and internal boundaries
    adapt: bool
        ​​if true, then FreeFem creates an adaptmesh
    plot: bool
        if true, plots the mesh


    Returns
    -------
    pyFreeFem object: Th object for FreeFem
    List of Th: list of boundaries retrieved from Th
    

    Note
    ----
    in Th x and y are stored according to segment number whereas in triangle they are stored according to node number
    Boundary segments are stored in a dictionnary whose keys are the border labels and vals are the segments.

    """

    # try to create a pyff mesh from the triangle mesh
    Th = pyff.TriMesh(xv, yv, mesh["triangles"])
    if borders:
        for b in borders:
            # add the dictionnary
            Th.add_boundary_edges(boundary_edges=b.nodes, label=b.rname)

    # else:
    #     Th.add_boundary_edges(border_dict)

    if adapt:
        Th = pyff.adaptmesh(Th, iso=1.)

    # label conversion
    print("============================================")
    print("Label conversion")
    print(Th.get_boundary_label_conversion()[0])
    print("Retrieving FreeFem mesh boundaries...")
    
    # on récupère la bordure selon freefem
    # print(Th.get_boundaries())
    Th_boundaries = []  # conteneur des limites

    segs = Th.get_boundaries()
    for key, val in segs.items():
        print("============================================")
        print(key)
        print(val)

        seg = segs[key][0]
        xTh = []
        yTh = []
        for s in seg:
            xTh.append(Th.x[s])
            yTh.append(Th.y[s])

        Th_river = river(key, np.array(xTh), np.array(yTh))
        Th_river.nodes = seg
        Th_boundaries.append(Th_river)

    if plot:
        plt.figure()
        Th.plot_triangles(color='k', alpha=.2, lw=.5)
        Th.plot_boundaries()
        # Th.plot_edges(color='r', labels='label')
        # Th.plot_triangles( labels = 'index' )
        # Th.plot_nodes(labels='index', color='tab:blue')
        plt.legend(title='Boundary label')

        plt.axis('square')
        # plt.savefig('FreeFemMesh.svg',bbox_inches='tight')
        plt.savefig('FreeFemMesh.pdf',bbox_inches='tight')

    return Th, Th_boundaries


def update_node_list(borders, rivers, nodes, mesh):
    """loop to add the boundary nodes created during mesh generation
    in the node lists


    Parameters
    ----------
        borders: list 
            list of border classes (outer borders)
        rivers: list 
            list of rivers (inner borders)
        nodes: list 
            list of nodes
        mesh: dictionnary 
            mesh dictionnary (be clearer...) 

    Returns
    -------
        list: list of nodes
        list: list of new_nodes
        list of list: list of border names


    Note
    ----
    pb d'affectation des noeuds.
    Il faut empêcher la double affectation sinon on est dans le caca !


    """

    print("Updating nodes lists...")

    bname_list = []
    for b in borders:
        # bname_list.append("contour")
        bname_list.append(b.rname)
    if rivers:  # si les rivières ajoute les listes de noeuds de rivières au noeuds de contour
        for r in rivers:
            nodes.append(r.nodes)
            bname_list.append(r.rname)

    new_nodes = []
    for n in nodes:
        new_nodes.append([])

    remaining = mesh["segments"]
    loop_number = 0

    while len(remaining) > 0:
        still_remaining = []
        for b in remaining:  # for each border segment
            # print(b)
            sn = b[0]  # starting  node
            en = b[1]  # ending  node
            continuer = True
            affected = False
            i = 0  # indice des noeuds
            while continuer:
                if is_inside(sn, nodes[i]) or is_inside(en, nodes[i]):
                    if not is_inside(sn, nodes[i]):
                        nodes[i].append(sn)
                    if not is_inside(en, nodes[i]):
                        nodes[i].append(en)
                    new_nodes[i].append(b)
                    continuer = False
                    affected = True
                else:
                    i += 1
                    if i >= len(nodes):
                        continuer = False
            if not affected:
                still_remaining.append(b)

        remaining = still_remaining
        loop_number += 1
        print(loop_number, len(mesh["segments"]), len(still_remaining))

    return nodes, new_nodes, bname_list


def build_border_dictionnary(mesh, nodes, bname_list):
    """obsolete: Creates a full border dictionary
    from the triangle-generated mesh
    """

    nodes_order = [[0, 1], [1, 2], [-1, 0]]
    edge_order = [0, 1, 2]
    border_dict = {}

    for b in mesh["segments"]:  # for each border segment
        # print(b)
        sn = b[0]  # starting  node
        en = b[1]  # ending  node
        if [en, sn] in mesh["segments"].tolist():
            print(mesh["segments"].tolist().index([en, sn]))
        for j in range(len(nodes)):
            if is_inside(sn, nodes[j]) or is_inside(en, nodes[j]):
                for i in range(len(nodes_order)):
                    try:
                        # try i sn.en order
                        bt = mesh["triangles"][:, nodes_order[i]].tolist().index([
                                                                         sn, en])
                        border_dict[(bt, edge_order[i])] = bname_list[j]
                    except:
                        try:
                            bt = mesh["triangles"][:, nodes_order[i]].tolist().index([
                                                                             en, sn])
                            border_dict[(bt, edge_order[i])
                                        ] = bname_list[j]
                        except:
                            pass

    for key, val in border_dict.items():
        print(key, val)
    print("Done with borders")
    return border_dict


def is_on(xv, yv, seg1, seg2, pval=False, epsilon=1e-5):
    """ Tests if a segment seg2 belongs to another segment seg1
    
    Used to retrieve additional edge nodes created by triangles when asked for an adapted 
    triangulation (equivalent to a FreeFem adaptmesh) under constraints

    Parameters
    ----------
    xv, yv: arrays 
        real coordinate vectors
    seg1: list 
        containing segment
    seg2: list 
        segment to test
    pval: bool
        if true, prints the slope coefficient values ​​of the two segments
    epsilon: real 
        margin of error tolerated in comparing the two segments. epsilon must be non-zero, due to numerical approximations, but very small to avoid confusion (double membership, for example). The default value seems to be suitable...

    Returns
    -------
    bool: true or false
     """
    if np.diff(xv[seg1]) == 0 and np.diff(xv[seg2]) == 0:
        if max(xv[seg2]) <= max(xv[seg1]) and min(xv[seg2]) >= min(xv[seg1]):
            if max(yv[seg2]) <= max(yv[seg1]) and min(yv[seg2]) >= min(yv[seg1]):
                return True
    elif np.diff(xv[seg1]) != 0 and np.diff(xv[seg2]) != 0:
        a = np.diff(yv[seg1])/np.diff(xv[seg1])
        b = np.diff(yv[seg2])/np.diff(xv[seg2])
        if pval:
            print(a, b)

        if b[0] - epsilon <= a[0] <= b[0] + epsilon:
            if pval:
                print('ok')
            if max(xv[seg2]) <= max(xv[seg1]) and min(xv[seg2]) >= min(xv[seg1]):
                if max(yv[seg2]) <= max(yv[seg1]) and min(yv[seg2]) >= min(yv[seg1]):
                    return True
        else:
            return False
    else:
        return False


def order_nodes(xv, yv, mesh, borders, plot=False):
    """Classifies the contour nodes.


    Parameters
    ----------
    xv, yv: arrays
        vectors of the triangle mesh node coordinates
    mesh: dictionary of vertices, segments, and nodes
        the triangle mesh 
    borders; list of rivers
        the boundaries
    plot: bool 
        if true, represents the graph of the oriented nodes

    Returns
    -------
    list or rivers: borders whose nodes have been reclassified

    Note
    ---- 
    This operation is very important because it is a necessary condition
    for successfully converting a triangle mesh to a FreeFem mesh.
    """

    segcount = 0
    segok = []
    for b in borders:
        nodes = b.nodes.copy()
        for j in range(len(b.nodes)-1):
            s = [b.nodes[j], b.nodes[j+1]]
            for ns in mesh["segments"]:
                # if ns[0] == 34 or ns[1] == 34:
                #     print(is_on(xv, yv, s, ns, pval=True))
                if is_on(xv, yv, s, ns):
                    b.segments.append(list(ns))
                    segcount += 1
                    segok.append(list(ns))
                    for n in ns:
                        if n not in nodes:
                            nodes.append(n)
                else:
                    if is_on(xv, yv, s, ns[::-1]):
                        b.segments.append(list(ns[::-1]))
                        segcount += 1
                        segok.append(list(ns))
                        for n in ns[::-1]:
                            if n not in nodes:
                                nodes.append(n)

        tmp_nodes = nodes.copy()
        compare = []

        sn = 0
        en = 1
        continuer = True
        print("=================================")
        while continuer:
            s = [tmp_nodes[sn], tmp_nodes[en]]
            cn = 0
            for j in range(en+1, len(tmp_nodes)):
                n = tmp_nodes[j]
                if is_on(xv, yv, s, [s[0], n]):
                    cn += 1
                    tmp_nodes.pop(j)
                    tmp_nodes.insert(en, n)

            print(cn, " ", end="")
            # print(s, nodes, tmp_nodes, cn)
            if cn == 0 and en < len(tmp_nodes)-1:
                sn += 1
                en += 1
            elif en == len(tmp_nodes)-1:
                continuer = False

        b.nodes = tmp_nodes.copy()
        print(b.rname, b.segments, b.nodes)
    print(segcount, len(mesh["segments"]))
    for ns in mesh["segments"]:
        if list(ns) not in segok:
            print(list(ns))
    if plot:
        fig, ax = plt.subplots(1)
        for b in borders:
            ax.plot(xv[b.nodes], yv[b.nodes], 'o-')
            for n in b.nodes:
                ax.text(xv[n], yv[n], n)

    return borders


def prepare_boundaries(borders,  lakes = None, rivers = None, for_FF=False):
    """Create the border contour from borders
    
    Parameters
    ----------
    borders: list of rivers
        external boundaries
    rivers: list of rivers
        internal rivers
    lakes: list of rivers
        island's lakes
    for_FF: bool
        if true, arranges segments and end nodes to match FreeFem mesh requirements

    Returns
    -------
    list of rivers: complete external and internal borders
    list of pairs of real numbers: vertices of border nodes
    list of pairs of integer numbers: border segments
    list of pairs of intergers: holes coordinates
    
    Note
    ----
    this supposes that the borders are oriented in a correct continuous fashion

    """
    print("Creating border contour...")

    #################################################
    #
    # Create vertices, borders (segments) and nodes
    #
    #################################################

    print("creating vertices, segments and nodes")
    vertices, segments, nodes = create_outer_vbn(borders, for_FF)

    ###########################################################
    #
    # Add river segments as internal boundaries
    # in progress
    #
    ###########################################################

    holes = []
    if lakes:
        print("Adding lakes as boundaries")
        len_xpt = len(borders[0].x)
        vertices, segments, lnodes, holes = add_lakes_vbn(vertices, segments, lakes, len_xpt)
    
    if rivers:
        print("Adding rivers as boundaries")
        vertices, segments, rnodes = add_inner_vbn(vertices, segments, rivers)

    if rivers:
        for r in rivers:
            borders.append(r)
    if lakes:
        for l in lakes:
            borders.append(l)
    for b in borders:
        print(b.rname, b.nodes[0], b.nodes[-1])

    # print(segments)

    return borders, vertices, segments, holes


def create_triangle_mesh(borders, vertices, segments, holes=None, plot=False, ttype='p'):
    """Create the mesh using the triangle library


    Parameters
    ----------
    borders: lits of rivers 
        the external and internal borders
    vertices: the boundary vertices
    segments: the boundary segments
    plot: bool
        if true, represents the resulting mesh
    ttype: str
        triangulation type, by default p is the simplest triangulation without adaptation and creation of new nodes. This is typically the one I use before a FreeFem adaptmesh.

    Returns
    -------
    array of real numbers: xv, yv, vertex coordinates
    mesh (dictionary): mesh, mesh produced by triangle
    list of rivers: borders, recalculated boundaries with reordered nodes    
    
    Note
    ----
    triangle allows constrained triangulation with internal borders.
    Personally, I find that the triangle meshes after adaptation are more balanced than those of FreeFem. 
    O.D. says I just don't know how to use adaptmesh. He's probably right.
    """
    print("Creating mesh...")

    if holes is not None:
        init = dict(vertices=vertices, segments=segments, holes=holes)
    else:
        init = dict(vertices=vertices, segments=segments)

    mesh = tr.triangulate(init, ttype)

    for key, val in mesh.items():
        print("=========================")
        print(key)
        print(val)

    xv = np.array(mesh["vertices"]).T[0]
    yv = np.array(mesh["vertices"]).T[1]

    #################################################
    #
    # Create ordered nodes list for Freefem
    #
    #################################################
    print("Ordering nodes...")
    borders = order_nodes(xv, yv, mesh, borders, plot = plot)

    #################################################
    #
    # plot
    #
    #################################################

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.axis('equal')

        # plot mesh
        for tri in mesh["triangles"]:
            for j in range(len(tri)):
                # x coordinate of first node
                XA = mesh["vertices"][tri[j-1]][0]
                XB = mesh["vertices"][tri[j]][0]  # x coordinate of second node
                # y coordinate of first node
                YA = mesh["vertices"][tri[j-1]][1]
                # y coordinate of seconde node
                YB = mesh["vertices"][tri[j]][1]
                ax.plot([XA, XB], [YA, YB], c="lightgray", alpha=0.2)

        # plot segments
        for s in mesh["segments"]:
            ax.plot(xv[s], yv[s], 'k--')
        for i in range(len(borders)):
            col = 'C%i' % (i)
            for s in borders[i].segments:
                ax.plot(xv[s], yv[s], '-', color=col)

        for b in borders:
            for n in b.nodes:
                ax.plot(xv[n], yv[n], 'o', color='C0', ms=12)
                ax.text(xv[n], yv[n], n, color='white',
                        fontsize=7, ha='center', va='center')

    print("Done with mesh creation")
    if plot:
        plt.savefig('triangleMesh.svg', bbox_inches='tight')
        plt.savefig('triangleMesh.pdf', bbox_inches='tight')

    return xv, yv, mesh, borders

def solve_fem(Th, bcs=None, fi=1e-4):
    """Solve Poisson equation

    Parameters
    ----------
    Th: Th object, 
    bcs: dic, 
        boundary conditions dictionary,
    fi: float, optional, 
        right-hand side of the Poisson equation :math:`(\Delta\\rho/\\rho) (P/K)`,

    Returns
    -------
    array: u, the solution for the potential (here, h^2)
    object: Th, the Th object that contains all the mesh information

    Note
    ----
    For now, :math:`fi = (\Delta\\rho/\\rho) (P/K)` is a constant, but nothing prevents non-uniform precipitation distributions
    
    """

    script = pyff.InputScript(Th=Th)

    script += '''
    fespace Vh( Th, P1 );
    Vh u, v;
    '''

    script += pyff.OutputScript(Th='mesh')

    db = Th.get_boundary_label_conversion()[0]
    print(db)

    matrice = {'stiffness': 'int2d(Th)( dx(u)*dx(v) +  dy(u)*dy(v) )',
                'Gramian': '-int2d(Th)( u*v)'}
    
    for key, val in db.items():
        matrice[key] = '-int1d(Th, %i)( u*v )' % (val)
        print("adding", matrice[key])

    script += pyff.VarfScript(**matrice)

    #####################################################
    # Get PDE Integrals and mesh
    #####################################################

    ff_output = script.get_output()
    Th = ff_output['Th']
    S = ff_output['stiffness']
    G = ff_output['Gramian']

    #####################################################
    # Boundary conditions
    #####################################################

    f = np.array([fi]*len(Th.x))

    # BC boundary
    bc_vals = {}
    print("=====================================")
    if bcs is not None:  # if boundary conditions dic passed as argument
        print("bcs given")
        for key, val in db.items():  # for each boundary
            # if condition is a list or is None
            if bcs[key] is not None:
                print(key, bcs[key])
                wanted_value = bcs[key]
                bc_vals[key] = boundary_values(val, wanted_value, Th)
            else:
                print("null flux on ", key)
                bc_vals[key] = None
    else:  # bcs = None 0 on every border
        print("No bcs given assuming 0 on every border")
        for key, val in db.items():
            wanted_value = [0]*len(Th.get_boundaries()[val][0])
            bc_vals[key] = boundary_values(val, wanted_value, Th)

    #####################################################
    # Solve PDE
    #####################################################

    epsilon = 1e-4
    e = 1./epsilon

    M = - S
    Source = G*f
    
    
    for key, val in db.items():
        if bc_vals[key] is not None:
            M += e*ff_output[key]
            Source += e*ff_output[key]*bc_vals[key]

    # Solve
    u = spsolve(M, Source)

    return u, Th


def boundary_values(label_number, wanted_value, Th):
    """Defines a set of value at a specific boundary

    * First step is to pick out the corresponding nodes in the mesh
    * Second step is to assign the values

    Parameters
    ----------
    label_number: Label of the boundary to modify (integer defined by FreeFem)
    wanted_values: List of boundary conditions
    Th: Th object containing the FreeFem mesh

    
    Returns
    -------
    array: list of boundary values corresponding to each nodes (boundary grammian)
    

    Note
    ----
    Source: Olivier Devauchelle & Céleste Romon

    
    """
    label = label_number
    segments = Th.get_boundaries()[label]  # Access the correct border
    if len(segments) > 1:  # si ff_output[Th] a foutu la merde !!!
        boundary_segments = []
        for seg in segments:
            for s in seg:
                if s not in boundary_segments:
                    boundary_segments.append(s)
    else:
        boundary_segments = segments[0]
    print(boundary_segments)
    L = np.array(Th.x)  # define vector
    c = 0
    for nodes in boundary_segments:  # node order is in increasing x value
        L[nodes] = wanted_value[c]  # add wanted value to specific nodes
        c = c+1

    print(c)
    return L




###########################################################
#
# 4- Fonctions de représentation graphique
#
###########################################################


def plot_u(Th, u, ax):
    """Tricontour plot of u the solution of poisson equaiton obtained with solve_fem

    Parameters
    ----------
     Th: Th object
     u: solution for the potential

    """
    if ax == None:
        fig, ax = plt.subplots(1)
    # plt.scatter(Th.x.Th.y, c=u)
    ax.tricontourf(Th, u, levels=20, cmap='Blues')
    Th.plot_boundaries(color='black')

def plot_psi(Th, psi, ax):
    """Tricontour plot de u la solution de poisson obtenue par solve_fem
    

    Parameters
    ----------
     Th: objet Th
     u: potential solution to be plotted

    Note
    ----
    The difference with plot_u is that the solution was obtained by solving for the symmetric boundary conditions hence 
    for the flow lines
    """
    if ax == None:
        fig.ax=plt.subplots(1)
    # plt.scatter(Th.x.Th.y, c=u)
    ax.tricontour(Th, psi, levels=10, cmap='Reds')
    Th.plot_boundaries(color='black')

def grid_plot_u(fig, ax, Th, u, stream_contour = 12, pstep=80 , fval = 0, dupuit=True, dxy = 10):
    """Graphical representation of :math:`\Phi = \sqrt{u}` on an interpolated grid
    with streamlines and contour mask

    Parameters
    ----------
    fig: figure on which to plot
    ax: axis on which to plot
    Th: Th object
    u: Poisson solution
    stream_contour: int, default 12
        contour number used to represent the streamlines
    pstep: int, default 80
        steps for the streamlines
    fval: int, default 0
        fill value for gridding
    dupuit: bool, default True
        is the solution of dupuit-boussinesq type
    dxy: int, default 10
        step in meters for the grid 

    Returns
    -------
    figure: fig 
    axis: ax 
    arrays: X, Y, Z
    arrays: dx, dy 
    object: itp regular grid interpolator for Z

    """
    #bordures de la zone de calcul
    xmin = min(Th.x)
    xmax = max(Th.x)
    ymin = min(Th.y)
    ymax = max(Th.y)

    x = np.arange(xmin, xmax, step = dxy)
    y = np.arange(ymin, ymax, step = dxy)

    dx = np.diff(x)[0] # still keep it to check and return
    dy = np.diff(y)[0]

    # création de la grille régulières
    X, Y = np.meshgrid(x, y)
    Z = griddata((Th.x, Th.y), u, (X, Y), fill_value = fval)
    Z = np.nan_to_num(Z) # replace nans by 0

    if dupuit:
        Z = np.sqrt(Z) 

    itp = RegularGridInterpolator( (x, y), Z.T, method='linear')

    cs = ax.contourf(X, Y, Z, levels=20, cmap='Blues_r')
    cbar = fig.colorbar(cs, label="Piezometric Height (m)")

    cs2 = ax.contour(X, Y, Z, levels=20, colors='k', linewidths=0.1)
    V, U = np.gradient(Z, y, x)
    # ax.clabel(cs, inline=True, fontsize=10)
    cbar.add_lines(cs2)

    #
    # Pour les streamlines on choisi de les créer à partir d 'une équipotentielle
    # en espaçant les points de façon ~ régulière
    #
    pos = list(cs.allsegs[stream_contour][0])  # le choix du contour est arbitraire et à tester !
    start = []
    for p in pos[::pstep]:  # idem le choix de l'echantillonnage est arbitraire et à tester
        # plt.plot(p[0], p[1], 'o', color='C3')
        for i in range(len(x)-1):
            if x[i] <= p[0] < x[i+1]:
                cx = i
        for j in range(len(y)-1):
            if y[j] <= p[1] < y[j+1]:
                cy = j

        start.append([x[cx], y[cy]])
    print("starting points", start)

    ax.streamplot(X, Y, -U, -V, color='C1', start_points = start, 
                  linewidth=1, density=2)

    Th.plot_boundaries(ax=ax, color='black')


    plt.axis("square")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    return fig, ax, X, Y, Z, dx, dy, itp

def add_mask(ax, px, py, outer=False, extent=[]):
    """Mask to remove the very ugly and especially false extrapolated part

    Parameters
    ----------
    ax: axis of the graphical representation
    px.py: outline for the mask
    outer: if true, use extent to draw an outer mask between the polygons
    formed by extent and outline. This is mainly used to remove the extrapolations when transitioning
    from a triangular mesh to a square grid by interpolation.
    extent: extension of the grid.

    Returns
    -------
    p (polygon) of the mask

    Note
    ----
    Adapted from https://www.matecdev.com/posts/shapely-point-polygon.html
    
    """

    p = shp.geometry.Polygon([[px[i], py[i]]
                              for i in range(len(px))])


    if outer:
        xmin = extent[0]
        xmax = extent[1]
        ymin = extent[2]
        ymax = extent[3]
        po = shp.geometry.Polygon(
            [[xmin, ymin], [xmin, ymax], [xmax, ymax], [xmax, ymin]])
        pdiff = po.difference(p)

        mask = geopandas.GeoSeries([pdiff])
    else:
        mask = geopandas.GeoSeries([p])

    mask.plot(ax=ax, color='white', ec='white', zorder=2)

    return p

def mask_data(X,Y,Z,p):
    """Masks the Z array using X,Y and p.
    this is useful because the interpolation of u, Th on a regular grids always leads to some smear outside the island contour

    Parameters
    ----------
    X : 2D array
        X positions
    Y : 2D array
        Y positions
    Z : 2D array
        Unmasked water table
    p : Shapely polygon
        island contour used as a mask

    Returns
    -------
    2D array
        masked water table
    """

    n,m = np.shape(Z)
    for i in range(m):
        for j in range(n):
            if not p.contains(shp.geometry.Point(X[j,i],Y[j,i])):
                Z[j,i] = 0

    return Z
###########################################################
#
# 6-Divers
#
###########################################################


def list_errors():
    """Small function to get the description of python errors.
    """
    import os
    import errno
    from pprint import pprint

    pprint({i: os.strerror(i) for i in sorted(errno.errorcode)})

