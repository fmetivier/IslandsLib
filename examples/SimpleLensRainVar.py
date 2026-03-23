
#############################################
#
# Mystery Island with varying precipitations
#
#############################################

# Libraries
import sys 
sys.path.append("./..")

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

import IslandsLib as il

#####################
# Set arguments
#####################

# Name of Island
islands = 'Mystery island'

# Path to filename containing contour
fname = "../data/Examples/MysteryIsland.txt"

#sub sampling
sub_sampling = 2

#clockwise
clockwise = False

# Triangulation settings
ttype = 'pq33a1000'

# Parameters
# Infiltration
R  = 0.19 # m/year
R = R / 365.25 # infiltration m/d
# Conductivity
K = 2.5e-5*86400 # conductivity m/d

#Poisson coefficient calculated for constant rainfall
fi = 2 * R * 25 / K / 1025

#######################################################
# Create variable fi matrix and use it as interpolator
#######################################################

f = open(fname,"r")

lines = f.readlines()
xi, yi = [], []
for line in lines:
     if line[0] != "#": # skip header
          l = line.strip('\n').split(',')
          xi.append( float(l[0]) )
          yi.append( float(l[1]) )

x = np.linspace(min(xi), max(xi),100)
y = np.linspace(min(yi), max(yi),100)

X,Y = np.meshgrid(x,y, indexing='ij')
Z = X*0 + 1

nx,ny = np.shape(Z)
for i in range(nx):
     for j in range(ny):
          if Y[i,j] > 0: Z[i,j] = fi
          else: Z[i,j] = 0.5*fi

ri = RegularGridInterpolator((x,y),Z) # interpolator to be passed to IslandLens


#####################
# Solve the problem
#####################

u, Th, X, Y, Zm, dx, dy, itp = il.IslandLens( islands = islands, fname = fname,\
     ttype = ttype, fi = ri, sub_sampling = sub_sampling , clockwise = clockwise,\
     lakes= [], plot=True)



###############################################################
# look at what our precipitation pattern looks loke on the mesh
###############################################################

fig, ax = plt.subplots(1)

boundaries = Th.get_boundaries()
for key,val in boundaries.items():
     boundary_nodes = val[0]
     ax.plot(Th.x[boundary_nodes], Th.y[boundary_nodes], '-', color='C1')

f = np.zeros(len(Th.x))
for i in range(len(f)):
     try:
          f[i] = ri((Th.x[i], Th.y[i]))
     except:
          print(Th.x[i], Th.y[i])

sc = ax.scatter(Th.x, Th.y, c=f )
cbar = plt.colorbar(sc)
cbar.set_label('R/K')
ax.set_xlabel("x-distance (m)")
ax.set_ylabel("y-distance (m)")
plt.savefig('Mystery_rain.svg', bbox_inches='tight')
plt.show()