
# Libraries
# import sys 
# sys.path.append("./..")

import matplotlib.pyplot as plt

import IslandsLib as il

#####################
# Set arguments
#####################

# Name of Island
islands = 'Mystery island'

# Path to filenmae containing contour
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

fi = 2 * R * 25 / K / 1025

#####################
# Solve the problem
#####################

u, Th, X, Y, Zm, dx, dy, itp = il.IslandLens( islands = islands, fname = fname,\
     ttype = ttype, fi = fi, sub_sampling = sub_sampling , clockwise = clockwise,\
     lakes= [], plot=True)

####################################
# Output global balance to csv file
####################################
il.IslandBalance("Mystery", Zm, dx, dy, R, K)