"""
    
    Solution of the Poisson equation for the Island of Petite Terre in Mayotte 
    Petite Terre ahas a large volcanic and hyper saline lake, the Dziani lake which elevation is at sea-level

"""
# Libraries
# import sys 
# sys.path.append("./..")

import matplotlib.pyplot as plt

import IslandsLib as il

#####################
# Set arguments
#####################

# Name of Island and path to filename containing contour
islands = 'Mystery island'
fname = "../data/Examples/MysteryIsland.txt"


# Name of lake and path to contour
lake = "Mystery Lake"
lake_fname = "../data/Examples/MysteryLake.txt"

lakes = [[lake, lake_fname, 0]] # list of lakes with file name and elevation of lake (masl)

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
     lakes= lakes, plot=False)

####################################
# Output global balance to csv file
####################################
il.IslandBalance("Mystery", Zm, dx, dy, R, K)


