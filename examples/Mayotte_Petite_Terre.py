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

# Name of Island and path to contour
islands = 'Petite Terre'
island_fname = "../data/Contours/Indian/Mayotte/Mayotte_Petite_Terre.txt"

# Name of lake and path to contour
lake = "Dziani"
lake_fname = "../data/Contours/Indian/Mayotte/Mayotte_Petite_Terre_dziani.txt"

lakes = [[lake, lake_fname, 0]] # list of lakes with file name and elevation of lake (masl)

#sub sampling
sub_sampling = 9

#clockwise
clockwise = True

# Triangulation settings
ttype = 'pq33a10000'

# Parameters
# Infiltration
R = 1. / (365.25) * (1-0.8) #m/d 
#Conductivity
K = 4e-5*86400 # m/d
# Porosity
por = 0.3


fi = 2 * R * 25 / K / 1025

#####################
# Solve the problem
#####################

u, Th, X, Y, Zm, dx, dy, itp = il.IslandLens( islands = islands, fname = island_fname, lakes = lakes, ttype = ttype,  fi = fi , sub_sampling = sub_sampling, clockwise = clockwise, plot = True)

il.IslandBalance("Petite Terre", Zm, dx, dy, R, K, por)