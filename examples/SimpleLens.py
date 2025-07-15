
# Libraries
import sys 
sys.path.append("./..")

import matplotlib.pyplot as plt

import IslandsLib as il

#####################
# Set arguments
#####################

# Name of Island
islands = 'Desirade'

# Path to filenmae containing contour
fname = "../data/Contours/Guadeloupe/Desirade.txt"

#sub sampling
sub_sampling = 50

#clockwise
clockwise = False

# Triangulation settings
ttype = 'pq33a10000'

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

u, Th, X, Y, Zm, dx, dy, itp = il.IslandLens( islands, fname, ttype, fi, sub_sampling , clockwise)

####################################
# Output global balance to csv file
####################################
il.IslandBalance("Desirade", Zm, dx, dy, R, K)