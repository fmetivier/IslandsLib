
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

# Poisson coefficient
# 2 Delta\rho/\rho R/K
# Densities
rho = 1000 # freshwater
rhos = 1025 # saltwater
# Recharge
R = 0.19 # m/year
R = R / 365.25 / 86400 # recharge m/s
# Conductivity
K = 2.5e-5 # conductivity m/s

fi = 2*R*(rhos-rho)/K/rhos

#####################
# Solve the problem
#####################

u, Th, X, Y, = il.IslandLens( islands, fname, ttype, fi, sub_sampling , clockwise)


########################
# Visualize the results
########################
plt.show()