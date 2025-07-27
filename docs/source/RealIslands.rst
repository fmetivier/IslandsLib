Real Islands
************

Let us now turn to real islands. 
This page is regularly updated when new islands are included and also when new data enable 
to strengthen the modeling and predictions

.. |date| date::

*Last modification* |date|

Summary table
=============

.. list-table:: Parameters used and water balance obtained
  :header-rows: 1

  * - Island
    - R (:math:`m/d`)
    - K (:math:`m/d`)
    - por
    - Volume (:math:`m^3`)
    - Surface (:math:`m^2`)
    - Recharge (:math:`m^3/yr`)
  

  * - Avatoru
    - 0.000828
    - 14.6
    - 0.25
    - 1.38e6 
    - 8.73e5
    - 2.64e5

  
  * - Petite Terre
    - 0.000548
    - 3.456000
    - 0.30
    - 15.11e6
    - 10.48e6
    - 2.97e6



Avatoru Motu in Rangiroa (Polynésie Française)
===============================================


Parameters
----------

.. code:: python

    #IslandLens Param
    sub_sampling = 3
    clockwise = True
    ttype = 'pq33a1000'
    R = 0.000827698# m/j
    por=0.25
    motu = 864
    K= 14.6
    por=0.25
    fi = 2*R*25/K/1025

    islands="Avatoru" 

    fname = "./IslandsLib/data/Contours/Pacific/Polynesie/Rangiroa/Rangirao_%s.txt" % (motu)

    u, Th, X, Y, Z, dx, dy, itp = il.IslandLens( islands, fname, ttype, fi, sub_sampling , clockwise, lakes=None, plot=False)

    il.IslandBalance(islands,Z,dx,dy,R,K,por)

Results
-------

* Values for :math:`R` and :math:`K` were obtained from the meteorologic station of the airport and from the BRGM study of Rangiroa. 
* One piezometer exists on the islands and the average height of the water table coincides with the reconstruction obtained. 

.. figure:: ./figures/Avatoru.svg



Petite Terre and Lake Dziani (Mayotte)
======================================

Lake Dziani is a hyper saline 
volcanic lake located north of the island
its level on average is 0. 

Parameters
----------

.. code:: python 

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


Results
-------

The form of the lens is given  in figure :ref:`fig-mayotte`. 
For a discussion on the shape of the water table and comparison with existing measurements see :cite:t:`metivier2024bilan` (spoiler: it seems to work :))

.. _fig-mayotte:

.. figure:: ./figures/Petite_Terre.svg

    Modeled water table of Petite Terre island in Mayotte


