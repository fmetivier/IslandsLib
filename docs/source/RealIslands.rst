Real Islands
************

Let us now turn to real islands. 


The Desirade island (Guadeloupe)
================================

Procedure
---------
The simplest way to model an island's freshwater lens using IslandsLib is to use the :code:`IslandLens()`
function that is meant to simplify FEM calculations in simple cases of Islands without river networks.

For example suppose we want to model the form of the the Desirade island fresh water lens in the French Antilles Guadeloupe Department.
there are two steps 

* We start by setting some values that will be passed to the main function as arguments namely the name of Island and the path to the file where the island contour is stored 

.. code:: python

    # Name of Island
    islands = 'Desirade'

    # Path to filenmae containing contour
    fname = "../data/Desirade.txt"

* Second, we set the sub sampling rate of the contour. Often the contours are far too precise for FEM computation at the scale of an entire island. 
by default the sub sampling rate is 10 but in the case of Desirade a sampling rate of 50 is used. 

.. code:: python

    #sub sampling
    sub_sampling = 50


.. warning::

    If the contour is complex and the sub sampling rate too low the triangulation will probably fail 


* Third, we define the parameter of the triangulation. 

.. code:: python

    # Triangulation settings
    ttype = 'pq33a10000'

The string ``ttype``  is composed of three caracters p, q, and a followed by numbers the meaning of which is as follows

   *  p: we are making a planar straight line graph (a collection of segments and vertices);
   *  q: we are imposing a constraint of 33° angle for the triangles (equilateral triangles);
   *  a: the area of the triangle is less or equal to 10000 :math:`m^2`. 

 There are other options and for more detail please refer to https://rufat.be/triangle/index.html


* Fourth, we set the recharge :math:`R` and Hydraulic conductivity :math:`K` that together make the poisson coefficient :math:`fi`.

.. code:: python

    # Poisson coefficient
    # 2 Delta\rho/\rho R/K
    # densities
    rho = 1000 # freshwater
    rhos = 1025 # saltwater
    R = 0.19 # recharge m/year
    R = R / 365.25  # recharge m/d
    K = 2.5e-5 * 86400 # conductivity m/d

    fi = 2*R*(rhos-rho)/K/rhos


.. figure:: ./figures/PhiPlot.svg


Comments
--------


.. literalinclude:: ../../examples/Desirade.py
    :language: python 
    :linenos:



Petite Terre and Lake Dziani (Mayotte)
======================================

Lake Dziani is a hyper saline 
volcanic lake located north of the island
its level on average is 0. 

As for the simple lens we start by defining the parameters before calling :func:`IslandLens`. 
First we define the contour names and file locations 

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

Then we create a lake list as :func:`IslandLens` can handle multiple lakes 

.. code:: python

    lakes = [[lake, lake_fname, 0]] # list of lakes with file name and elevation of lake (masl)

We proceed with the parameters need to solve the Poisson equation hence 

* contour sub_sampling, 
* Recharge and Conductvity 
* Poisson coefficient 

.. code:: python 

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
    K = 9e-5*86400 # m/d


    fi = 2 * R * 25 / K / 1025

With all this set we can call :func:`IslandLens`

.. code:: python

    u, Th, X, Y, Zm, dx, dy, itp = il.IslandLens( islands = islands, fname = island_fname, lakes = lakes, ttype = ttype,  fi = fi , sub_sampling = sub_sampling, clockwise = clockwise, plot = True)


The form of the lens is given  in figure :ref:`fig-mayotte`. 
For a discussion on the shape of the water table and comparison with existing measurements see :cite:t:`metivier2024bilan` (spoiler: it seems to work :))

.. _fig-mayotte:

.. figure:: ./figures/Petite_Terre.svg

    Modeled water table of Petite Terre island in Mayotte



The full script is given below

.. literalinclude:: ../../examples/IslandWithLake.py
    :language: python 
    :linenos:
