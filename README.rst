====================================================
IslandsLib: 2D modeling of small island water tables
====================================================

Introduction
============

IslandsLib provides a simple function `IslandLens()` to model freshwater lenses under the Dupuits-Boussinesq approximation 
without a prior knowledge of Finite Element Modeling. 
Given a set of island contour data and hydrogeological constraints if will model the freshwater lens of a *small* island. 
A *small* island is here defined as an island without a river network. Lakes can be included (see examples).

IslandsLib also provides a set of functions and classes to help 

* solve the poisson equation to model simple freshwater lenses  islands using the pyfreefem library and FreeFem++ software; 
* Create single coastline contours from shapefile for this purpose.

River networks will be included in  later versions of `IslendLens`. If you want to solve problems with river networks 
you need to use the set of functions provided or to use pyFreeFem directly (the documentation provides an example).

IslandsLib also provides a data repository of Island contours. Please feel free to contribute !


Installation
============

Requirements
------------

To use IslandsLib you need to Install

* FreeFem++ an opensource finite element partial differential equation solver (https://freefem.org/)
* pyFreeFem a python wrapper for FreeFem (https://github.com/odevauchelle/pyFreeFem)

classical libraries such as 

* numpy
* scipy
* matplotlib

are also needed.

Installation
------------

download the IslandsLib library from the github site, open a terminal in the main repository and enter

.. code:: bash 

    python setup.py install

This will install the python Dependancies but not FreeFem++ nor pyFreeFem.

Documentation
=============

Examples and code detail can be found in the docs/build/html and docs/build/pdf folder or on readthedocs