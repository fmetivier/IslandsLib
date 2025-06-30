################
About IslandsLib
################


Modeling an island's freshwater lens 
====================================

Under certain conditions, the water table of an island can be modeled using the following form of the *Poisson* equation,
named after the French Mathematician Simeon Denis Poisson (1781-1840):

.. math::
    \Delta z_d^2 = \frac{2R(\rho_s-\rho_d)}{K\rho_s}.

where :math:`\Delta z_d` is the Laplacian of the water table elevation :math:`z_d` above sea level, :math:`R` is the recharge (the water that infiltrates), 
:math:`K` is the average hydraulic conductivity, and :math:`\rho_s,\rho_d` are the densities of seawater and freshwater respectively.


The resulting stationnary water table corresponds to an average level. This model assumes that

#. the lens is fully developped, hence there is salwater everywhere beneath the freshwater;
#. the vertical component of velocity in the lens is neglected (Dupuit-Boussinesq approximation);
#. the flow velocity in the salwater is negligible and pressure balance at the saltwater-freshwater interface is hydrostatic;
#. the interface between salt and freshwater is thin.
  
Under theses assumptions the depth of the Freshwater-saltwater interface :math:`z_s` can be deduced from the water table by

.. math::
    z_s = \left(\frac{\rho_d}{\rho_s-\rho_d}\right)z_d

For a complete derivation and discussion see for example Metivier et al. (2024) https://hal.science/hal-04632890v1)


What Islandslib is
==================

IslandsLib provides a simple function `IslandLens()` to model freshwater lenses without a prior knowledge of Finite Element Modeling. 
Given a set of data and constraints if will model the freshwater lens of a *small* island. 
A *small* island is here defined as an island without a river network. Lakes can be included (see examples).

IslandsLib also provides a set of functions and classes to help 

* model the poisson equation to model simple freshwater lenses  islands using the pyfreefem library and FreeFem++ software; 
* Create single coastline contours from shapefile for this purpose.


River networks will be included in  later versions of `IslendLens`. If you want to solve problems with river networks 
you need to use the set of functions provided or to use pyFreeFem directly.

IslandsLib also provides a data repository of Island contours. Please feel free to contribute !

What IslandsLib is not
======================


IslandsLib is not a Finite Element Model Solver
----------------------------------------------- 
It relies on the pyFreeFem library developped by Olivier Devauchelle (https://github.com/odevauchelle/pyFreeFem) a python wrapper 
around the FreeFem++ Solver (https://freefem.org/). 
As such **both** pyFreeFem and FreeFem++ must be installed on your computer in order to use IslandsLib

IslandsLib is not a plug and play IGN converter
-----------------------------------------------
functions are provided to help transform contours, especially those disclosed by the French Geographic Institute (IGN). The functions and procedures 
can be used with  basically any set of contours but theyr require a bit of work. 
An example of single contour création from in IGN shapefile is given


Installation
============

* Download the library from github
* Open a terminal in the library directory
* Run python setup.py install

**Dependancies**

* FreeFem++
* pyFreeFem