==========================================================
IslandsLib: 2D Modeling of Small Islands Freshwater Lenses
==========================================================


IslandsLib is a library that provides  a set of functions and classes to help 

* solve the Poisson equation to model simple freshwater lenses  islands using the pyfreefem library and FreeFem++ software; 
* process IGN (French National Geography Institute) products, for example to create single coastline contours from shapefile for the purpose of modeling.

For now the aquifer has a constant average conductivity. 
The source term (Precipitation or withdrawal) can be constant or vary in space

IslandsLib also provides a data repository of Island contours. **Contributions welcome** !

For more information (installation procedure, sample scripts, motivations) please see the documentation on readthedocs : https://islandslib.readthedocs.io/en/latest/index.html
