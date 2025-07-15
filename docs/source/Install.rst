############
Installation
############


Install FreeFem++
=================
FreeFem++ is a free and open source finite element PDE solver.
FreeFem++ can be downloaded from https://freefem.org/.
Running versions exist for all platforms. 

If you're using linux just use the package repository. Open a terminal and type

.. code:: bash

    sudo apt-get install freefem++

This will install FreeFem properly.



Install pyFreeFem
=================

At present there is no automated install (via pip for example). You must download the pyFreeFem from github as a zip
of by cloning 

.. code:: bash

    git clone https://github.com/odevauchelle/pyFreeFem.git




to use pyFreeFem you must then add its path in your python file:

.. code:: python

    import sys
    sys.path.append('/Your/Path/To/pyFreeFem-master/')
    
    import pyFreeFem as pyff

Upon installation IslandsLib you will have to include the path manually. 


Install IslandsLib
==================

* Download the library from github
* *régler le pb du chemin vers pyfreefem*
* Open a terminal (in admin mode) in the library directory and run
  
.. code:: bash

    pip install .
  
For linux users :code:`sudo` recommended !


