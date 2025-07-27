Installation
************

To work with IslandsLib you need first to install FreeFem++ *then* pyFreeFem and *eventually* IslandsLib

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

to install pyFreeFem via pip, you must download the pyFreeFem library from the follwing fork on Github

.. code:: bash

    git clone https://github.com/fmetivier/pyFreeFem.git


You must then open a terminal window in the main directory and enter

.. code:: bash

    pip install .

Depending on your system and the way you work you may need administrator priviledges


Install IslandsLib
==================

Once FreeFem++ and pyFreeFem are installed you can download IslandsLib from github 
.. code:: bash

    git clone https://github.com/fmetivier/IslandsLib.git


open a terminal in the main directory and type 

.. code:: bash

    pip install .
  
Again depending on your environment you may need admin priviledges.


