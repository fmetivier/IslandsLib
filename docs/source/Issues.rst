Known issues and Todo list
**************************

Issues
------

* If nodes on the borders are two close to each other the mesh generation will fail.
* Contours must **all** be either clockwise or anticlockwise not a mix
* When drawing standard plots the streamlines can be a bit *moche*  :)
* Install
   * Windows 11, dell optiplex. Installation runs ok but pyFreeFem crashes because of :code:`tempfile` library issues
   * Mac OS catalina 10.15 on imac with i5 processors. :code:`FreeFem++` and :code:`pyFreeFem` ok. :code:`IslandsLib` installation failed because pyproj library could not be installed easily.
   

Todo
----

* Allow for different the choice of types of graphic output 
* Include rivers 
  