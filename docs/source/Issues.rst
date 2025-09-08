Known issues and Todo list
**************************

Issues
------

* If nodes on the borders are two close to each other the mesh generation will fail.
* Contours must **all** be either clockwise or anticlockwise not a mix
* When drawing standard plots the streamlines can be a bit *moche*  :)
* Install
   * Mac OS catalina 10.15 on imac with i5 processors. FreeFemm++ and pyFreeFem ok. IslandsLib failed because gdal library could not be installed easily.
   * Windows 11, dell optiplex. Installation ok but pyFreeFem crashes because of :code:`tempfile` library issues

Todo
----

* Allow for different the choice of types of graphic output 
* Include rivers 
  