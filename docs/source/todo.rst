Add variable source term 
========================

fi: source term, at present constant

* create polygons with fi values or give a function with fi(x,y)
* check whether a node Th.x, Th.y is inside the polygon or interpolate the function fi to nodes ? seems much more simple to me !
* apply fi*G (simple)

En fait donner une série de points x,y,z et 
* soit créer un interpolateur (plus proche voisin ou linéaire) sur le polygone entourant l'ile (précipitations) et interpoler la source en chaque point
* soit affecter les points aux noeuds les plus proches (pompages) et ajouter une valeur de base à tous les points (précipitations)
    