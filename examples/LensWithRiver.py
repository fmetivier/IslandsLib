
# Libraries
import sys 
sys.path.append("./..")

import matplotlib.pyplot as plt

import IslandsLib as il

def interp_bc(river,boundaries):


    for i in range(len(river.x)-1):
        s = shapely.LineString([(river.x[i],river.y[i],(rivier.x[i+1],river.y[x+1]))])

def river_at_sl():
    """Poisson test in a simple square (!!) island with a river
    here the river is at see level as the coast
    """
    multi = True
    b = []

    if multi:
        b.append(il.river("sud", x=[0, 9], y=[0, 0]))
        b.append(il.river("est", x=[10, 10], y=[0, 9]))
        b.append(il.river("nord", x=[10, 1], y=[10, 10]))
        b.append(il.river("ouest", x=[0, 0], y=[10, 1]))
    else:
        b.append(il.river("total", x=[0, 10, 10, 0], y=[0, 0, 10, 10]))

    r = []
    r.append(il.river("milieu", x=[5, 4, 3, 7], y=[0, 2, 3, 7]))

    borders, vertices, segments = il.prepare_boundaries(
        b, r, for_FF=True)

    ##############################################################
    #
    # effectue la triangulation
    # 1) avec triangle
    # 2) avec Freefem
    #
    ##############################################################

    ttype = 'pq33a0.1'
    # ttype = 'p'
    xv, yv, mesh, borders = il.create_triangle_mesh(
        borders = borders, vertices = vertices, segments = segments, plot = True, ttype = ttype )

    print("Creating FreeFem mesh...")
    Th, Th_boundaries = il.create_FreeFem_mesh(
        xv, yv, mesh, adapt = False, borders = borders, plot = True)

    # for edge in Th.get_boundary_edges():
    #     print(edge)

    bcs = {}
    for T in Th_boundaries:
        if T.rname in ['milieu','sud','est','nord','ouest']:
            bcs[T.rname] = 0 * T.x

    for key, val in bcs.items():
        print(key, val)

    u, Th = il.solve_fem( Th = Th, bcs = bcs)

    fig, ax = plt.subplots(1)
    il.plot_u(Th, u, ax)


if __name__ == "__main__":

    river_at_sl()
    plt.savefig("river.svg", bbox_inches='tight')
    plt.show()    