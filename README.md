# ACOSA
A compilation of spherical algorithms.

This code collection contains C++ classes to calculate the Voronoi-tesselation
and Delaunay-triangulation (**VDTesselation**), and the convex hull (**ConvexHull**)
of a set of nodes on a sphere. The algorithms used to calculate the tesselation
and the convex hull are O(N\*log(N)).

The geometricgraph.hpp header contains an algorithm for the determination of links of a geometric
graph, a spatially embedded network with links between all pairs of nodes closer than
a threshold. The algorithm's complexity is O(N\*log(N) + M\*N + M^2) where N is the
number of nodes and M the mean degree of the resulting graph.

Additionally, it can be used in Python using the Cython module. To install it,
change to the ACOSA directory and type (as root):

    pip install .

Note that there are still things left to do, e.g. some unhandled degeneracies that
currently (hopefully all) purposefully crash the code. So be sure to check the results.
For randomly distributed node sets, however, things should run smoothly.
