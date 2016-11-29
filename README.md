# ACOSA
A compilation of spherical algorithms.

This code collection contains C++ classes to calculate the Voronoi-tesselation
and Delaunay-triangulation (**VDTesselation**), and the convex hull (**ConvexHull**)
of a set of nodes on a sphere. The algorithms used to calculate the tesselation
and the convex hull are O(N*log(N)).

Additionally, it can be used in Python using the Cython module. To install it,
change to the ACOSA directory and type (as root):

    pip install .

Note that there are still things left to do, e.g. some unhandled degeneracies that
currently (hopefully all) purposefully crash the code. So be sure to check the results.
For randomly distributed node sets, however, things should run smoothly.
