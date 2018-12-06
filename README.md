# ACOSA
A Collection Of Spherical Algorithms.

This code collection contains C++ classes and a Python interface to calculate
- the Voronoi-tesselation and Delaunay-triangulation (**VDTesselation**),
- the convex hull (**ConvexHull**),
- and the negative alpha spectrum and alpha shapes (**AlphaSpectrum**)

of a set of nodes on a sphere. The algorithms used to calculate the tesselation and the
convex hull are O(N\*log(N)), alpha spectra and shapes work on the such
generated data and take additional O(N) time.

The algorithms should run without errors between N=10^2 and 10^6. Note that while some
tests have been done without errors, testing of the latest revision most likely has
not been extensive. Make sure to do some checks agains your results if you plan to use
the code. Most notably in two cases:
- The code has not been tested for very small node sets (N\<20)
- For very big node sets (N >> 10^6), double precision accuracy may become significant.
  Chances are good, though, that such errors are noted for the triangulation (aka 
  VDTesselation), as it empirically reacts chaotically to significant rounding errors.
  Adjusting the tolerance parameter may help in such cases.

Please message me if you encounter any errors.

## Python
The code is available under the module **acosa** (not yet on PyPI). To install it, change to
the ACOSA directory and install via
```bash
pip install .
```


## C++
The geometricgraph.hpp header contains an algorithm for the determination of links of a geometric
graph, a spatially embedded network with links between all pairs of nodes closer than
a threshold. The algorithm's complexity is O(N\*log(N) + M\*N + M^2) where N is the
number of nodes and M the mean degree of the resulting graph.

The Fortune's algorithm used to determine the Delaunay triangulation is implemented
from:  
[Xiaoyu Zheng et al.: A Plane Sweep Algorithm for the Voronoi Tesselation of the
 Sphere,  
 in: electronic-Liquid Crystal Communications, 2011-12-13]
(http://www.e-lc.org/docs/2011_12_05_14_35_11)

The alpha sphere and shape algorithms are implemented from:  
[Herbert Edelsbrunner et al.: On the Shape of a Set of Points in the
 Plane,  
 in: IEEE Transactions on Information Theory, Vol. 29, No. 4, July 1983]
