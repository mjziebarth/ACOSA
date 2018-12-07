# -*- coding: utf-8 -*-
#
# Python interface for the ACOSA C++ classes
# Copyright (C) 2016 Malte Ziebarth
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Interface to python world:
import numpy as np

# Cython imports:
import cython
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference, preincrement
np.import_array()

####################################################################################
#                                    C++ types                                     #
####################################################################################

# Basic types:
cdef extern from "basic_types.hpp" namespace "ACOSA":
	cdef struct Triangle:
		size_t i
		size_t j
		size_t k
	
	cdef struct Node:
		double lon
		double lat
	
	cdef struct Link:
		size_t i
		size_t j

# The Voronoi-/Delaunay-tesselation class:
cdef extern from "vdtesselation.hpp" namespace "ACOSA":
	cdef cppclass VDTesselation :
		const size_t N
		
		size_t size() const
		
		VDTesselation(const vector[Node]& nodes, double tolerance) except +
		
		void delaunay_triangulation(vector[Link]& links) const
		
		const vector[Triangle]& delaunay_triangles() const
		
		void voronoi_tesselation(vector[Node]& voronoi_nodes,
		                         vector[Link]& voronoi_links) const
		
		void voronoi_cell_areas(vector[double]& areas) const
		
		void associated_nodes(const vector[size_t]& voronoi_nodes,
			vector[size_t]& associated) const;

# The ConvexHull class:
cdef extern from "convexhull.hpp" namespace "ACOSA":
	cdef cppclass CppConvexHull "ACOSA::ConvexHull":
		CppConvexHull(const vector[Node]& nodes, const Node& inside, double tolerance) except +
		
		vector[size_t].const_iterator begin() const
		
		vector[size_t].const_iterator end() const
		
		size_t size() const
		
		bool empty() const
		
		bool is_contained(const Node& node) const
		
		void distance_to_border(const vector[Node]& nodes,\
                                vector[double]& distances) const

# The AlphaSpectrum and AlphaShape classes:
cdef extern from "alphaspectrum.hpp" namespace "ACOSA":
	cdef cppclass AlphaShape :
		AlphaShape()
		
		AlphaShape(const vector[size_t]& nodes,
		           const vector[Link]& links)
		
		const vector[size_t]& nodes() const
		
		const vector[Link]& links() const
		
		vector[Link] source_network_links() const
	
	
	cdef cppclass AlphaSpectrum :
		AlphaSpectrum(const vector[Node]& nodes,
		              const VDTesselation& tesselation)
		
		AlphaShape operator()(double alpha) const



####################################################################################
#                              The Python classes                                  #
####################################################################################
cdef class VoronoiDelaunayTesselation:
	"""
	A wrapper class around the VDTesselation c++ class, computing the
	Voronoi tesselation and Delaunay triangulation of a set of points on
	the unit sphere. The algorithm used is a O(N*log(N)) sweepline
	algorithm [1]. After initialization, the class may be queried about
	the tesselations' network structure and the Voronoi cells spatial
	layout.

	Keeping an instance of this class alive allows caching of computed values.

	[1] Xiaoyu Zheng et al.: A Plane Sweep Algorithm for the Voronoi
	    Tesselation of the Sphere, in: electronic-Liquid Crystal
	    Communications, 2011-12-13.
	    http://www.e-lc.org/docs/2011_12_05_14_35_11

	Note: A similar algorithm, which I was not aware of at the time
	of implementation, was derived earlier by Dinis & Marmede (2010)
	(doi:10.1109/ISVD.2010.32).
	"""
	cdef VDTesselation* tesselation

	# Constructor:
	def __cinit__(self, lon, lat, tolerance = 1e-10):
		"""
		Initialize a VoronoiDelaunayTesselation instance, running the
		sweepline algorithm and creating data structures for both the
		Voronoi tesselation and the Delanay triangulation.

		Required arguments:
		   lon : Array of longitude coordinates in degrees.
		   lat : Array of latitude coordinates in degrees.

		Optional keyword argument:
		   tolerance : The absolute tolerance to use for various
		               floating point comparisons. Must be
		               convertible to float. Experimenting with
		               this parameter may sometimes solve problems
		               in numerical instable configurations.
		               (Default: 1e-10)

		The input coordinates need to be convertible to numpy arrays
		and are flattened before executing the algorithm. The flattened
		arrays are required to be of same length.
		Will throw a runtime error if VDTesselation failed.
		"""
		# Input safety:
		cdef np.ndarray[double, ndim=1, cast=True] lon_fix \
		   = np.atleast_1d(np.deg2rad(lon).astype(float).flatten())
		cdef np.ndarray[double, ndim=1, cast=True] lat_fix \
		   = np.atleast_1d(np.deg2rad(lat).astype(float).flatten())
		cdef double _tolerance = float(tolerance)

		# Sanity checks:
		cdef size_t N
		N = len(lon_fix)
	
		if (len(lat_fix) != N):
			raise Exception("VoronoiDelaunayTesselation() :\nLength of longitude "
				"and latitude arrays not equal!")

		# Copy numpy arrays to c++ vectors:
		cdef size_t i
		cdef vector[Node] nodes
		nodes.resize(N)
		for i in range(N):
			nodes[i].lon = lon_fix[i]
			nodes[i].lat = lat_fix[i]

		# Create VDTesselation object:
		# TODO put this into a smart pointer.
		self.tesselation = new VDTesselation(nodes,_tolerance)

		if not self.tesselation:
			raise Exception("VoronoiDelaunayTesselation() :\nCould not allocate "
				"VDTesselation object!")

	# Destructor:
	def __dealloc__(self):
		if self.tesselation:
			del self.tesselation

	# Voronoi cell areas:
	def voronoi_cell_areas(self):
		# Sanity check:
		if not self.tesselation:
			raise Exception("VoronoiDelaunayTesselation() :\nVDTesselation "
				"was not initialized!\n")

		# Obtain vector of areas:
		cdef vector[double] areas
		dereference(self.tesselation).voronoi_cell_areas(areas)

		# Copy to numpy array:
		cdef np.ndarray[float, ndim=1] np_areas = np.zeros(areas.size(), 
		                                                   dtype=np.float32)

		cdef size_t i
		for i in range(areas.size()):
			np_areas[i] = areas[i]

		# Return numpy array:
		return np_areas


	def voronoi_tesselation(self):
		"""
		Obtain the Voronoi tesselation.

		Returns (lon, lat, links):
		   lon   : One-dimensional numpy array of longitude
		           coordinates of the Voronoi tesselation vertices.
		   lat   : One-dimensional numpy array of latitude
		           coordinates of the Voronoi tesselation vertices.
		   links : Two-dimensional Mx2 numpy array representing the
		           Voronoi tesselation link list.

		Both lon and lat are length N, where N is the number of
		vertices of the tesselation. The pair (i,j)=links[k,:]
		indicates a link between the Voronoi tesselation vertices
		(lon[i],lat[i]) and (lon[j],lat[j]).
		"""
		# Sanity check:
		if not self.tesselation:
			raise Exception("VoronoiDelaunayTesselation() :\nVDTesselation "
				"was not initialized!\n")

		# Obtain vectors of nodes and links:
		cdef vector[Node] nodes
		cdef vector[Link] links
		dereference(self.tesselation).voronoi_tesselation(nodes, links)

		# Copy to numpy arrays:
		cdef np.ndarray[long, ndim=2] np_links = np.zeros((links.size(),2), 
		                                                  dtype=long)
		cdef np.ndarray[float, ndim=1] lon = np.zeros(nodes.size(), 
		                                              dtype=np.float32)
		cdef np.ndarray[float, ndim=1] lat = np.zeros(nodes.size(), 
		                                              dtype=np.float32)

		cdef size_t i
		cdef double r2d = 180.0/np.pi
		for i in range(nodes.size()):
			lon[i] = r2d*nodes[i].lon
			lat[i] = r2d*nodes[i].lat

		for i in range(links.size()):
			np_links[i,0] = links[i].i
			np_links[i,1] = links[i].j

		# Return numpy arrays:
		return lon, lat, np_links


	def delaunay_triangulation(self):
		"""
		Obtain the Delaunay triangulation as a link list.

		Returns:
		   links : Two-dimensional Mx2 numpy array representing
		           the Delauany triangulation link list.

		The pair (i,j)=links[k,:] indicates a link between the
		nodes i and j of the generator node set. Both indices
		refer to the position in the flattened generating
		longitude/latitude arrays.
		"""
		# Sanity check:
		if not self.tesselation:
			raise Exception("VoronoiDelaunayTesselation() :\nVDTesselation "
				"was not initialized!\n")

		# Obtain vector of links:
		cdef vector[Link] links
		dereference(self.tesselation).delaunay_triangulation(links)

		# Copy to numpy array:
		cdef np.ndarray[long, ndim=2] np_links = np.zeros((links.size(),2), 
		                                                  dtype=long)

		cdef size_t i
		for i in range(links.size()):
			np_links[i,0] = links[i].i
			np_links[i,1] = links[i].j

		# Return numpy arrays:
		return np_links


	def delaunay_triangles(self):
		"""
		Obtain the Delaunay triangulation as a set of triangles.

		Returns:
		   triangles : A Tx3 numpy array representing the triangles
		               of the Delauany triangulation.

		The triple (i,j,k)=triangles[l,:] indicates a triangle
		consisting of the nodes i, j, and k of the generator
		node set. Node numbering refers to flattened indices
		of the generating longitude/latitude arrays.
		"""
		# Sanity check:
		if not self.tesselation:
			raise Exception("VoronoiDelaunayTesselation() :\nVDTesselation "
				"was not initialized!\n")

		# Copy triangles to numpy array:
		cdef size_t N = dereference(self.tesselation).delaunay_triangles().size()
		cdef np.ndarray[long, ndim=2] triangles = np.zeros((N,3),dtype=long)

		for i in range(N):
			triangles[i,0] = dereference(self.tesselation).delaunay_triangles()[i].i
			triangles[i,1] = dereference(self.tesselation).delaunay_triangles()[i].j
			triangles[i,2] = dereference(self.tesselation).delaunay_triangles()[i].k

		return triangles


	def associated_nodes(self, np.ndarray[long, ndim=1] voronoi_nodes):
		"""
		Obtain the combined set of generator nodes associated
		to a set of nodes of the Voronoi tesselation, i.e.
		whose Voronoi cells have at least one of the given
		Voronoi node set as their vertices.

		Required arguments:
		   voronoi_nodes : The set of nodes of the Voronoi
		                   tesselation for each of which the
		                   neighbourhood shall be calculated.
		                   Needs to be convertible to a flat
		                   numpy array.

		Returns:
		   A flat numpy array of node indices of generator
		   nodes that are associated in the described sense
		   to the given Voronoi tesselation vertices.
		   The node indices refer to the linear index in
		   the flat generating longitude/latitude arrays.
		"""
		# Sanity checks:
		if not self.tesselation:
			raise Exception("VoronoiDelaunayTesselation.associated_nodes():\n"
				"VDTesselation was not initialized!\n")

		cdef size_t VN = dereference(self.tesselation).size()
		if np.any(voronoi_nodes >= VN) or np.any(voronoi_nodes < 0):
			raise ValueError("VoronoiDelaunayTesselation.associated_nodes() :\n"
				"At least one Voronoi node index out of bounds!")

		# Copy Voronoi index to c++ container:
		cdef vector[size_t] vnodes
		vnodes.resize(len(voronoi_nodes))
		cdef size_t i
		for i in range(len(voronoi_nodes)):
			vnodes[i] = voronoi_nodes[i]

		# Obtain vector of associated nodes:
		cdef vector[size_t] associated
		dereference(self.tesselation).associated_nodes(vnodes, associated)
		vnodes.clear()

		# Copy to numpy output:
		cdef np.ndarray[long, ndim=1] out = np.zeros(associated.size(),
		                                             dtype=long)
		for i in range(associated.size()):
			out[i] = associated[i]

		return out



################################################################################
cdef class ConvexHull:
	"""
	ConvexHull(lon, lat, lon_inside, lat_inside, tolerance=1e-12)

	A wrapper class around the ConvexHull C++ class, computing the
	on-sphere convex hull of a set of points on the unit sphere in
	O(N*log(N)) time using Graham's scan.

	This implementation may show numerical instabilities if points
	lie too dense. If a non-valid convex hull was detectet, a
	RuntimeError is thrown. In that case, using the AlphaShape for
	alpha=-1/pi may be worthwhile to test.

	Required arguments:
	   lon        : Set of longitude coordinates of the point set.
	   lat        : Set of latitude coordinates of the point set.
	                Both lon and lat need to be flattenable, of
	                equal size, and their flat indexing needs to
	                be in same order.
	   lon_inside : Longitude of a point inside the convex hull,
	                which is used as the Graham's scan anchor.
	   lat_inside : Latitude of that point.

	Optional keyword arguments:
	   tolerance  : A numerical tolerance parameter used in the
	                geometrical comparisons. Tuning this
	                parameter may be worthwhile if the algorithm
	                fails.
	                (Default: 1e-12)
	"""
	cdef CppConvexHull* hull


	def __cinit__(self, lon, lat, lon_inside, lat_inside, tolerance=1e-12):
		"""
		Constructor.
		"""
		# Input safety:
		cdef np.ndarray[double, ndim=1, cast=True] lon_fix \
		   = np.atleast_1d(np.deg2rad(lon).astype(float).flatten())
		cdef np.ndarray[double, ndim=1, cast=True] lat_fix \
		   = np.atleast_1d(np.deg2rad(lat).astype(float).flatten())
		cdef double _tolerance = float(tolerance)

		# Sanity check of input array dimensions:
		cdef size_t N
		N = len(lon)

		if (len(lat) != N):
			raise Exception("ConvexHull() :\nLength of longitude "
				"and latitude arrays not equal!")

		# Copy numpy arrays to c++ vectors:
		cdef size_t i
		cdef double d2r = np.pi/180.0
		cdef double lon_i = lon_inside
		cdef double lat_i = lat_inside
		cdef Node inside = Node(d2r*lon_i, d2r*lat_i)

		cdef vector[Node] nodes
		nodes.resize(N)
		for i in range(N):
			nodes[i].lon = lon_fix[i]
			nodes[i].lat = lat_fix[i]

		# Create VDTesselation object:
		self.hull= new CppConvexHull(nodes, inside, tolerance=_tolerance)

		if not self.hull:
			raise Exception("ConvexHull() :\nCould not allocate "
				"C++ ConvexHull object!")


	def __dealloc__(self):
		"""
		The destructor.
		"""
		if self.hull:
			del self.hull


	def nodes(self):
		"""
		Return the set of nodes that form the convex hull.

		Returns:
		   A one-dimensional numpy array containing the indices
		   of the convex hull nodes. The indices refer to the
		   flattened generator set.
		"""
		if not self.hull:
			return None

		cdef size_t N = dereference(self.hull).size()
		cdef vector[size_t].const_iterator it = dereference(self.hull).begin()
		cdef size_t i

		cdef np.ndarray[ndim=1, dtype=size_t] nodes = np.zeros(N, dtype='uint64')

		for i in range(N):
			nodes[i] = dereference(it)
			preincrement(it)

		return nodes


	def contains(self, lon, lat):
		"""
		Query whether the convex hull contains a set of
		longitude and latitude coordinates.

		Required arguments:
		   lon : Set of longitudes in degrees.
		   lat : Set of latitudes in degrees.

		The coordinates need to be convertible to one-dimensional
		numpy ndarrays and the size of longitude and latitude
		arrays need to be equal.

		Returns:
		   An array of booleans, each indicating whether the
		   corresponding coordinate pair is contained in the
		   convex hull.
		"""
		# Input safety:
		cdef np.ndarray[double, ndim=1, cast=True] lon_fix \
		   = np.atleast_1d(np.deg2rad(lon).astype(float).flatten())
		cdef np.ndarray[double, ndim=1, cast=True] lat_fix \
		   = np.atleast_1d(np.deg2rad(lat).astype(float).flatten())

		# Handle shaped input arrays:
		if isinstance(lon,np.ndarray):
			shape = lon.shape
			assert isinstance(lat,np.ndarray)
			assert np.array_equal(shape,lat.shape)
		else:
			shape = None

		# Sanity checks:
		cdef size_t N
		N = len(lon_fix)

		if (len(lat_fix) != N):
			raise Exception("ConvexHull.contains() :\nLength of longitude "
				"and latitude arrays not equal!")

		if not self.hull:
			raise Exception("ConvexHull.contains() :\nNo hull object was "
				"allocated!")

		cdef np.ndarray[np.uint8_t, ndim=1, cast=True] contained \
			= np.zeros(N, dtype=np.bool_)

		cdef size_t i
		cdef Node node
		for i in range(N):
			# When writing ".is_contained(Node(lon[i],lat[i]))", Cython creates
			# a CONSTANT global struct that it tries to assign the lon and lat
			# values and the pass to .is_contained. Passing to the constant struct
			# then obviously does not work...
			# Workaround:
			node.lon = lon_fix[i]
			node.lat = lat_fix[i]
			contained[i] = dereference(self.hull).is_contained(node)

		# Handle shaped input arrays:
		if shape is None:
			return contained
		else:
			return contained.reshape(shape)


	def distance_to_border(self, lon, lat):
		"""
		Calculates the distance of nodes inside this convex
		hull to the hull's borders.

		Required arguments:
		   lon, lat : Numpy arrays of longitude and latitude
		              coordinates in degrees.

		Raises an exception if any of the given coordinates
		are outside the hull.
		"""
		# Convert latlon to known types:
		cdef np.ndarray[double, ndim=1, cast=True] lon_fix \
		   = np.atleast_1d(np.deg2rad(lon).astype(float).flatten())
		cdef np.ndarray[double, ndim=1, cast=True] lat_fix \
		   = np.atleast_1d(np.deg2rad(lat).astype(float).flatten())

		# Sanity checks:
		cdef size_t N
		N = len(lon_fix)

		if (len(lat_fix) != N):
			raise Exception("ConvexHull.distance_to_border() :\nLength of "
			                "longitude and latitude arrays not equal!")

		if not self.hull:
			raise Exception("ConvexHull.distance_to_border() :\nNo hull "
			                "object was allocated!")

		# Create node vector:
		cdef vector[Node] nodes
		cdef size_t i
		nodes.resize(N)
		for i in range(N):
			nodes[i].lon = lon_fix[i]
			nodes[i].lat = lat_fix[i]

		# Calculate distance:
		cdef vector[double] dist_vec
		try:
			dereference(self.hull).distance_to_border(nodes, dist_vec)
		except:
			raise Exception("ConvexHull empty!")

		# Copy to numpy array:
		cdef np.ndarray[double, ndim=1] distance = np.zeros(N, dtype=float)
		for i in range(N):
			distance[i] = dist_vec[i]

		return distance



################################################################################
cdef class PyAlphaSpectrum:
	"""
	A cython interface to the ACOSA C++ implementation of the AlphaSpectrum
	class that implements the alpha sphere algorithm in [1] for negative
	alpha on a spherical topology.
	
	It is initialized passing numpy arrays of longitude and latitude coordinates
	of the point set to calculate the alpha shape of.
	
	Afterwards, a linear-time algorithm can be used calling the PyAlphaSpectrum
	instance on a specific value of alpha to obtain the alpha shape.
	
	Implemented methods:
	
	__cinit__
	__dealloc__
	__call__
	
	
	Example:
	
	> lon = 120.0 * np.random.random(200)
	> lat = 90.0 * np.random.random(200)
	> alpha = -12.0
	> spectrum = AlphaSpectrum(lon, lat)
	> shape = spectrum(alpha)
	
	
	Bibliography:
	
	[1] Herbert Edelsbrunner et al.: On the Shape of a Set of Points in the
	    Plane, in: IEEE Transactions on Information Theory, Vol. 29, No. 4,
	    July 1983
	"""
	cdef AlphaSpectrum* spectrum

	# Constructor:
	def __cinit__(self, lon, lat, vdtesselation=None):
		"""
		Initialize an AlphaSpectrum.

		Required arguments:
		   lon : Longitude coordinates in degrees.
		   lat : Latitude coordinates in degrees.

		Both coordinate arrays need to be convertible to
		one-dimensional numpy arrays. These one-dimensional
		arrays define the indexing that is used later on.

		Optional keyword arguments:
		   vdtesselation : An instance of a VoronoiDelauanyTesselation of
		                     the point set if already calculated. Optional: If
		                     omitted, the tesselations will be calculated
		                     internally.
		"""
		# Convert latlon to known types:
		cdef np.ndarray[double, ndim=1, cast=True] lon_fix \
		   = np.atleast_1d(np.deg2rad(lon).astype(float).flatten())
		cdef np.ndarray[double, ndim=1, cast=True] lat_fix \
		   = np.atleast_1d(np.deg2rad(lat).astype(float).flatten())
		assert vdtesselation is None or \
		   isinstance(vdtesselation,VoronoiDelaunayTesselation)

		# Sanity checks:
		cdef size_t N
		N = len(lon)

		if (len(lat) != N):
			raise Exception("PyAlphaSpectrum() :\nLength of "
			                "longitude and latitude arrays not equal!")

		# Create node vector:
		cdef vector[Node] nodes
		cdef size_t i
		nodes.resize(N)
		for i in range(N):
			nodes[i].lon = lon_fix[i]
			nodes[i].lat = lat_fix[i]

		# If no Voronoi-Delaunay-Tesselation is given, compute it:
		cdef VoronoiDelaunayTesselation vdt
		if vdtesselation is None:
			vdt = VoronoiDelaunayTesselation(lon, lat)
		else:
			vdt = vdtesselation

		# Calculate Spectrum:
		cdef VDTesselation* tess = vdt.tesselation
		if not tess:
			raise Exception("PyAlphaSpectrum() :\nTesselation not initialized!")

		self.spectrum = new AlphaSpectrum(nodes,dereference(tess))

		if not self.spectrum:
			raise Exception("PyAlphaSpectrum() :\nCould not allocate "
			                "AlphaSpectrum object!")


	# Destructor:
	def __dealloc__(self):
		if self.spectrum:
			del self.spectrum
	
	
	# Create alpha shape:
	def __call__(self, float alpha):
		"""
		Compute the alpha shape for a specific alpha.
		
		:type alpha: float
		:arg  alpha: Alpha value to calculate the shape for.
		
		:return: Tuple (nodes,links) of node indices and links between them that
		         are alpha-extreme (see [1]).
		         The indices in nodes refer to the original node set, the
		         indices in links to nodes.
		"""
		# Sanity check:
		if alpha > -1.0/np.pi:
			raise Exception("PyAlphaSpectrum.__call__() : Only implemented "
			                "for alpha <= -1.0/pi !")
		
		if not self.spectrum:
			raise Exception("PyAlphaSpectrum.__call__() : Spectrum was not "
			                "created!")
		
		# Create AlphaShape object:
		cdef AlphaShape shape = dereference(self.spectrum)(alpha)
		
		# Obtain vectors from alpha shape:
		cdef np.ndarray[ndim=1,dtype=long] nodes \
		    = np.zeros(shape.nodes().size(),dtype=int)
		cdef np.ndarray[ndim=2,dtype=long] links \
		    = np.zeros((shape.links().size(),2),dtype=int)
		
		cdef size_t i
		for i in range(shape.nodes().size()):
			nodes[i] = shape.nodes()[i]
		
		for i in range(shape.links().size()):
			links[i,0] = shape.links()[i].i
			links[i,1] = shape.links()[i].j
		
		# Return numpy arrays:
		return nodes, links
