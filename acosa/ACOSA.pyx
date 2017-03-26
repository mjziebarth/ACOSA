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
	cdef cppclass ConvexHull :
		ConvexHull(const vector[Node]& nodes, const Node& inside)
		
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
	This class is a wrapper class around the VDTesselation c++ class.
	Keeping an instance of this class alive allows caching of computed values.
	"""
	cdef VDTesselation* tesselation
	
	# Constructor:
	def __cinit__(self, np.ndarray[float, ndim=1] lon,
	              np.ndarray[float, ndim=1] lat,
	              double tolerance = 1e-10
	    ):
		# Sanity checks:
		cdef size_t N
		N = len(lon)
	
		if (len(lat) != N):
			raise Exception("VoronoiDelaunayTesselation() :\nLength of longitude "
				"and latitude arrays not equal!")
		
		# Copy numpy arrays to c++ vectors:
		cdef size_t i
		cdef double d2r = np.pi/180.0
	
		cdef vector[Node] nodes
		nodes.resize(N)
		for i in range(N):
			nodes[i].lon = d2r*lon[i]
			nodes[i].lat = d2r*lat[i]
		
		# Create VDTesselation object:
		self.tesselation = new VDTesselation(nodes,tolerance)
		
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
	
	# Voronoi tesselation:
	def voronoi_tesselation(self):
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
	
	# Delaunay triangulation:
	def delaunay_triangulation(self):
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
	
	
	# Delaunay triangles:
	def delaunay_triangles(self):
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
	
	# Nodes of the network associated to a set of Voronoi nodes:
	def associated_nodes(self, np.ndarray[long, ndim=1] voronoi_nodes):
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
cdef class PyConvexHull:
	cdef ConvexHull* hull
	
	# Constructor:
	def __cinit__(self, np.ndarray[float, ndim=1] lon,
	              np.ndarray[float, ndim=1] lat,
	              float lon_inside, float lat_inside
	    ):
		# Sanity checks:
		cdef size_t N
		N = len(lon)
	
		if (len(lat) != N):
			raise Exception("PyConvexHull() :\nLength of longitude "
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
			nodes[i].lon = d2r*lon[i]
			nodes[i].lat = d2r*lat[i]
		
		# Create VDTesselation object:
		self.hull= new ConvexHull(nodes, inside)
		
		if not self.hull:
			raise Exception("PyConvexHull() :\nCould not allocate "
				"ConvexHull object!")
	
	# Destructor:
	def __dealloc__(self):
		if self.hull:
			del self.hull
	
	# Obtain list of hull nodes:
	def nodes(self):
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
	
	# Check if nodes are contained in the convex hull:
	def contains(self, np.ndarray[float, ndim=1] lon,
	             np.ndarray[float, ndim=1] lat
	    ):
		
		# Sanity checks:
		cdef size_t N
		N = len(lon)
	
		if (len(lat) != N):
			raise Exception("PyConvexHull.contains() :\nLength of longitude "
				"and latitude arrays not equal!")
		
		if not self.hull:
			raise Exception("PyConvexHull.contains() :\nNo hull object was "
				"allocated!")
		
		
		# Convert to radians:
		cdef double d2r = np.pi/180.0
		lon = d2r*lon
		lat = d2r*lat
		
		
		cdef np.ndarray[np.uint8_t, ndim=1, cast=True] contained \
			= np.zeros(N, dtype=np.bool_)
		
		cdef size_t i
		for i in range(N):
			contained[i] = dereference(self.hull).is_contained(Node(lon[i],lat[i]))
		
		return contained
	
	# Calculate distances to border of convex hull:
	def distance_to_border(self, np.ndarray[float, ndim=1] lon,
	                       np.ndarray[float, ndim=1] lat
	    ):
	    
		"""
		Calculates the distance of nodes inside this convex hull to the
		hull's borders.
		
		Raises exception if given any of the given coordinates are outside the hull.
		"""
		
		# Sanity checks:
		cdef size_t N
		N = len(lon)
	
		if (len(lat) != N):
			raise Exception("PyConvexHull.distance_to_border() :\nLength of "
			                "longitude and latitude arrays not equal!")
		
		if not self.hull:
			raise Exception("PyConvexHull.distance_to_border() :\nNo hull "
			                "object was allocated!")
		
		# Create node vector:
		cdef vector[Node] nodes
		cdef double d2r = np.pi/180.0
		cdef size_t i
		nodes.resize(N)
		for i in range(N):
			nodes[i].lon = d2r*lon[i]
			nodes[i].lat = d2r*lat[i]
		
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
	
	[2] Herbert Edelsbrunner et al.: On the Shape of a Set of Points in the
	    Plane, in: IEEE Transactions on Information Theory, Vol. 29, No. 4,
	    July 1983
	"""
	cdef AlphaSpectrum* spectrum
	
	# Constructor:
	def __cinit__(self, np.ndarray[float, ndim=1] lon,
	              np.ndarray[float, ndim=1] lat,
	              VoronoiDelaunayTesselation vdtesselation=None
	    ):
		
		"""
		Initialize an AlphaSpectrum.
		
		Parameters:
		:type lon: 1d numpy array (float)
		:arg  lon: Longitude coordinates of point set.
		
		:type lat: 1d numpy array (float)
		:arg  lat: Latitude coordinates of point set.
		
		:type vdtesselation: VoronoiDelaunayTesselation instance or None
		:arg  vdtesselation: An instance of a Voronoi-Delaunay-tesselation of
		                     the point set if already calculated. Optional: If
		                     omitted, the tesselations will be calculated
		                     internally.
		"""
		
		# Sanity checks:
		cdef size_t N
		N = len(lon)
	
		if (len(lat) != N):
			raise Exception("PyAlphaSpectrum() :\nLength of "
			                "longitude and latitude arrays not equal!")
		
		# Create node vector:
		cdef vector[Node] nodes
		cdef double d2r = np.pi/180.0
		cdef size_t i
		nodes.resize(N)
		for i in range(N):
			nodes[i].lon = d2r*lon[i]
			nodes[i].lat = d2r*lat[i]
		
		
		# If no Voronoi-Delaunay-Tesselation is given, compute it:
		if vdtesselation is None:
			vdtesselation = VoronoiDelaunayTesselation(lon, lat)
		
		
		# Calculate Spectrum:
		cdef VDTesselation* tess = vdtesselation.tesselation
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
		if alpha > 0:
			raise Exception("PyAlphaSpectrum.__call__() : Only implemented "
			                "for alpha<0 !")
		
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
