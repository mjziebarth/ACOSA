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
cdef extern from "voronoi2/basic_types.hpp" namespace "ACOSA":
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
cdef extern from "voronoi2/vdtesselation.hpp" namespace "ACOSA":
	cdef cppclass VDTesselation :
		const size_t N
		
		size_t size() const
		
		VDTesselation(const vector[Node]& nodes) except +
	
		void delaunay_triangulation(vector[Link]& links) const
	
		void voronoi_tesselation(vector[Node]& voronoi_nodes,
		                         vector[Link]& voronoi_links) const
	
		void voronoi_cell_areas(vector[double]& areas) const
		
		void associated_nodes(const vector[size_t]& voronoi_nodes,
			vector[size_t]& associated) const;

# The ConvexHull class:
cdef extern from "voronoi2/convexhull.hpp" namespace "ACOSA":
	cdef cppclass ConvexHull :
		ConvexHull(const vector[Node]& nodes, const Node& inside)
		
		vector[size_t].const_iterator begin() const
		
		vector[size_t].const_iterator end() const
		
		size_t size() const
		
		bool is_contained(const Node& node) const



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
	              np.ndarray[float, ndim=1] lat
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
		self.tesselation = new VDTesselation(nodes)
		
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
		cdef np.ndarray[long, ndim=1] out = np.zeros(associated.size(), dtype=long)
		for i in range(associated.size()):
			out[i] = associated[i]
		
		return out



####################################################################################
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
