#!/bin/python
#
# Test script for the various parts of ACOSA.
# Copyright (C) 2016 Malte Ziebarth

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
# 

import numpy as np
from acosa import VoronoiDelaunayTesselation as VDT

######################################################################################
#                                Configuration:                                      #
######################################################################################
random_network_sizes = [20, 2000, 200000]
random_runs          = [10000, 1000, 10]

regular_lattice_1_sizes = np.arange(20,20000,20,dtype=int)
regular_lattice_2_sizes = regular_lattice_1_sizes



######################################################################################
#                               Helper functions:                                    #
######################################################################################
def longitude_grid_points(N):
	# Start a guess such that guess*(guess/2) = N
	guess = int(np.sqrt(2*N))

	# Now find biggest smaller factor of N:
	while guess != 0 and (N % guess) != 0:
		guess -= 1

	if guess == 0:
		raise Exception("ERROR : N does not have any factors smaller than sqrt(N/2)!")

	return guess



######################################################################################
#                                     Script:                                        #
######################################################################################

# 1) Random lattices:
print("1) Random lattices.")

for i in range(len(random_network_sizes)):
	N = random_network_sizes[i]
	print("  N="+str(N))
	for j in range(random_runs[i]):
		# Create random nodes:
		lon = 360.0*np.random.random(N)
		lat = 180.0/np.pi*np.arcsin(2*np.random.random(N)-1)
	
		# Create Delaunay triangulation, Voronoi tesselation and Voronoi cell
		# areas:
		vdt =VDT(np.array(lon,dtype=np.float32), np.array(lat,dtype=np.float32))
		vdt.delaunay_triangulation()
		vdt.voronoi_tesselation()


# 2) Regular lattices 1:
print("2) Regular lattices (rectangular lon/lat grid).")

for N in regular_lattice_1_sizes:
	# Create grid:
	n_lon = longitude_grid_points(N)
	n_lat = int(N / n_lon)

	i,j = np.meshgrid(np.arange(n_lon), np.arange(n_lat))

	lon = (360.0*i/float(n_lon)).flatten()
	lat = (-180.0*(0.5-(j+1.0)/(n_lat+1.0))).flatten()
	
	# Create Delaunay triangulation, Voronoi tesselation and Voronoi cell
	# areas:
	vdt =VDT(np.array(lon,dtype=np.float32), np.array(lat,dtype=np.float32))
	vdt.delaunay_triangulation()
	vdt.voronoi_tesselation()


# 3) Regular lattices 2:
print("3) Regular lattices (hexagonalish lon/lat grid).")

for N in regular_lattice_2_sizes:
	# Create grid:
	n_lon = longitude_grid_points(N)
	n_lat = int(N / n_lon)

	i,j = np.meshgrid(np.arange(n_lon), np.arange(n_lat))

	lon = (360.0*(i + 0.5*(j % 2))/n_lon).flatten()
	lat = (-180.0*(0.5-(j+1.0)/(n_lat+1.0))).flatten()
	
	# Create Delaunay triangulation, Voronoi tesselation and Voronoi cell
	# areas:
	vdt =VDT(np.array(lon,dtype=np.float32), np.array(lat,dtype=np.float32))
	vdt.delaunay_triangulation()
	vdt.voronoi_tesselation()



print("finished!")
