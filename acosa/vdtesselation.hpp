/* ACOSA class computing the Voronoi tesselation and Delaunay
 * triangulation of a set of nodes on a sphere, using the Fortune's
 * variant in [1].
 * Copyright (C) 2016 Malte Ziebarth
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * Bibliography:
 * [1] Xiaoyu Zheng et al.: A Plane Sweep Algorithm for the Voronoi 
 *     Tesselation of the Sphere, in: electronic-Liquid Crystal
 *     Communications, 2011-12-13
 *     http://www.e-lc.org/docs/2011_12_05_14_35_11
 */


#ifndef ACOSA_VDTESSELATION_H
#define ACOSA_VDTESSELATION_H

#include <vector>
#include <basic_types.hpp>

namespace ACOSA {

/*!
 * \brief A class representing both the Delaunay- and Voronoi-
 *        tesselation on a unit sphere.
 * 
 * This class represents the tesselation of a fixed set of nodes and is,
 * as such, immutably tied to the originally given set.
 * 
 * This class is not (yet?) thread safe.
 */
class VDTesselation {
	public:
		/*! \brief An enumeration of algorithms available to calculate
		 *         the Delaunay tesselation.*/
		enum delaunay_algorithm_t {
			/*! \brief The fortune's sweepline algorithms adapted for
			 *         spherical geometries from [1].
			 *         Its complexity is O(N*log(N))
			 */
			FORTUNES,
			/*! \brief An explicit algorithms that tests triples of
			 *         nodes for the Delaunay condition.
			 *         Its complexity is O(N^4).
			 *         Though not useful for practical purposes, it
			 *         can be used for debugging / testint purposes.
			 *
			 * Note that this algorithm is untested and most likely broken for
			 * lattices that have more than 3 nodes on any circumcircle (e.g.
			 * regular lattices).
			 */
			BRUTE_FORCE
		};


		/* The number of nodes of the original network: */
		const size_t N;
		
		/* The number of nodes in the Voronoi tesselation: */
		size_t size() const;
		
		/*!
		 * \brief Constructs a Vorono-Delaunay-tesselation object from
		 *        a set of nodes.
		 * \param nodes Set of input nodes of length N. May not contain
		 *              duplicates within tolerance. Otherwise, an
		 *              std::domain_error is thrown.
		 * \param tolerance A tolerance to account for numerical errors
		 *                  when calculating the geometric points. This
		 *                  parameter is crucial especially for regular
		 *                  grids, where rounding error may lead to
		 *                  chaotic behaviour of the Fortune's algorithm
		 *                  if tolerance is too low.
		 *                  Random node sets are, empirically more
		 *                  resistant to lower tolerance
		 *
		 * This method executes the O(N*log(N)) sweepline algorithm
		 * from [1].
		 *
		 * If errors are detected, an std::runtime_error is thrown.
		 */
		VDTesselation(const std::vector<Node>& nodes,
		              double tolerance = 1e-8,
					  delaunay_algorithm_t algorithm = FORTUNES);
		
		/*!
		 * \brief Obtain the set of links of the Delaunay triangulation.
		 * \param links Output array of links of the triangulation.
		 * 
		 * Each link contains the indices of two nodes of the original
		 * set of nodes, ordered as given to the constructor, that are
		 * connected in the Delaunay triangulation.
		 * 
		 * Since the links are undirected, they will be return only once
		 * with the smaller index being the first index.
		 */
		void delaunay_triangulation(std::vector<Link>& links) const;
		
		/*!
		 * \brief Obtain the set of spatially embedded nodes and their
		 *        links of the Voronoi tesselation of the source set of
		 *        nodes.
		 * \param voronoi_nodes Output array of spatially embedded
		 *                      Voronoi nodes
		 * \param voronoi_links Output array of links between the
		 *                      Voronoi links.
		 * 
		 * Link indices refer to position in voronoi_links.
		 * 
		 * Since the links are undirected, they will be return only once
		 * with the smaller index being the first index.
		 */
		void voronoi_tesselation(std::vector<Node>& voronoi_nodes,
		                        std::vector<Link>& voronoi_links) const;
		
		/*!
		 * \brief Obtain the areas of the Voronoi cells of the original
		 *        set of nodes.
		 * \param areas Output array of Voronoi cell areas. Its order
		 *              corresponds to the order of the source node set
		 *              given in the constructor.
		 */
		void voronoi_cell_areas(std::vector<double>& areas) const;
		
		/*!
		 * \brief Obtains all nodes of the original network that are
		 * associated with a set of Voronoi nodes.
		 * \param voronoi_nodes Set of indices of the set of Voronoi
		 *                      nodes.
		 * \param associated Output array of indices of the associated
		 *                   nodes in the Voronoi tesselation.
		 * 
		 * The returned indices correspond to the order as given by
		 * voronoi_tesselation().
		 */
		void associated_nodes(const std::vector<size_t>& voronoi_nodes,
			std::vector<size_t>& associated) const;

		/*!
		 * \brief Print a debug output of the current state to standard
		 *        output.
		 */
		void print_debug(bool sort_triangles = true) const;
	
	private:
		/* This variable holds the tolerance that has been set: */
		const double tolerance;

		/* This variable holds the initial delaunay triangulation
		 * in form of a list of triangles. */
		mutable std::vector<Triangle> delaunay_triangles;

		/* A cluster merge might have been done. In that case,
		 * the map of Delaunay triangles to Voronoi nodes is surjective
		 * but not injective: Multiple Delaunay triangles may belong
		 * to the same Voronoi nodes (which is the case if more than
		 * 3 nodes of the original network lie on a circle). */
		mutable std::vector<size_t>   delaunay2voronoi;
		
		/* Cached variables: */
		mutable unsigned char cache_state;
		
		/* The original nodes (deleted in case everything else has been
		 * calculated): */
		mutable std::vector<Node> nodes;
		
		/* Delaunay triangulation: */
		mutable std::vector<Link> delaunay_links;
		
		/* Voronoi tesselation: */
		mutable std::vector<Node> voronoi_nodes;
		mutable std::vector<Link> voronoi_links;
		mutable std::vector<double> voronoi_areas;
		
		void calculate_delaunay_links() const;
		
		void calculate_voronoi_nodes() const;
		
		void calculate_voronoi_network() const;
		
		void calculate_voronoi_cell_areas() const;

		void merge_clusters() const;
		
		void tidy_up_cache() const;
};


} // NAMESPACE ACOSA

#endif // ACOSA_VDTESSELATION_H
