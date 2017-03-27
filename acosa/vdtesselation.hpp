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
#include <forward_list>
#include <basic_types.hpp>

namespace ACOSA {

class AlphaSpectrum;

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

	friend class AlphaSpectrum;

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


		/*!
		 * \brief Do no consisteny check.
		 */
		constexpr static int CHECK_NOTHING = 0;

		/*! \brief Calculate the Voronoi tesselation and check if it's
		 *         possible to calculate the dual link for each link of the
		 *         Delaunay tesselation.
		 *
		 * Complexity is N*log(N)
		 */
		constexpr static int CHECK_DUAL_LINKS = 1;

		    /*! \brief Calculate the Voronoi tesselation and the Voronoi cells'
			 *         areas, and check if their sum matches up to 4pi within
			 *         a magnitude of the provided tolerance.
			 */
		constexpr static int CHECK_VORONOI_CELL_AREAS = 2;



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
		 * \param checks Determines which checks to do after the Delaunay
		 *               triangulation has been created.
		 *               Default: CHECK_DUAL_LINKS | CHECK_VORONOI_CELL_AREAS
		 * \param on_error_display_nodes If check fails, the node coordinates
		 *                               are written to std::cerr if true.
		 *                               This can be useful for debugging on
		 *                               randomly generated networks.
		 *
		 * This method executes the O(N*log(N)) sweepline algorithm
		 * from [1].
		 *
		 * If errors are detected, an std::runtime_error is thrown.
		 */
		VDTesselation(const std::vector<Node>& nodes,
		              double tolerance = 1e-10,
		              delaunay_algorithm_t algorithm = FORTUNES,
					  int checks = CHECK_DUAL_LINKS | CHECK_VORONOI_CELL_AREAS,
					  bool on_error_display_nodes = true);
		
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
		 * \brief Obtain the set of Delaunay triangles.
		 * \return Constant reference to vector of Delaunay triangles.
		 *
		 * Each triangle contains indices of three nodes of the original set
		 * of nodes that form a Delaunay triangle.
		 */
		const std::vector<Triangle>& delaunay_triangles() const;
		
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
		mutable std::vector<Triangle> delaunay_triangles_;

		/* A cluster merge might have been done. In that case,
		 * the map of Delaunay triangles to Voronoi nodes is surjective
		 * but not injective: Multiple Delaunay triangles may belong
		 * to the same Voronoi nodes (which is the case if more than
		 * 3 nodes of the original network lie on a circle). */
		mutable std::vector<size_t>   delaunay2voronoi;

		/* Also because of the cluster merge, there may be more than one
		 * Delaunay triangle associated with a Voronoi node.
		 * This is a map that stores the Delaunay triangles that contribute
		 * to a Voronoi node. */
		mutable std::vector<std::forward_list<size_t>> voronoi2delaunay;
		
		/* Cached variables: */
		mutable unsigned char cache_state;
		
		/* The original nodes (deleted in case everything else has been
		 * calculated): */
		mutable std::vector<Node> nodes;
		
		/* Delaunay triangulation: */
		mutable std::vector<Link> delaunay_links;

		/* Mapping links of the Delaunay triangulation to links of the
		 * Voronoi tesselation: */
		mutable std::vector<size_t> dual_link_delaunay2voronoi;
		
		/* Voronoi tesselation: */
		mutable std::vector<Node> voronoi_nodes;
		mutable std::vector<Link> voronoi_links;
		mutable std::vector<double> voronoi_areas;
		
		void calculate_delaunay_links() const;
		
		void calculate_voronoi_nodes() const;
		
		void calculate_voronoi_network() const;
		
		void calculate_voronoi_cell_areas() const;

		void calculate_dual_links() const;

		/*!
		 * \brief For a link of the Delaunay triangulation, return the dual
		 *        link of the Voronoi tesselation.
		 * \param link A link of the Delaunay triangulation.
		 * \param node2delaunay A vector of lists that for each node of the
		 *                      original node set contain the Delaunay triangles
		 *                      it is part of.
		 * \return The dual link of the Voronoi tesselation, with indices
		 *         therein.
		 */
		Link dual_link_d2v(const Link& link, const
		                    std::vector<std::forward_list<size_t>>&
		                        node2delaunay) const;

		void merge_clusters() const;
		
		void tidy_up_cache() const;
};


} // NAMESPACE ACOSA

#endif // ACOSA_VDTESSELATION_H
