/* ACOSA class computing the (negative) alpha shapes of a set of points on a
 * sphere given the closest-distance Voronoi- & Delaunay-tesselation of the
 * set using the algorithm described in [2].
 * Copyright (C) 2017 Malte Ziebarth
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
 * [2] Herbert Edelsbrunner et al.: On the Shape of a Set of Points in the
 *     Plane, in: IEEE Transactions on Information Theory, Vol. 29, No. 4,
 *     July 1983
 */

#ifndef ALPHASHAPE_HPP
#define ALPHASHAPE_HPP

#include <vdtesselation.hpp>

namespace ACOSA {


/*!
 * \brief This class holds an alpha shape of a set of points.
 *
 * Basically, it consists of a set of node indices, which mark the nodes of the
 * original point set that belong to the alpha shape, and a list of links
 * between those nodes.
 *
 * The links are given in internal coordinates: (i,j), where i and j refer to
 * indices in the node array. Thus, refering to the original network, the
 * nodes (nodes[i],nodes[j]) are linked.
 */
class AlphaShape {
	public:
		/*!
		 * \brief Create an empty alpha shape.
		 */
		AlphaShape();

		AlphaShape(const std::vector<size_t> nodes,
		           const std::vector<Link>& links);

		const std::vector<size_t>& nodes() const;

		const std::vector<Link>& links() const;

		/*!
		 * \brief Returns the links of the alpha shape such that the indices
		 *        refer to nodes of the source point set.
		 * \return
		 */
		std::vector<Link> source_network_links() const;

	private:
		std::vector<size_t> nodes_;
		std::vector<Link>   links_;

};


/*!
 * \brief This class is an implementation of the alpha spectrum used to
 *        compute alpha shapes of sets of points on a sphere.
 *
 * This is an implementation of the algorithm described in [2] for set of points
 * in a plane. It is based on the Voronoi- & Delaunay-tesselations calculated in
 * the VDTesselation class. Since these tesselations are closest-distance, only
 * cases for alpha <= 0 can be handled (alpha=0 resulting in the convex hull
 * calculated by the ConvexHull class. Mind: What should happen here is likely
 * not well-defined or at least unlike the planar case because of the
 * finiteness of the spherical embedding space. Returning the convex hull is
 * thus the programmer's choice.).
 *
 * Because of the spherical topology, some adjustments have to be
 * made compared to the planar case, stemming from the following differences:
 *
 * 1) On a sphere, a convex hull of a set of points may not exist if the
 *    minimum spanning circle covers more than a half-sphere.
 *
 * 2) On a sphere, there are no semi-infinite lines of the Voronoi diagram: It
 *    is a well-defined closed network.
 *
 * The adjustments made to accomodate for these facts are as follows:
 * - Whenever there are cases handling nodes of the convex hull, these cases may
 *   not occur at all in the spherical case. No adjustment needs to be done here
 *   however.
 *
 * - Lemma 3 proof, case c) considers the semi-infinite Voronoi edges of nodes
 *   on the convex hull. Since the Voronoi edges on a sphere are finite, this
 *   case can instead be handled like case a)
 */
class AlphaSpectrum
{
	public:
		/*!
		 * \brief Constructs an AlphaSpectrum object that can be used to
		 *        construct the shapes for different (negative) alpha.
		 * \param nodes The point set for which to compute the shape.
		 * \param tesselation The existing Voronoi- & Delaunay-Tesselation of
		 *                    the point set.
		 *
		 * This methods calculates the alpha bounds for each node of the point
		 * set, see [2], so that no more references to the nodes or the
		 * VDTesselation are required.
		 *
		 * As this method is backed by the calculated closest-distance Voronoi
		 * tesselation, only negative alpha shapes can be calculated.
		 */
		AlphaSpectrum(const std::vector<Node>& nodes,
					  const VDTesselation& tesselation);

		/*!
		 * \brief Obtain the (negative) alpha shape of the original set of
		 *        nodes.
		 * \param alpha The negative inverse radius of the disks that compose
		 *              the alpha shape. Only alpha < 0 accepted.
		 *
		 * Note: This method is not optimal, since it's O(N). Using more
		 *       complicated data structures (see [2]), it can be accelerated
		 *       to O(log(N)+M) where M is the number of elements of the
		 *       shape.
		 */
		AlphaShape operator()(double alpha) const;

	private:
		/* In this vector, we save the upper bounds of alpha, for which each
		 * node begins to be part of the alpha shape.
		 * Note: To increase execution speed, use RTree as described in [2]. */
		std::vector<double> node_max_alpha;

		/* In another vector, we save for each link of the Delaunay the alpha
		 * interval in which it is part of the shape.
		 * Note: To increase execution speed, use RTree as described in [2]. */
		struct interval_t {
			double min;
			double max;
		};

		std::vector<interval_t> alpha_intervals;

		std::vector<Link> delaunay_links;
};

} // NAMESPACE ACOSA

#endif // ALPHASHAPE_HPP
