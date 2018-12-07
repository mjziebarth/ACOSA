/* ACOSA class computing the convex hull of a set of nodes on a sphere.
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
 */

#ifndef ACOSA_CONVEXHULL_HPP
#define ACOSA_CONVEXHULL_HPP

#include <basic_types.hpp>
#include <spherics.hpp>

#include <vector>
#include <map>

namespace ACOSA {

class ConvexHull {
	public:
		/*!
		 * \brief Constructs the convex hull of a set of nodes on a unit
		 *        sphere using Graham's scan.
		 * \param nodes Set of nodes whose convex hull is requested
		 * \param inside A vector pointing to the inside of the set.
		 *        The point furthest away from this vector needs to be
		 *        outside the convex hull.
		 * \param tolerance A tolerance parameter for geometric comparisons.
		 * \param mode The mode how to apply the tolerance parameter.
		 *        Default: ToleranceMode::INCLUSIVE.
		 * \param euclidean_backend Whether to use the Euclidean backend
		 *        for convexity assessment (true) or the spherical one (false).
		 *        Default: false
		 * \param sanity_check If true, it is checked whether all nodes are
		 *        contained in the convex hull within tolerance using
		 *        is_contained. If one or more nodes are not contained, a
		 *        std::runtime_error is thrown.
		 * \param throw_on_fail If true, will throw a runtime_exception if the
		 *        hull could not be constructed. Otherwise, will return an
		 *        empty hull in that case.
		 * 
		 * 'inside' is used to sort the nodes for Graham's scan.
		 * Complexity of algorithm is O(N*log(N)).
		 */
		ConvexHull(const std::vector<Node>& nodes, const Node& inside,
		           double tolerance = 1e-12,
		           ToleranceMode mode = ToleranceMode::INCLUSIVE,
		           const bool euclidean_backend=true,
		           bool sanity_check=true,
		           bool throw_on_fail=true);

		std::vector<size_t>::const_iterator begin() const;

		std::vector<size_t>::const_iterator end() const;

		size_t size() const;

		bool empty() const;

		bool is_contained(const Node& node, ToleranceMode mode = ToleranceMode::INCLUSIVE) const;

		/*!
		 * \brief For each node, calculate the shortest distance to the hull's
		 *        border.
		 * \param nodes Vector of nodes to calculate distances to border.
		 * \param distances Target vector of distances. Order equals the order
		 *                  of nodes passed to the constructor.
		 *
		 * Complexity is O(N*M) where N is the number of nodes in the set and
		 * M the number of segments of the hull.
		 */
		void distance_to_border(const std::vector<Node>& nodes,
		                        std::vector<double>& distances) const;

	private:
		/* Hull size: */
		size_t size_;

		/* Point inside the convex hull: */
		SphereVector inside_point;

		/* The indices of the nodes that form the convex hull: */
		std::vector<size_t> hull_node_ids;

		/* Normal vectors onto the great circle planes that form
		 * the hull segments.
		 * This is used in the Euclidean backend. */
		std::vector<SphereVectorEuclid>   hull_segment_normals;

		/* The vectors of the nodes that form the convex hull.
		 * This is used in the spherical backend. */
		std::map<double,SphereVector> hull_vertices;
		double azimuth_0;

		/* A tolerance parameter for determining whether a node is inside
		 * the hull: */
		double tolerance;
		ToleranceMode tolerance_mode;

		/* A switch to toggle between the backends: */
		bool euclidean_backend;


		void init_euclidean(const std::vector<Node>& nodes);

		void init_spherical(const std::vector<Node>& nodes);
};

} // NAMESPACE ACOSA

#endif // ACOSA_CONVEXHULL_HPP
