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
		 * 
		 * 'inside' is used to sort the nodes for Graham's scan.
		 * Complexity of algorithm is O(N*log(N)).
		 */
		ConvexHull(const std::vector<Node>& nodes, const Node& inside);
		
		std::vector<size_t>::const_iterator begin() const;
		
		std::vector<size_t>::const_iterator end() const;
		
		size_t size() const;
		
		bool is_contained(const Node& node) const;
	
	private:
		/* The indices of the nodes that form the convex hull: */
		std::vector<size_t> hull_node_ids;
		
		/* Normal vectors onto the great circle planes that form
		 * the hull segments.
		 * They are not necessarily normalized. */
		std::vector<SphereVectorEuclid>   hull_segment_normals;
};

} // NAMESPACE ACOSA

#endif // ACOSA_CONVEXHULL_HPP
