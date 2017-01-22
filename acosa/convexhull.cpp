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

#include <convexhull.hpp>

#include <cmath>
#include <queue>
#include <math.h>
#include <iostream>

namespace ACOSA {

static bool is_convex(const SphereVectorEuclid& l,
	const SphereVectorEuclid& m, const SphereVectorEuclid& r,
	double lon, const SphereVectorEuclid& x0,
	const SphereVectorEuclid& y0)
{
	/* Calculate the normal vector to the plane that corresponds to the
	 * great circle defined by the segment {l,r}.
	 * Because of the ordering of {l, m, r}, the positive direction of
	 * normal points to the inside of the hull.
	 *       n = (r) x (l-r)
	 * If the dot product of m and the normal vector is negative, then
	 * m lies outside of the hypothetical hull if {l,r} were part of
	 * it. If the dot product is zero, it lies on the hull.
	 * In both cases, the segment {l,m,r} is part of the convex hull. */
	return m * r.cross(l-r) <= 0.0;
}


ConvexHull::ConvexHull(const std::vector<Node>& nodes,
                       const Node& inside)
{
	/* This uses Graham's scan (-like algorithm?). */
	
	/* Convert to different data structure: */
	const size_t N = nodes.size();
	
	
	struct internal_node {
		size_t id;
		double lon; // This comes later.
		SphereVectorEuclid vec;
		
		bool operator>(const internal_node& other) const {
			return lon > other.lon;
		}
	};
	
	std::vector<internal_node> points(N);
	for (size_t i=0; i<N; ++i){
		points[i] = {i, 0.0, SphereVectorEuclid(nodes[i])};
	}
	
	/* First, we look for the node that is furthest away from inside: */
	SphereVectorEuclid z(inside);
	double furthest_dot_product = SphereVectorEuclid(nodes[0]) * z;
	size_t furthest_id = 0;
	for (size_t i=1; i<N; ++i){
		double dotp = points[i].vec * z;
		if (dotp < furthest_dot_product){
			furthest_dot_product = dotp;
			furthest_id = i;
		}
	}
	
	/* A unit vector orthogonal to z with direction x0: */
	SphereVectorEuclid x0(nodes[furthest_id]);
	x0 -= (z * x0)*z;
	x0 *= 1.0 / x0.norm();
	
	/* Now we set z to be the north pole, and the farthest node to be
	 * lon=0. The lon=90° coordinate is defined by the vector
	 *      y0 = (x0) x (z)
	 * Then we can sort nodes by their longitude coordinate. */
	SphereVectorEuclid y0 = x0.cross(z);
	y0 *= 1.0 / y0.norm();
	
	std::priority_queue<internal_node, std::vector<internal_node>,
		std::greater<internal_node>>     point_queue;
		
	for (size_t i=0; i<N; ++i){
		if (i == furthest_id){
			/* Make sure the furthest is actually the first element of
			 * the queue: */
			points[i].lon = -1.0;
		} else {
			double x = x0 * points[i].vec;
			double y = y0 * points[i].vec;
			double lon = std::atan2(y, x);
			if (lon < 0)
				lon += 2*M_PI;
			points[i].lon = lon;
		}
		point_queue.push(points[i]);
	}
	points.clear();
	
	/* Now we can do the scan. */
	std::vector<internal_node> hull;
	
	/* Start with first 2 elements: */
	hull.push_back(point_queue.top());
	point_queue.pop();
	while (!point_queue.empty() && point_queue.top().lon == 0.0){
		hull.push_back(point_queue.top());
		point_queue.pop();
	}
	if (point_queue.empty()){
		/* Degeneracy: All nodes are in a straight line of same longitude.
		 * Thus, the convex hull consists of the two extreme nodes in latitude
		 * direction. */
		if (hull.size() > 1){
			hull[1] = hull.back();
			hull.resize(2);
		}
	} else {
		if (hull.size() > 1){
			hull.resize(1);
		}
		hull.push_back(point_queue.top());
		point_queue.pop();

		/* Now we iterate through the queue. For each element a, we check
		 * the convexity requirement for the nodes {i, i+1, a}, where
		 * i and i+1 are the last elements on the hull stack.
		 * If they are not convex, successively remove the hull stack's
		 * last element until the condition is fulfilled. */
		while (!point_queue.empty()){
			internal_node node = point_queue.top();
			point_queue.pop();

			if (node.lon == hull.back().lon){
				/* If the longitudes are equal, only the southest node
				 * remains: */
				if (node.vec * z < hull.back().vec * z){
					/* The hull's last node was rejected.
					 * We can still reject further nodes. */
					hull.pop_back();
				} else{
					/* The new node was rejected. Continue: */
					continue;
				}
			}

			/* Now remove all nodes from the hull that are not convex
			 * anymore: */
			while (hull.size() > 1 &&
			       !is_convex((hull.end()-2)->vec, hull.back().vec,
			                  node.vec, hull.back().lon, x0, y0))
			{
				hull.pop_back();
			}

			/* Add new node: */
			hull.push_back(node);

		}

		/* Finally, act like we would add the first node again to ensure
		 * periodicity: */
		while (hull.size() > 1 &&
		       !is_convex((hull.end()-2)->vec, hull.back().vec,
		                  SphereVectorEuclid(nodes[furthest_id]),
		                  hull.back().lon, x0, y0))
		{
			hull.pop_back();
		}
	}
	
	/* Now, hull contains the ordered set of the convex hull. */
	hull_segment_normals.resize(hull.size());
	hull_node_ids.resize(hull.size());
	for (size_t i=0; i<hull.size(); ++i){
		hull_node_ids[i] = hull[i].id;
		size_t next = (i+1) % hull.size();
		SphereVectorEuclid v1(nodes[hull[i].id]);
		SphereVectorEuclid v2(nodes[hull[next].id]);
		/* Usually not normalized normal vectors to plane defined by
		 * the great circle through points v1 and v2.
		 * Its positive direction points to the inside. */
		hull_segment_normals[i] = v2.cross(v1-v2); 
	}
}


std::vector<size_t>::const_iterator ConvexHull::begin() const
{
	return hull_node_ids.cbegin();
}
		

std::vector<size_t>::const_iterator ConvexHull::end() const
{
	return hull_node_ids.cend();
}

size_t ConvexHull::size() const
{
	return hull_node_ids.size();
}

		
bool ConvexHull::is_contained(const Node& node) const
{
	SphereVectorEuclid vec(node);
	for (const SphereVectorEuclid& segment : hull_segment_normals){
		if (segment * vec < 0.0)
			return false;
	}
	return true;
}

} // NAMESPACE ACOSA
