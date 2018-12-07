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

//#define _USE_MATH_DEFINES

#include <cmath>
#include <queue>
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <utility>
#include <iomanip>

namespace ACOSA {


constexpr double TWO_PI = 2*M_PI;

/* Euclidean code: */

static bool is_convex(const SphereVectorEuclid& l,
    const SphereVectorEuclid& m, const SphereVectorEuclid& r,
    double tolerance, ToleranceMode mode)
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
	const SphereVectorEuclid normal = r.cross(l-r);
	return (m * normal) / normal.norm() <= tolerance;
}



void ConvexHull::init_euclidean(const std::vector<Node>& nodes)
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
	SphereVectorEuclid z(inside_point);
	double furthest_dot_product = points[0].vec * z;
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
				lon += TWO_PI;
			points[i].lon = lon;
		}
		point_queue.push(points[i]);
	}
	points.clear();

	/* Now we can do the scan. */
	std::vector<internal_node> hull;

	/* Start with first 2 elements: */
	hull.push_back(point_queue.top());
	hull[0].lon = 0.0;
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
				/* If the longitudes are equal, only the southernmost node
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
			                  node.vec, tolerance, tolerance_mode))
			{
				hull.pop_back();
			}

			/* Add new node: */
			hull.push_back(node);
		}

		/* Finally, act like we would add the first node again to ensure
		 * periodicity: */
		while (hull.size() > 2 &&
		       !is_convex((hull.end()-2)->vec, hull.back().vec, hull[0].vec,
		                  tolerance, tolerance_mode))
		{
			hull.pop_back();
		}
	}

	/* Now, hull contains the ordered set of the convex hull. */
	size_ = hull.size();
	hull_segment_normals.resize(size_);
	hull_node_ids.resize(size_);
	for (size_t i=0; i<size_; ++i){
		hull_node_ids[i] = hull[i].id;
		size_t next = (i+1) % size_;
		SphereVectorEuclid v1(nodes[hull[i].id]);
		SphereVectorEuclid v2(nodes[hull[next].id]);
		/* Usually not normalized normal vectors to plane defined by
		 * the great circle through points v1 and v2.
		 * Its positive direction points to the inside. */
		SphereVectorEuclid v3 = v2.cross(v1-v2);

		/* Norm: */
		v3 *= 1.0/v3.norm();

		hull_segment_normals[i] = v3;
	}
}






/*************************************************************************
 *                                                                       *
 *                      SPHERICAL INITIALIZATION                         *
 *                                                                       *
 *************************************************************************/

/* Stereographic projection: */
struct Vector2d {
	Vector2d(double x, double y) : x(x),y(y)
	{};

	double x;
	double y;
};

static Vector2d oblique_stereographic(const SphereVector& x, const SphereVector& x0)
{
	const double dlon = x.lon()-x0.lon();
	const double cos_lat = cos(x.lat());
	return Vector2d(cos_lat * sin(dlon),
	                cos(x0.lat()) * sin(x.lat()) - sin(x0.lat()) * cos_lat * cos(dlon));
}

static bool is_convex_stereographic(const SphereVector& l,
    const SphereVector& m, const SphereVector& r, const SphereVector& x0,
    double tolerance, ToleranceMode mode)
{
	/* Determine whether the set {l,m,r} is convex or collinear,
	 * seen from x0, or not (i.e. concave).
	 * We do this in a stereographic projection centered on m in
	 * this way:
	 *  0) Project to stereographic projection
	 *  1) Calculate azimuths of l, r, and x0.
	 *  2) Order the azimuths starting from l, assign azi(l)=0.
	 *     Then the collinear case of {l,m,r} is when azi(r)=180°.
	 *     This case is also convex. Other cases are convex if
	 *     azi(r) is closer to 180° than azi(x0) and the differences
	 *     to 180° have the same sign.
	 *
	 * First project to oblique stereographic projection at point m,
	 * following Snyder (1987): */
	auto l_2d = oblique_stereographic(l, m);
	auto r_2d = oblique_stereographic(r, m);
	auto x0_2d = oblique_stereographic(x0, m);

	/* Determine the azimuths in these projections: */
	const double a_l  = atan2(l_2d.y, l_2d.x);
	double a_r  = atan2(r_2d.y, r_2d.x);
	double a_x0 = atan2(x0_2d.y, x0_2d.x);

	/* Gauge them to the azimuth of point l: */
	a_r = fmod(a_r - a_l + TWO_PI, TWO_PI);
	a_x0 = fmod(a_x0 - a_l + TWO_PI, TWO_PI);

	/* Compute the difference to 180°: */
	const double d_r  = M_PI - a_r;
	const double d_x0 = M_PI - a_x0;

	/* Now we have convexity if |d_r| <= |d_x0| and the signs
	 * are equal. Continue depending on tolerance mode: */
	if (mode == ToleranceMode::INCLUSIVE){
		return std::abs(d_r) <= std::abs(d_x0) + tolerance
		       && ((d_r < 0.0) == (d_x0 < 0.0)
		           || std::abs(d_r) < tolerance || std::abs(d_x0) < tolerance);

	} else if (mode == ToleranceMode::EXCLUSIVE) {
		return std::abs(d_r) <= std::abs(d_x0) - tolerance
		       && (d_r < 0.0) == (d_x0 < 0.0);

	} else if (mode == ToleranceMode::EXACT) {
		return std::abs(d_r) <= std::abs(d_x0)
		       && (d_r < 0.0) == (d_x0 < 0.0);
	}

	/* This should never be reached. */
	throw std::runtime_error("Unreacheable code reached!");
}



void ConvexHull::init_spherical(const std::vector<Node>& nodes)
{
	/* This uses Graham's scan (-like algorithm?). */

	/* Convert to different data structure: */
	const size_t N = nodes.size();

	struct internal_node {
		size_t id;
		double azimuth; // This comes later.
		SphereVector vec;

		bool operator>(const internal_node& other) const {
			return azimuth > other.azimuth;
		}
	};

	std::vector<internal_node> points(N);
	for (size_t i=0; i<N; ++i){
		points[i] = {i, inside_point.azimuth_to(nodes[i]), SphereVector(nodes[i])};
	}

	/* First, we look for the node that is furthest away from inside.
	 * That extreme point is part of the convex hull: */
	double furthest_distance = inside_point.distance(points[0].vec);
	size_t furthest_id = 0;
	for (size_t i=1; i<N; ++i){
		double dist = inside_point.distance(points[i].vec);
		if (dist > furthest_distance){
			furthest_distance = dist;
			furthest_id = i;
		}
	}

	/* Calculate the azimuth to that node: */
	azimuth_0 = points[furthest_id].azimuth;

	/* Calculate all nodes' azimuths relative to that azimuth: */
	for (internal_node& n : points){
		n.azimuth = std::max(std::fmod(n.azimuth - azimuth_0 + TWO_PI, TWO_PI), 0.0);
	}

	/* Make sure the furthest is actually the first element of
	 * the queue: */
	points[furthest_id].azimuth = -1.0;

	/* Order all points by increasing relative azimuth: */
	std::priority_queue<internal_node,
	                    std::vector<internal_node>,
	                    std::greater<internal_node>>
	    point_queue(std::greater<internal_node>(), std::move(points));

	/* Now we can do the scan. */
	std::vector<internal_node> hull;

	/* Start with first 2 elements: */
	hull.push_back(point_queue.top());
	point_queue.pop();
	while (!point_queue.empty() && point_queue.top().azimuth == 0.0){
		hull.push_back(point_queue.top());
		point_queue.pop();
	}
	if (point_queue.empty()){
		/* Degeneracy: All nodes are in a straight line of same azimuth from center.
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

			if (node.azimuth == hull.back().azimuth){
				/* If the azimuths are equal, only the southernmost node
				 * remains: */
				if (node.vec.distance(inside_point) > hull.back().vec.distance(inside_point))
				{
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
				   !is_convex_stereographic((hull.end()-2)->vec, hull.back().vec,
				              node.vec, inside_point, tolerance, tolerance_mode))
			{
				hull.pop_back();
			}

			/* Add new node: */
			hull.push_back(node);

		}

		/* Finally, act like we would add the first node again to ensure
		 * periodicity: */
		/* Use the spherical backend: */
		while (hull.size() > 2 &&
		       !is_convex_stereographic((hull.end()-2)->vec, hull.back().vec,
		                  hull[0].vec, inside_point, tolerance,
		                  tolerance_mode))
		{
			hull.pop_back();
		}
	}

	/* Set first element proper azimuth again: */
	if (hull[0].azimuth < 0.0)
		hull[0].azimuth = 0.0;

	/* Now, hull contains the ordered set of the convex hull. */
	size_ = hull.size();
	hull_node_ids.resize(size_);
	for (size_t i=0; i<size_; ++i){
		/* Save ID and vector: */
		hull_node_ids[i] = hull[i].id;
		hull_vertices[hull[i].azimuth] = hull[i].vec;
	}
}


/*************************************************************************
 *                                                                       *
 *                        GENERAL CLASS METHODS                          *
 *                                                                       *
 *************************************************************************/




ConvexHull::ConvexHull(const std::vector<Node>& nodes,
                       const Node& inside, double tolerance,
                       ToleranceMode mode, const bool euclid_backend,
                       bool sanity_check, bool throw_on_fail)
    : tolerance(tolerance), euclidean_backend(euclid_backend),
      tolerance_mode(mode), inside_point(inside)
{
	/* Initialize: */
	if (euclidean_backend){
		init_euclidean(nodes);
	} else {
		init_spherical(nodes);
	}

	/* Sanity checks: */
	if (sanity_check){
		/* Sanity check 1: Throw if hull is empty or smaller than plausible: */
		if (size_ < 3 && size_ < nodes.size()){
			/* TODO : Test this part of the code. This is written and never
			 *        tested. */
			throw std::runtime_error("Degenerate case not handled.");
			std::cerr << "Warning: You are entering uncharted territory in the "
			             "source code! (convexhull.cpp, lines " << __LINE__ << "+)\n";
			if (size_ == 0){
				/* If nodes exist, hull must contain at least one element. */
				if (throw_on_fail)
					throw std::runtime_error("Convex hull failed: Node set not "
					                         "empty but hull is.");
			} else if (size_ == 1){
				/* If hull is just one element, all nodes need to be equal. */
				SphereVectorEuclid vc(nodes[0]);
				if (!std::all_of(nodes.cbegin(),nodes.cend(),
				                 [tolerance,&vc](const SphereVectorEuclid& v)
				                   {return v.distance(vc) <= tolerance;}))
				{
					if (throw_on_fail){
						throw std::runtime_error("Convex hull failed: Hull contains "
						                         "only one element but node set "
						                         "is not singular.");
					} else {
						hull_node_ids.clear();
						hull_segment_normals.clear();
					}
				}
			} else if (size_ == 2){
				/* If hull is just two elements, all elements need to be
				 * collinear, i.e. all the hull segment normals need to
				 * be equal. */
				SphereVectorEuclid vc(hull_segment_normals[0]);
				if (!std::all_of(hull_segment_normals.cbegin(),
				                 hull_segment_normals.cend(),
				                 [tolerance,&vc](const SphereVectorEuclid& v)
				                   {return std::min(v.distance(vc),v.distance(-vc))
				                           <= tolerance;}
				                ))
				{
					if (throw_on_fail){
						throw std::runtime_error("Convex hull failed: Hull contains "
						                         "only two elements but is not "
						                         "collinear.");
					} else {
						hull_node_ids.clear();
						hull_segment_normals.clear();
					}
				}
			}
			/* End of TODO */
		}

		/* Sanity check 2: Throw if not all nodes contained. */
		size_t contained = 0;
		for (const Node& n: nodes){
			if (is_contained(n))
				++contained;
		}
		if (contained != nodes.size()){
			if (throw_on_fail){
				/* Throw if hull could not be constructed: */
				std::cerr << "#nodes:     " << nodes.size() << "\n";
				std::cerr << "#contained: " << contained << "\n";
				std::cerr << "#hull:      " << hull_node_ids.size() << "\n";
				throw std::runtime_error("Convex hull failed: Not all "
					                     "nodes contained.\n");
			} else {
				/* This is the more silent way: Simply continue with
				 * empty hull. */
				hull_node_ids.clear();
				hull_segment_normals.clear();
				hull_vertices.clear();
			}
		}
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
	return size_;
}

bool ConvexHull::empty() const
{
	return hull_node_ids.empty();
}

bool ConvexHull::is_contained(const Node& node, ToleranceMode mode) const
{
	/* Sanity: */
	if (size_ == 0){
		return false;
	}

	/* Switch between backends: */
	if (euclidean_backend){
		SphereVectorEuclid vec(node);
		/* The direction of the hill segment normals is the inside of the
		 * hull. Thus, a vector is inside the hull (we include the border) if
		 * its dot product with all segment normals is positive or zero.
		 * Thus, if any dot product is negative, it lies outside the hull.
		 * This is the exlusion criterion. To implement the tolerance,
		 * proceed as follows:
		 *    - INCLUSIVE:  Allow small negative dot products with a maximum
		 *                  value of -tolerance.
		 *    - EXCLUSIVE:  Be strict about containment, exclude also dot
		 *                  product values smaller than tolerance.
		 *    - EXACT:      Exclude all dot product values smaller zero,
		 *                  no tolerance in either direction.
		 */
		const double compare = (mode == ToleranceMode::INCLUSIVE) ? -tolerance :
			                   (mode == ToleranceMode::EXCLUSIVE) ? tolerance :
			                   0.0;

		for (const SphereVectorEuclid& segment : hull_segment_normals){
			if (segment * vec < compare){
				return false;
			}
		}
	} else {
		/* One more sanity check: */
		if (size_ == 1){
			return hull_vertices.begin()->second == node;
		}
		/* Calculate the azimuth of the looked-up node relative to the
		 * reference point: */
		double azi = std::max(std::fmod(inside_point.azimuth_to(node) - azimuth_0 + TWO_PI,
		                                TWO_PI), 0.0);

		/* Find the hull segment in that direction: */
		auto it = hull_vertices.lower_bound(azi);
		if (it == hull_vertices.end()){
			/* Azimuth in wraparound interval: */
			return !is_convex_stereographic((--it)->second, node, hull_vertices.begin()->second,
			                                inside_point, tolerance, ToleranceMode::EXCLUSIVE);
		} else if (it->first == azi){
			/* The node's azimuth coincides with one hull vertex, reducing
			 * the problem to a distance comparison: */
			return inside_point.distance(node) <= inside_point.distance(it->second);
		} else {
			/* We know that it-1 and it are proper elements of hull_vertices
			 * and that they have strictly smaller and larger azimuth,
			 * respectively. Thus, call convexity check: */
			auto right = it->second;
			auto left = (--it)->second;
			return !is_convex_stereographic(left, node, right, inside_point,
			                                tolerance, ToleranceMode::EXCLUSIVE);
		}
		std::cout << "Spherical backend not implemented.\n";
		exit(-1);
	}

	return true;
}

void ConvexHull::distance_to_border(const std::vector<Node>& nodes,
                                    std::vector<double>& distances) const
{
	/* Sanity check: */
	if (hull_segment_normals.empty()){
		distances.resize(nodes.size(), 0.0);
		return;
	}

	const size_t N = nodes.size();
	distances.resize(N);

	for (size_t i=0; i<N; ++i){
		SphereVectorEuclid vec(nodes[i]);

		/* Iterate over all border segments and find the one with smallest
		 * distance: */
		double closest_dotp = std::abs(hull_segment_normals.front() * vec);
		for (size_t j=1; j<hull_segment_normals.size(); ++j){
			double dotp = vec*hull_segment_normals[j];
			if (dotp < -tolerance){
				throw std::domain_error("ConvexHull::distance_to_border():\n"
				                        "Node not inside hull.\n");
			} else if (std::abs(dotp) < closest_dotp){
				closest_dotp = std::abs(dotp);
			}
		}

		/* Convert that dot product to a distance.
		 * Since hull_segment_normals are normed to 1.0, their dot product
		 * with vec equals the z coordinate of vec in a coordinate system
		 * with z axis parallel to the respective segment normal vector.
		 * We can set x=sqrt(1-z^2) and the resulting vector's latitude will
		 * be vecs distance to the closest hull segment. */
		distances[i] = SphereVectorEuclid(std::sqrt(closest_dotp), 0.0,
		                                  closest_dotp).lat();
	}
}

} // NAMESPACE ACOSA
