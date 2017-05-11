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

#include <alphaspectrum.hpp>
#include <convexhull.hpp>
#include <unordered_set>
#include <cmath>

#include <iostream>
#include <limits>

using ACOSA::AlphaSpectrum;
using ACOSA::AlphaShape;

//------------------------------------------------------------------------------
static void calculate_convex_hull_members(const std::vector<ACOSA::Node>& nodes,
                                          std::vector<bool>& is_ch_node,
                                          const std::vector<ACOSA::Link>& links,
                                          std::vector<bool>& is_ch_link)
{
	const size_t N = nodes.size();
	const size_t M = links.size();

	/* Resize vectors first: */
	is_ch_node.resize(N, false);
	is_ch_link.resize(M, false);

	/* Calculate convex hull.
	 * A vector inside the convex hull will be calculated as the node sets
	 * mean (replace this by smalles circle). */
	ACOSA::SphereVectorEuclid mean_vec(0,0,0);
	for (size_t i=0; i<N; ++i){
		mean_vec += ACOSA::SphereVectorEuclid(nodes[i].lon, nodes[i].lat);
	}
	mean_vec /= mean_vec.norm();

	ACOSA::ConvexHull hull(nodes, ACOSA::Node(mean_vec.lon(), mean_vec.lat()));

	/* Mark all nodes inside the hull and create ordered set of hull links: */
	std::unordered_set<ACOSA::Link> hull_links;
	auto it = hull.begin();
	if (it != hull.end()){
		/* Mark first node: */
		size_t i = *it;
		is_ch_node[i] = true;

		/* For all following nodes, we mark the node and create a link to
		 * the previous node: */
		for (++it; it != hull.end(); ++it)
		{
			/* Mark node: */
			size_t j = *it;
			is_ch_node[j] = true;

			/* Create link: */
			if (i < j){
				hull_links.emplace(i,j);
			} else {
				hull_links.emplace(j,i);
			}

			/* Next link starts from this node: */
			i=j;
		}

		/* Create last link that completes the circle: */
		if (hull.size() > 1){
			size_t j=*hull.begin();
			if (i < j){
				hull_links.emplace(i,j);
			} else {
				hull_links.emplace(j,i);
			}
		}
	}

	/* Mark all links inside the hull: */
	for (size_t i=0; i<M; ++i){
		if (hull_links.count(links[i]))
			is_ch_link[i] = true;
	}
}


//------------------------------------------------------------------------------
AlphaSpectrum::AlphaSpectrum(const std::vector<ACOSA::Node>& nodes,
                             const ACOSA::VDTesselation& tesselation)
    : node_max_alpha(nodes.size())
{
	const size_t N = nodes.size();

	/* Step 1: Ensure tesselation has all necessary variables cached.
	 *         Then initialize the vectors now that we're sure we've got the
	 *         right dimensions: */
	tesselation.calculate_voronoi_nodes();
	tesselation.calculate_delaunay_links();
	tesselation.calculate_dual_links();


	/* Step 2: Calculate the convex hull of the node set: */
	std::vector<bool> is_convex_hull_node;
	std::vector<bool> is_convex_hull_link;
	calculate_convex_hull_members(nodes, is_convex_hull_node,
	                              tesselation.delaunay_links,
	                              is_convex_hull_link);


	/* Step 3: Create a map from node indices to associated Voronoi cells.
	 *         The map may contain duplicates as some Delaunay triangles
	 *         may have been merged to the same Voronoi node, but that
	 *         should not make a big performance difference in most cases,
	 *         alas we do not check for duplicates. */
	std::vector<std::forward_list<size_t>> node2delaunay(N);
	for (size_t i=0; i<tesselation.delaunay_triangles_.size(); ++i){
		const Triangle& t = tesselation.delaunay_triangles_[i];
		node2delaunay[t.i].push_front(i);
		node2delaunay[t.j].push_front(i);
		node2delaunay[t.k].push_front(i);
	}


	alpha_intervals.resize(tesselation.delaunay_links.size());
	delaunay_links = tesselation.delaunay_links;


	/* Step 4: For each node, calculate the maximum alpha where it is
	 * alpha-extreme (read: part of the alpha-shape): */
	for (size_t i=0; i<N; ++i){
		if (is_convex_hull_node[i]){
			/* Case a) in [2], Lemma 2:
			 * Convex hull points are trivially alpha-extreme for non-positive
			 * alpha, which we're limited to here: */
			node_max_alpha[i] = 0.0;
		} else {
			/* Case b) in [2], Lemma 2:
			 * alph_max = -1/d_max, where d_max is the maximum distance of
			 * node i to a point x of its Voronoi cell. This point is one
			 * of the Voronoi cell's vertices, thus we need to check all
			 * associated Voronoi vertices and take the distance of the
			 * farthest away: */
			ACOSA::SphereVector vec(nodes[i].lon, nodes[i].lat);
			double max_dist = 0.0;
			for (size_t j : node2delaunay[i]){
				const Node& n2 = tesselation.voronoi_nodes[tesselation
				                    .delaunay2voronoi[j]];
				double d = vec.distance(ACOSA::SphereVector(n2.lon, n2.lat));
				if (d > max_dist)
					max_dist = d;
			}
			node_max_alpha[i] = -1.0/max_dist;

		}
	}


	/* Step 5: For each link of the Delaunay-tesselation, calculate the alpha
	 *         bounds inside which it is part of the alpha shape of the point
	 *         set: */
	for (size_t i=0; i<tesselation.delaunay_links.size(); ++i){
		/* Obtain the id (in terms of Voronoi link array) of the i'th Delaunay
		 * tesselation link. In cases where there are more than 3 cocircular
		 * points, the Delaunay tesselation is not unique and, more importantly,
		 * two or more Voronoi nodes are equal, meaning they have a link of
		 * spatial length 0.
		 * Convention here is that between these nodes (they are merged) there
		 * is no link. Thus, the corresponding Delaunay tesselation has no dual
		 * link, which we have to check.*/
		size_t dual_link = tesselation.dual_link_delaunay2voronoi[i];
		if (dual_link == ACOSA::NO_LINK)
			continue;

		/* If we have a dual link, obtain both the Delaunay and the dual Voronoi
		 * link: */
		ACOSA::Link dl = tesselation.delaunay_links[i];
		ACOSA::Link vl = tesselation.voronoi_links[dual_link];

		/* Case a) in [2], Lemma 3:
		 *		alpha \in [alpha_min, alpha_max]
		 * where alpha_min = -1/a(j,i), alpha_max = -1/b(j,i)
		 * and a(j,i) and b(j,i) are the closest and farthest points from
		 * node j (one of the nodes j1 and j2 belonging to the Voronoi edge) to
		 * points of the selected Voronoi edge which belongs to the Delaunay
		 * link i.
		 *
		 * We know the following things:
		 * - a, the distance to the closest point on the Voronoi edge, is either
		 *    - half the distance between the points j1 and j2   or
		 *    - the shortest of the distances to the Voronoi nodes.
		 * - b, the distance to the farthest point on the Voronoi edge, is
		 *   the distance to one of the end points of the Voronoi edge (-->
		 *   one of the Voronoi nodes).
		 *
		 * We can thus calculate a and b. */
		ACOSA::SphereVector n1(nodes[dl.i].lon, nodes[dl.i].lat);
		ACOSA::SphereVector n2(nodes[dl.j].lon, nodes[dl.j].lat);
		ACOSA::SphereVector vn1(tesselation.voronoi_nodes[vl.i].lon,
		                        tesselation.voronoi_nodes[vl.i].lat);
		ACOSA::SphereVector vn2(tesselation.voronoi_nodes[vl.j].lon,
		                        tesselation.voronoi_nodes[vl.j].lat);

		/* Calculating b is straightforward: */
		double d_n1_v1 = n1.distance(vn1);
		double d_n1_v2 = n1.distance(vn2);
		double b = std::max(d_n1_v1,d_n1_v2);

		/* To determine a, we need to determine whether the line between the two
		 * Voronoi nodes crosses the line between the two nodes.
		 * That's the case if the distance from vn1 to vn2 is longer than the
		 * distance from vn1 to the middle between n1 and n2.
		 * Using spherical geometry, these distances are: */
		double halfdist = 0.5*n1.distance(n2);
		double d_cmp_1 = vn1.distance(vn2);
		double d_cmp_2 = std::acos(std::cos(b)/std::cos(halfdist));

		double a = (d_cmp_1 > d_cmp_2) ? halfdist : std::min(d_n1_v1, d_n1_v2);

		alpha_intervals[i].min = -1.0/a;
		alpha_intervals[i].max = -1.0/b;

	}
}

//------------------------------------------------------------------------------
AlphaShape AlphaSpectrum::operator()(double alpha) const
{
	/* If (alpha > -1/pi), the alpha shape is empty since pi is the maximum
	 * radius on the sphere: */
	std::vector<size_t> shape_nodes;
	std::vector<ACOSA::Link> shape_links;
	if (alpha > -1.0/M_PI)
		return AlphaShape(shape_nodes, shape_links);

	/* Determine all nodes of the shape at current alpha. Also generate a
	 * map to map these nodes' indices in the original point set to indices in
	 * the alpha shape point set: */
	size_t N = 0;
	for (double a : node_max_alpha){
		if (alpha < a)
			++N;
	}
	shape_nodes.reserve(N);
	for (size_t i=0; i<node_max_alpha.size(); ++i){
		if (alpha < node_max_alpha[i]){
			shape_nodes.push_back(i);
		}
	}

	/* Determine all links of the shape at current alpha: */
	size_t M = 0;
	for (const interval_t& I : alpha_intervals){
		if (alpha <= I.max && alpha >= I.min)
			++M;
	}
	shape_links.reserve(M);
	for (size_t i=0; i<alpha_intervals.size(); ++i){
		if (alpha <= alpha_intervals[i].max &&
		    alpha >= alpha_intervals[i].min)
		{
			shape_links.emplace_back(delaunay_links[i].i,
			                         delaunay_links[i].j);
		}
	}

	/* Create alpha shape: */
	return AlphaShape(shape_nodes, shape_links);
}

//------------------------------------------------------------------------------
AlphaShape::AlphaShape()
{
}

//------------------------------------------------------------------------------
AlphaShape::AlphaShape(const std::vector<size_t> nodes,
                       const std::vector<ACOSA::Link>& links)
    : nodes_(nodes)
{
	/* Determine maximum node index: */
	size_t max_i = 0;
	for (size_t i : nodes){
		if (i > max_i)
			max_i = i;
	}

	/* Create a map from external node indices (which refer to the original
	 * point set) to internal indices (which refer to the nodes of the alpha
	 * shape): */
	std::vector<size_t> ext2int(max_i+1);
	for (size_t i=0; i<nodes.size(); ++i){
		ext2int[nodes[i]] = i;
	}

	/* Copy links but adjust indices so that they refer to the indices in this
	 * object's vector: */
	links_.resize(links.size());
	for (size_t i=0; i<links.size(); ++i){
		links_[i] = ACOSA::Link(ext2int[links[i].i], ext2int[links[i].j]);
	}
}

//------------------------------------------------------------------------------
const std::vector<size_t>&AlphaShape::nodes() const
{
	return nodes_;
}

//------------------------------------------------------------------------------
const std::vector<ACOSA::Link>&AlphaShape::links() const
{
	return links_;
}

//------------------------------------------------------------------------------
std::vector<ACOSA::Link> AlphaShape::source_network_links() const
{
	std::vector<ACOSA::Link> source_links(links_.size());
	for (size_t i=0; i<links_.size(); ++i){
		source_links[i] = ACOSA::Link(nodes_[links_[i].i], nodes_[links_[i].j]);
	}
	return source_links;
}
