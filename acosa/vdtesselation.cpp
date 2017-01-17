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

#include <vdtesselation.hpp>
#include <spherics.hpp>
#include <fortunes_sphere.hpp>
#include <geometricgraph.hpp>

#include <map>
#include <set>
#include <forward_list>
#include <iostream>
#include <math.h>
#include <algorithm>

namespace ACOSA {

/* Cache states: */
static constexpr unsigned char DELAUNAY_LINKS_CACHED = 1;
static constexpr unsigned char VORONOI_NODES_CACHED =  2;
static constexpr unsigned char VORONOI_LINKS_CACHED =  4;
static constexpr unsigned char VORONOI_CELLS_CACHED =  8;


static void
delaunay_triangulation_brute_force(const std::vector<Node>& nodes,
                                   std::vector<Triangle>& triangles,
                                   double tolerance)
{
	const size_t N = nodes.size();

	/* Create internal data structure for nodes: */
	std::vector<SphereVectorEuclid> v(N);
	for (size_t i=0; i<N; ++i){
		v[i] = SphereVectorEuclid(nodes[i].lon, nodes[i].lat);
	}

	for (size_t i=0; i<N; ++i){
		for (size_t j=i+1; j<N; ++j){
			for (size_t k=j+1; k<N; ++k){
				/* Step 1: Check triangle (i,j,k) */
				SphereVectorEuclid cc =
				    SphereVectorEuclid::circumcenter(v[i], v[j], v[k]);
				double max_distance =
				        (cc.distance(v[i]) + cc.distance(v[j])
				         + cc.distance(v[k])) / 3.0
				        - tolerance;
				bool is_delaunay = true;
				for (size_t m=0; m<N; ++m){
					if (m != i && m != j && m != k &&
					    cc.distance(v[m]) <= max_distance)
					{
						is_delaunay = false;
					}
				}
				if (is_delaunay){
					triangles.emplace_back(i,j,k);
				}
				/* Step 2: Check triangle (i,k,j)
				 *         (the antipodal one) */
				cc = SphereVectorEuclid::circumcenter(v[i], v[k],
				                                      v[j]);
				max_distance =
				        (cc.distance(v[i]) + cc.distance(v[j])
				         + cc.distance(v[k])) / 3.0
				        - tolerance;
				is_delaunay = true;
				for (size_t m=0; m<N; ++m){
					if (m != i && m != j && m != k &&
					    cc.distance(v[m]) <= max_distance)
					{
						is_delaunay = false;
					}
				}
				if (is_delaunay){
					triangles.emplace_back(i,k,j);
				}
			}
		}
	}
}



//----------------------------------------------------------------------
VDTesselation::VDTesselation(const std::vector<Node>& nodes,
                             double tolerance,
                             delaunay_algorithm_t algorithm)
    : nodes(nodes), cache_state(0), N(nodes.size()), tolerance(tolerance)
{
	if (algorithm == FORTUNES){
		/* Do Fortune's algorithm: */
		delaunay_triangulation_sphere(nodes, delaunay_triangles,
		                              tolerance);
	} else if (algorithm == BRUTE_FORCE) {
		/* Do a brute force algorithm: */
		delaunay_triangulation_brute_force(nodes, delaunay_triangles,
		                                   tolerance);
	}
}

//----------------------------------------------------------------------
size_t VDTesselation::size() const
{
	return delaunay_triangles.size();
}


//----------------------------------------------------------------------
void VDTesselation::delaunay_triangulation(std::vector<Link>& links)
	const
{
	/* Make sure Delaunay links are cached: */
	calculate_delaunay_links();
	
	/* Copy cache: */
	links = delaunay_links;
}


//----------------------------------------------------------------------
void VDTesselation::voronoi_tesselation(std::vector<Node>& nodes,
	std::vector<Link>& links) const
{
	/* Make sure Voronoi network is cached: */
	calculate_voronoi_network();
	
	/* Copy cache: */
	nodes = voronoi_nodes;
	links = voronoi_links;
}


//----------------------------------------------------------------------
void VDTesselation::voronoi_cell_areas(std::vector<double>& areas) const
{
	/* Make sure Voronoi areas are cached: */
	calculate_voronoi_cell_areas();
	
	/* Copy cache: */
	areas = voronoi_areas;
}



/* ************************** Caching ******************************* */


	
/* Refurbished Link-Class that orders its indices.
 * All links are undirected, so we don't need to count them twice
 * in different order (i,j), (j,i).
 * Instead, order them so that i<j. */
struct link_t {
	/* Constructor: */
	link_t(size_t i, size_t j){
		if (i < j){
			this->i = i;
			this->j = j;
		} else {
			this->i = j;
			this->j = i;
		}
	}
	
	/* Comparator: */
	bool operator<(const link_t& other) const {
		if (i == other.i){
			return j < other.j;
		}
		return i < other.i;
	}
	
	/* Members: */
	size_t i, j;
};



//----------------------------------------------------------------------
void VDTesselation::calculate_delaunay_links() const
{
	/* Check if we've previously calculated the Delaunay links: */
	if (cache_state & DELAUNAY_LINKS_CACHED)
		return;
	
	/* We use a set to order the Delaunay links by their indices.
	 * This will also make sure that links are added only once (all
	 * links of the Delaunay triangulation are undirected, so we need
	 * only one of each pair (i,j) and (j,i)) */
	std::set<link_t> links;
	
	for (const Triangle& t : delaunay_triangles){
		links.insert(link_t(t.i, t.j));
		links.insert(link_t(t.i, t.k));
		links.insert(link_t(t.j, t.k));
	}

	delaunay_links.reserve(links.size());
	for (const link_t& l : links){
		delaunay_links.emplace_back(l.i, l.j);
	}
	
	/* Update cache state: */
	cache_state |= DELAUNAY_LINKS_CACHED;
	tidy_up_cache();
}



//----------------------------------------------------------------------
void VDTesselation::merge_clusters() const
{
	/* First, an injective map: */
	delaunay2voronoi.resize(delaunay_triangles.size());
	for (size_t i=0; i<delaunay_triangles.size(); ++i){
		delaunay2voronoi[i] = i;
	}

	/* Now determine clustered Voronoi nodes (that actually belong to
	 * the same node): */
}




void VDTesselation::calculate_voronoi_nodes() const
{
	/* Check if we've previously calculated the Voronoi nodes: */
	if (cache_state & VORONOI_NODES_CACHED)
		return;
	
	/* Voronoi nodes are at the circumcenter of the three nodes of the
	 * Delaunay triangles: */
	voronoi_nodes.reserve(delaunay_triangles.size());
//	std::cout << "voronoi nodes:\n[";
	for (const Triangle& t : delaunay_triangles){
		SphereVector vec = SphereVector::circumcenter(
						SphereVector(nodes[t.i].lon, nodes[t.i].lat),
						SphereVector(nodes[t.j].lon, nodes[t.j].lat),
						SphereVector(nodes[t.k].lon, nodes[t.k].lat));
		voronoi_nodes.emplace_back(vec.lon(), vec.lat());
//		std::cout << "  [" << 180.0/M_PI*vec.lon() << "," << 180.0/M_PI*vec.lat() << "],\n";
	}
//	std::cout << "]\n";

	/* Merge clusters if needed: */
	merge_clusters();
	
	/* Cache state: */
	cache_state |= VORONOI_NODES_CACHED;
	tidy_up_cache();
}


//----------------------------------------------------------------------
void VDTesselation::calculate_voronoi_network() const
{
	/* Check if we've previously calculated the Voronoi network: */
	if (cache_state & VORONOI_LINKS_CACHED)
		return;

	/* First make sure that Voronoi nodes are calculated: */
	calculate_voronoi_nodes();

	/* Init Voronoi cell area vector: */
	voronoi_areas.resize(nodes.size(), 0.0);

	/* Each Delaunay triangle corresponds to a Voronoi node. In general,
	 * there will be a constant number of Voronoi nodes per Voronoi cell,
	 * so that we use indexing by node (of the input network) to generate
	 * the links of the Delaunay triangulation in linear time.
	 *
	 * Here, we reference, by node index, all triangles that node is
	 * a vertex of.
	 */
	std::vector<std::forward_list<size_t>> node2triangle(nodes.size());
	for (size_t i=0; i<delaunay_triangles.size(); ++i){
		node2triangle[delaunay_triangles[i].i].push_front(i);
		node2triangle[delaunay_triangles[i].j].push_front(i);
		node2triangle[delaunay_triangles[i].k].push_front(i);
	}

	/* Now we iterate over each node index (of the input network) and
	 * create a closed path, the Voronoi cell: */
	std::vector<size_t> path;

	struct helper_t {
		bool   available;
		size_t id;
		Triangle t;
	};
	std::vector<helper_t> local_triangles;
	for (size_t i=0; i<nodes.size(); ++i){
		/* Create a local copy of the set of triangles belonging to node
		 * i: */
		local_triangles.resize(0);
		for (size_t t : node2triangle[i]){
			local_triangles.push_back({true, t, delaunay_triangles[t]});
		}
		node2triangle[i].clear();

		/* Find a closed path: */
		path.resize(0);
		path.push_back(local_triangles[0].id);
		Triangle last = local_triangles[0].t;
		for (size_t j=1; j<local_triangles.size(); ++j){
			bool found = false;
			for (size_t k=1; k<local_triangles.size(); ++k){
				if (!local_triangles[k].available)
					continue;
				if (last.common_border(local_triangles[k].t)){
					path.push_back(local_triangles[k].id);
					last = local_triangles[k].t;
					local_triangles[k].available = false;
					found = true;
					break;
				}
			}
			if (!found){
				std::cerr << "FOUND NO CLOSED PATH!\n";
			}
		}
		if (!delaunay_triangles[path.front()].common_border(delaunay_triangles[path.back()])){
			std::cerr << "PATH FOUND BUT NOT CLOSED! (node=" << i << ") ("
			          << 180.0/M_PI*nodes[i].lon << "," << 180.0/M_PI*nodes[i].lat << ")\n";
		}

		/* Along this closed path, we can now calculate the cell areas
		 * of the Voronoi cells using adjacent triangles.
		 * Also, we can obtain links of the Voronoi tesselation.
		 * Since each link of the Voronoi tesselation is part of the Voronoi
		 * cell of two nodes, we will end up with each link being duplicate.*/
		double area = 0.0;
		size_t last_i = path[0];
		SphereVectorEuclid last_vec(voronoi_nodes[last_i]);
		SphereVectorEuclid v_i(nodes[i]);
		for (size_t j=1; j<path.size(); ++j){
			size_t next_i = path[j];
			size_t l1 = delaunay2voronoi[next_i];
			size_t l2 = delaunay2voronoi[last_i];
			if (l1 != l2){
				if (l1 < l2){
					voronoi_links.emplace_back(l1, l2);
				} else {
					voronoi_links.emplace_back(l2, l1);
				}
				SphereVectorEuclid next(voronoi_nodes[l1]);
				area += SphereVectorEuclid::triangle_area(last_vec, next,
				                                          v_i);
				last_vec = next;
			}
			last_i = next_i;
		}
		if (delaunay2voronoi[last_i] != delaunay2voronoi[path[0]]){
			size_t l1 = delaunay2voronoi[path[0]];
			size_t l2 = delaunay2voronoi[last_i];
			if (l1 < l2){
				voronoi_links.emplace_back(l1, l2);
			} else {
				voronoi_links.emplace_back(l2, l1);
			}
			area += SphereVectorEuclid::triangle_area(last_vec, v_i,
			                    SphereVectorEuclid(voronoi_nodes[l1]));

		}

		voronoi_areas[i] = area;
	}

	/* Now we have successfully created the Voronoi tesselation and
	 * calculated the Voronoi weights. Clean up:
	 */
	node2triangle.clear();
	local_triangles.clear();
	path.clear();

	/* Finally, we have to sort the Voronoi links... */
	auto cmp = [&](const Link& l, const Link& r)->bool
	    {
		    if (l.i != r.i){
				return l.i < r.i;
			}
			return l.j < r.j;
	    };

	std::sort<std::vector<Link>::iterator, decltype(cmp)>
	        (voronoi_links.begin(), voronoi_links.end(), cmp);


	/* ... and remove duplicates! */
	for (size_t i=1; i<voronoi_links.size()/2; ++i){
		if (voronoi_links[2*i] != voronoi_links[2*i+1]){
			std::cerr << "ERROR : NO DUPLICATES! (2*i=" << 2*i
			          << " -> (" << voronoi_links[2*i].i << ","
			          << voronoi_links[2*i].j << "), ("
			          << voronoi_links[2*i+1].i << ","
			          << voronoi_links[2*i+1].j << ")\n\n[";
			for (size_t i=0; i<voronoi_links.size(); ++i){
				std::cerr << "(" << voronoi_links[i].i << ","
				          << voronoi_links[i].j << "), ";
			}
			std::cerr << "]\n";
			std::cerr << "TRIANGULATION:\n[\n";
			for (size_t i=0; i<delaunay_triangles.size(); ++i){
				Triangle t = delaunay_triangles[i];
				std::cout << "   " << i << ": [" << t.i << "," << t.j << "," << t.k << "] "
				          << " (" << voronoi_nodes[i].lon << "," << voronoi_nodes[i].lat << "),\n";
			}
			std::cout << "]\n";
			throw 0;
		}
		voronoi_links[i] = voronoi_links[2*i];
	}


	/* Cache state: */
	cache_state |= (VORONOI_LINKS_CACHED | VORONOI_CELLS_CACHED);
	tidy_up_cache();
}



//void VDTesselation::calculate_voronoi_network() const
//{
//	/* Check if we've previously calculated the Voronoi network: */
//	if (cache_state & VORONOI_LINKS_CACHED)
//		return;
	
//	/* First make sure that Voronoi nodes are calculated: */
//	calculate_voronoi_nodes();
	
//	/* Temporary class for saving voronoi edges as values: */
//	struct vedge_t {
//		/* Constructor: */
//		vedge_t(){
//			N = 0;
//			ids[0] = -1;
//			ids[1] = -1;
//		}
		
//		/* Append one id: */
//		void append(size_t id){
//			if (N == 2){
//				std::cerr << "Voronoi edges connect exactly two "
//				             "Delaunay triangles!\n"
//				             "   id[0] = " << ids[0] << "\n   id[1] = "
//				          << ids[1] << "\n";
//				throw 0;
//			}
//			ids[N++] = id;
//		}
		
		
//		unsigned char N;
//		size_t        ids[2];
		
//	};
	
//	/* In a map, we save triangle references by their links: */
//	std::map<link_t, vedge_t> references;
	
//	for (size_t i=0; i<delaunay_triangles.size(); ++i){
//		references[link_t(delaunay_triangles[i].i,
//		                  delaunay_triangles[i].j)].append(i);
//		references[link_t(delaunay_triangles[i].i,
//		                  delaunay_triangles[i].k)].append(i);
//		references[link_t(delaunay_triangles[i].j,
//		                  delaunay_triangles[i].k)].append(i);
//	}
	
//	/* Extract links from map: */
//	std::set<link_t> links;
//	for (const std::pair<link_t,vedge_t>& p : references){
//		if (p.second.N != 2){
//			std::cerr << "Voronoi edges connect exactly two "
//			             "Delaunay triangles!\n"
//			             "   N = " << (int)p.second.N << "\n";
//			throw 0;
//		}
//		links.insert(link_t(p.second.ids[0], p.second.ids[1]));
//	}
//	references.clear();

//	/* Insert links to return array: */
//	voronoi_links.reserve(links.size());
//	for (auto it = links.begin(); it != links.end(); ++it){
//		voronoi_links.emplace_back(it->i, it->j);
//	}
	
////	/* Merge clusters: */
////	merge_clusters(voronoi_nodes, voronoi_links);

//	/* Cache state: */
//	cache_state |= VORONOI_LINKS_CACHED;
//	tidy_up_cache();
	
//}

//----------------------------------------------------------------------
void VDTesselation::calculate_voronoi_cell_areas() const
{
	/* Check if we've previously calculated the Voronoi cells: */
	if (cache_state & VORONOI_CELLS_CACHED)
		return;


	/* Calculation of Voronoi links makes calculating the weight easy: */
	calculate_voronoi_network();
	
//	/* We need the Voronoi nodes: */
//	calculate_voronoi_nodes();
	
//	/* For each node, we need all surrounding Voronoi nodes.
//	 * Each Voronoi node corresponds to a Delaunay triangle, so we
//	 * first add references to all Delaunay triangles a node is part
//	 * of: */
//	std::vector<std::forward_list<size_t>> node2triangle(nodes.size());
//	for (size_t i=0; i<delaunay_triangles.size(); ++i){
//		node2triangle[delaunay_triangles[i].i].push_front(i);
//		node2triangle[delaunay_triangles[i].j].push_front(i);
//		node2triangle[delaunay_triangles[i].k].push_front(i);
//	}
	
//	/* Init destination vector: */
//	voronoi_areas.resize(nodes.size(), 0.0);
	
//	/* Now we iterate over all nodes and calculate the Voronoi cell
//	 * areas: */
//	std::vector<size_t> path;
	
//	struct helper_t {
//		bool   available;
//		size_t id;
//		Triangle t;
//	};
//	std::vector<helper_t> local_triangles;
//	for (size_t i=0; i<nodes.size(); ++i){
//		/* First, calculate a closed path to traverse the edge of the
//		 * Voronoi cell.
//		 * Create a local copy of the set of triangles belonging to node
//		 * i: */
//		local_triangles.resize(0);
//		for (size_t t : node2triangle[i]){
//			local_triangles.push_back({true, t, delaunay_triangles[t]});
//		}
//		node2triangle[i].clear();
		
//		/* Find a closed path: */
//		path.resize(0);
//		path.push_back(local_triangles[0].id);
//		Triangle last = local_triangles[0].t;
//		for (size_t j=1; j<local_triangles.size(); ++j){
//			bool found = false;
//			for (size_t k=1; k<local_triangles.size(); ++k){
//				if (!local_triangles[k].available)
//					continue;
//				if (last.common_border(local_triangles[k].t)){
//					path.push_back(local_triangles[k].id);
//					last = local_triangles[k].t;
//					local_triangles[k].available = false;
//					found = true;
//					break;
//				}
//			}
//		}
		
//		/* Along this closed path, we can now calculate the cell areas
//		 * of the Voronoi cells using adjacent triangles: */
//		double area = 0.0;
//		SphereVectorEuclid last_vec(voronoi_nodes[path[0]]);
//		SphereVectorEuclid v_i(nodes[i]);
//		for (size_t j=1; j<path.size(); ++j){
//			SphereVectorEuclid next(voronoi_nodes[path[j]]);
//			area += SphereVectorEuclid::triangle_area(last_vec, next,
//			                                          v_i);
//			last_vec = next;
//		}
//		area += SphereVectorEuclid::triangle_area(last_vec, v_i,
//		                    SphereVectorEuclid(voronoi_nodes[path[0]]));
//		voronoi_areas[i] = area;
//	}
	
//	/* Cache state: */
//	cache_state |= VORONOI_CELLS_CACHED;
//	tidy_up_cache();
}


//----------------------------------------------------------------------
void VDTesselation::associated_nodes(
	const std::vector<size_t>& voronoi_nodes,
	std::vector<size_t>& associated) const
{
	/* We iterate over all Voronoi nodes in the given set and mark
	 * the associated nodes of the original network.
	 * In the end, we collect all marked nodes. */
	std::vector<bool> marked(N, false);
	
	for (size_t node : voronoi_nodes){
		marked[delaunay_triangles[node].i] = true;
		marked[delaunay_triangles[node].j] = true;
		marked[delaunay_triangles[node].k] = true;
	}
	
	size_t count = 0;
	for (bool b : marked){
		if (b)
			++count;
	}
	
	associated.resize(count);
	size_t j=0;
	for (size_t i=0; i<N; ++i){
		if (marked[i])
			associated[j++] = i;
	}
}

//----------------------------------------------------------------------
void VDTesselation::tidy_up_cache() const
{
	/* Check if nodes still need to be stored: */
	if ((cache_state & VORONOI_NODES_CACHED) &&
	    (cache_state & VORONOI_CELLS_CACHED))
	{
		nodes.clear();
	}
}

//----------------------------------------------------------------------
void VDTesselation::print_debug(bool sort_triangles) const
{
	/* First print Delaunay tesselation: */
	std::cout << "--- VDTesselation debug output ---\n\nDelaunay tesselation:";
	size_t column = 0;

	/* Order triangles by smallest index first, then by j: */
	const size_t M = delaunay_triangles.size();
	std::vector<Triangle> triangles_copy(M);
	for (size_t i=0; i<M; ++i){
		Triangle t = delaunay_triangles[i];
		if (t.k < t.i && t.k < t.j){
			triangles_copy[i] = Triangle(t.k, t.i, t.j);
		} else if (t.j < t.i && t.j < t.k){
			triangles_copy[i] = Triangle(t.j, t.k, t.i);
		} else {
			triangles_copy[i] = t;
		}
	}
	std::vector<size_t> index_map(M);
	for (size_t i=0; i<M; ++i){
		index_map[i] = i;
	}
	if (sort_triangles){
		auto cmp = [&](size_t i, size_t j){
			    if (triangles_copy[i].i < triangles_copy[j].i){
					return true;
				} else if (triangles_copy[i].i > triangles_copy[j].i){
					return false;
				}
				if (triangles_copy[i].j < triangles_copy[j].j){
					return true;
				} else if (triangles_copy[i].j > triangles_copy[j].j){
					return false;
				}
				return triangles_copy[i].k < triangles_copy[j].k;
		    };

		std::sort<std::vector<size_t>::iterator,decltype(cmp)>(
		    index_map.begin(), index_map.end(), cmp);
	}



	std::cout << "[";
	for (size_t i=0; i<M; ++i){
		if (column == 0){
			std::cout << "\n  ";
		}
		Triangle t = triangles_copy[index_map[i]];
		std::cout << " [" << t.i << "," << t.j << "," << t.k << "]";
		if (i != M-1){
			std::cout << ",";
		}
		column = (column+1) % 1;
	}
	std::cout << "]\n";
	if (cache_state & VORONOI_NODES_CACHED){
		std::cout << "\nVoronoi nodes (in °):\n[";
		column = 0;
		for (size_t i=0; i<M; ++i){
			if (column == 0){
				std::cout << "\n  ";
			}

			std::cout << " [" << 180.0/M_PI*voronoi_nodes[index_map[i]].lon << ","
			          << 180.0/M_PI*voronoi_nodes[index_map[i]].lat << "]";
			if (i != M-1){
				std::cout << ",";
			}
			column = (column+1) % 3;
		}
		std::cout << "]\n";
	}
}


} // NAMESPACE ACOSA
