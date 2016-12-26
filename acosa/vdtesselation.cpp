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

#include <map>
#include <set>
#include <forward_list>
#include <iostream>

namespace ACOSA {

/* Cache states: */
static constexpr unsigned char DELAUNAY_LINKS_CACHED = 1;
static constexpr unsigned char VORONOI_NODES_CACHED =  2;
static constexpr unsigned char VORONOI_LINKS_CACHED =  4;
static constexpr unsigned char VORONOI_CELLS_CACHED =  8;



//----------------------------------------------------------------------
VDTesselation::VDTesselation(const std::vector<Node>& nodes) 
	: nodes(nodes), cache_state(0), N(nodes.size())
{
	/* Do Fortune's algorithm: */
	delaunay_triangulation_sphere(nodes, delaunay_triangles);
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
	
	/* Update cache state: */
	cache_state |= DELAUNAY_LINKS_CACHED;
	tidy_up_cache();
}



//----------------------------------------------------------------------
void VDTesselation::calculate_voronoi_nodes() const
{
	/* Check if we've previously calculated the Voronoi nodes: */
	if (cache_state & VORONOI_NODES_CACHED)
		return;
	
	/* Voronoi nodes are at the circumcenter of the three nodes of the
	 * Delaunay triangles: */
	voronoi_nodes.reserve(delaunay_triangles.size());
	for (const Triangle& t : delaunay_triangles){
		SphereVector vec = SphereVector::circumcenter(
						SphereVector(nodes[t.i].lon, nodes[t.i].lat),
						SphereVector(nodes[t.j].lon, nodes[t.j].lat),
						SphereVector(nodes[t.k].lon, nodes[t.k].lat));
		voronoi_nodes.emplace_back(vec.lon(), vec.lat());
	}
	
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
	
	/* Temporary class for saving voronoi edges as values: */
	struct vedge_t {
		/* Constructor: */
		vedge_t(){
			N = 0;
			ids[0] = -1;
			ids[1] = -1;
		}
		
		/* Append one id: */
		void append(size_t id){
			if (N == 2){
				std::cerr << "Voronoi edges connect exactly two "
						     "Delaunay triangles!\n"
						     "   id[0] = " << ids[0] << "\n   id[1] = "
						  << ids[1] << "\n";
				throw 0;
			}
			ids[N++] = id;
		}
		
		
		unsigned char N;
		size_t        ids[2];
		
	};
	
	/* In a map, we save triangle references by their links: */
	std::map<link_t, vedge_t> references;
	
	for (size_t i=0; i<delaunay_triangles.size(); ++i){
		references[link_t(delaunay_triangles[i].i,
						  delaunay_triangles[i].j)].append(i);
		references[link_t(delaunay_triangles[i].i, 
						  delaunay_triangles[i].k)].append(i);
		references[link_t(delaunay_triangles[i].j,
						  delaunay_triangles[i].k)].append(i);
	}
	
	/* Extract links from map: */
	std::set<link_t> links;
	for (const std::pair<link_t,vedge_t>& p : references){
		links.insert(link_t(p.second.ids[0], p.second.ids[1]));
	}
	references.clear();
	
	/* Insert links to return array: */
	voronoi_links.reserve(links.size());
	for (auto it = links.begin(); it != links.end(); ++it){
		voronoi_links.emplace_back(it->i, it->j);
	}
	
	/* Cache state: */
	cache_state |= VORONOI_LINKS_CACHED;
	tidy_up_cache();
	
}

//----------------------------------------------------------------------
void VDTesselation::calculate_voronoi_cell_areas() const
{
	/* Check if we've previously calculated the Voronoi cells: */
	if (cache_state & VORONOI_CELLS_CACHED)
		return;
	
	/* We need the Voronoi nodes: */
	calculate_voronoi_nodes();
	
	/* For each Voronoi node, we need all surrounding Voronoi edges.
	 * Each Voronoi edge corresponds to a Delaunay triangle, so we
	 * first add references to all Delaunay triangles a node is part
	 * of: */
	std::vector<std::forward_list<size_t>> node2triangle(nodes.size());
	for (size_t i=0; i<delaunay_triangles.size(); ++i){
		node2triangle[delaunay_triangles[i].i].push_front(i);
		node2triangle[delaunay_triangles[i].j].push_front(i);
		node2triangle[delaunay_triangles[i].k].push_front(i);
	}
	
	/* Init destination vector: */
	voronoi_areas.resize(nodes.size(), 0.0);
	
	/* Now we iterate over all nodes and calculate the Voronoi cell
	 * areas: */
	std::vector<size_t> path;
	
	struct helper_t {
		bool   available;
		size_t id;
		Triangle t;
	};
	std::vector<helper_t> local_triangles;
	std::vector<bool> local_triangles_available;
	for (size_t i=0; i<nodes.size(); ++i){
		/* First, calculate a closed path to traverse the edge of the
		 * Voronoi cell.
		 * Create a local copy of the set of triangles belonging to node
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
		}
		
		/* Along this closed path, we can now calculate the cell areas
		 * of the Voronoi cells using adjacent triangles: */
		double area = 0.0;
		SphereVectorEuclid last_vec(voronoi_nodes[path[0]]);
		SphereVectorEuclid v_i(nodes[i]);
		for (size_t j=1; j<path.size(); ++j){
			SphereVectorEuclid next(voronoi_nodes[path[j]]);
			area += SphereVectorEuclid::triangle_area(last_vec, next,
			                                          v_i);
			last_vec = next;
		}
		area += SphereVectorEuclid::triangle_area(last_vec, v_i,
							SphereVectorEuclid(voronoi_nodes[path[0]]));
		voronoi_areas[i] = area;
	}
	
	/* Cache state: */
	cache_state |= VORONOI_CELLS_CACHED;
	tidy_up_cache();
}


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


} // NAMESPACE ACOSA
