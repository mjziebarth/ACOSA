/* A fast algorithm for the determination of links of a geometric
 * graph. Part of ACOSA.
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
 * [1] Xiaoyu Zheng et al.: A Plane Sweep Algorithm for the Voronoi
 *     Tesselation of the Sphere, in: electronic-Liquid Crystal
 *     Communications, 2011-12-13
 *     http://www.e-lc.org/docs/2011_12_05_14_35_11
 */


#include <geometricgraph.hpp>

#include <queue>
#include <set>
#include <cmath>
#include <algorithm>


namespace ACOSA {


struct node_t {
	size_t id;
	double lon;
	double lat;
};




/*************** REFACTOR THIS *************************************/

/* TODO : Remove this section and rewrite the other code, so that
 * functionality from spherics.cpp is used! */
static double greatcircle_distance(double sin_lat_1, double cos_lat_1,
    double sin_lon_1, double cos_lon_1, double sin_lat_2,
    double cos_lat_2, double sin_lon_2, double cos_lon_2)
{
	/* Dot product of vectors i and j: */
	double d = sin_lat_1*sin_lat_2 + cos_lat_1*cos_lat_2*(
	    sin_lon_1*sin_lon_2 + cos_lon_1*cos_lon_2);

	/* Arccos to obtain distance: */
	if (d >= 1.0){
		d = 0.0;
	} else if (d <= -1.0){
		d = M_PI;
	} else {
		d = acos(d);
	}

	return d;
}

/* Saving the trigonometric function values of a (lon,lat) pair: */
struct coord_trigf_t{
	double sin_lon;
	double cos_lon;
	double sin_lat;
	double cos_lat;
};


/*************** REFACTOR THIS *************************************/

static void geometric_graph_links_pairwise(
    const std::vector<node_t>::const_iterator& begin,
    const std::vector<node_t>::const_iterator& end,
    std::vector<Link>& links,
    double sigma_0)
{
	for (auto it=begin; it!=end; ++it){
		double slon_1 = std::sin(it->lon);
		double clon_1 = std::cos(it->lon);
		double slat_1 = std::sin(it->lat);
		double clat_1 = std::cos(it->lat);
		for (auto it2=it+1; it2 != end; ++it2){
			if (greatcircle_distance(slat_1, clat_1, slon_1, clon_1,
			        std::sin(it2->lat), std::cos(it2->lat),
			        std::sin(it2->lon), std::cos(it2->lon))
			    < sigma_0)
			{
				links.push_back({it->id, it2->id});
				links.push_back({it2->id, it->id});
			}
		}
	}

}



void geometric_graph_links(
    const std::vector<Node>& coordinates,
    std::vector<Link>& links,
    double sigma_0)
{
	const size_t N = coordinates.size();

	/* If N is small or sigma_0 is big, there is no advantage of the
	 * more complicated algorithm, so we can take a simple pairwise
	 * check.
	 * These parameters should most likely be subject to
	 * machine-dependent tweaking.
	 * Note that sigma_0 < PI/2 is expected lateron (it is also required
	 * to limit the longitude interval in which to search). */
	constexpr double sigma_limit = std::min(M_PI_2, 17.0/180.0 * M_PI);
	if (N <= 100 || sigma_0 > sigma_limit){
		std::vector<coord_trigf_t> vec(N);
		for (size_t i=0; i<N; ++i){
			vec[i].sin_lon = std::sin(coordinates[i].lon);
			vec[i].cos_lon = std::cos(coordinates[i].lon);
			vec[i].sin_lat = std::sin(coordinates[i].lat);
			vec[i].cos_lat = std::cos(coordinates[i].lat);
		}

		for (size_t i=0; i<N; ++i){
			for (size_t j=i+1; j<N; ++j){
				if (greatcircle_distance(vec[i].sin_lat, vec[i].cos_lat,
				        vec[i].sin_lon, vec[i].cos_lon,
				        vec[j].sin_lat, vec[j].cos_lat, vec[j].sin_lon,
				        vec[j].cos_lon)
				    < sigma_0)
				{
					links.push_back({i, j});
					links.push_back({j, i});
				}
			}
		}

		return;
	}


	/* We order the nodes by increasing latitude, so we need to keep
	 * track of their indices: */

	auto cmp_lat = [](const node_t& l, const node_t& r)->bool
	    {
		    return l.lat < r.lat;
	    };

	std::vector<node_t> nodes(N);

	for (size_t i=0; i<N; ++i){
		nodes[i] = {i, coordinates[i].lon, coordinates[i].lat};
	}

	std::sort<std::vector<node_t>::iterator,decltype(cmp_lat)>(
	    nodes.begin(), nodes.end(), cmp_lat);


	/* The belt: */
	auto cmp_lon = [](const node_t* l, const node_t* r)->bool
	    {
		    if (l->lon == r->lon){
				return l->id < r->id;
			}
			return l->lon < r->lon;
	    };

	std::set<node_t*, decltype(cmp_lon)> belt(cmp_lon);

	std::queue<decltype(belt)::iterator> leave_events;


	/* Step 1: Interconnect all nodes inside the southpole's sigma_0
	 *         circle. Also add them to belt for step 2. */
	auto enter = nodes.begin();
	for (;enter != nodes.end(); ++enter){
		if (enter->lat > sigma_0-M_PI_2){
			break;
		}
		leave_events.push(belt.insert(&*enter).first);
	}

	geometric_graph_links_pairwise(nodes.begin(),
	    enter, links, sigma_0);


	/* Step 2: Use belt for most of the other nodes: */
	std::vector<node_t>::iterator arctic_begin = nodes.end();
	bool arctic_entered = false;
	const double sin_s0 = std::sin(sigma_0);
	while (enter != nodes.end()){
		if (leave_events.empty() ||
		    enter->lat <= (*leave_events.front())->lat+sigma_0)
		{
			/* Connect to nodes in belt.
			 * First we determine the longitude bounds in which we have
			 * to search the belt: */

			/* Half of the symmetric maximum longitude interval in which
			 * we have to search can be calculated by constructing a
			 * spherical triangle that touches the sigma_0 circle around
			 * 'enter' in one of its points, where it thus a right angle.
			 * The other two points of the triangle are given by the
			 * southpole (the corresponding angle being delta_lon) and
			 * by enter's position. Two sides are then known: enter's
			 * latitude and the circle diameter sigma_0.
			 * Using spherical rectangular triangle formula, we arrive
			 * at the following equation: */
			double delta_lon = std::asin(sin_s0 / std::cos(enter->lat));
			double lon_left  = enter->lon - delta_lon;
			double lon_right = enter->lon + delta_lon;



			/* Because sigma_0 < pi/2 and all events in this iteration
			 * have latitude bigger than sigma_0, we can be sure that
			 * the longitude bounds do not overlap. */



			/* Find insertion position: */
			auto insert_pos = belt.lower_bound(&*enter);
			auto left_bound = insert_pos;

			/* Search the belt for all nodes we have to link the new
			 * node to: */
			if (!belt.empty()){
				/* Now iterate to left and right from insertion position,
				 * checking all nodes inside the longitude bounds: */
				double slon_1 = std::sin(enter->lon);
				double clon_1 = std::cos(enter->lon);
				double slat_1 = std::sin(enter->lat);
				double clat_1 = std::cos(enter->lat);

				/* Iterate to left: */
				auto it = insert_pos;
				if (lon_left < 0.0){
					/* The left longitude border wraps around the
					 * longitude periodic border. We can thus iterate
					 * over all nodes up to that border, wrap around and
					 * then continue as usual: */
					while (it != belt.begin()){
						--it;
						if (greatcircle_distance(slat_1, clat_1, slon_1,
						        clon_1,
						        std::sin((*it)->lat),
						        std::cos((*it)->lat),
						        std::sin((*it)->lon),
						        std::cos((*it)->lon))
						    < sigma_0)
						{
							links.push_back({(*it)->id, enter->id});
							links.push_back({enter->id, (*it)->id});
						}
					}
					/* We're now at belt.begin() and have checked it
					 * already.
					 * Wrap iterator around border. Since iterator will
					 * be decremented before being dereferenced, we have
					 * to set it to end: */
					it = belt.end();
					lon_left += M_2_PI;
				}
				/* Iterate until we reach the longitude border: */
				while (it != belt.begin() && (*(--it))->lon >= lon_left)
				{
					if (greatcircle_distance(slat_1, clat_1, slon_1,
					        clon_1,
					        std::sin((*it)->lat), std::cos((*it)->lat),
					        std::sin((*it)->lon), std::cos((*it)->lon))
					    < sigma_0)
					{
						links.push_back({(*it)->id, enter->id});
						links.push_back({enter->id, (*it)->id});
					}
				}


				/* Iterate right.
				 * In this iteration, we also handle insert_pos. */
				it = insert_pos;

				if (lon_right > M_2_PI){
					/* The right longitude border wraps around the
					 * longitude periodic border. We can thus iterate
					 * over all nodes up to that border, wrap around and
					 * then continue as usual: */
					while (it != belt.end()){
						if (greatcircle_distance(slat_1, clat_1, slon_1,
						        clon_1,
						        std::sin((*it)->lat),
						        std::cos((*it)->lat),
						        std::sin((*it)->lon),
						        std::cos((*it)->lon))
						    < sigma_0)
						{
							links.push_back({(*it)->id, enter->id});
							links.push_back({enter->id, (*it)->id});
						}
						++it;
					}
					it = belt.begin();
					lon_right -= M_2_PI;

				}
				/* Iterate until we reach the longitude border: */
				while (it != belt.end() && (*it)->lon <= lon_right){
					if (greatcircle_distance(slat_1, clat_1, slon_1,
					        clon_1,
					        std::sin((*it)->lat), std::cos((*it)->lat),
					        std::sin((*it)->lon), std::cos((*it)->lon))
					    < sigma_0)
					{
						links.push_back({(*it)->id, enter->id});
						links.push_back({enter->id, (*it)->id});
					}
					++it;
				}
			}

			/* Insert node into belt, if it is not part of the arctic: */
			if (enter->lat < M_PI_2 - sigma_0){
				/* Inset to belt: */
				leave_events.push(belt.insert(&*enter).first);
			} else {
				/* We entered the arctic. Cool stuff!
				 * Now for science and stuff, we mark the arctic's
				 * border: */
				arctic_begin = enter;
				break;
			}

			++enter;
		} else {
			/* Leave event. Simply remove node from belt: */
			belt.erase(leave_events.front());
			leave_events.pop();
		}
	}


	/* Step 3: Connect the remaining nodes of the belt to all its
	 *         remaining neighbours: */
	for (const node_t* node : belt){
		double slon_1 = std::sin(node->lon);
		double clon_1 = std::cos(node->lon);
		double slat_1 = std::sin(node->lat);
		double clat_1 = std::cos(node->lat);

		for (auto it = arctic_begin; it != nodes.end(); ++it){
			if (greatcircle_distance(slat_1, clat_1, slon_1, clon_1,
			        std::sin(it->lat), std::cos(it->lat),
			        std::sin(it->lon), std::cos(it->lon))
			    < sigma_0)
			{
				links.push_back({it->id, node->id});
				links.push_back({node->id, it->id});
			}
		}
	}


	/* Step 4: Connect all nodes inside the arctic (the north pole's
	 *         sigma_0 circle): */
	geometric_graph_links_pairwise(arctic_begin,
	    nodes.end(), links, sigma_0);

	/* Et voila, we're done! */
}

} // NAMESPACE ACOSA
