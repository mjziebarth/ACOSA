/* Implementation of the spherical Fortune's algorithm of [1]. Part of
 * ACOSA.
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

#include <fortunes_sphere.hpp>
#include <spherics.hpp>
#include <beach.hpp>
#include <circleevent.hpp>


#include <queue>
#include <memory>
#include <algorithm>

#include <iostream>

#include <limits>


namespace ACOSA {


//######################################################################

struct site_event_t {
	site_event_t(double lon, double lat, size_t id)
	    : vec(lon, lat), id(id)
	{}

	bool operator<(const site_event_t& other) const {
		if (vec.lat() == other.vec.lat()){
			return vec.lon() < other.vec.lon();
		}
		return vec.lat() > other.vec.lat();
	}

	SphereVector vec;
	size_t id;
};



//----------------------------------------------------------------------
static std::shared_ptr<bool> add_circle_event(
    eventqueue_t& queue, const SphereVector& v1,
    const SphereVector& v2, const SphereVector& v3, double tide,
    const BeachIterator& beach_site, double tolerance)
{
	CircleEvent event(v1, v2, v3, beach_site, tolerance);
	if (event.lat()+tolerance >= tide){
		event.set_valid();
		queue.push(event);
		return event.valid_;
	}
	return nullptr;
}

//----------------------------------------------------------------------
template<typename site_queue_t>
static inline Beach init_beach(site_queue_t& site_events,
                            eventqueue_t& circle_events,
                            std::vector<Triangle>& delaunay_triangles)
{
	/* Step 1): Obtain the first two site events: */
	site_event_t e1 = site_events.top();
	site_events.pop();
	site_event_t e2 = site_events.top();
	site_events.pop();

	/* Step 2): Continue depending on whether we have a degeneracy
	 *          where the first three sites have the same latitude
	 *          (see [1])
	 */
	if (e1.vec.lat() == e2.vec.lat() &&
	    site_events.top().vec.lat() == e1.vec.lat())
	{
		/* Degeneracy.
		 * Collect all sites with the same latitude: */
		std::vector<site_event_t> sites;
		sites.push_back(e1);
		sites.push_back(e2);

		while (site_events.top().vec.lat() == e1.vec.lat()){
			sites.push_back(site_events.top());
			site_events.pop();
		}

		/* Sort them by longitude: */
		auto cmp = [](const site_event_t& l, const site_event_t& r)
		{
			return l.vec.lon() < r.vec.lon();
		};

		std::sort<std::vector<site_event_t>::iterator,
		          decltype(cmp)>
		    (sites.begin(), sites.end(), cmp);


		/* Create one of the possible Delaunay triangulations
		 * of the area inside the circumcircle: */
		for (size_t i=2; i<sites.size(); ++i){
			delaunay_triangles.emplace_back(sites[i].id, sites[i-1].id,
			    sites[0].id);
		}

		/* Fill vectors for beach constructor: */
		std::vector<SphereVector> vecs(sites.size());
		std::vector<size_t> ids(sites.size());

		for (size_t i=0; i<sites.size(); ++i){
			vecs[i] = sites[i].vec;
			ids[i] = sites[i].id;
		}

		/* Create beach: */
		return Beach(ids, vecs, circle_events);

	} else {
		/* No degeneracy: */
		return Beach(e1.id, e1.vec, e2.id, e2.vec);
	}
}


//----------------------------------------------------------------------
void delaunay_triangulation_sphere(const std::vector<Node>& nodes,
    std::vector<Triangle>& delaunay_triangles, double tolerance)
{
	const size_t N = nodes.size();
	
	/* 1) Priority queue of site events. */
	/*    Site event comparison lambda: */
	auto cmp_site = 
	    [](const site_event_t& l, const site_event_t& r)->bool
		{
			if (l.vec.lat() == r.vec.lat()){
				return l.vec.lon() < r.vec.lon();
			}
			return l.vec.lat() > r.vec.lat();
		};
	/*    1.3: The queue: */
	std::priority_queue<site_event_t, std::vector<site_event_t>,
	                    decltype(cmp_site)>
	    site_events(cmp_site);
	
	/*    1.4: Fill the queue: */
	for (size_t i=0; i<N; ++i){
		site_events.emplace(nodes[i].lon, nodes[i].lat, i);
	}
	
	
	/* 2) Priority queue of circle events: */
	eventqueue_t circle_events;
	
	
	/* 3) Beach line (ordered): */
	Beach beach = init_beach(site_events, circle_events,
	                         delaunay_triangles);
	
	
	
	/* Now we're ready for the algorithm! */
	while (!site_events.empty() || !circle_events.empty())
	{
		/* Decide whether the next event is a circle event or a site
		 * event: */
		if (!site_events.empty() && (circle_events.empty()
		     || site_events.top().vec.lat() <
			    circle_events.top().lat()))
		{
			/* Site event comes first: */
			SphereVector vec = site_events.top().vec;
			size_t       id  = site_events.top().id;
			site_events.pop();
			double tide = vec.lat();
			
			/* Find iterator where to insert the site event into the
			 * beach.
			 * This returns the site which fulfills
			 *    site->lon() <= it->lon_left()
			 * In terms of [1] this means that if:
			 *    site->lon() <  it->lon_left()  :  it == p_3
			 *    site->lon() == it->lon_left()  :  it == p_j
			 */
			BeachIterator it = 
			    beach.find_insert_position(vec.lon(), vec.lat(),
			                               tolerance);
			
			
			/* Check for degenerate case 
			 * 			vec.lon() == it->lon_left()
			 */
			if (it.lon_left_equal(vec.lon(), tolerance))
			{
				/* Degenerate case. The closest-point-line originating
				 * from the inserted node (see [1]) meets the arc
				 * intersect of 'it' and its left neighbour. Thus, we
				 * have instantly found a voronoi node.
				 * 
				 * This is conceptually equivalent to
				 * splitting the arc 'it' and instantly removing the
				 * split-off arc (of length 0).
				 * Thus, we have to perform a combination of a site and
				 * circle event. */

				/* Obtain beach site attributes.
				 * Since the site is inserted between two arcs, we can
				 * sort all relevant arcs into left and right arcs: */
				SphereVector v_r1 = it.vec();
				size_t       i_r1 = it.id();
				++it; // p_r2
				SphereVector v_r2 = it.vec();
				size_t       i_r2 = it.id();
				--it; // p_r1
				--it; // p_l1
				SphereVector v_l1 = it.vec();
				size_t       i_l1 = it.id();
				--it; // p_l2
				SphereVector v_l2 = it.vec();
				size_t       i_l2 = it.id();
				++it; // p_l1

				/* Invalidate the current beach event for node p_l1: */
				beach.invalidate_circle_event(it);

				/* Invalidate the current beach event for node p_r1: */
				++it; // p_r1
				beach.invalidate_circle_event(it);

				/* Emplace the new beach site. This
				 * automatically updates the beach site before it: */
				it = beach.insert_before(it, id, vec);

				/* Check for circle event [p_l2, p_l1, p_i]:  */
				--it; // p_l1
				auto ce = add_circle_event(circle_events, v_l2, v_l1, vec,
				                           tide, it, tolerance);
				if (ce){
					beach.register_circle_event(it, ce);
				}

				/* Check for circle event [p_i, p_r1, p_r2]: */
				++it; // p_i
				++it; // p_r1
				ce = add_circle_event(circle_events, vec, v_r1, v_r2,
				                      tide, it, tolerance);
				if (ce){
					beach.register_circle_event(it, ce);
				}

				/* Finally, create Delaunay triangle: */
				delaunay_triangles.emplace_back(id, i_l1, i_r1);

			} else {
				/* Non-degenerate case. We split a beach site and
				 * obtain two finite-length arcs. */
				
				/* Obtain beach site attributes. Notation follows [1].
				 * After the search, we have it == p_3: */
				SphereVector v_3 = it.vec();
				size_t       i_3 = it.id();
				--it; // p_j
				SphereVector v_j = it.vec();
				size_t       i_j = it.id();
				--it; // p_2
				SphereVector v_2 = it.vec();
				size_t       i_2 = it.id();
				++it; // p_j          
				
				/* Invalidate the current beach event for node p_j: */
				beach.invalidate_circle_event(it);
				
				/* Now we can emplace the new beach sites. This
				 * automatically updates the beach site before it,
				 * so we do not need to care about the current
				 * beach site for p_j here: */
				it = beach.insert_before(it, id, vec);
				it = beach.insert_before(it, i_j, v_j);
				     
				/* Update the circle event for the split-off node v_j
				 * if it exists ( [p_2, p_j, pi] ). */
				auto ce = add_circle_event(circle_events, v_2, v_j, vec,
				                           tide, it, tolerance);
				if (ce){
					beach.register_circle_event(it, ce);
				}
				
				/* Update the circle event for the current node v_j
				 * if it exists ( [p_i, p_j, p3] ): */
				++it;
				++it;
				ce = add_circle_event(circle_events, vec, v_j, v_3,
				                           tide, it, tolerance);
				if (ce){
					beach.register_circle_event(it, ce);
				}
			}
			
			
			
		} else {
			/* Circle event comes first: */
			CircleEvent event = circle_events.top();
			circle_events.pop();
			if (!event.valid()){
				continue;
			}
			
			/* Valid event. Remove it from beach.
			 * This will also invalidate circle events of the
			 * neighbouring nodes that contain it and update 
			 * its right neighbour's left-vector: */
			BeachIterator it = event.beach_site;
			size_t        id = it.id();
			double      tide = event.lat();

			beach.set_tide(it.tide());
			it = beach.erase(it);
			
			/* Left and right neighbours of removed node: */
			ArcIntersect r = it->first;
			--it;
			ArcIntersect l = it->first;
			
			if (beach.size() > 2){
				/* Check left circle event: */
				auto ce = add_circle_event(circle_events, l.left(), 
				                           l.vec(),
				                           r.vec(), tide, it,
				                           tolerance);
				if (ce){
					beach.register_circle_event(it, ce);
				}
				
				/* Check right circle event: */
				++it;
				++it;
				ArcIntersect rr = it->first;
				--it;
				ce = add_circle_event(circle_events, r.left(), 
				                      r.vec(),
				                      rr.vec(), tide, it,
				                      tolerance);
				if (ce){
					beach.register_circle_event(it, ce);
				}
			}
			
			/* Valid event. Create Delaunay triangle: */
			delaunay_triangles.emplace_back(l.id(), id, r.id());
			
		}
	}
}


} // NAMESPACE ACOSA
