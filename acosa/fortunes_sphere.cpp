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


#include <queue>
#include <memory>

#include <iostream>
#include <beach.hpp>

namespace ACOSA {

//######################################################################

class CircleEvent;

typedef std::priority_queue<CircleEvent, std::vector<CircleEvent>,
	                        std::greater<CircleEvent>>
	    eventqueue_t;




//######################################################################
class CircleEvent {
	public:
		CircleEvent(const SphereVector& v1, const SphereVector& v2,
					const SphereVector& v3, const BeachIterator& it);
		
		double lat() const;
		
		bool operator<(const CircleEvent& other) const;
		
		bool operator>(const CircleEvent& other) const;
		
		void invalidate();
		void set_valid();
		
		bool valid() const;
		
		void set_beach_site_iterator(const BeachIterator& iterator);
		
	public:
		std::shared_ptr<bool> valid_;
		double                lat_;
		BeachIterator     beach_site;
};



//----------------------------------------------------------------------
CircleEvent::CircleEvent(const SphereVector& v1, const SphereVector& v2,
	const SphereVector& v3, const BeachIterator& it)
	: beach_site(it)
{
	/* Calculate circumcenter: */
	SphereVectorEuclid e1(v1), e2(v2), e3(v3);
	SphereVectorEuclid cc = SphereVectorEuclid::circumcenter(e1,e2,e3);
	
	/* Calculate maximum latitude of circumcircle: */
	double d = cc.distance(e1);
	
	lat_ = cc.lat() + d;
	
	/* Copy iterator and set new latitude: */

	beach_site.tide_ = lat_;
}



//----------------------------------------------------------------------
double CircleEvent::lat() const
{
	return lat_;
}


//----------------------------------------------------------------------
bool CircleEvent::operator<(const CircleEvent& other) const
{
	return lat_ < other.lat_;
}

//----------------------------------------------------------------------
bool CircleEvent::operator>(const CircleEvent& other) const
{
	return lat_ > other.lat_;
}

//----------------------------------------------------------------------
void CircleEvent::set_beach_site_iterator(const BeachIterator&
	iterator)
{
	beach_site = iterator;
}

//----------------------------------------------------------------------
void CircleEvent::invalidate(){
	if (valid_){
		*valid_ = false;
		valid_.reset();
	}
}

//----------------------------------------------------------------------
void CircleEvent::set_valid() {
	valid_ = std::make_shared<bool>(true);
}

//----------------------------------------------------------------------
bool CircleEvent::valid() const
{
	return valid_ && *valid_;
}


//######################################################################


//----------------------------------------------------------------------
void voronoi_tesselation_sphere(const std::vector<Node>& nodes,
	std::vector<Node>& voronoi_nodes, std::vector<Link>& voronoi_links)
{
	/* Call Delaunay triangulation: */
	std::vector<Triangle> delaunay_triangles;
	delaunay_triangulation_sphere(nodes, delaunay_triangles);
}

void voronoi_weights_sphere(const std::vector<Node>& nodes,
	std::vector<double>& weights);



//----------------------------------------------------------------------
static std::shared_ptr<bool> add_circle_event(
	eventqueue_t& queue, const SphereVector& v1,
    const SphereVector& v2, const SphereVector& v3, double tide,
    const BeachIterator& beach_site)
{
	CircleEvent event(v1, v2, v3, beach_site);
	if (event.lat() > tide){
		event.set_valid();
		queue.push(event);
		return event.valid_;
	} 
	return nullptr;
}


//----------------------------------------------------------------------
void delaunay_triangulation_sphere(const std::vector<Node>& nodes,
	std::vector<Triangle>& delaunay_triangles)
{
	const size_t N = nodes.size();
	
	/* 1) Priority queue of site events.
	 *    1.1: Site event type: */
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
	
	/* TODO: Remove this and replace by operator! */
	/*    1.2: Site event comparison lambda: */
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
	SphereVector v1 = site_events.top().vec;
	size_t       i1 = site_events.top().id;
	site_events.pop();
	SphereVector v2 = site_events.top().vec;
	size_t       i2 = site_events.top().id;
	site_events.pop();
	Beach beach(i1, v1, i2, v2);
	
	
	
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
				beach.find_insert_position(vec.lon(), vec.lat());
			
			
			/* Check for degenerate case 
			 * 			vec.lon() == it->lon_left()
			 */
			if (it.lon_left_equal(vec.lon()))
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
				
				
				std::cerr << "ERROR : Degenerate case not handled.\n";
				throw;
				
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
									  tide, it);          
				if (ce){
					beach.register_circle_event(it, ce);
				}
				
				/* Update the circle event for the current node v_j
				 * if it exists ( [p_i, p_j, p3] ): */
				++it;
				++it;
				ce = add_circle_event(circle_events, vec, v_j, v_3,
					                       tide, it);
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
										   r.vec(), event.lat(), it);
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
										 rr.vec(), event.lat(), it);
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
