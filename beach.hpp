/* ACOSA classes representing the beach line in the spherical Fortune's
 * alrogithm in [1].
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

#ifndef ACOSA_BEACH_HPP
#define ACOSA_BEACH_HPP

#include <map>
#include <stddef.h>
#include <math.h>
#include <spherics.hpp>
#include <memory>

namespace ACOSA {

class Beach;
class CircleEvent;

class ArcIntersect {
	friend struct ArcIntersectComparator;
	friend class Beach;
	
	public:
		ArcIntersect(double position, size_t id, const SphereVector& vec,
		          const SphereVector& left);
	
	
		double lon_left(double tide, double anchor) const;
		
		
		size_t id() const;
		
		const SphereVector& vec() const;
		
		const SphereVector& left() const;
		
	
	private:
		ArcIntersect();
		
		size_t       id_;
		SphereVector vec_;
		SphereVector left_;

		
		/* This is used for ordering: */
		double pos;
};






class BeachSiteData {
	public:
		BeachSiteData(const std::shared_ptr<bool>& valid = 
					  std::shared_ptr<bool>());
	
		void register_circle_event_ptr(const std::shared_ptr<bool>& 
		                               valid);
		
		std::shared_ptr<bool> valid() const;
		
		void invalidate();
		
	private:
		/* A reference to the circle event: */
		std::shared_ptr<bool> valid_;
	
	
};






class BeachIterator;


class Beach {
	friend class BeachIterator;
	
	public:
		Beach(size_t id1, const SphereVector& v1, size_t id2,
		      const SphereVector& v2);
		
		BeachIterator begin(double tide);
		
		BeachIterator find_insert_position(double d, double tide);
		
		/* Warning: The iterator 'pos' is invalidated! */
		BeachIterator insert_before(const BeachIterator& pos, size_t id,
		                            const SphereVector& vec);
		
		void register_circle_event(BeachIterator& pos,
		                           const std::shared_ptr<bool>& valid);
		
		void invalidate_circle_event(BeachIterator& pos);
		
		BeachIterator erase(BeachIterator& pos);
		
		void set_tide(double d);
		
		size_t size() const;
		
		void print_debug() const;
		
		bool check_consistency(double tide) const;
	
	private:
		constexpr static bool check_increasing_tide = true;
		
		mutable double tide   = 0.0;
		mutable double anchor = 0.0;
	
		struct ArcIntersectComparator
		{
			using is_transparent = std::true_type;
			
			ArcIntersectComparator(double* tide, double* anchor)
				: tide(tide), anchor(anchor)
			{}
			
			
			/* Comparison data: */
			double* tide = nullptr;
			double* anchor = nullptr;
			
			/* Comparison operators: */
			
			/* Compare beach sites by their ordering parameter: */
			bool operator()(const ArcIntersect& l, const ArcIntersect& r) const
			{
				return l.pos < r.pos;
			}
			
			/* Compare a beach site with a double by comparing the double to
			 * the left border of the arc: */
			bool operator()(const ArcIntersect& l, double lon) const
			{
				lon -= *anchor;
				if (lon < 0.0)
					lon += 2*M_PI;
				
				return l.lon_left(*tide, *anchor) < lon;
			}
			
			bool operator()(double lon, const ArcIntersect& r) const
			{
				lon -= *anchor;
				if (lon < 0.0)
					lon += 2*M_PI;
				
				return lon < r.lon_left(*tide, *anchor);
			}
		};

	
		ArcIntersectComparator compare;
	
		std::map<ArcIntersect,BeachSiteData,ArcIntersectComparator>  data;
	
};

class CircleEvent;

/* An iterator  */
class BeachIterator {
	friend class Beach;
	friend class CircleEvent;
	
	
	public:
		typedef
		std::map<ArcIntersect, BeachSiteData, 
			Beach::ArcIntersectComparator>::iterator
		internal_iterator;
	
		BeachIterator(Beach* beach, internal_iterator it,
		              double tide);
	
		std::pair<const ArcIntersect, BeachSiteData>& operator*();
		
		std::pair<const ArcIntersect, BeachSiteData>* operator->();
		
		const std::pair<const ArcIntersect, BeachSiteData>* operator->()
			const;
		
		BeachIterator& operator++();
		
		BeachIterator& operator--();
		
		const SphereVector& vec() const;
		
		size_t id() const;
	
		double tide() const;
		
		bool lon_left_equal(double lon) const;
	
		/* For debugging purposes: */
		bool is_valid() const;
	
	private:
		double tide_;
	
		internal_iterator it;
		
		Beach* beach;
};





typedef std::pair<const ArcIntersect, BeachSiteData> BeachSite;

} // NAMESPACE ACOSA

#endif // ACOSA_BEACH_HPP
