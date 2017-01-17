/* Implementation of a class representing the circle events of [1].
 * Part of ACOSA.
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

#include <circleevent.hpp>

namespace ACOSA {

//----------------------------------------------------------------------
CircleEvent::CircleEvent(const SphereVector& v1, const SphereVector& v2,
    const SphereVector& v3, const BeachIterator& it, double tolerance)
    : beach_site(it)
{
	/* Calculate circumcenter: */
	SphereVectorEuclid e1(v1), e2(v2), e3(v3);
	SphereVectorEuclid cc = SphereVectorEuclid::circumcenter(e1,e2,e3);

	/* Calculate maximum latitude of circumcircle: */
	double d = cc.distance(e1);

	lat_ = cc.lat() + d;

	/* Be certain to avoid numerical errors: */
	double vec_lat = std::max(v1.lat(), std::max(v2.lat(),
	                          v3.lat()));
	double dist = lat_ - vec_lat;
	if (dist > -tolerance && dist < tolerance){
		lat_ = vec_lat;
	}

	/* Copy iterator and set new latitude: */
	beach_site.tide_ = lat_;
}


//----------------------------------------------------------------------
CircleEvent::CircleEvent(double lat, const BeachIterator& it)
    : beach_site(it), lat_(lat)
{
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

} // NAMESPACE ACOSA
