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


#include <beach.hpp>
#include <memory>

namespace ACOSA {

class BeachIterator;

//######################################################################
class CircleEvent {
    public:
	    CircleEvent(const SphereVector& v1, const SphereVector& v2,
		            const SphereVector& v3, const BeachIterator& it,
		            double tolerance);

		/* This constructor has been created to allow the creation of a
		 * circle event for the north pole. */
		CircleEvent(double lat, const BeachIterator& it);

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



} // NAMESPACE ACOSA
