/* The spherical implementation of the versatile-voronoi collection.
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

#ifndef ACOSA_SPHERICS_HPP
#define ACOSA_SPHERICS_HPP

#include <basic_types.hpp>

namespace ACOSA {

class SphereVector {
	public:
		SphereVector();
	
		SphereVector(double lon, double lat);
		
		bool null() const;
		
		double lon() const;
		
		double lat() const;
		
		static SphereVector circumcenter(
			const SphereVector& v1, const SphereVector& v2,
			const SphereVector& v3);
	
	private:
		double lon_;
		double lat_;
		
};

//######################################################################


/* Use threedimensional vectors internally: */
class SphereVectorEuclid{
	friend class VConfSphere;
	friend class ExternalVector;
		
	public:
		SphereVectorEuclid();
	
		SphereVectorEuclid(double x, double y, double z);
		
		SphereVectorEuclid(double lon, double lat);
		
		SphereVectorEuclid(const SphereVector& v);
		
		SphereVectorEuclid(const Node& node);
	
		/* Arithmetic operations: */
		SphereVectorEuclid operator-(const SphereVectorEuclid& other)
			const;
			
		SphereVectorEuclid operator-() const;
			
		SphereVectorEuclid operator+(const SphereVectorEuclid& other)
			const;
			
		SphereVectorEuclid cross(const SphereVectorEuclid& other) const;
		
		double operator*(const SphereVectorEuclid& other) const;
		
		SphereVectorEuclid operator*(double d) const;
		
		SphereVectorEuclid& operator*(double d);
		
		void operator*=(double d);
		
		void operator-=(const SphereVectorEuclid& other);
		
		/* Attributes: */
		double norm() const;
		
		/* Squared norm */
		double norm2() const;
		
		double lon() const;
		
		double lat() const;
		
		bool is_null() const;
		
		/* Sphere operations: */
		double distance(const SphereVectorEuclid& other) const;
		
		double angle(const SphereVectorEuclid& to1,
			const SphereVectorEuclid& to2) const;
		
		
		/* Static: */
		static SphereVectorEuclid geodesic_center(
			const SphereVectorEuclid& v1, const SphereVectorEuclid& v2);
			
		static SphereVectorEuclid circumcenter(
			const SphereVectorEuclid& v1, const SphereVectorEuclid& v2,
			const SphereVectorEuclid& v3);
		
		static double triangle_area(const SphereVectorEuclid& v1,
			const SphereVectorEuclid& v2, const SphereVectorEuclid& v3);
		
	private:
		double x;
		double y;
		double z;
	
};

SphereVectorEuclid operator*(double d, const SphereVectorEuclid& v);

//######################################################################

} // NAMESPACE ACOSA

#endif // SPHERICS_HPP
