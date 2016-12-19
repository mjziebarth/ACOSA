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

#include <spherics.hpp>

#include <cmath>
#include <math.h>

namespace ACOSA {

/* ****************************************************************** */
/*                           SphereVector                             */
/* ****************************************************************** */
SphereVector::SphereVector()
	: lon_(0), lat_(M_PI)
{
}

//----------------------------------------------------------------------
SphereVector::SphereVector(double lon, double lat)
	: lon_(lon), lat_(lat)
{
}

//----------------------------------------------------------------------
bool SphereVector::null() const {
	return lat_ == M_PI;
}

//----------------------------------------------------------------------
double SphereVector::lon() const
{
	return lon_;
}

//----------------------------------------------------------------------
double SphereVector::lat() const
{
	return lat_;
}


//----------------------------------------------------------------------
SphereVector SphereVector::circumcenter(const SphereVector& v1,
	const SphereVector& v2, const SphereVector& v3)
{
	SphereVectorEuclid ve1(v1);
	SphereVectorEuclid ve2(v2);
	SphereVectorEuclid ve3(v3);
	SphereVectorEuclid cce = SphereVectorEuclid::circumcenter(ve1, ve2,
	                                                          ve3);
	
	return SphereVector(cce.lon(), cce.lat());
}




/* ****************************************************************** */
/*                        SphereVectorEuclid                          */
/* ****************************************************************** */

SphereVectorEuclid::SphereVectorEuclid() : x(0), y(0), z(0)
{
}

//----------------------------------------------------------------------
SphereVectorEuclid::SphereVectorEuclid(double x, double y, double z)
	: x(x), y(y), z(z)
{
}

//----------------------------------------------------------------------
SphereVectorEuclid::SphereVectorEuclid(const SphereVector& vec)
{
	double clat = std::cos(vec.lat());
	x = std::cos(vec.lon())*clat;
	y = std::sin(vec.lon())*clat;
	z = std::sin(vec.lat());
}

//----------------------------------------------------------------------
SphereVectorEuclid::SphereVectorEuclid(const Node& node)
{
	double clat = std::cos(node.lat);
	x = std::cos(node.lon)*clat;
	y = std::sin(node.lon)*clat;
	z = std::sin(node.lat);
}

//----------------------------------------------------------------------
double SphereVectorEuclid::lat() const 
{
	if (z > 1.0){
		return 0.5*M_PI;
	} else if (z < -1.0) {
		return -0.5*M_PI;
	}
	return std::asin(z);
}

//----------------------------------------------------------------------
double SphereVectorEuclid::lon() const 
{
	double lon_ = std::atan2(y, x);
	if (lon_ < 0){
		lon_ += 2*M_PI;
	}
	return lon_;
}

//----------------------------------------------------------------------
double SphereVectorEuclid::operator*(const SphereVectorEuclid& other)
	const
{
	return x*other.x + y*other.y + z*other.z;
}

//----------------------------------------------------------------------
SphereVectorEuclid
SphereVectorEuclid::operator+(const SphereVectorEuclid& other) const
{
	return SphereVectorEuclid(x+other.x, y+other.y, z+other.z);
}

//----------------------------------------------------------------------
double SphereVectorEuclid::distance(const SphereVectorEuclid& other)
	const
{
	/* Snippet adapted from pik-copan/pyunicorn.
	 * On a unit sphere, the native distance metric is the greatcircle
	 * distance, which equals angles between the point vectors.
	 * Angles between vectors v1 and v2 can be calculated using the
	 * dot product relation
	 *     v1 * v2 = |v1| |v2| cos(angle(v1,v2))
	 */
	double d = operator*(other);
	if (d > 1.0){
		return 0.0;
	} else if (d < -1.0){
		return M_PI;
	}
	return std::acos(d);
}

//----------------------------------------------------------------------
SphereVectorEuclid SphereVectorEuclid::cross(
	const SphereVectorEuclid& other) const
{
	return SphereVectorEuclid(y * other.z  -  other.y * z,
	                          z * other.x  -  other.z * x,
	                          x * other.y  -  other.x * y);
}

//----------------------------------------------------------------------
SphereVectorEuclid SphereVectorEuclid::operator-(
	const SphereVectorEuclid& other) const
{
	return SphereVectorEuclid(x-other.x, y-other.y, z-other.z);
}

//----------------------------------------------------------------------
SphereVectorEuclid SphereVectorEuclid::operator-() const
{
	return SphereVectorEuclid(-x, -y, -z);
}

//----------------------------------------------------------------------
void SphereVectorEuclid::operator-=(const SphereVectorEuclid& other)
{
	x -= other.x;
	y -= other.y;
	z -= other.z;
}

//----------------------------------------------------------------------
void SphereVectorEuclid::operator*=(double d)
{
	x *= d;
	y *= d;
	z *= d;
}

//----------------------------------------------------------------------
SphereVectorEuclid& SphereVectorEuclid::operator*(double d)
{
	operator*=(d);
	return *this;
}

SphereVectorEuclid SphereVectorEuclid::operator*(double d) const
{
	return SphereVectorEuclid(d*x, d*y, d*z);
}

//----------------------------------------------------------------------
SphereVectorEuclid operator*(double d, const SphereVectorEuclid& v)
{
	return v * d;
}

//----------------------------------------------------------------------
double SphereVectorEuclid::norm() const
{
	return std::sqrt(x*x + y*y + z*z);
}

//----------------------------------------------------------------------
double SphereVectorEuclid::norm2() const
{
	return x*x + y*y + z*z;
}

//----------------------------------------------------------------------
double SphereVectorEuclid::angle(const SphereVectorEuclid& to1,
	const SphereVectorEuclid& to2) const
{
	/* Calculate vectors (from --> to1) and (from --> to2).
	 * They are inside the  */
	SphereVectorEuclid v1 = to1 - *this;
	SphereVectorEuclid v2 = to2 - *this;
	
	/* Now calculate the cross product between those vectors and from.
	 * That way it is guaranteed that the vectors in question are 
	 * tangential to the sphere at point from: */
	v1 = cross(v1);
	v2 = cross(v2);
	
	/* Now the angle between both is (take care since they're not
	 * normalized): */
	double c = v1 * v2 / (v1.norm()*v2.norm());
	double angle;
	if (c >= 1.0)
		angle = 0.0;
	else if (c <= -1.0)
		angle = M_PI;
	else
		angle = std::acos(c);
	
	return angle;
}


//----------------------------------------------------------------------
double SphereVectorEuclid::triangle_area(const SphereVectorEuclid& v1,
	const SphereVectorEuclid& v2, const SphereVectorEuclid& v3)
{
	double angle1 = v1.angle(v2, v3),
		angle2 = v2.angle(v1, v3),
		angle3 = v3.angle(v1, v2);
		
	return angle1 + angle2 + angle3 - M_PI;
}


//----------------------------------------------------------------------
SphereVectorEuclid
SphereVectorEuclid::circumcenter(const SphereVectorEuclid& v1,
	const SphereVectorEuclid& v2, const SphereVectorEuclid& v3)
{
	/* I have not checked the details why, but the minus is needed.
	 * Rest of formula taken from [1].
	 */
	SphereVectorEuclid v = -(v1 - v2).cross(v3-v2);
	v *= 1.0/ v.norm();
	return v;
}

} // NAMESPACE ACOSA
