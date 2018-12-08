/* Euclidean vectors for ACOSA.
 * Copyright (C) 2018 Malte Ziebarth
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
 */


#include <euclid.hpp>
#include <utility>

#include <cmath>
//#include <math.h>

namespace ACOSA {

/* Vector2d: */
Vector2d::Vector2d(double x, double y) : x_(x), y_(y)
{
}

double Vector2d::x() const {
	return x_;
}

double Vector2d::y() const {
	return y_;
}

/* Arithmetic operators: */
Vector2d Vector2d::operator+(const Vector2d& other) const
{
	return Vector2d(x_ + other.x_, y_ + other.y_);
}

Vector2d Vector2d::operator-(const Vector2d& other) const
{
	return Vector2d(x_ - other.x_, y_ - other.y_);
}

Vector2d Vector2d::operator*(double d) const
{
	return Vector2d(x_*d, y_*d);
}

Vector2d Vector2d::operator/(double d) const
{
	return operator*(1.0/d);
}

double Vector2d::operator*(const Vector2d& other) const
{
	return x_*other.x_ + y_*other.y_;
}


/* Convenience methods: */

double Vector2d::norm() const
{
	return std::sqrt(x_*x_ + y_*y_);
}

/* Return a normed normal vector: */
Vector2d Vector2d::normal() const
{
	Vector2d nrm(-y_, x_);
	return nrm / nrm.norm();
}


/* Static operators: */
Vector2d operator*(double d, const Vector2d& vec)
{
	return vec * d;
}


} // NAMESPACE ACOSA
