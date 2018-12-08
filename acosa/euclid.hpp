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

#ifndef EUCLID_HPP
#define EUCLID_HPP

#include <basic_types.hpp>

namespace ACOSA {

class Vector2d {
	public:
		Vector2d(double x, double y);

		/* Coordinate access: */
		double x() const;

		double y() const;

		/* Arithmetic operators: */
		Vector2d operator+(const Vector2d& other) const;

		Vector2d operator-(const Vector2d& other) const;

		Vector2d operator*(double d) const;

		Vector2d operator/(double d) const;

		double operator*(const Vector2d& other) const;


		/* Convenience methods: */
		double norm() const;

		/* Return a normed normal vector: */
		Vector2d normal() const;

	private:
		double x_;
		double y_;
};

/* Operators: */
Vector2d operator*(double d, const Vector2d& vec);


} // NAMESPACE ACOSA

#endif
