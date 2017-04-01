/* Basic data types used in ACOSA
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
 */

#ifndef ACOSA_BASIC_TYPES_HPP
#define ACOSA_BASIC_TYPES_HPP

#include <stddef.h>
#include <functional>
#include <limits>

namespace ACOSA {


static constexpr size_t NO_LINK = std::numeric_limits<size_t>::max();

struct Triangle {
	public:
		Triangle();
		Triangle(size_t i, size_t j, size_t k);
		
		bool common_border(const Triangle& other) const;
	
	size_t i;
	size_t j;
	size_t k;
};


struct Link {
	public:
		Link();
		Link(size_t i, size_t j);

		bool operator==(const Link& other) const;
		bool operator!=(const Link& other) const;

		bool operator<(const Link& other) const;
	
	size_t i;
	size_t j;
};


struct Node {
	public:
		Node();
		Node(double lon, double lat);
	
	double lon;
	double lat;
};

} // NAMESPACE ACOSA

namespace std
{
    template<> struct hash<ACOSA::Link>
    {
		size_t operator()(const ACOSA::Link& link) const;
    };
}


#endif // ACOSA_BASIC_TYPES_HPP
