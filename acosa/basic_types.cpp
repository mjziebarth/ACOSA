/* Basic data types used in the ACOSA
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


#include <basic_types.hpp>
#include <iostream>
#include <stdexcept>
#include <math.h>

namespace ACOSA {

//######################################################################

Triangle::Triangle(size_t i, size_t j, size_t k)
	: i(i), j(j), k(k)
{
	if (i==j || i==k || j==k)
	{
		std::string str("Created a Delaunay triangle with two equal "
		    "borders!\n\ti=");
		str.append(std::to_string(i)).append("\n\tj=").append(std::to_string(j))
		   .append("\n\tk=").append(std::to_string(k)).append("\n");
		throw std::runtime_error(str);
	}
}

Triangle::Triangle(): i(0), j(0), k(0)
{
}

bool Triangle::common_border(const Triangle& other) const
{
	int common_nodes = 0;
	if (i == other.i || i == other.j || i == other.k){
		++common_nodes;
	}
	if (j == other.i || j == other.j || j == other.k){
		++common_nodes;
	}
	if (k == other.i || k == other.j || k == other.k){
		++common_nodes;
	}
	return common_nodes == 2;
}

//######################################################################

Link::Link(size_t i, size_t j) : i(i), j(j)
{
}

Link::Link() : i(0), j(0)
{
}

bool Link::operator==(const Link& other) const
{
	return i == other.i && j == other.j;
}

bool Link::operator!=(const Link& other) const
{
	return i != other.i || j != other.j;
}

bool Link::operator<(const Link& other) const
{
	if (i != other.i){
		return i < other.i;
	}
	return j < other.j;
}

//######################################################################

Node::Node(double lon, double lat) : lon(lon), lat(lat)
{
	if (lon > 2.0*M_PI || lon < 0 || lat > M_PI_2 || lat < -M_PI_2)
	{
		throw std::domain_error("ERROR : ACOSA::Node() :\nlon or lat out of"
		                        "bounds. Expected ranges are [0, 2pi] for "
		                        "longitude and [-pi/2, pi/2] for latitude.\n");
	}
}

Node::Node() : lon(0), lat(0)
{
}

} // NAMESPACE ACOSA
