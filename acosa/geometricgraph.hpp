/* A fast algorithm for the determination of links of a geometric
 * graph. Part of ACOSA.
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

#include <basic_types.hpp>
#include <vector>


namespace ACOSA {


/*!
 * \brief Calculate the links of a geometric graph on a sphere.
 * \param coordinates Coordinates of the nodes of the graph.
 * \param links       Target vector for the links of the graph.
 * \param sigma_0     The geometric graph's connection threshold.
 *                    All node pairs closer than sigma_0 will be
 *                    connected.
 *
 * A geometric graph is a spatially embedded graph in which
 * nodes are connected pairwise iff they are distanced equal to
 * or less than a threshold distance (sigma_0).
 *
 * The algorithm used is based off the spherically adapted
 * Fortune's sweepline algorithm from [1]. It has a complexity
 * of O(N*log(N) + N*M + M^2) where N is the number of nodes and
 * M the mean degree of the geometric graph. For uniform spatial
 * node distributions, M~sigma_0^2.
 */
void geometric_graph_links(const std::vector<Node>& coordinates,
    std::vector<Link>& links, double sigma_0);

}
