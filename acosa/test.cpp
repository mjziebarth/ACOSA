/* ACOSA command line test suite.
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

#include <vdtesselation.hpp>
#include <convexhull.hpp>
#include <order_parameter.hpp>

#include <random>
#include <math.h>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <sstream>
#include <limits>
#include <unistd.h>
#include <cstdlib>
#include <list>
#include <chrono>


const size_t N = 1000000;

/* This method obtains a random seed for testing. Adjust this method
 * if you want to reproduce a test run. */
static long random_seed(){
//	return 5186709571096577860; // An example test run that produced an error.
	
	std::random_device rd;
    std::uniform_int_distribution<long>
		dist(std::numeric_limits<long>::min(),
		     std::numeric_limits<long>::max());
    
    
    return dist(rd);
}



struct configuration {
	size_t N;
	size_t runs;
	bool   test_order_param;
	bool   regular_grid;
	char   grid_type;
	size_t r;
	bool r_selected;
	bool    scaled_dbg_output;
	std::string testfile;
};


static configuration get_config(int argc, char **argv){
	configuration conf = {0,  1, false, false, 0, 0, false, false};
	
	char *Nvalue = nullptr;
	char *Rvalue = nullptr;
	char *rvalue = nullptr;
	char *gridtype = nullptr;
	char *file = nullptr;
	int index;
	int c;

	opterr = 0;

	while ((c = getopt (argc, argv, "R:ON:r:G:Df:")) != -1){
		switch (c)
		{
			case 'r':
				rvalue = optarg;
				std::cout << "rvalue: '" << rvalue << "'\n";
				break;
			case 'R':
				Rvalue = optarg;
				std::cout << "Rvalue: '" << Rvalue << "'\n";
				break;
			case 'O':
				conf.test_order_param = true;
				std::cout << "Testing order parameter!\n";
				break;
		    case 'G':
			    gridtype = optarg;
			    conf.regular_grid = true;
				std::cout << "Using a regular grid.\n";
			    break;
		    case 'D':
			    conf.scaled_dbg_output = true;
				std::cout << "Using extra debug output.\n";
			    break;
			case 'N':
				Nvalue = optarg;
				std::cout << "Nvalue: '" << Nvalue << "'\n";
				break;
			case 'f':
				file = optarg;
				std::cout << "Using test data file '" << file << "'\n";
				conf.testfile = std::string(file);
				conf.runs = 1;
				break;
			case '?':
				if (optopt == 'c')
					std::cerr << "Option -" << optopt
					          << " requires an argument.\n";
				else if (isprint(optopt))
					std::cerr << "Unknown option '-" << (char)optopt << "\'.\n";
				else
					std::cerr << "Unknown option character '"
							  << (char)optopt  << "'\n";
			default:
				return {0,  1, false, false, 0, 0, false, false};
		}
	}
	if (Nvalue){
		conf.N = std::atol(Nvalue);
	}
	if (Rvalue){
		conf.runs = std::atol(Rvalue);
	}
	if (rvalue){
		conf.r_selected = true;
		conf.r = std::atol(rvalue);
	}
	if (gridtype){
		conf.grid_type = std::atoi(gridtype);
	}
	return conf;
}


/*!
 * This method tests the OrderParameter class. It starts by creating a
 * list of two OrderParameters: max() and min()
 * Then, it successively inserts new OrderParameter in the ordered list
 * by randomly selecting two existing neighbours of the list,
 * calculating their mean, and inserting the mean into the list.
 * 
 * At each step, the order left < mean < right is checked.
 */
void test_order_parameter(size_t N) {
	/* RNG to calculate insert positions: */
	std::mt19937_64 engine(random_seed());
	std::uniform_real_distribution<double> generator;
	
	/* Ordered list of parameters: */
	std::list<ACOSA::OrderParameter> order_params;
	order_params.push_back(ACOSA::OrderParameter::max());
	order_params.push_front(ACOSA::OrderParameter::min());
	for (size_t i=2; i<N; ++i){
		/* Calculate the insert position, alway between the first and
		 * last element: */
		size_t pos = generator(engine) * i; 
		pos = std::min(std::max(pos, (size_t)1), i-1);
		auto it = order_params.begin();
		for (size_t j=0; j<pos; ++j){
			++it;
		}
		--it;
		ACOSA::OrderParameter left = *it;
		++it;
		ACOSA::OrderParameter right = *it;
		ACOSA::OrderParameter middle 
			= ACOSA::OrderParameter::between(left, right);
		if (!(middle > left && middle < right)){
			std::cerr <<   "ERROR in ordering!\n\tleft:  "
			          << left.to_string() << "\n\tright: "
			          << right.to_string()<< "\n\tmid:   "
			          << middle.to_string() << "\n";
			throw 0;
		}
		order_params.insert(it, middle);
	}
}


/*!
 * \brief longitude_grid_points
 * \param N
 * \return
 */
size_t longitude_grid_points(size_t N)
{
	/* Start a guess such that guess*(guess/2) = N */
	size_t guess = sqrt(2*N);

	/* Now find biggest smaller factor of N: */
	while (guess && N % guess){
		--guess;
	}

	if (!guess){
		std::cerr << "ERROR : N does not have any factors smaller than "
		             "sqrt(N/2)!\n";
		exit(-1);
	}

	return guess;
}


/*!
 * \brief Read a number of comma seperated double values from a character
 *        stream.
 * \param stream The character stream. Should contain only comma seperated
 *               floating point values. No error checking is executed.
 *               Everything between commas is directly passed to std::stod.
 * \param vec Target vector that will be filled.
 */
static void read_double_vector(std::ifstream& stream, std::vector<double>& vec)
{
	char c;
	/* Skip white spaces: */
	while (!stream.eof() && (c = stream.get()) != '-' &&
	       (c < '0' || c > '9'))
	{
		if (c == '}'){
			std::cerr << "ERROR : Empty vector braces initialization given in "
			             "data file. Aborting!\n";
			exit(-1);
		}
	}
	/* Read whole braced block: */
	std::stringstream sstream;
	stream.get(*sstream.rdbuf(), '}');
	std::string numbers(sstream.str());
	size_t pos = 0;
	while (pos < numbers.size())
	{
		/* Find next comma: */
		size_t next_pos = numbers.find(',', pos);

		/* Convert substring to number: */
		if (next_pos == std::string::npos){
			vec.push_back(std::stod(numbers.substr(pos, std::string::npos)));
			break;
		} else {
			vec.push_back(std::stod(numbers.substr(pos, next_pos-pos)));
		}

		/* Continue beyond comma: */
		pos = next_pos+1;
	}

}

/*!
 * \brief Reads a data file containing a list of each latitude and longitude
 *        coordinates.
 * \param filename A path to the file.
 * \return A vector of nodes: (lon,lat)-pairs.
 *
 * The file should contain two lists of comma seperated floating point values
 * marked by 'lon={' and 'lat={' and each delimited by '}'.
 *
 * Everything else is ignored. Only the first occurence of the block begin
 * keywords are evaluated. Not much processing is done.
 */
static std::vector<ACOSA::Node> read_data_file(const std::string& filename)
{
	std::ifstream instream;
	instream.open(filename);
	if (!instream.good()){
		std::cerr << "Failed to open file '" << filename << "'. Aborting.\n";
		exit(-1);
	}
	if (instream.eof()){
		std::cerr << "Empty file '" << filename << "'. Aborting.\n";
		exit(-1);
	}

	char tag[5] = {0, 0, 0, 0, 0};
	bool has_lon = false;
	bool has_lat = false;

	std::vector<double> lon;
	std::vector<double> lat;

	while ((!has_lon || !has_lat) && !instream.eof())
	{
		/* Read one char: */
		while (!instream.eof() && instream.get() != 'l')
		{
		}

		/* Read 5 chars: */
		instream.get(tag, 5);

		if (tag[0] == 'a' && tag[1] == 't' && tag[2] == '=' && tag[3] == '{'){
			/* Read latitutde: */
			read_double_vector(instream, lat);
			has_lat = true;
		}
		else if (tag[0] == 'o' && tag[1] == 'n' && tag[2] == '=' &&
		         tag[3] == '{')
		{
			/* Read longitutde: */
			read_double_vector(instream, lon);
			has_lon = true;
		}
	}

	/* Sanity checks: */
	if (!has_lon || !has_lat)
	{
		std::cerr << "ERROR : At least one of keywords \"lat={\" or \"lon={\" "
		             "not found! Aborting.\n";
		exit(-1);
	}
	if (lon.size() != lat.size()){
		std::cerr << "ERROR : Longitude and latitude vectors are not equally "
		            "sized in data file! Aborting.\n";
		exit(-1);
	}

	/* Now fill return vector: */
	std::vector<ACOSA::Node> vec(lon.size());
	for (size_t i=0; i<lon.size(); ++i){
		vec[i] = ACOSA::Node(M_PI/180.0*lon[i], M_PI/180.0*lat[i]);
	}

	return vec;
}





/*!
 * This compiles into a small testing application. The application
 * creates a number of sets of randomly distributed points on a unit
 * sphere and calculates their VDTesselation, the Voronoi network,
 * the Voronoi cell areas, the Delaunay network, and a convex hull using
 * a random inside direction (this is not yet really sensible since
 * for a randomly distributed network it's hard to define the
 * "outside").
 * 
 * The application can be called from the command line using the
 * following parameters:
 * "-N x" : Defines N := x, the size of each of the random set of nodes.
 *          Needs to be supplied as there is no default.
 * "-R x" : Defines R := x, the number of runs or number of different
 *          sets generated.
 * "-r x" : Selects r := x, the number of the first executed run (all
 *          other runs are empty generated nodes. This is a feature to
 *          select a known bad run for a certain seed for debugging).
 * "-G x" : Nodes are distributed on a regular grid instead of being
 *          randomly distributed. The grid type is tweaked using x:
 *          x=0: Grid is rectangular in longitude/latitude space.
 *          x=1: Grid is (near) hexagonal in lon/lat space (every
 *               second line is shifted by half a grid distance
 *               in longitude).
 * "-D"   : Print debug output that scales with N.
 * "-O"   : A different test mode is chosen where the OrderParameter
 *          class is tested.
 */
int main(int argc, char **argv){
	/* Parse command line argument: */
	configuration c = get_config(argc, argv);
	
    #ifdef ACOSA_HIST
	ACOSA::OrderParameter::hist.clear();
    #endif

	if (c.test_order_param){
		test_order_parameter(c.N);
		return 0;
	}
	
	
	size_t N = c.N;
	if (!N && c.testfile.empty()){
		std::cerr << "N==0, returning!\n";
		return -1;
	}
	
	/* Initialize random number generator: */
	long seed = random_seed();
	if (!c.regular_grid){
		std::cout << "Create random nodes. (seed=" << seed << ")\n";
	}
	std::mt19937_64 engine(random_seed());
	std::uniform_real_distribution<double> generator;
	
	for (size_t r=0; r<c.runs; ++r){
		std::cout << "run " << r << "/" << c.runs << "\n";

		/* Create node positions:: */
		std::vector<ACOSA::Node> nodes(N);

		if (!c.testfile.empty()){
			std::cout << "Using grid in file '" << c.testfile << "'.\n";
			nodes = read_data_file(c.testfile);
			N = nodes.size();
		} else if (c.regular_grid){
			/* Create regular grid: */
			size_t n_lon = longitude_grid_points(N);
			size_t n_lat = N/n_lon;
			std::cout << "Create regular grid.\n\tn_lon=" << n_lon <<
			             "\n\tn_lat=" << n_lat << "\n";
			for (size_t j=0; j<n_lat; ++j){
				for (size_t i=0; i<n_lon; ++i){
					double lon;
					if (c.grid_type == 0){
						lon = 2*M_PI*((double)i)/n_lon;
					} else if (c.grid_type == 1){
						lon = 2*M_PI*((double)i + 0.5*(j % 2))/n_lon;
					}
					double lat = -M_PI*(0.5-((double)j+1)/(n_lat+1));
					nodes[n_lon*j+i] = ACOSA::Node(lon, lat);
				}
			}
		} else {
			/* Create random nodes: */
			for (size_t i=0; i<N; ++i){
				nodes[i] = ACOSA::Node(2*M_PI*generator(engine),
				                       M_PI*(0.5-generator(engine)));
			}
		}

		if (c.r_selected && r != c.r){
			/* If we want a certain run, fast forward: */
			std::cout << "  --> skipping.\n";
			continue;
		}
		
		/* Create tesselation: */
		std::cout << "Create tesselation.\n";
		auto t1 = std::chrono::high_resolution_clock::now();
		ACOSA::VDTesselation tesselation(nodes);
		auto t2 = std::chrono::high_resolution_clock::now();
		std::cout << "  --> elapsed: "
		          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1)
		             .count()
		          << "s\n";
		if (c.scaled_dbg_output){
			tesselation.print_debug(false);
		}
		
		/* Obtain network: */
		std::cout << "Obtain Voronoi network.\n";
		std::vector<ACOSA::Node> voronoi_nodes;
		std::vector<ACOSA::Link> voronoi_links;
		t1 = std::chrono::high_resolution_clock::now();
		tesselation.voronoi_tesselation(voronoi_nodes, voronoi_links);
		t2 = std::chrono::high_resolution_clock::now();
		std::cout << "  --> elapsed: "
		          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1)
		             .count()
		          << "s\n";

		voronoi_nodes.clear();
		voronoi_links.clear();
		
		/* Obtain areas: */
		std::cout << "Obtain Voronoi areas.\n";
		std::vector<double> weights;
		tesselation.voronoi_cell_areas(weights);
		double total_weight = 0.0;
		for (double w : weights){
			total_weight += w;
		}
		std::cout << "  --> total weight: " << total_weight << "\n"
		             "      (that's " << total_weight-4.0*M_PI
		          << " difference to 4pi)\n";
		
		/* Obtain Delaunay tesselation: */
		std::cout << "Obtain Delaunay tesselation.\n";
		std::vector<ACOSA::Link> delaunay_links;
		tesselation.delaunay_triangulation(delaunay_links);
		
		/* Obtain hull: */
		std::cout << "Obtain hull.\n";
		ACOSA::Node inside = ACOSA::Node(2*M_PI*generator(engine),
										 M_PI*(0.5-generator(engine)));
		ACOSA::ConvexHull hull(nodes, inside);
		std::cout << "hull size = " << hull.size();
	
	}
		
	/* Informative output: */
    #ifdef ACOSA_HIST
	std::cout << "Histogram of OrderParameter lengths:\n";
	for (size_t i=0; i<ACOSA::OrderParameter::hist.size(); ++i){
		std::cout << "\t[ " << i << ": " << ACOSA::OrderParameter::hist[i] 
				  << " ]\n";
	}
	std::cout << "\n";
    #endif
	
	std::cout << "Finished.\n";
	return 0;
}
