/* ACOSA classes representing the beach line in the spherical Fortune's
 * alrogithm in [1].
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

#include <beach.hpp>

#include <cmath>
#include <iostream>

namespace ACOSA {
	
//typedef ArcIntersect::pos_t pos_t;

/* ****************************************************************** */
/*                        class ArcIntersect                          */
/* ****************************************************************** */
ArcIntersect::ArcIntersect()
{
}

ArcIntersect::ArcIntersect(const OrderParameter& position, size_t id,
	const SphereVector& vec, const SphereVector& left)
	: vec_(vec), id_(id), left_(left), pos(position)
{
}


//----------------------------------------------------------------------
size_t ArcIntersect::id() const
{
	return id_;
}

//----------------------------------------------------------------------		
const SphereVector& ArcIntersect::vec() const
{
	return vec_;
}

//----------------------------------------------------------------------		
const SphereVector& ArcIntersect::left() const
{
	return left_;
}

//----------------------------------------------------------------------	

//----------------------------------------------------------------------
static double asin_bounded(double s){
	if (s > 1.0){
		return 0.5*M_PI;
	} else if (s < -1.0){
		return -0.5*M_PI;
	}
	return std::asin(s);
}

//----------------------------------------------------------------------
static double acos_bounded(double c){
	if (c > 1.0){
		return 0.0;
	} else if (c < -1.0){
		return M_PI;
	}
	return std::acos(c);
}

//----------------------------------------------------------------------
double ArcIntersect::lon_left(double tide, double anchor) const
{	
	/* Calculation following [1].*/
	double stide = std::sin(tide);
	double clat1 = std::cos(vec_.lat());
	double slat1 = std::sin(vec_.lat());
	double clat2 = std::cos(left_.lat());
	double slat2 = std::sin(left_.lat());
	double lon1 = vec_.lon();
	double lon2 = left_.lon();
	
	/* We have theta = pi/2-lat, so
	 *    cos(theta) --> sin(lat)   ;    sin(theta) --> cos(lat) */
	double e = (slat1 - slat2) * std::cos(tide);
	double a =   (stide - slat2) * clat1 * std::cos(lon1)
	           - (stide - slat1) * clat2 * std::cos(lon2);
	double b =   (stide - slat2) * clat1 * std::sin(lon1)
	           - (stide - slat1) * clat2 * std::sin(lon2);
	
	double inv_norm = 1.0 / std::sqrt(a*a+b*b);
	
	/* Calculate gamma using the cosine relation: */
	double gamma = acos_bounded(b * inv_norm);
	{
		/* The cosine relation still leaves gamma and -gamma as possible
		 * solutions. Decide between those by evaluating the sine
		 * relation for both solutions and deciding which of both has
		 * the least error (one should be in order of machine precision)
		 */
		double sgamma = std::sin(gamma);
		double err1 = std::abs(sgamma - a * inv_norm);
		double err2 = std::abs(sgamma + a * inv_norm);
		if (err1 > err2){
			gamma = -gamma;
		}
	}
	
	/* Now we can calculate the longitude: */
	double s = e*inv_norm;
	double lon = -gamma + asin_bounded(e * inv_norm);
	
	/* We need to make sure the longitude is inside the bounds: */
	if (lon < 0.0){
		lon += 2*M_PI;
	}
	
	/* Start from anchor: */
	lon -= anchor;
	if (lon < 0.0){
		lon += 2*M_PI;
	}
	
	return lon;
}





/* ****************************************************************** */
/*                         class BeachSiteData                        */
/* ****************************************************************** */
BeachSiteData::BeachSiteData(const std::shared_ptr<bool>& valid)
	: valid_(valid)
{
}

//----------------------------------------------------------------------	
void BeachSiteData::register_circle_event_ptr(
	const std::shared_ptr<bool>& valid)
{
	valid_ = valid;
}

//----------------------------------------------------------------------	
void BeachSiteData::invalidate()
{
	if (valid_){
		*valid_ = false;
		valid_.reset();
	}
}

//----------------------------------------------------------------------	
std::shared_ptr<bool> BeachSiteData::valid() const
{
	return valid_;
}




/* ****************************************************************** */
/*                         class BeachIterator                        */
/* ****************************************************************** */
BeachIterator::BeachIterator(Beach* beach,
                             BeachIterator::internal_iterator it_,
                             double tide)
    : beach(beach), it(it_), tide_(tide)
{
	/* Periodic wrapping: */
	if (it == beach->data.end()){
		it = beach->data.begin();
	}
}


//----------------------------------------------------------------------	
std::pair<const ArcIntersect, BeachSiteData>& BeachIterator::operator*() 
{
	return *it;
}

//----------------------------------------------------------------------	
const std::pair<const ArcIntersect, BeachSiteData>*
	BeachIterator::operator->() const 
{
	const std::pair<const ArcIntersect, BeachSiteData>* dat
		 = it.operator->();
	return dat;
}

//----------------------------------------------------------------------	
std::pair<const ArcIntersect, BeachSiteData>* 
	BeachIterator::operator->()
{
	std::pair<const ArcIntersect, BeachSiteData>* dat
		 = it.operator->();
	return dat;
}

//----------------------------------------------------------------------	
BeachIterator& BeachIterator::operator++()
{
	++it;
	if (it == beach->data.end()){
		it = beach->data.begin();
	}
	return *this;
}

//----------------------------------------------------------------------	
BeachIterator& BeachIterator::operator--()
{
	if (it == beach->data.begin()){
		it = beach->data.end();
	}
	--it;
	return *this;
}

//----------------------------------------------------------------------	
double BeachIterator::tide() const
{
	return tide_;
}

//----------------------------------------------------------------------
const SphereVector& BeachIterator::vec() const
{
	return it->first.vec();
}

//----------------------------------------------------------------------
size_t BeachIterator::id() const
{
	return it->first.id();
}


//----------------------------------------------------------------------	
bool BeachIterator::lon_left_equal(double lon) const
{
	/* Anchor from beach's first element: */
	double anchor = beach->data.begin()->first.lon_left(tide_, 0.0);
	
	/* Transform longitude using anchor: */
	lon -= anchor;
	if (lon < 0.0)
		lon += 2*M_PI;
	
	return it->first.lon_left(tide_, anchor) == lon;
}



//----------------------------------------------------------------------	
bool BeachIterator::is_valid() const
{
	for (auto it2=beach->data.begin(); it2 != beach->data.end(); ++it2){
		if (&*it == &*it2){
			return true;
		}
	}
	return false;
}


/* ****************************************************************** */
/*                            class Beach                             */
/* ****************************************************************** */
Beach::Beach(size_t id1, const SphereVector& v1, size_t id2,
	         const SphereVector& v2)
	: compare(&tide, &anchor), data(compare)
{
	if (v1.lat() == v2.lat()){
		/* TODO! */
		std::cerr << "Beach::Beach() :\nUnhandled degeneracy detected. "
			"Aborting!\n";
		throw;
	}
	tide =   v2.lat();
	anchor = v1.lon();
	OrderParameter middle 
		= OrderParameter::between(OrderParameter::max(),
		                          OrderParameter::min());
	OrderParameter left
		= OrderParameter::between(OrderParameter::min(), middle);
	OrderParameter right 
		= OrderParameter::between(OrderParameter::max(), middle);
//	std::cout << "middle=" << middle.to_string() << "\n";
//	std::cout << "left=" << left.to_string() << "\n";
//	std::cout << "right=" << right.to_string() << "\n";
	ArcIntersect b1(left, id1, v1, v2);
	ArcIntersect b2(right, id2, v2, v1);
	data.insert(BeachSite(b1, BeachSiteData()));
	data.insert(BeachSite(b2, BeachSiteData()));
}

//----------------------------------------------------------------------
BeachIterator Beach::begin(double tide)
{
	return BeachIterator(this, data.begin(), tide);
}

//----------------------------------------------------------------------
BeachIterator Beach::find_insert_position(double d, double tide)
{
	/* A check if tide increases (or at least behaves monotonely): */
	if (check_increasing_tide && tide < this->tide){
		std::cerr << "ERROR : Beach::find_insert_position() :\n"
			"Received decreasing tide!\nOld tide: " << this->tide
			<< "\nNew tide: " << tide << "\n";
		throw;
	}
	
	/* Adjust compare: */
	this->tide   = tide;
	anchor = data.begin()->first.lon_left(tide, 0.0);
	
	/* Find iterator: */
	BeachIterator::internal_iterator it = data.lower_bound<double>(d);
		
	/* Construct iterator: */
	return BeachIterator(this, it, tide);
	
}

//----------------------------------------------------------------------
BeachIterator Beach::insert_before(const BeachIterator& pos, size_t id,
		                           const SphereVector& vec)
{
	/* Get information about the node pointed to by it.
	 * This is necessary because we have to replace it (its 'left'
	 * vector needs to be replaced). */
	ArcIntersect r = pos->first;
	r.left_ = vec;
	
	
	/* Replace r: */
	BeachSiteData site_data = pos->second;
	auto it = data.erase(pos.it);
	it = data.insert(it, BeachSite(r, site_data));
	
	/* Get information about left node: */
	ArcIntersect l;
	if (it == data.begin()){
		it = data.end();
		--it;
		l = it->first;
		it = data.begin();
	} else {
		--it;
		l = it->first;
		++it;
	}
	
	
	/* Construct beach site to be inserted: */
	OrderParameter p;
	if (l.pos < r.pos){
		/* Set position to be middle between both: */
		p = OrderParameter::between(l.pos, r.pos);
	} else if (l.pos > r.pos){
		/* Set position to either be in between left and 2.0
		 * or right and 0.0 - whichever is the bigger interval.
		 * This should cause the position interval [0.0, 2.0] to
		 * stay somewhat uniformly sampled. */
		double dl = OrderParameter::max() - l.pos;
		double dr = OrderParameter::min() - r.pos;
		if (dr > dl){
			p = OrderParameter::between(r.pos, OrderParameter::min());
		} else {
			p = OrderParameter::between(l.pos, OrderParameter::max());
		}
	} else {
		std::cerr << "ERROR : Beach::insert_before() :\nTwo positions "
			"equal!\n";
		std::cerr << "\tleft = " << l.pos.to_string() << "\n"
		             "\tright= " << r.pos.to_string() << "\n";
		throw;
	}
	ArcIntersect m(p, id, vec, l.vec_);
	
	
	/* Insert beach site: */
	it = data.insert(it, BeachSite(m, BeachSiteData()));
	
	return BeachIterator(this, it, pos.tide());
}


//----------------------------------------------------------------------
void Beach::register_circle_event(BeachIterator& pos,
		                          const std::shared_ptr<bool>& valid)
{
	/* This is easy: */
	pos->second.register_circle_event_ptr(valid);
}


//----------------------------------------------------------------------
void Beach::invalidate_circle_event(BeachIterator& pos)
{
	/* This is easy: */
	pos->second.invalidate();
}


//----------------------------------------------------------------------
BeachIterator Beach::erase(BeachIterator& it)
{	
	/* Invalidate left neighbour: */
	--it;
	it->second.invalidate();
	++it;
	
	/* Invalidate circle event of erased iterator: */
	BeachSite b = *it;
	b.second.invalidate();
	double tide = it.tide();
	
	/* Obtain new iterator: */
	auto it_new = data.erase(it.it);
	if (it_new == data.end())
		it_new = data.begin();
	
	/* For this new iterator (which is the erase element's right
	 * neighbour), we have to update the left neighbour.
	 * This can be done by replacing the beach site: */
	std::pair<ArcIntersect, BeachSiteData> right = *it_new;
	right.second.invalidate();
	right.first.left_ = b.first.left_;
	it_new = data.insert(data.erase(it_new), right);
	
	/* Create and return new BeachIterator: */
	return BeachIterator(this, it_new, tide);
}


//----------------------------------------------------------------------
void Beach::set_tide(double d)
{
	tide = d;
}


//----------------------------------------------------------------------
size_t Beach::size() const
{
	return data.size();
}


//----------------------------------------------------------------------
void Beach::print_debug() const
{
	anchor = data.begin()->first.lon_left(tide, 0.0);
	std::cout << "\nTIDE: " << tide << ", ANCHOR: " << anchor;
	std::cout << "\nID and sort parameter:\n{";
	for (auto it = data.begin(); it != data.end(); ++it){
		std::cout << it->first.id_ << "[" << it->first.pos.to_string()
		          <<  "],"; 
	}
	std::cout.precision(3);
	std::cout << "}\ncoordinates:\n{";
	for (auto it = data.begin(); it != data.end(); ++it){
		std::cout << "(" << it->first.vec().lon() << "," 
		          << it->first.vec().lat() <<  "), "; 
	}
	std::cout << "}\nleft coordinates:\n{";
	for (auto it = data.begin(); it != data.end(); ++it){
		std::cout << "(" << it->first.left().lon() << "," 
		          << it->first.left().lat() <<  "), "; 
	}
	std::cout.precision(17);
	std::cout << "}\nparabola intersects:\n{";
	for (auto it = data.begin(); it != data.end(); ++it){
		std::cout << it->first.lon_left(tide, anchor) 
		          << ", "; 
	}
	std::cout << "}\n\n";
	
}

//----------------------------------------------------------------------
bool Beach::check_consistency(double tide) const
{
	double anchor = data.begin()->first.lon_left(tide, 0.0);
	double last_lon = 0.0;
	for (auto it = data.begin(); it != data.end(); ++it){
		double lon = it->first.lon_left(tide, anchor);
		if (lon + 1e-6 < last_lon){
			return false;
		}
		last_lon = lon;
	}
	
	return true;
}


} // NAMESPACE ACOSA
