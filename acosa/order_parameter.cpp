

#include <order_parameter.hpp>

#include <limits>
#include <iostream>
#include <stdio.h>

using ACOSA::OrderParameter;

#ifdef ACOSA_HIST
std::vector<unsigned long long> OrderParameter::hist;
#endif


OrderParameter OrderParameter::max() {
	return OrderParameter(std::numeric_limits<unsigned long>::max());
}



OrderParameter OrderParameter::min() {
	return OrderParameter(std::numeric_limits<unsigned long>::min());
}


	
OrderParameter::OrderParameter(unsigned long l) : first_digits(l)
{	
}
		
OrderParameter::OrderParameter() : first_digits(0)
{
}


OrderParameter OrderParameter::between_sub2(const OrderParameter& op1,
	const OrderParameter& op2, bool verbose)
{
	OrderParameter p(op1.first_digits);
	
	/* Find position where op1 and op2 differ: */
	size_t i=0;
	size_t op1_md_size = op1.more_digits.size();
	size_t op2_md_size = op2.more_digits.size();
	size_t max_size = std::max(op1_md_size, op2_md_size);
	while (i<max_size){
		if (op1.digit_at(i) != op2.digit_at(i))
			break;
		++i;
	}
	
	/* We know that the fraction vector of p must have at least
	 * i+1 elements, in some cases 2. To keep it simple, always
	 * use i+2 elements.
	 */ 
	p.more_digits.resize(i+2, 0);
	
	/* The first i elements are equal: */
	for (size_t j=0; j<i; ++j){
		p.more_digits[j] = op1.digit_at(j);
	}
	
	/* Debugging information: */
	#ifdef ACOSA_DEBUG
	if (verbose)
		std::cout << "OrderParameter::between_sub2() : i==" << i
		          << "\n";
	#endif
	
	/* Now mean between those two.
	 * We save in a long the last two ints of both arrays where
	 * l1 is the first int in which both differ: */
	unsigned long l1, l2;
	if (i >= op1_md_size && i >= op2_md_size){
		/* Both numbers were equal to the maximum number of digits! */
		std::string str("two identical operators!\n\t");
		str.append(op1.to_string()).append("\n\t").append(op2.to_string())
		   .append("\n");
		throw std::runtime_error(str);
	}
	
	/* Obtain the digits i and i+1 as a single unsigned long for op1: */
	if (i >= op1_md_size){
		l1 = 0;
	} else if (i+1 == op1_md_size) {
		l1 = ((unsigned long) op1.more_digits.back()) << 32;
	} else {
		l1 = (((unsigned long) op1.more_digits[i]) << 32) 
			 + op1.more_digits[i+1];
	}
	
	/* Obtain the digits i and i+1 as a single unsigned long for op1: */
	if (i >= op2_md_size){
		l2 = 0;
	} else if (i+1 == op2_md_size) {
		l2 = ((unsigned long) op2.more_digits.back()) << 32;
	} else {
		l2 = (((unsigned long) op2.more_digits[i]) << 32) 
			 + op2.more_digits[i+1];
	}
	
	unsigned long diff;
	unsigned long smaller;
	if (l1 > l2){
		diff = l1 - l2;
		smaller = l2;
	} else {
		diff = l2 - l1;
		smaller = l1;
	}
	unsigned long m = smaller + diff/2;
	unsigned int i1 = m >> 32;
	unsigned int i0 = m;
	p.more_digits[i]   = i1;
	p.more_digits[i+1] = i0;
	
	if (diff == 1){
		/* We need to append yet another digit because at the current
		 * resolution we cannot fit another number between the two:
		 * We would just end up with m == smaller.
		 * We need to find the first digit k were we can insert a number
		 * between bigger.digit_at(k) and digit_t::max(): */
		size_t k=i+1;
		constexpr digit_t max = std::numeric_limits<digit_t>::max();
		digit_t last;
		if (l1 > l2){
			while (op2.digit_at(++k) == max);
			last = op2.digit_at(k) + (max - op2.digit_at(k))/2;
		} else {
			while (op1.digit_at(++k) == max);
			last = op1.digit_at(k) + (max - op1.digit_at(k))/2;
		}
		p.more_digits.resize(k+1);
		for (size_t j=i+2; j<k; ++j){
			p.more_digits[j] = max;
		}
		p.more_digits[k] = last;
	}
	
	/* Optional debug output: */
	#ifdef ACOSA_DEBUG
	if (verbose){
		printf("in hex:\n");
		printf("l1 = %016llx\n",l1);
		printf("l2 = %016llx\n",l2);
		printf("m  = %016llx\n", m);
		printf("i0 = %016x\n",i0);
		printf("i1 = %016x\n",i1);
		std::cerr << "l1 = " << l1 << "\nl2 = " << l2 << "\ni1 = "
			<< i1 << "\ni0 = " << i0 << "\n";
	}
	#endif
	
	return p;
}


//----------------------------------------------------------------------
OrderParameter OrderParameter::between_sub1(const OrderParameter& 
	smaller)
{
	/* Search for an index where not all bits are set, e.g. where
	 * we can insert a number between the current number and the maximum
	 * unsigned int: */
	size_t i=0;
	OrderParameter p(smaller.first_digits);
	size_t md_size = smaller.more_digits.size();
	while (i<md_size){
		if (smaller.more_digits[i]<std::numeric_limits<digit_t>::max())
			break;
		++i;
	}
	/* We need to resize in either case: */
	p.more_digits.resize(i+1);
	
	/* Copy accordingly: */
	if (i < md_size){
		/* Copy digits until we reached element i: */
		for (size_t j=0; j<i; ++j){
			p.more_digits[j] = smaller.more_digits[j];
		}
		
		/* Instead of i, we copy the mean between smaller.more_digits[i]
		 * and the maximum unsigned int (rounded up): */
		p.more_digits[i] = smaller.more_digits[i] +
			std::max((std::numeric_limits<digit_t>::max() -
			          smaller.more_digits[i]) / 2, 1u);
	} else {
		/* We need to append: */
		for (size_t j=0; j<md_size; ++j){
			p.more_digits[j] = smaller.more_digits[j];
		}
		/* At last digit, we set 0.5: */
		p.more_digits.back() = std::numeric_limits<digit_t>::max() / 2;
	}
	
	return p;
}


//----------------------------------------------------------------------
bool is_between(const OrderParameter& m, const OrderParameter& b1,
	const OrderParameter& b2)
{
	if (b1 < b2){
		return b1 < m && m < b2;
	} else if (b2 < b1){
		return b2 < m && m < b1;
	} else {
		throw std::runtime_error("ERROR : is_between(): b1 == b2\n");
	}
}


//----------------------------------------------------------------------
OrderParameter OrderParameter::between(const OrderParameter& o1,
	const OrderParameter& o2)
{
	OrderParameter p;
	if (o1.first_digits > o2.first_digits){
		unsigned long delta = o1.first_digits - o2.first_digits;
		if (delta > 1){
			/* Easy case. Just use the first digits: */
			p.first_digits = o2.first_digits + delta/2;
		} else {
			/* Harder case. Search index of digit where we can still fit
			 * another number: */
			p = OrderParameter::between_sub1(o2);
			
			/* Debugging checks: */
			#ifdef ACOSA_DEBUG
			if (!is_between(p, o1, o2)){
				std::cerr << "OrderParameter::between_sub1  failed 1\n";
				throw 0;
			}
			#endif
		}
	} else if (o1.first_digits < o2.first_digits){
		unsigned long delta = o2.first_digits - o1.first_digits;
		if (delta > 1){
			/* Easy case. Just use the first digits: */
			p.first_digits = o1.first_digits + delta/2;
		} else {
			/* Harder case. Search index of digit where we can still fit
			 * another number: */
			p = OrderParameter::between_sub1(o1);
			
			/* Debugging checks: */
			#ifdef ACOSA_DEBUG
			if (!is_between(p, o1, o2)){
				std::cerr << "OrderParameter::between_sub1  failed 2\n";
				throw 0;
			}
			#endif
		}
	} else {
		/* Calculate the mean for two order parameters that are equal
		 * in the first digits: */
		p = OrderParameter::between_sub2(o1, o2);
		
		/* Debugging checks: */
		#ifdef ACOSA_DEBUG
		if (!is_between(p, o1, o2)){
			std::cerr.flush();
			std::cerr << "OrderParameter::between_sub2  failed\n";
			p = OrderParameter::between_sub2(o1, o2, true);
			std::cerr << "p: " << p.to_string() << "\n"
					  << "l: " << this->to_string() << "\n"
					  << "r: " << other.to_string() << "\n";
			throw 0;
		}
		#endif
	}
	
	#ifdef ACOSA_HIST
	if (p.more_digits.size() >= hist.size()){
		hist.resize(p.more_digits.size()+1, 0);
	}
	++hist[p.more_digits.size()];
	#endif
	
	return p;
}


//----------------------------------------------------------------------	
double OrderParameter::operator-(const OrderParameter& other) const
{
	/* Make it fast: Use only the first digits. */
	if (other < *this){
		return first_digits - other.first_digits;
	} else {
		return -(double)(other.first_digits - first_digits);
	}
}


//----------------------------------------------------------------------
bool OrderParameter::operator<(const OrderParameter& other) const
{
	/* Fast comparison of first digits: */
	if (first_digits != other.first_digits){
		return first_digits < other.first_digits;
	}
	
	/* Otherwise iterate over more digits: */
	size_t max_size = std::max(more_digits.size(),
	                           other.more_digits.size());
	for (size_t i=0; i<max_size; ++i){
		if (digit_at(i) != other.digit_at(i))
			return digit_at(i) < other.digit_at(i);
	}
	
	/* If we have come this far, both OrderParameters are equal.
	 * Thus, the relational comparison returns false: */
	return false;
}

//----------------------------------------------------------------------
bool OrderParameter::operator>(const OrderParameter& other) const
{
	return other < *this;
}


//----------------------------------------------------------------------
bool OrderParameter::operator==(const OrderParameter& other) const
{
	/* Fast comparison of first digits: */
	if (first_digits != other.first_digits){
		return false;
	}
	
	/* Otherwise iterate over more digits: */
	size_t max_size = std::max(more_digits.size(),
	                           other.more_digits.size());
	for (size_t i=0; i<max_size; ++i){
		if (digit_at(i) != other.digit_at(i))
			return false;
	}
	
	/* If we have come this far, both OrderParameters are equal */
	return true;
}


//----------------------------------------------------------------------
std::string OrderParameter::to_string() const
{
	std::string str = std::to_string(first_digits);
	if (more_digits.size()){
		str.append(".");
		for (size_t i=0; i<more_digits.size(); ++i){
			str.append("[");
			str.append(std::to_string(more_digits[i]));
			str.append("]");
		}
	}
	return str;
}
		
//----------------------------------------------------------------------
OrderParameter::operator double() const
{
	double d = first_digits;
	double divisor = std::numeric_limits<unsigned int>::max();
	size_t md_size = more_digits.size();
	for (size_t i=0; i<md_size; ++i){
		d += more_digits[i] / divisor;
		divisor *= std::numeric_limits<unsigned int>::max();
	}
	return d;
}

//----------------------------------------------------------------------
OrderParameter::digit_t OrderParameter::digit_at(size_t i) const
{
	return (i < more_digits.size()) ? more_digits[i] : 0;
}
