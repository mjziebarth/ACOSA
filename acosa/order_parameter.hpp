#ifndef ACOSA_ORDER_PARAMETER_H
#define ACOSA_ORDER_PARAMETER_H

#include <string>
#include <vector>
#include <memory>

namespace ACOSA {


/*!
 * This class is a simple arbitrary precision fixed-point-number-like
 * class. It is used to immutably sort elements in the Beach class.
 * 
 * Ordering is achieved by selecting two adjacend OrderParameters o1<o2
 * between which an element should be inserted and then calling
 * between(o1,o2) or between(o2,o1) which will create an OrderParameter
 * o3 satisfying o1 < o3 < o2.
 * 
 * The OrderParameter class will try to calculate a new element o3
 * by using an unsigned long first_digits between the respective
 * member of o1 and o2. If that is not possible, more "digits" will be
 * appended to o3's arbitrary length digit vector more_digits until
 * a fixed-point-like number between o1 and o2 is created.
 * 
 * There are two preprocessor defines that enable debugging / statistics
 * features that are disabled by default:
 * - ACOSA_DEBUG: Enables further consistency checks and debug output.
 * - ACOSA_HIST:  Enables histogram statistics of OrderParameter digit
 *                length.
 */
class OrderParameter {
	public:
        #ifdef ACOSA_HIST
		static std::vector<unsigned long long> hist;
        #endif
	
		static OrderParameter max();
		static OrderParameter min();
	
		OrderParameter(unsigned long);
		
		OrderParameter();
		
		static OrderParameter between(const OrderParameter& o1,
			const OrderParameter& o2);
		
		double operator-(const OrderParameter& other) const;
		
		bool operator<(const OrderParameter& other) const;
		
		bool operator>(const OrderParameter& other) const;
		
		bool operator==(const OrderParameter& other) const;
		
		std::string to_string() const;
		
		operator double() const;
	
	private:
		typedef unsigned int digit_t;
	
	
	
	
		unsigned long first_digits;
		std::vector<digit_t> more_digits;
		
		
		digit_t digit_at(size_t i) const;
		
		
		
		
		static OrderParameter between_sub1(const OrderParameter& 
		                                   smaller);
		
		static OrderParameter between_sub2(const OrderParameter& op1,
			const OrderParameter& op2, bool verbose=false);
};

} // NAMESPACE ACOSA

#endif // ACOSA_ORDER_PARAMETER_H
