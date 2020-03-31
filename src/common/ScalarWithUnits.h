#ifndef NCPA_SCALARWITHUNITS_H_DEFINED
#define NCPA_SCALARWITHUNITS_H_DEFINED

#include <map>
#include <stack>
#include <iostream>
#include "units.h"


namespace NCPA {

	class ScalarWithUnits {
	protected:
		double value_;
		std::stack< NCPA::units_t > units_;
		
		void do_units_conversion_( NCPA::units_t fromUnits, NCPA::units_t toUnits );

	public:
		ScalarWithUnits( double value, units_t property_units );
		~ScalarWithUnits();

		double get() const;

		void convert_units( units_t new_units );
		units_t get_units() const;
		void revert_units();

		
	};
	std::ostream &operator<<( std::ostream &output, const ScalarWithUnits &D );

}

#endif