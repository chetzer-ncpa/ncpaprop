#ifndef NCPA_VECTORWITHUNITS_H_INCLUDED
#define NCPA_VECTORWITHUNITS_H_INCLUDED

#include <map>
#include <stack>
#include "units.h"


namespace NCPA {

	class VectorWithUnits {
	protected:
		double *values_;
		std::stack< NCPA::units_t > units_;
		size_t n_;

		void do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits );

	public:
		VectorWithUnits();
		VectorWithUnits( size_t n_points, double *values, units_t units );
		VectorWithUnits( const VectorWithUnits &source );
		~VectorWithUnits();

		void convert_units( units_t new_units );
		units_t get_units() const;
		void revert_units();

		size_t size() const;
		void get_vector( double *buffer, units_t *buffer_units ) const;
	};

}

#endif