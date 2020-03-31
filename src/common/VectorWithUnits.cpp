#include "VectorWithUnits.h"
#include "units.h"
#include "util.h"
#include <cstring>
#include <stdexcept>
#include <sstream>




NCPA::VectorWithUnits::VectorWithUnits( size_t n_points, double *property_values, units_t property_units ) {

	values_ = new double[ n_points ];
	//units_ = property_units;
	//units_last_ = property_units;
	units_.push( property_units );
	std::memcpy( values_, property_values, n_points*sizeof(double) );
	n_ = n_points;
}

NCPA::VectorWithUnits::VectorWithUnits( const NCPA::VectorWithUnits &source ) {
	n_ = source.size();
	NCPA::units_t u;
	values_ = new double[ n_ ];
	source.get( values_, &u );
	units_.push( u );
}

NCPA::VectorWithUnits::~VectorWithUnits() {
	delete [] values_;
}

NCPA::units_t NCPA::VectorWithUnits::get_units() const {
	return units_.top();
}

void NCPA::VectorWithUnits::convert_units( NCPA::units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	if (new_units != units_.top()) {
		do_units_conversion_( n_, values_, units_.top(), new_units );
	}
	units_.push( new_units );
}

void NCPA::VectorWithUnits::revert_units() {
	if (units_.size() < 2) {
		return;
	}

	NCPA::units_t current_units, last_units;
	current_units = units_.top();
	units_.pop();
	last_units = units_.top();
	if (current_units != last_units) {
		try {
			do_units_conversion_( n_, values_, current_units, last_units );
		} catch (std::out_of_range &oor) {
			units_.push( current_units );
			throw;
		}
	}
}

void NCPA::VectorWithUnits::do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double *units_buffer = new double[ n_points ];
	std::memset( units_buffer, 0, n_points * sizeof( double ) );
	
	// throws out_of_range if conversion is undefined
	NCPA::Units::convert( inplace, n_points, fromUnits, toUnits, units_buffer );

	// successful, so record the units change
	std::memcpy( inplace, units_buffer, n_points * sizeof( double ) );
	delete [] units_buffer;
}


size_t NCPA::VectorWithUnits::size() const {
	return n_;
}

void NCPA::VectorWithUnits::get( double *buffer, units_t *buffer_units ) const {
	*buffer_units = units_.top();
	std::memcpy( buffer, values_, n_ * sizeof( double ) );
}

