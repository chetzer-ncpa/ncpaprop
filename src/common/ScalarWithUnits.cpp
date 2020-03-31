#include "ScalarWithUnits.h"
#include "units.h"
#include "util.h"
#include <cstring>
#include <stdexcept>
#include <sstream>
#include <iostream>




NCPA::ScalarWithUnits::ScalarWithUnits( double value, units_t units ) {
	value_ = value;
	units_.push( units );
}

NCPA::ScalarWithUnits::~ScalarWithUnits() { }


NCPA::units_t NCPA::ScalarWithUnits::get_units() const {
	return units_.top();
}

void NCPA::ScalarWithUnits::convert_units( NCPA::units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	if (new_units != units_.top()) {
		do_units_conversion_( units_.top(), new_units );
	}
	units_.push( new_units );
}

void NCPA::ScalarWithUnits::revert_units() {
	if (units_.size() < 2) {
		return;
	}

	NCPA::units_t current_units, last_units;
	current_units = units_.top();
	units_.pop();
	last_units = units_.top();
	if (current_units != last_units) {
		try {
			do_units_conversion_( current_units, last_units );
		} catch (std::out_of_range &oor) {
			units_.push( current_units );
			throw;
		}
	}
}

void NCPA::ScalarWithUnits::do_units_conversion_( NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double units_buffer = 0.0;
	
	// throws out_of_range if conversion is undefined
	units_buffer = NCPA::Units::convert( value_, fromUnits, toUnits );

	// successful, so record the units change
	value_ = units_buffer;
}

double NCPA::ScalarWithUnits::get() const {
	return value_;
}

std::ostream &NCPA::operator<<( std::ostream &output, const NCPA::ScalarWithUnits &D ) { 
	output << D.get() << " " << NCPA::Units::toStr( D.get_units() );
	return output;            
}