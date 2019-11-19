#include "units.h"
#include <map>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <iostream>


/* 
Constructor for the UnitConverter class.
Internally the class uses a std::map to associate a 
std::pair< UNITS_TYPE, UNITS_TYPE > key with a function pointer that performs the
indicated conversion.  The constructor populates this map with all defined
conversions.  No memory is dynamically allocated within the constructor so no
explicit destructor is required.
*/
NCPA::UnitConverter::UnitConverter() {
	
	// set up valid unit pairs: temperature conversion
	_map[ get_unit_pair_( UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_FAHRENHEIT ) ] 
		= convert_temperature_c_to_f_;
	_map[ get_unit_pair_( UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_KELVIN ) ] 
		= convert_temperature_c_to_k_;
	_map[ get_unit_pair_( UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_FAHRENHEIT ) ] 
		= convert_temperature_k_to_f_;
	_map[ get_unit_pair_( UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_CELSIUS ) ]
		= convert_temperature_k_to_c_;
	_map[ get_unit_pair_( UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_KELVIN ) ] 
		= convert_temperature_f_to_k_;
	_map[ get_unit_pair_( UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_CELSIUS ) ]
		= convert_temperature_f_to_c_;
	
	// length conversion
	_map[ get_unit_pair_( UNITS_DISTANCE_METERS, UNITS_DISTANCE_KILOMETERS ) ] 
		= convert_distance_m_to_km_;
	_map[ get_unit_pair_( UNITS_DISTANCE_KILOMETERS, UNITS_DISTANCE_METERS ) ] 
		= convert_distance_km_to_m_;
	
	// speed conversion
	_map[ get_unit_pair_( UNITS_SPEED_METERS_PER_SECOND, UNITS_SPEED_KILOMETERS_PER_SECOND ) ] 
		= convert_speed_mps_to_kmps_;
	_map[ get_unit_pair_( UNITS_SPEED_KILOMETERS_PER_SECOND, UNITS_SPEED_METERS_PER_SECOND ) ] 
		= convert_speed_kmps_to_mps_;
	
	// pressure conversion
	_map[ get_unit_pair_( UNITS_PRESSURE_PASCALS, UNITS_PRESSURE_MILLIBARS ) ] 
		= convert_pressure_pa_to_mbar_;
	_map[ get_unit_pair_( UNITS_PRESSURE_MILLIBARS, UNITS_PRESSURE_PASCALS ) ] 
		= convert_pressure_mbar_to_pa_;
	
	// density conversion
	_map[ get_unit_pair_( UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER, UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER ) ] 
		= convert_density_kgpm3_to_gpcm3_;
	_map[ get_unit_pair_( UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER, UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER ) ] 
		= convert_density_gpcm3_to_kgpm3_;
}


double NCPA::UnitConverter::convert( double in, NCPA::UNITS_TYPE type_in, NCPA::UNITS_TYPE type_out ) {
	double out = 0.0;
	try {
		// Call the vector conversion with one-element vectors
		this->convert( &in, 1, type_in, type_out, &out );
	} catch (const std::out_of_range& oor) {
		// didn't find the requested conversion, kick it upstairs, user
		// will have been notified in the called function
		throw;
	}
	return out;
}

/*
Converts one or more double values from one unit to another.
*/
void NCPA::UnitConverter::convert( const double *in, unsigned int nSamples, 
	NCPA::UNITS_TYPE type_in, NCPA::UNITS_TYPE type_out, double *out ) {
		
	// Is the conversion from one type to the same type?
	if (type_in == type_out) {
		
		// If the identity conversion is in place, just return
		if (in == out) {
			return;
		}
		
		// Copy the old values to the new values with no conversion
		std::memcpy( out, in, nSamples*sizeof(double) );
		return;
	}
	
	// Create a pair with the requested in and out units
	conversion_pair cpair( type_in, type_out );
	
	// function will go here if found
	conversion_function fPtr;
	try {
		// see if the requested conversion exists in the map.  If not, it throws
		// an out_of_range exception that is caught later
		fPtr = _map.at( cpair );
		
		// perform the requested conversion
		for (unsigned int i = 0; i < nSamples; i++) {
			out[ i ] = fPtr( in[ i ] );
		}
		
	} catch (const std::out_of_range& oor) {
		// didn't find the conversion, notify the user and kick it upstairs
		std::cerr << "Invalid conversion requested: " << oor.what() << std::endl;
		throw;
	}
	
}


/*
Generates and returns a std::pair with the two unit types, for use as a map key.
*/
conversion_pair NCPA::UnitConverter::get_unit_pair_( NCPA::UNITS_TYPE t1, NCPA::UNITS_TYPE t2 ) {
	return std::make_pair( t1, t2 );
}