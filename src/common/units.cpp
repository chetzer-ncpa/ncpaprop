#include "units.h"
#include <map>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <iostream>

conversion_map NCPA::Units::_map = conversion_map();

std::string NCPA::toString( UNITS_TYPE u ) {
	switch (u) {
		case UNITS_NONE:
			return "No units";
		case UNITS_TEMPERATURE_KELVIN:
			return "degrees Kelvin";
		case UNITS_TEMPERATURE_CELSIUS:
			return "degrees Celsius";
		case UNITS_TEMPERATURE_FAHRENHEIT:
			return "degrees Fahrenheit";
		case UNITS_DISTANCE_METERS:
			return "meters";
		case UNITS_DISTANCE_KILOMETERS:
			return "kilometers";
		case UNITS_SPEED_METERS_PER_SECOND:
			return "meters per second";
		case UNITS_SPEED_KILOMETERS_PER_SECOND:
			return "kilometers per second";
		case UNITS_PRESSURE_PASCALS:
			return "Pascals";
		case UNITS_PRESSURE_MILLIBARS:
			return "millibars";
		case UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER:
			return "kilograms per cubic meter";
		case UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER:
			return "grams per cubic centimeter";
		case UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH:
			return "degrees clockwise from North";
		case UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST:
			return "degrees counterclockwise from East";
		case UNITS_ANGLE_DEGREES:
			return "degrees";
		case UNITS_ANGLE_RADIANS:
			return "radians";
		default:
			throw std::out_of_range( "Unrecognized units type" );
	}
}


std::string NCPA::toStr( UNITS_TYPE u ) {
	switch (u) {
		case UNITS_NONE:
			return "N/A";
		case UNITS_TEMPERATURE_KELVIN:
			return "K";
		case UNITS_TEMPERATURE_CELSIUS:
			return "C";
		case UNITS_TEMPERATURE_FAHRENHEIT:
			return "F";
		case UNITS_DISTANCE_METERS:
			return "m";
		case UNITS_DISTANCE_KILOMETERS:
			return "km";
		case UNITS_SPEED_METERS_PER_SECOND:
			return "m/s";
		case UNITS_SPEED_KILOMETERS_PER_SECOND:
			return "km/s";
		case UNITS_PRESSURE_PASCALS:
			return "Pa";
		case UNITS_PRESSURE_MILLIBARS:
			return "mbar";
		case UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER:
			return "kg/m3";
		case UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER:
			return "g/cm3";
		case UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH:
			return "deg CW from N";
		case UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST:
			return "deg CCW from E";
		case UNITS_ANGLE_DEGREES:
			return "deg";
		case UNITS_ANGLE_RADIANS:
			return "rad";
		default:
			throw std::out_of_range( "Unrecognized units type" );
	}
}


/* 
Constructor for the UnitConverter class.
Internally the class uses a std::map to associate a 
std::pair< UNITS_TYPE, UNITS_TYPE > key with a function pointer that performs the
indicated conversion.  The constructor populates this map with all defined
conversions.  No memory is dynamically allocated within the constructor so no
explicit destructor is required.
*/
void NCPA::Units::initialize_() {
	
	_map.clear();
	
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
	
	// direction conversion
	_map[ get_unit_pair_( UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST ) ] 
		= convert_direction_geo_to_math_;
	_map[ get_unit_pair_( UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST, UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH ) ] 
		= convert_direction_math_to_geo_;
}

bool NCPA::Units::ready_() {
	return !(_map.empty());
}

double NCPA::Units::convert( double in, NCPA::UNITS_TYPE type_in, NCPA::UNITS_TYPE type_out ) {
	double out = 0.0;
	try {
		// Call the vector conversion with one-element vectors
		NCPA::Units::convert( &in, 1, type_in, type_out, &out );
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
void NCPA::Units::convert( const double *in, unsigned int nSamples, 
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
	
	if ( ! NCPA::Units::ready_() ) {
		NCPA::Units::initialize_();
	}
	
	// Create a pair with the requested in and out units
	conversion_pair cpair( type_in, type_out );
	
	// function will go here if found
	conversion_function fPtr;
	try {
		// see if the requested conversion exists in the map.  If not, it throws
		// an out_of_range exception that is caught later
		fPtr = NCPA::Units::_map.at( cpair );
		
		// perform the requested conversion
		for (unsigned int i = 0; i < nSamples; i++) {
			out[ i ] = fPtr( in[ i ] );
		}
		
	} catch (const std::out_of_range& oor) {
		// didn't find the conversion, notify the user and kick it upstairs
		throw std::out_of_range( "Undefined conversion requested from " 
			+ NCPA::toString( type_in ) + " to " + NCPA::toString( type_out ) );
	}
	
}


/*
Generates and returns a std::pair with the two unit types, for use as a map key.
*/
conversion_pair NCPA::Units::get_unit_pair_( NCPA::UNITS_TYPE t1, NCPA::UNITS_TYPE t2 ) {
	return std::make_pair( t1, t2 );
}