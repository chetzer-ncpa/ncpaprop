#include "units.h"
#include <map>
#include <utility>
#include <stdexcept>
#include <cstring>
#include <iostream>



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
	this->convert( &in, 1, type_in, type_out, &out );
	return out;
}


void NCPA::UnitConverter::convert( const double *in, unsigned int nSamples, 
	NCPA::UNITS_TYPE type_in, NCPA::UNITS_TYPE type_out, double *out ) {
		
	if (type_in == type_out) {
		if (in == out) {
			return;
		}
		std::memcpy( out, in, nSamples*sizeof(double) );
		return;
	}
	
	conversion_pair cpair( type_in, type_out );
	conversion_function fPtr;
	try {
		fPtr = _map.at( cpair );
		for (unsigned int i = 0; i < nSamples; i++) {
			out[ i ] = fPtr( in[ i ] );
		}
	} catch (const std::out_of_range& oor) {
		std::cerr << "Invalid conversion requested: " << oor.what() << std::endl;
		throw;
	}
	
}



conversion_pair NCPA::UnitConverter::get_unit_pair_( NCPA::UNITS_TYPE t1, NCPA::UNITS_TYPE t2 ) {
	return std::make_pair( t1, t2 );
}