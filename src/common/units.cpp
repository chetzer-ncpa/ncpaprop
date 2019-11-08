#include "units.h"
#include <stdexcept>


// convert temperature units.
// Note: to expand for new units, add case statements for input->degC and degC->output
double NCPA::convert_units( double input, NCPA::TEMPERATURE_UNITS units_in, NCPA::TEMPERATURE_UNITS units_out ) {
	
	// check for identity conversion
	if (units_in == units_out) {
		return input;
	}
	
	double intermediate, output;
	
	// intermediate temperature units: degrees C
	switch (units_in) {
		case NCPA::TEMPERATURE_UNITS_C:
			intermediate = input;
			break;
		case NCPA::TEMPERATURE_UNITS_K:
			intermediate = input - 273.15;
			break;
		case NCPA::TEMPERATURE_UNITS_F:
			intermediate = ( input - 32.0 ) * 5.0 / 9.0;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented input temperature units" );
	}
	
	// now convert to proper output units
	switch (units_out) {
		case NCPA::TEMPERATURE_UNITS_C:
			output = intermediate;
			break;
		case NCPA::TEMPERATURE_UNITS_K:
			output = intermediate + 273.15;
			break;
		case NCPA::TEMPERATURE_UNITS_F:
			output = ( intermediate * 9.0 / 5.0 ) + 32.0;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented output temperature units" );
	}
	
	return output;
}



// convert distance units.
// Note: to expand for new units, add case statements for input->m and m->output
double NCPA::convert_units( double input, NCPA::DISTANCE_UNITS units_in, NCPA::DISTANCE_UNITS units_out ) {
	
	// check for identity conversion
	if (units_in == units_out) {
		return input;
	}
	
	double intermediate, output;
	
	// intermediate distance units: meters
	switch (units_in) {
		case NCPA::DISTANCE_UNITS_M:
			intermediate = input;
			break;
		case NCPA::DISTANCE_UNITS_KM:
			intermediate = input * 1000.0;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented input distance units" );
	}
	
	// now convert to proper output units
	switch (units_out) {
		case NCPA::DISTANCE_UNITS_M:
			output = intermediate;
			break;
		case NCPA::DISTANCE_UNITS_KM:
			output = intermediate * 0.001;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented output distance units" );
	}
	
	return output;
}


// convert speed units.
// Note: to expand for new units, add case statements for input->m/s and m/s->output
double NCPA::convert_units( double input, NCPA::SPEED_UNITS units_in, NCPA::SPEED_UNITS units_out ) {
	
	// check for identity conversion
	if (units_in == units_out) {
		return input;
	}
	
	double intermediate, output;
	
	// intermediate speed units: meters/sec
	switch (units_in) {
		case NCPA::SPEED_UNITS_MPS:
			intermediate = input;
			break;
		case NCPA::SPEED_UNITS_KMPS:
			intermediate = input * 1000.0;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented input speed units" );
	}
	
	// now convert to proper output units
	switch (units_out) {
		case NCPA::SPEED_UNITS_MPS:
			output = intermediate;
			break;
		case NCPA::SPEED_UNITS_KMPS:
			output = intermediate * 0.001;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented output speed units" );
	}
	
	return output;
}


// convert pressure units.
// Note: to expand for new units, add case statements for input->Pa and Pa->output
double NCPA::convert_units( double input, NCPA::PRESSURE_UNITS units_in, NCPA::PRESSURE_UNITS units_out ) {
	
	// check for identity conversion
	if (units_in == units_out) {
		return input;
	}
	
	double intermediate, output;
	
	// intermediate pressure units: Pascals
	switch (units_in) {
		case NCPA::PRESSURE_UNITS_PA:
			intermediate = input;
			break;
		case NCPA::PRESSURE_UNITS_MBAR:
			intermediate = input * 100.0;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented input pressure units" );
	}
	
	// now convert to proper output units
	switch (units_out) {
		case NCPA::PRESSURE_UNITS_PA:
			output = intermediate;
			break;
		case NCPA::PRESSURE_UNITS_MBAR:
			output = intermediate * 0.01;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented output pressure units" );
	}
	
	return output;
}



// convert density units.
// Note: to expand for new units, add case statements for input->kg/m3 and kg/m3->output
double NCPA::convert_units( double input, NCPA::DENSITY_UNITS units_in, NCPA::DENSITY_UNITS units_out ) {
	
	// check for identity conversion
	if (units_in == units_out) {
		return input;
	}
	
	double intermediate, output;
	
	// intermediate density units: kg/m3
	switch (units_in) {
		case NCPA::DENSITY_UNITS_KGPM3:
			intermediate = input;
			break;
		case NCPA::DENSITY_UNITS_GPCM3:
			intermediate = input * 1000.0;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented input density units" );
	}
	
	// now convert to proper output units
	switch (units_out) {
		case NCPA::DENSITY_UNITS_KGPM3:
			output = intermediate;
			break;
		case NCPA::DENSITY_UNITS_GPCM3:
			output = intermediate * 0.001;
			break;
		default:
			throw std::invalid_argument( "Unknown or unimplemented output density units" );
	}
	
	return output;
}
