/*
Functions and constants for unit conversion.
*/

#ifndef NCPA_UNITS_H__
#define NCPA_UNITS_H__

namespace NCPA {
	
	enum UNITS_TEMPERATURE : unsigned int {
		UNITS_TEMPERATURE_NONE = 0,
		UNITS_TEMPERATURE_K,
		UNITS_TEMPERATURE_C,
		UNITS_TEMPERATURE_F
	};
	
	enum UNITS_DISTANCE : unsigned int {
		UNITS_DISTANCE_NONE = 0,
		UNITS_DISTANCE_M,
		UNITS_DISTANCE_KM
	};
	
	enum UNITS_SPEED : unsigned int {
		UNITS_SPEED_NONE = 0,
		UNITS_SPEED_MPS,
		UNITS_SPEED_KMPS
	};
	
	enum UNITS_PRESSURE : unsigned int {
		UNITS_PRESSURE_NONE = 0,
		UNITS_PRESSURE_PA,
		UNITS_PRESSURE_MBAR
	};
	
	enum UNITS_DENSITY : unsigned int {
		UNITS_DENSITY_NONE = 0,
		UNITS_DENSITY_KGPM3,
		UNITS_DENSITY_GPCM3
	};
	
	/*
	 * Convert a number from one unit to another.  Specifying the same units for both
	 * input and output will return the original number with no processing, otherwise
	 * numbers are internally converted to an intermediate SI unit and then reconverted
	 * to the specified output unit.
	 */
	double convert_units( double input, UNITS_TEMPERATURE units_in, UNITS_TEMPERATURE units_out );
	double convert_units( double input, UNITS_DISTANCE units_in, UNITS_DISTANCE units_out );
	double convert_units( double input, UNITS_SPEED units_in, UNITS_SPEED units_out );
	double convert_units( double input, UNITS_PRESSURE units_in, UNITS_PRESSURE units_out );
	double convert_units( double input, UNITS_DENSITY units_in, UNITS_DENSITY units_out );
	
}



#endif