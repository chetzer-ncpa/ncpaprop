#ifndef NCPA_UNITS_H__
#define NCPA_UNITS_H__

namespace NCPA {
	
	enum TEMPERATURE_UNITS : unsigned int {
		TEMPERATURE_UNITS_K,
		TEMPERATURE_UNITS_C,
		TEMPERATURE_UNITS_F
	};
	
	enum DISTANCE_UNITS : unsigned int {
		DISTANCE_UNITS_M,
		DISTANCE_UNITS_KM
	};
	
	enum SPEED_UNITS : unsigned int {
		SPEED_UNITS_MPS,
		SPEED_UNITS_KMPS
	};
	
	enum PRESSURE_UNITS : unsigned int {
		PRESSURE_UNITS_PA,
		PRESSURE_UNITS_MBAR
	};
	
	enum DENSITY_UNITS : unsigned int {
		DENSITY_UNITS_KGPM3,
		DENSITY_UNITS_GPCM3
	};
	
	double convert_units( double input, TEMPERATURE_UNITS units_in, TEMPERATURE_UNITS units_out );
	double convert_units( double input, DISTANCE_UNITS units_in, DISTANCE_UNITS units_out );
	double convert_units( double input, SPEED_UNITS units_in, SPEED_UNITS units_out );
	double convert_units( double input, PRESSURE_UNITS units_in, PRESSURE_UNITS units_out );
	double convert_units( double input, DENSITY_UNITS units_in, DENSITY_UNITS units_out );
	
}



#endif