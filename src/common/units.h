/*
Classes, functions and constants for unit conversion.

To add a unit and its associated conversions, the following should be done:
1. Add a symbol to the UNITS_TYPE enum
2. Add all appropriate inline static double conversion functions to the 
   protected section of the UnitConverter class
3. Map the conversion functions to the appropriate unit pairs in the body 
   of the UnitConverter constructor
*/

#ifndef NCPA_UNITS_H__
#define NCPA_UNITS_H__

#include <map>
#include <utility>


namespace NCPA {
	
	/**
	 * An enum of unit constants.
	 * Constants that can be used to identify or specify units.
	 */
	enum UNITS_TYPE : unsigned int {
		UNITS_NONE = 0,					/**< Indicates no units */
		
		UNITS_TEMPERATURE_KELVIN,			/**< Temperature in Kelvin */
		UNITS_TEMPERATURE_CELSIUS,			/**< Temperature in Celsius */
		UNITS_TEMPERATURE_FAHRENHEIT,			/**< Temperature in Fahrenheit */
		
		UNITS_DISTANCE_METERS,				/**< Distance in meters */
		UNITS_DISTANCE_KILOMETERS,			/**< Distance in kilometers */
		
		UNITS_SPEED_METERS_PER_SECOND,			/**< Speed in m/s */
		UNITS_SPEED_KILOMETERS_PER_SECOND,		/**< Speed in km/s */
		
		UNITS_PRESSURE_PASCALS,				/**< Pressure in Pa */
		UNITS_PRESSURE_MILLIBARS,			/**< Pressure in mbar */
		
		UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER,	/**< Density in kg/m^3 */
		UNITS_DENSITY_GRAMS_PER_CUBIC_CENTIMETER	/**< Density in g/cm^3 */
	};
}


typedef std::pair< NCPA::UNITS_TYPE, NCPA::UNITS_TYPE > conversion_pair;
typedef double (*conversion_function)(double);
typedef std::map< conversion_pair, conversion_function > conversion_map;

namespace NCPA {
	
	/**
	 * A class for converting units.
	 */
	class UnitConverter {
	
	public:
		/**
		 * Default constructor.
		 * Creates a new UnitConverter object and populates the conversions.
		 */
		UnitConverter();
		
		/**
		 * Convert an array of numbers from one unit to another.
		 * @param in 		A pointer to an array of double values
		 * @param nSamples 	The number of consecutive samples to convert
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @param out		A pointer to a preallocated array to place the converted values in
		 * @see UNITS_TYPE
		 */
		void convert( const double *in, unsigned int nSamples, 
			UNITS_TYPE type_in, UNITS_TYPE type_out, double *out );
	
		/**
		 * Convert a single double value from one unit to another.
		 * @param in		A double value to convert.
		 * @param type_in	The units to convert from
		 * @param type_out	The units to convert to
		 * @return 		The converted value
		 */
		double convert( double in, UNITS_TYPE type_in, UNITS_TYPE type_out );
		
	
	protected:
		
		conversion_map _map; 

		conversion_pair get_unit_pair_( UNITS_TYPE t1, UNITS_TYPE t2 );

		// These functions are for specific unit conversions 
		inline static double convert_temperature_f_to_c_( double in ) { return ( in - 32.0 ) * 5.0 / 9.0; }
		inline static double convert_temperature_f_to_k_( double in ) { 
			return convert_temperature_c_to_k_( convert_temperature_f_to_c_( in ) ); 
		}
		inline static double convert_temperature_c_to_f_( double in ) { return ( in * 9.0 / 5.0 ) + 32.0; }
		inline static double convert_temperature_c_to_k_( double in ) { return in + 273.15; }
		inline static double convert_temperature_k_to_c_( double in ) { return in - 273.15; }
		inline static double convert_temperature_k_to_f_( double in ) {
			return convert_temperature_c_to_f_( convert_temperature_k_to_c_( in ) );
		}
	
		inline static double convert_distance_m_to_km_( double in ) { return in * 0.001; }
		inline static double convert_distance_km_to_m_( double in ) { return in * 1000.0; }
	
		inline static double convert_speed_mps_to_kmps_( double in ) { return in * 0.001; }
		inline static double convert_speed_kmps_to_mps_( double in ) { return in * 1000.0; }
	
		inline static double convert_pressure_pa_to_mbar_( double in ) { return in * 100.0; }
		inline static double convert_pressure_mbar_to_pa_( double in ) { return in * 0.01; }
	
		inline static double convert_density_kgpm3_to_gpcm3_( double in ) { return in * 0.001; }
		inline static double convert_density_gpcm3_to_kgpm3_( double in ) { return in * 1000.0; }
	
		inline static double convert_no_conversion( double in ) { return in; }
	};

}



#endif