#ifndef NCPAPROP_ATMOSPHERE2D_H_INCLUDED
#define NCPAPROP_ATMOSPHERE2D_H_INCLUDED

#include "Atmosphere1D.h"
#include <vector>
#include <climits>

namespace NCPA {

	bool sort_profiles_by_range_( Atmosphere1D *p1, Atmosphere1D *p2 );

	class Atmosphere2D {

	public:
		Atmosphere2D();
		~Atmosphere2D();

		// setup of profiles
		void insert_profile( const Atmosphere1D *profile, double range );
		void set_insert_range_units( units_t u );
		void sort_profiles();
		void convert_range_units( NCPA::units_t new_units );

		// data retrieval, single values
		double get( double range, std::string key );    // retrieve scalar quantity from profile
		double get( double range, std::string key, double altitude );
		double get_first_derivative( double range, std::string key, double altitude );
		double get_second_derivative( double range, std::string key, double altitude );

		// data retrieval, arrays
		size_t get_profile_index( double range );
		size_t nz( double range );
		void get_altitude_vector( double range, double *buffer, units_t *buffer_units );
		void get_property_vector( double range, std::string key, double *buffer, units_t *buffer_units );
		void get_altitude_vector( double range, double *buffer );
		void get_property_vector( double range, std::string key, double *buffer );
		units_t get_altitude_units( double range );
		units_t get_property_units( double range, std::string key );
		bool contains_scalar( double range, std::string key );
		bool contains_vector( double range, std::string key );
		bool contains_key( double range, std::string key );

		// metadata
		double get_minimum_altitude( double range );
		double get_overall_minimum_altitude();
		double get_maximum_altitude( double range );
		double get_overall_maximum_altitude();
		//double get_overall_maximum_altitude() const;

		// bulk calculations
		void calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key, 
			units_t wind_units );
		void calculate_sound_speed_from_pressure_and_density( std::string new_key, std::string pressure_key, 
			std::string density_key, units_t wind_units );
		void calculate_wind_speed( std::string new_key, std::string we_wind_speed_key, std::string sn_wind_speed_key );
		void calculate_wind_direction( std::string new_key, std::string we_wind_speed_key, std::string sn_wind_speed_key, 
			units_t direction_units = NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );
		void calculate_attenuation( std::string new_key, std::string temperature_key, std::string pressure_key, 
			std::string density_key, double freq, double tweak_factor = 1.0 );
		void calculate_wind_component( std::string new_key, std::string wind_speed_key, std::string wind_direction_key, 
			double azimuth );
		void calculate_effective_sound_speed( std::string new_key, std::string sound_speed_key, std::string wind_component_key );
		void convert_altitude_units( units_t new_units );
		void convert_property_units( std::string key, units_t new_units );
		

	protected:
		void clear_last_index_();
		void set_last_index_( size_t ind );
		void calculate_midpoints_();
		

		// data storage
		std::vector< Atmosphere1D * > profiles_;
		std::vector< double > midpoints_;

		// for keeping track of the last index calculated, speeds things up a little
		size_t last_index_;
		double last_index_min_range_, last_index_max_range_;
		units_t range_units_;
		bool sorted_;

	};


}







#endif