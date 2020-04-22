#ifndef NCPAPROP_ATMOSPHERE1D_H_INCLUDED
#define NCPAPROP_ATMOSPHERE1D_H_INCLUDED

#include "AtmosphericModel.h"
#include "AtmosphericProperty1D.h"
#include "geographic.h"
#include "units.h"
#include <string>
#include <map>
#include <stack>
#include <fstream>

namespace NCPA {

	class Atmosphere1D : public AtmosphericModel {

	public:
		Atmosphere1D( size_t n_altitude_points, double *altitude_points, units_t altitude_units );
		Atmosphere1D( std::istream& in );
		Atmosphere1D( std::string filename );
		~Atmosphere1D();

		void read_from_stream( std::istream& in );

		void add_property( std::string key, size_t n_points, double *quantity_points, units_t quantity_units );    // vector quantity
		void add_property( std::string key, double value, units_t units );    // scalar quantity
		void remove_property( std::string key );

		double get_minimum_altitude() const;
		double get_maximum_altitude() const;

		// derived quantities
		void calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key, units_t wind_units );
		void calculate_sound_speed_from_pressure_and_density( std::string new_key, std::string pressure_key, std::string density_key, units_t wind_units );
		void calculate_wind_speed( std::string new_key, std::string we_wind_speed_key, std::string sn_wind_speed_key );
		void calculate_wind_direction( std::string new_key, std::string we_wind_speed_key, std::string sn_wind_speed_key, 
			units_t direction_units = NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH );

		void calculate_wind_component( std::string new_key, std::string wind_speed_key, std::string wind_direction_key, 
			double azimuth );
		void calculate_effective_sound_speed( std::string new_key, std::string sound_speed_key, std::string wind_component_key );

		double get( std::string key ) const;    // scalars
		double get( std::string key, double altitude ) const;
		double get_first_derivative( std::string key, double altitude ) const;
		double get_second_derivative( std::string key, double altitude ) const;

		void convert_altitude_units( units_t new_units );
		//void revert_altitude_units();
		void convert_property_units( std::string key, units_t new_units );
		//void revert_property_units( std::string key );
		units_t get_property_units( std::string key );

		size_t get_basis_length() const;
		size_t nz() const;
		void get_altitude_vector( double *buffer, units_t *buffer_units ) const;
		void get_property_vector( std::string key, double *buffer, units_t *buffer_units ) const;
		void get_altitude_vector( double *buffer ) const;
		void get_property_vector( std::string key, double *buffer ) const;

		void resample( double new_dz );

		units_t get_altitude_units() const;
		units_t get_property_units( std::string key ) const;

		std::vector< std::string > get_keys() const;
		bool contains_scalar( std::string key ) const;
		bool contains_vector( std::string key ) const;

		void print_atmosphere( const std::vector< std::string >& columnorder, std::string altitude_key = "Z", std::ostream& os = std::cout );
		void print_atmosphere( std::string altitude_key = "Z", std::ostream& os = std::cout );

	protected:
		// internal storage
		std::map< std::string, NCPA::AtmosphericProperty1D * > contents_;
		std::map< std::string, NCPA::ScalarWithUnits * > scalar_contents_;
		
		//double *z_;
		//size_t nz_;
		//std::stack< units_t > z_units_;
		NCPA::VectorWithUnits *z_;

		void do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits );

		// @todo: location
	};

}

#endif
