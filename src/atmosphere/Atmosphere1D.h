#ifndef NCPAPROP_ATMOSPHERE1D_H_INCLUDED
#define NCPAPROP_ATMOSPHERE1D_H_INCLUDED

#include "Atmosphere.h"
#include "AtmosphericProperty1D.h"
#include "geographic.h"
#include "units.h"
#include <string>
#include <map>
#include <stack>

namespace NCPA {

	class Atmosphere1D : public Atmosphere {

	public:
		Atmosphere1D( size_t n_altitude_points, double *altitude_points, units_t altitude_units );
		~Atmosphere1D();

		void add_quantity( std::string key, size_t n_points, double *quantity_points, units_t quantity_units );

		double get_minimum_altitude() const;
		double get_maximum_altitude() const;

		// double get_minimum_altitude( units_t altitude_units ) const;
		// double get_maximum_altitude( units_t altitude_units ) const;

		void calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key );
		void calculate_sound_speed_from_pressure_and_density( std::string new_key, std::string pressure_key, std::string density_key );

		double get( std::string key, double altitude ) const;
		double get_first_derivative( std::string key, double altitude ) const;
		double get_second_derivative( std::string key, double altitude ) const;

		void convert_altitude_units( units_t new_units );
		void revert_altitude_units();
		void convert_property_units( std::string key, units_t new_units );
		void revert_property_units( std::string key );

		// double get( std::string key, units_t quantity_units, double altitude, units_t altitude_units ) const;
		// double get_first_derivative( std::string key, units_t quantity_units, double altitude, units_t altitude_units ) const;
		// double get_second_derivative( std::string key, units_t quantity_units, double altitude, units_t altitude_units ) const;

		size_t get_basis_length() const;
		void get_altitude_basis( double *buffer, units_t *buffer_units ) const;
		void get_property_basis( std::string key, double *buffer, units_t *buffer_units ) const;

		units_t get_altitude_units() const;
		units_t get_property_units( std::string key ) const;

		std::vector< std::string > get_keys() const;

	protected:
		// internal storage
		std::map< std::string, NCPA::AtmosphericProperty1D * > contents_;
		double *z_;
		size_t nz_;
		std::stack< units_t > z_units_;

		void do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits );

		// @todo: location

	};

}




#endif
