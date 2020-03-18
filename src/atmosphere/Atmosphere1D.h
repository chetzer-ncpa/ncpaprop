#ifndef NCPAPROP_ATMOSPHERE1D_H_INCLUDED
#define NCPAPROP_ATMOSPHERE1D_H_INCLUDED

#include "geographic.h"
#include "units.h"
#include <string>
#include <map>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

namespace NCPA {

	class AtmosphericProperty1D {
	private:
		double *values_, *z_;
		gsl_interp_accel *accel_;
		gsl_spline *spline_;
		units_t units_, z_units_;
		size_t nz_;

	public:
		AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units );
		~AtmosphericProperty1D();

		get( double altitude, units_t altitude_units, units_t quantity_units );
		get_first_derivative( double altitude, units_t altitude_units, units_t quantity_units );
		get_second_derivative( double altitude, units_t altitude_units, units_t quantity_units );
	};

	class Atmosphere1D : public Atmosphere {

	public:
		Atmosphere1D( size_t n_altitude_points, double *altitude_points, units_t altitude_units );
		~Atmosphere1D();

		void add_quantity( std::string key, size_t n_points, double *quantity_points, units_t quantity_units );

		double get_minimum_altitude( units_t altitude_units );
		double get_maximum_altitude( units_t altitude_units );

		void calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key );
		void calculate_sound_speed_from_pressure_and_density( std::string new_key, std::string pressure_key, std::string density_key );

		double get( std::string key, units_t quantity_units, double altitude, units_t altitude_units );
		double get_first_derivative( std::string key, units_t quantity_units, double altitude, units_t altitude_units );
		double get_second_derivative( std::string key, units_t quantity_units, double altitude, units_t altitude_units );

	protected:
		// internal storage
		std::map< std::string, NCPA::AtmosphericProperty1D * > contents_;
		double *z_;
		size_t nz_;
		units_t z_units_;

	};




}




#endif
