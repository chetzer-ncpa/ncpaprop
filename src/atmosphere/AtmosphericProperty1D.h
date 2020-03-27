#ifndef NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED
#define NCPAPROP_ATMOSPHERICPROPERTY1D_H_DEFINED

#include <map>
#include <stack>
#include "units.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"


namespace NCPA {

	class AtmosphericProperty1D {
	protected:
		double *values_, *z_;
		gsl_interp_accel *accel_ = NULL;
		gsl_spline *spline_ = NULL;
		std::stack< NCPA::units_t > units_, z_units_;
		//units_t units_, z_units_, units_last_, z_units_last;
		size_t nz_;

		void check_altitude_( double z_req ) const;
		void build_splines_();
		void delete_splines_();

		void do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits );

	public:
		AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units );
		~AtmosphericProperty1D();

		double get( double altitude ) const;
		double get_first_derivative( double altitude ) const;
		double get_second_derivative( double altitude ) const;

		void convert_altitude_units( units_t new_units );
		units_t get_altitude_units() const;
		void revert_altitude_units();
		void convert_units( units_t new_units );
		units_t get_units() const;
		void revert_units();

		//double get( double altitude, units_t altitude_units, units_t quantity_units );
		//double get_first_derivative( double altitude, units_t altitude_units, units_t quantity_units );
		//double get_second_derivative( double altitude, units_t altitude_units, units_t quantity_units );

		size_t get_basis_length() const;
		void get_altitude_basis( double *buffer, units_t *buffer_units ) const;
		void get_property_basis( double *buffer, units_t *buffer_units ) const;
	};

}

#endif