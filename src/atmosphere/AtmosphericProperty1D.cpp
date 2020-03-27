#include "AtmosphericProperty1D.h"
#include "units.h"
#include "util.h"
#include <cstring>
#include <stdexcept>
#include <sstream>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"




NCPA::AtmosphericProperty1D::AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units ) {

	z_ = new double[ n_points ];
	//z_units_ = altitude_units;
	//z_units_last_ = altitude_units;
	z_units_.push( altitude_units );
	std::memcpy( z_, altitude_points, n_points*sizeof(double) );
	values_ = new double[ n_points ];
	//units_ = property_units;
	//units_last_ = property_units;
	units_.push( property_units );
	std::memcpy( values_, property_values, n_points*sizeof(double) );
	nz_ = n_points;

	build_splines_();
}

NCPA::AtmosphericProperty1D::~AtmosphericProperty1D() {
	delete_splines_();
	delete [] z_;
	delete [] values_;
}

NCPA::units_t NCPA::AtmosphericProperty1D::get_altitude_units() const {
	return z_units_.top();
}

void NCPA::AtmosphericProperty1D::convert_altitude_units( NCPA::units_t new_units ) {

	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	if (new_units != z_units_.top()) {
		do_units_conversion_( nz_, z_, z_units_.top(), new_units );
	}
	z_units_.push( new_units );
}

void NCPA::AtmosphericProperty1D::revert_altitude_units() {
	if (z_units_.size() < 2) {
		return;
	}

	NCPA::units_t current_units, last_units;
	current_units = z_units_.top();
	z_units_.pop();
	last_units = z_units_.top();
	if (current_units != last_units) {
		try {
			do_units_conversion_( nz_, z_, current_units, last_units );
		} catch (std::out_of_range &oor) {
			z_units_.push( current_units );
			throw;
		}
	}
}

NCPA::units_t NCPA::AtmosphericProperty1D::get_units() const {
	return units_.top();
}

void NCPA::AtmosphericProperty1D::convert_units( NCPA::units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	if (new_units != z_units_.top()) {
		do_units_conversion_( nz_, values_, units_.top(), new_units );
	}
	units_.push( new_units );
}

void NCPA::AtmosphericProperty1D::revert_units() {
	if (units_.size() < 2) {
		return;
	}

	NCPA::units_t current_units, last_units;
	current_units = units_.top();
	units_.pop();
	last_units = units_.top();
	if (current_units != last_units) {
		try {
			do_units_conversion_( nz_, values_, current_units, last_units );
		} catch (std::out_of_range &oor) {
			units_.push( current_units );
			throw;
		}
	}
}

void NCPA::AtmosphericProperty1D::do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double *units_buffer = new double[ n_points ];
	std::memset( units_buffer, 0, n_points * sizeof( double ) );
	
	// throws out_of_range if conversion is undefined
	NCPA::Units::convert( inplace, n_points, fromUnits, toUnits, units_buffer );

	// successful, so record the units change
	std::memcpy( inplace, units_buffer, n_points * sizeof( double ) );
	build_splines_();
	delete [] units_buffer;
}


void NCPA::AtmosphericProperty1D::build_splines_() {
	delete_splines_();
	// construct spline
	accel_ = gsl_interp_accel_alloc();
	spline_ = gsl_spline_alloc( gsl_interp_cspline, nz_ );
	gsl_spline_init( spline_, z_, values_, nz_ );
}

void NCPA::AtmosphericProperty1D::delete_splines_() {
	if (spline_ != NULL) {
		delete spline_;
		spline_ = NULL;
	}
	if (accel_ != NULL) {
		delete accel_;
		accel_ = NULL;
	}
}

size_t NCPA::AtmosphericProperty1D::get_basis_length() const {
	return nz_;
}

void NCPA::AtmosphericProperty1D::get_altitude_basis( double *buffer, units_t *buffer_units ) const {
	*buffer_units = z_units_.top();
	std::memcpy( buffer, z_, nz_ * sizeof(double) );
}

void NCPA::AtmosphericProperty1D::get_property_basis( double *buffer, units_t *buffer_units ) const {
	*buffer_units = units_.top();
	std::memcpy( buffer, values_, nz_ * sizeof( double ) );
}

void NCPA::AtmosphericProperty1D::check_altitude_( double z_req ) const {
	if ( z_req < z_[0] || z_req > z_[ nz_ - 1 ] ) {
		std::ostringstream oss;
		oss << "Requested altitude " << z_req << " " << NCPA::Units::toStr( z_units_.top() ) << " outside profile bounds.";
		throw std::range_error( oss.str() );
	}
}

double NCPA::AtmosphericProperty1D::get( double z_req ) const {
	check_altitude_( z_req );
	return gsl_spline_eval( spline_, z_req, accel_ );
}

double NCPA::AtmosphericProperty1D::get_first_derivative( double z_req ) const {
	check_altitude_( z_req );
	return gsl_spline_eval_deriv( spline_, z_req, accel_ );
}

double NCPA::AtmosphericProperty1D::get_second_derivative( double z_req ) const {
	check_altitude_( z_req );
	return gsl_spline_eval_deriv2( spline_, z_req, accel_ );
}






// double NCPA::AtmosphericProperty1D::get( double altitude, units_t altitude_units, units_t quantity_units ) {

// 	// make sure units are consistent
// 	double z_req = NCPA::Units::convert( altitude, altitude_units, z_units_ );
// 	if ( z_req < z_[0] || z_req > z_[ nz_ - 1 ] ) {
// 		std::ostringstream oss;
// 		oss << "Requested altitude " << altitude << " " << NCPA::Units::toStr( altitude_units ) << " outside profile bounds.";
// 		throw std::range_error( oss.str() );
// 	}

// 	return NCPA::Units::convert( gsl_spline_eval( spline_, z_req, accel_ ), units_, quantity_units );
// }

// double NCPA::AtmosphericProperty1D::get_first_derivative( double altitude, units_t altitude_units, units_t quantity_units ) {

// 	// make sure units are consistent
// 	double z_req = NCPA::Units::convert( altitude, altitude_units, z_units_ );
// 	if ( z_req < z_[0] || z_req > z_[ nz_ - 1 ] ) {
// 		std::ostringstream oss;
// 		oss << "Requested altitude " << altitude << " " << NCPA::Units::toStr( altitude_units ) << " outside profile bounds.";
// 		throw std::range_error( oss.str() );
// 	}


// 	return NCPA::Units::convert_first_derivative( gsl_spline_eval_deriv( spline_, z_req, accel_ ), units_, quantity_units );
// }

// double NCPA::AtmosphericProperty1D::get_second_derivative( double altitude, units_t altitude_units, units_t quantity_units ) {

// 	// make sure units are consistent
// 	double z_req = NCPA::Units::convert( altitude, altitude_units, z_units_ );
// 	if ( z_req < z_[0] || z_req > z_[ nz_ - 1 ] ) {
// 		std::ostringstream oss;
// 		oss << "Requested altitude " << altitude << " " << NCPA::Units::toStr( altitude_units ) << " outside profile bounds.";
// 		throw std::range_error( oss.str() );
// 	}


// 	return NCPA::Units::convert_second_derivative( gsl_spline_eval_deriv2( spline_, z_req, accel_ ), units_, quantity_units );
// }

