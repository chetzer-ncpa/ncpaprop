#include "Atmosphere1D.h"
#include <map>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <sstream>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#define GAMMA_FOR_C 1.4
#define R_FOR_C 287.0



NCPA::Atmosphere1D::Atmosphere1D( size_t n_altitude_points, double *altitude_points, units_t altitude_units ) {
	contents_.clear();
	nz_ = n_altitude_points;
	z_ = new double[ nz_ ];
	std::memcpy( z_, altitude_points, nz_ );
	z_units_ = altitude_units;
}

NCPA::Atmosphere1D::~Atmosphere1D() {
	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it=contents_.begin(); it != contents_.end(); ++it) {
		delete it->second;
	}
	contents_.clear();
}


double NCPA::Atmosphere1D::get_minimum_altitude( units_t altitude_units ) {
	return NCPA::Units::convert( z_[0], z_units_, altitude_units );
}

double NCPA::Atmosphere1D::get_maximum_altitude( units_t altitude_units ) {
	return NCPA::Units::convert( z_[ nz_-1 ], z_units_, altitude_units );
}

void NCPA::Atmosphere1D::add_quantity( std::string key, size_t n_points, double *quantity_points, units_t quantity_units ) {

	// see if we already have one with that key
	NCPA::AtmosphericProperty1D *prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		prop = contents_.at( key );
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	if (n_points != nz_) {
		std::ostringstream oss;
		oss << "Array length " << n_points << " for key " << key << " does not match number of altitude points " << nz_;
		throw std::runtime_error( oss.str() );
	}

	prop = new AtmosphericProperty1D( n_points, z_, z_units_, quantity_points, quantity_units );
	contents_[ key ] = prop;

}

double NCPA::Atmosphere1D::get( std::string key, units_t quantity_units, double altitude, units_t altitude_units ) {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get( altitude, altitude_units, quantity_units );
}

double NCPA::Atmosphere1D::get_first_derivative( std::string key, units_t quantity_units, double altitude, units_t altitude_units ) {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_first_derivative( altitude, altitude_units, quantity_units );
}


double NCPA::Atmosphere1D::get_second_derivative( std::string key, units_t quantity_units, double altitude, units_t altitude_units ) {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_second_derivative( altitude, altitude_units, quantity_units );
}


void NCPA::Atmosphere1D::calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *t_prop = contents_.at( temperature_key );
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( GAMMA_FOR_C * R_FOR_C * t_prop->get( z_[i], z_units_, NCPA::UNITS_TEMPERATURE_KELVIN ) );
	}

	add_quantity( new_key, nz_, c, NCPA::UNITS_SPEED_METERS_PER_SECOND );
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_pressure_and_density( std::string new_key, 
			std::string pressure_key, std::string density_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *p_prop = contents_.at( pressure_key );
	NCPA::AtmosphericProperty1D *r_prop = contents_.at( density_key );
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( 
			GAMMA_FOR_C * p_prop->get( z_[i], z_units_, NCPA::UNITS_PRESSURE_PASCALS ) 
			/ r_prop->get( z_[i], z_units_, NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER ) 
			);
	}

	add_quantity( new_key, nz_, c, NCPA::UNITS_SPEED_METERS_PER_SECOND );
}





NCPA::AtmosphericProperty1D::AtmosphericProperty1D( size_t n_points, double *altitude_points, units_t altitude_units,
			double *property_values, units_t property_units ) {

	z_ = new double[ n_points ];
	z_units_ = altitude_units;
	std::memcpy( z_, altitude_points, n_points*sizeof(double) );
	values_ = new double[ n_points ];
	units_ = property_units;
	std::memcpy( values_, property_values, n_points*sizeof(double) );
	nz_ = n_points;

	// construct spline
	accel_ = gsl_interp_accel_alloc();
	spline_ = gsl_spline_alloc( gsl_interp_cspline, nz_ );
	gsl_spline_init( spline_, z_, values_, nz_ );

}

NCPA::AtmosphericProperty1D::~AtmosphericProperty1D() {
	delete spline_;
	delete accel_;
	delete [] z_;
	delete [] values_;
}

NCPA::AtmosphericProperty1D::get( double altitude, units_t altitude_units, units_t quantity_units ) {

	// make sure units are consistent
	double z_req = NCPA::Units::convert( altitude, altitude_units, z_units_ );
	if ( z_req < z_[0] || z_req > z[ nz_ - 1 ] ) {
		std::ostringstream oss;
		oss << "Requested altitude " << altitude << " " << NCPA::toStr( altitude_units ) << " outside profile bounds.";
		throw std::range_error( oss.str() );
	}

	return NCPA::Units::convert( gsl_spline_eval( spline_, z_req, accel_ ), units_, quantity_units );
}

NCPA::AtmosphericProperty1D::get_first_derivative( double altitude, units_t altitude_units, units_t quantity_units ) {

	// make sure units are consistent
	double z_req = NCPA::Units::convert( altitude, altitude_units, z_units_ );
	if ( z_req < z_[0] || z_req > z[ nz_ - 1 ] ) {
		std::ostringstream oss;
		oss << "Requested altitude " << altitude << " " << NCPA::toStr( altitude_units ) << " outside profile bounds.";
		throw std::range_error( oss.str() );
	}


	return NCPA::Units::convert_first_derivative( gsl_spline_eval_deriv( spline_, z_req, accel_ ), units_, quantity_units );
}

NCPA::AtmosphericProperty1D::get_second_derivative( double altitude, units_t altitude_units, units_t quantity_units ) {

	// make sure units are consistent
	double z_req = NCPA::Units::convert( altitude, altitude_units, z_units_ );
	if ( z_req < z_[0] || z_req > z[ nz_ - 1 ] ) {
		std::ostringstream oss;
		oss << "Requested altitude " << altitude << " " << NCPA::toStr( altitude_units ) << " outside profile bounds.";
		throw std::range_error( oss.str() );
	}


	return NCPA::Units::convert_second_derivative( gsl_spline_eval_deriv2( spline_, z_req, accel_ ), units_, quantity_units );
}

