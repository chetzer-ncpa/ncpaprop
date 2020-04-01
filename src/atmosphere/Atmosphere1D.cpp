#include "Atmosphere.h"
#include "Atmosphere1D.h"
#include "AtmosphericProperty1D.h"
#include <map>
#include <iostream>
#include <cstring>

#define GAMMA_FOR_C 1.4
#define R_FOR_C 287.0



NCPA::Atmosphere1D::Atmosphere1D( size_t n_altitude_points, double *altitude_points, units_t altitude_units ) {
	contents_.clear();
	nz_ = n_altitude_points;
	z_ = new double[ nz_ ];
	std::memcpy( z_, altitude_points, nz_ * sizeof(double) );
	z_units_.push( altitude_units );
}

NCPA::Atmosphere1D::~Atmosphere1D() {
	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it=contents_.begin(); it != contents_.end(); ++it) {
		delete it->second;
	}
	contents_.clear();
	while (! z_units_.empty()) {
		z_units_.pop();
	}
}

size_t NCPA::Atmosphere1D::get_basis_length() const {
	return nz_;
}

void NCPA::Atmosphere1D::get_altitude_vector( double *buffer, units_t *buffer_units ) const {
	*buffer_units = z_units_.top();
	std::memcpy( buffer, z_, nz_* sizeof( double ) );
}

void NCPA::Atmosphere1D::get_property_vector( std::string key, double *buffer, units_t *buffer_units ) const {
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );
	prop->get_vector( buffer, buffer_units );
}


double NCPA::Atmosphere1D::get_minimum_altitude() const {
	return z_[0];
	//return NCPA::Units::convert( z_[0], z_units_, altitude_units );
}

double NCPA::Atmosphere1D::get_maximum_altitude() const {
	return z_[ nz_ - 1 ];
	//return NCPA::Units::convert( z_[ nz_-1 ], z_units_, altitude_units );
}

void NCPA::Atmosphere1D::add_property( std::string key, size_t n_points, double *quantity_points, units_t quantity_units ) {

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

	prop = new AtmosphericProperty1D( n_points, z_, z_units_.top(), quantity_points, quantity_units );
	contents_[ key ] = prop;
}

void NCPA::Atmosphere1D::add_property( std::string key, double value, NCPA::units_t units ) {
	NCPA::ScalarWithUnits *scalar;
	try {
		scalar = scalar_contents_.at( key );
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	scalar = new ScalarWithUnits( value, units );
	scalar_contents_[ key ] = scalar;
}

double NCPA::Atmosphere1D::get( std::string key, double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop;
	try {
		prop = contents_.at( key );
	} catch (std::out_of_range& oor) {
		throw std::out_of_range( "No vector quantity \"" + key + "\" found" );
	}
	return prop->get( altitude );
}

double NCPA::Atmosphere1D::get( std::string key ) const {
	NCPA::ScalarWithUnits *prop;
	try {
		prop = scalar_contents_.at( key );
	} catch (std::out_of_range& oor) {
		throw std::out_of_range( "No scalar quantity \"" + key + "\" found" );
	}
	return prop->get();
}

bool NCPA::Atmosphere1D::contains_vector( std::string key ) const {
	NCPA::AtmosphericProperty1D *prop;
	try {
		prop = contents_.at( key );
	} catch (std::out_of_range& oor) {
		return false;
	}
	return true;
}

bool NCPA::Atmosphere1D::contains_scalar( std::string key ) const {
	NCPA::ScalarWithUnits *prop;
	try {
		prop = scalar_contents_.at( key );
	} catch (std::out_of_range& oor) {
		return false;
	}
	return true;
}

double NCPA::Atmosphere1D::get_first_derivative( std::string key, double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_first_derivative( altitude );
}


double NCPA::Atmosphere1D::get_second_derivative( std::string key, double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_second_derivative( altitude );
}


void NCPA::Atmosphere1D::calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *t_prop = contents_.at( temperature_key );
	t_prop->convert_units( NCPA::UNITS_TEMPERATURE_KELVIN );
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( GAMMA_FOR_C * R_FOR_C * t_prop->get( z_[i] ) );
	}
	t_prop->revert_units();

	add_property( new_key, nz_, c, NCPA::UNITS_SPEED_METERS_PER_SECOND );
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_pressure_and_density( std::string new_key, 
			std::string pressure_key, std::string density_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *p_prop = contents_.at( pressure_key );
	p_prop->convert_units( NCPA::UNITS_PRESSURE_PASCALS );
	NCPA::AtmosphericProperty1D *r_prop = contents_.at( density_key );
	r_prop->convert_units( NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER );
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( GAMMA_FOR_C * p_prop->get( z_[i] ) / r_prop->get( z_[i] ) );
	}

	p_prop->revert_units();
	r_prop->revert_units();
	add_property( new_key, nz_, c, NCPA::UNITS_SPEED_METERS_PER_SECOND );
}

void NCPA::Atmosphere1D::convert_altitude_units( units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	if (new_units != z_units_.top()) {
		do_units_conversion_( nz_, z_, z_units_.top(), new_units );
	}
	z_units_.push( new_units );

	// update all contents
	for ( std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it = contents_.begin();
						it != contents_.end(); ++it ) {
		(*it).second->convert_altitude_units( new_units );
	}
}

void NCPA::Atmosphere1D::revert_altitude_units() {
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

void NCPA::Atmosphere1D::convert_property_units( std::string key, units_t new_units ) {
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );
	prop->convert_units( new_units );
}

void NCPA::Atmosphere1D::revert_property_units( std::string key ) {
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );
	prop->revert_units();
}



void NCPA::Atmosphere1D::do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double *units_buffer = new double[ n_points ];
	std::memset( units_buffer, 0, n_points * sizeof( double ) );
	
	// throws out_of_range if conversion is undefined
	NCPA::Units::convert( inplace, n_points, fromUnits, toUnits, units_buffer );

	// successful, so record the units change
	std::memcpy( inplace, units_buffer, n_points * sizeof( double ) );
	delete [] units_buffer;
}

std::vector< std::string > NCPA::Atmosphere1D::get_keys() const {
	std::vector< std::string > keys;
	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::const_iterator it = contents_.begin();
			it != contents_.end(); ++it) {
		keys.push_back( (*it).first );
	}
	return keys;
}