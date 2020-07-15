#include "Atmosphere2D.h"
#include "Atmosphere1D.h"
#include <vector>
#include <climits>
#include <cfloat>
#include <stdexcept>
#include <algorithm>


NCPA::Atmosphere2D::Atmosphere2D() {
	profiles_.clear();
	midpoints_.clear();
	clear_last_index_();
	range_units_ = NCPA::UNITS_NONE;
}

NCPA::Atmosphere2D::~Atmosphere2D() {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator i = profiles_.begin();
		i != profiles_.end(); ++i ) {
		delete (*i);
	}
	profiles_.clear();
	midpoints_.clear();
}

bool NCPA::sort_profiles_by_range_( Atmosphere1D *p1, Atmosphere1D *p2 ) {
	// assumes consistent units
	return p1->get( "_RANGE_" ) < p2->get( "_RANGE_" );
}

void NCPA::Atmosphere2D::insert_profile( const NCPA::Atmosphere1D *profile, double range ) {

	// throw an exception if units haven't been specified yet
	if ( range_units_ == NCPA::UNITS_NONE ) {
		throw std::runtime_error( "Range units have not been specified!" );
	}
	NCPA::Atmosphere1D *newProfile = new NCPA::Atmosphere1D( *profile );
	newProfile->add_property( "_RANGE_", range, range_units_ );
	profiles_.push_back( newProfile );
	sorted_ = false;
	clear_last_index_();
}

void NCPA::Atmosphere2D::set_insert_range_units( units_t u ) {
	range_units_ = u;
}

void NCPA::Atmosphere2D::clear_last_index_() {
	last_index_ = INT_MAX;
	last_index_min_range_ = 0.0;
	last_index_max_range_ = 0.0;
}

void NCPA::Atmosphere2D::set_last_index_( size_t ind ) {
	if (profiles_.size() == 0) {
		// should never happen, but you never know
		throw std::runtime_error( "No profiles have been added to 2-D atmosphere!" );
	} else if (profiles_.size() == 1) {
		last_index_ = 0;
		last_index_min_range_ = 0.0;
		last_index_max_range_ = DBL_MAX;
	} else if (ind == 0) {
		last_index_ = 0;
		last_index_min_range_ = 0.0;
		last_index_max_range_ = midpoints_[ 0 ];
	} else if (ind == (profiles_.size() - 1) ) {
		last_index_ = ind;
		last_index_min_range_ = midpoints_[ ind-1 ];
		last_index_max_range_ = DBL_MAX;
	} else {
		last_index_ = ind;
		last_index_min_range_ = midpoints_[ ind-1 ];
		last_index_max_range_ = midpoints_[ ind ];
	}
}

void NCPA::Atmosphere2D::sort_profiles() {
	std::sort( profiles_.begin(), profiles_.end(), NCPA::sort_profiles_by_range_ );
	midpoints_.clear();
	clear_last_index_();
	sorted_ = true;
}

void NCPA::Atmosphere2D::calculate_midpoints_() {
	midpoints_.clear();
	if (profiles_.size() == 1) {
		return;
	}

	if (!sorted_) {
		sort_profiles();
	}

	for (size_t i = 1; i < profiles_.size(); i++) {
		double r0 = profiles_.at( i-1 )->get( "_RANGE_" );
		double r1 = profiles_.at( i )->get( "_RANGE_" );
		double mid = r0 + ( r1 - r0 ) * 0.5;
		midpoints_.push_back( mid );
	}
}

size_t NCPA::Atmosphere2D::get_profile_index_( double range ) {
	if (profiles_.size() == 0) {
		throw std::runtime_error( "No profiles have been added to 2-D atmosphere!" );
	}

	if (profiles_.size() == 1) {
		return 0;
	}

	if (midpoints_.size() == 0) {
		calculate_midpoints_();
		clear_last_index_();
	}

	// first check to see if we have a valid last range, to avoid going through the search again
	if (last_index_ < INT_MAX) {
		if (range >= last_index_min_range_ && range < last_index_max_range_) {
			return last_index_;
		}
	}

	size_t ind = 0;
	while ( ind < midpoints_.size() && midpoints_[ ind ] <= range ) {
		ind++;
	}
	set_last_index_( ind );
	return ind;
}

double NCPA::Atmosphere2D::get( std::string key, double range ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get( key );
}

double NCPA::Atmosphere2D::get( std::string key, double range, double altitude ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get( key, altitude );
}

double NCPA::Atmosphere2D::get_first_derivative( std::string key, double range, double altitude ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get_first_derivative( key, altitude );
}

double NCPA::Atmosphere2D::get_second_derivative( std::string key, double range, double altitude ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get_second_derivative( key, altitude );
}

size_t NCPA::Atmosphere2D::nz( double range ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->nz();
}

void NCPA::Atmosphere2D::get_altitude_vector( double range, double *buffer, units_t *buffer_units ) {
	size_t ind = get_profile_index_( range );
	profiles_.at( ind )->get_altitude_vector( buffer, buffer_units );
}

void NCPA::Atmosphere2D::get_property_vector( double range, std::string key, double *buffer, units_t *buffer_units ) {
	size_t ind = get_profile_index_( range );
	profiles_.at( ind )->get_property_vector( key, buffer, buffer_units );
}

void NCPA::Atmosphere2D::get_altitude_vector( double range, double *buffer ) {
	size_t ind = get_profile_index_( range );
	profiles_.at( ind )->get_altitude_vector( buffer );
}

void NCPA::Atmosphere2D::get_property_vector( double range, std::string key, double *buffer ) {
	size_t ind = get_profile_index_( range );
	profiles_.at( ind )->get_property_vector( key, buffer );
}

NCPA::units_t NCPA::Atmosphere2D::get_altitude_units( double range ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get_altitude_units();
}

NCPA::units_t NCPA::Atmosphere2D::get_property_units( double range, std::string key ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get_property_units( key );
}

double NCPA::Atmosphere2D::get_minimum_altitude( double range ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get_minimum_altitude();
}

double NCPA::Atmosphere2D::get_maximum_altitude( double range ) {
	size_t ind = get_profile_index_( range );
	return profiles_.at( ind )->get_maximum_altitude();
}

void NCPA::Atmosphere2D::calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key, 
		units_t wind_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_sound_speed_from_temperature( new_key, temperature_key, wind_units );
	}
}

void NCPA::Atmosphere2D::calculate_sound_speed_from_pressure_and_density( std::string new_key, std::string pressure_key, 
		std::string density_key, NCPA::units_t wind_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_sound_speed_from_pressure_and_density( new_key, pressure_key, density_key, wind_units );
	}
}

void NCPA::Atmosphere2D::calculate_wind_speed( std::string new_key, std::string we_wind_speed_key, std::string sn_wind_speed_key ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_wind_speed( new_key, we_wind_speed_key, sn_wind_speed_key );
	}
}

void NCPA::Atmosphere2D::calculate_wind_direction( std::string new_key, std::string we_wind_speed_key, std::string sn_wind_speed_key, 
		NCPA::units_t direction_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_wind_direction( new_key, we_wind_speed_key, sn_wind_speed_key, direction_units );
	}
}

void NCPA::Atmosphere2D::calculate_attenuation( std::string new_key, std::string temperature_key, std::string pressure_key, 
		std::string density_key, double freq, double tweak_factor ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_attenuation( new_key, temperature_key, pressure_key, density_key, freq, tweak_factor );
	}
}

void NCPA::Atmosphere2D::calculate_wind_component( std::string new_key, std::string wind_speed_key, std::string wind_direction_key, 
		double azimuth ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_wind_component( new_key, wind_speed_key, wind_direction_key, azimuth );
	}
}

void NCPA::Atmosphere2D::calculate_effective_sound_speed( std::string new_key, std::string sound_speed_key, 
	std::string wind_component_key ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->calculate_effective_sound_speed( new_key, sound_speed_key, wind_component_key );
	}
}

void NCPA::Atmosphere2D::convert_altitude_units( NCPA::units_t new_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->convert_altitude_units( new_units );
	}
}

void NCPA::Atmosphere2D::convert_property_units( std::string key, NCPA::units_t new_units ) {
	for ( std::vector< NCPA::Atmosphere1D * >::iterator it = profiles_.begin();
		  it != profiles_.end(); ++it ) {
		(*it)->convert_property_units( key, new_units );
	}
}