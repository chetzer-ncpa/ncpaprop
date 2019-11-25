#include <cmath>
#include "AtmosphericProfile.h"
#include "util.h"
#include "units.h"
#include <stdexcept>
#include <iostream>

#define GAM 1.4
#ifndef PI
#define PI 3.14159
#endif
#define R 287.0

NCPA::AtmosphericProfile::~AtmosphericProfile() { }

bool NCPA::AtmosphericProfile::good() {
	return good_;
}

void NCPA::AtmosphericProfile::setOrigin( double lat, double lon ) {
	delete origin_;
	origin_ = new Location( lat, lon );
}

double NCPA::AtmosphericProfile::lat() const {
	return origin_->lat();
}

double NCPA::AtmosphericProfile::lon() const {
	return origin_->lon();
}


void NCPA::AtmosphericProfile::setUnits( NCPA::ATMOSPHERIC_QUANTITY quantity, NCPA::UNITS_TYPE units ) {
	units_map_[ quantity ] = units;
}

NCPA::UNITS_TYPE NCPA::AtmosphericProfile::getUnits( ATMOSPHERIC_QUANTITY quantity ) {
	try {
		return units_map_.at( quantity );
	} catch (const std::out_of_bounds& oob) {
		throw std::out_of_bounds( "No units specified" );
	}
}



/* deprecated */
bool NCPA::AtmosphericProfile::hasW() {
	return this->has( NCPA::QUANTITY_WIND_SPEED_VERTICAL );
	//return hasW_;
}

/* deprecated */
bool NCPA::AtmosphericProfile::hasP() {
	return this->has( NCPA::QUANTITY_AIR_PRESSURE );
	//return hasP_;
}

/* deprecated */
bool NCPA::AtmosphericProfile::hasRho() {
	return this->has( NCPA::QUANTITY_AIR_DENSITY );
	//return hasRho_;
}

double NCPA::AtmosphericProfile::c0( double z ) {
	return this->calculate_c0_( z );
}


// Returns sound speed calculated from P and Rho if they are both present,
// from T otherwise.  Altitude z is assumed to be in the specified output units.
double NCPA::AtmosphericProfile::calculate_c0_( double z ) {
	//double zi = NCPA::convert_units( z, z_units_output_, z_units_internal_ );
        //if (this->hasRho() && this->hasP()) {
	if (this->has( NCPA::QUANTITY_AIR_DENSITY) && this->has( NCPA::QUANTITY_AIR_PRESSURE ) ) {
        	return this->calculate_c0_using_p_( z );
        } else if (this->has( NCPA::QUANTITY_TEMPERATURE ) ) {
        	return this->calculate_c0_using_t_( z );
        } else {
        	throw( std::out_of_bounds( "No atmospheric quantities available to calculate sound speed" ) );
        }
}

// Calculate sound speed using sqrt( GAM * R * t ) as a function of z
// Returns c in the specified output speed units
double NCPA::AtmosphericProfile::calculate_c0_using_t_( double z ) {
	
	// calculation uses t in Kelvin, c in m/s
	double t;
	try {
		t = this->get( NCPA::QUANTITY_TEMPERATURE, z, NCPA::UNITS_TEMPERATURE_KELVIN );
	} catch (const std::out_of_bounds &oob) {
		std::cerr << "Invalid temperature conversion attempted: " << oob.what() << std::endl;
		throw;
	}
		
	try {
		return NCPA::Units::convert( std::sqrt( GAM * R * t ), 
			NCPA::UNITS_SPEED_METERS_PER_SECOND, 
			this->getUnits( NCPA::QUANTITY_STATIC_SOUND_SPEED ) );
	} catch (const std::out_of_bounds &oob) {
		std::cerr << "Invalid sound speed conversion attempted: " << oob.what() << std::endl;
		throw;
	}
}

// Calculate sound speed using sqrt( GAM * P / Rho ) as a function of z, 
// assuming z has already been converted to the internal altitude units
// Returns c in the specified output speed units
double NCPA::AtmosphericProfile::calculate_c0_using_p_( double z ) {
	
	// calculation uses p in Pa, rho in kg/m3
	double p, rho;
	
	// get pressure in right units
	try {
		p = this->get( NCPA::QUANTITY_PRESSURE, z, NCPA::UNITS_PRESSURE_PASCALS );
	} catch (const std::out_of_bounds &oob) {
		std::cerr << "Invalid pressure conversion attempted: " << oob.what() << std::endl;
		throw;
	}
	try {
		rho = this->get( NCPA::QUANTITY_DENSITY, z, NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER );
	} catch (const std::out_of_bounds &oob) {
		std::cerr << "Invalid density conversion attempted: " << oob.what() << std::endl;
		throw;
	}
	
	return NCPA::Units::convert( std::sqrt( GAM * p / rho ),
		NCPA::UNITS_SPEED_METERS_PER_SECOND, 
		this->getUnits( NCPA::QUANTITY_STATIC_SOUND_SPEED ) );
}


double NCPA::AtmosphericProfile::calculate_total_wind_speed( double z ) {
	
	if ( this->has( NCPA::QUANTITY_WIND_SPEED_VERTICAL ) ) {
		return std::sqrt(
			std::pow( this->calculate_horizontal_wind_speed( z ), 2.0 ) + 
			std::pow( this->get( NCPA::QUANTITY_WIND_SPEED_VERTICAL ), 2.0 )
		);
	} else {
		return this->calculate_horizontal_wind_speed( z );
	}
}


double NCPA::AtmosphericProfile::calculate_horizontal_wind_speed( double z ) {
	return std::sqrt(
		std::pow( this->get( NCPA::QUANTITY_WIND_SPEED_WEST_TO_EAST, z ), 2.0 ) +
		std::pow( this->get( NCPA::QUANTITY_WIND_SPEED_SOUTH_TO_NORTH, z ), 2.0 )
	);
}

double NCPA::AtmosphericProfile::calculate_wind_direction_( double z ) {
	
	// calculate in radians
	double angle_math = std::atan2( 
		this->get( NCPA::QUANTITY_WIND_SPEED_SOUTH_TO_NORTH, z ),
		this->get( NCPA::QUANTITY_WIND_SPEED_WEST_TO_EAST, z ) );
	
	// convert to degrees
	NCPA::Units::convert( &angle_math, 1, NCPA::UNITS_ANGLE_RADIANS, NCPA::UNITS_ANGLE_DEGREES, &angle_math );
	
	// convert to internal direction convention
	return NCPA::Units::convert( angle_math, NCPA::UNITS_DIRECTION_DEGREES_COUNTERCLOCKWISE_FROM_EAST,
		this->getUnits( NCPA::QUANTITY_WIND_DIRECTION ) );
}

double NCPA::AtmosphericProfile::calculate_wind_component_( double z, double phi ) {
	
	return this->calculate_horizontal_wind_speed_( z ) * std::cos(
		NCPA::Units::convert( 
			phi - this->calculate_wind_direction_( z ), 
			NCPA::UNITS_ANGLE_DEGREES, NCPA::UNITS_ANGLE_RADIANS
		
	) );
}

double NCPA::AtmosphericProfile::calculate_effective_sound_speed_( double z, double phi ) {
	
	// convert wind speed units to match sound speed units if necessary
	return this->get( NCPA::QUANTITY_STATIC_SOUND_SPEED, z ) + 
		NCPA::Units::convert( this->wcomponent( z, phi ), 
			this->getUnits( QUANTITY_WIND_SPEED_HORIZONTAL ),
			this->getUnits( NCPA::QUANTITY_STATIC_SOUND_SPEED ) );
}









// Calculate the effective sound speed approximation for the specified
// propagation direction.  Z units are in the specified output units
double NCPA::AtmosphericProfile::ceff( double z, double phi ) {
	return this->calculate_effective_sound_speed_( z, phi );
	
}


// Vertical derivative of atmospheric quantities.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::get_dZ( ATMOSPHERIC_QUANTITY quantity, double z ) {
	
	double zminus = z - eps_z;
	double zplus  = z + eps_z;
	double q, qplus, qminus;
	try {
		if (zminus < minimum_altitude()) {
			qplus = this->get( quantity, zplus );
			q = this->get( quantity, z );
			return (qplus - q) / eps_z;
		} else if (zplus > maximum_altitude()) {
			q = this->get( quantity, z );
			qminus = this->get( quantity, zminus );
			return (q - qminus) / eps_z;
		} else {
			qplus = this->get( quantity, zplus );
			qminus = this->get( quantity, zminus );
			return (qplus - qminus) / (2.0 * eps_z);
		}
	} catch (const std::out_of_bounds& oob) {
		throw std::out_of_bounds( "Requested altitude outside of valid bounds" );
	}
}

double NCPA::AtmosphericProfile::get_ddZ( ATMOSPHERIC_QUANTITY quantity, double z ) {
	double zminus = z - eps_z;
	double zplus  = z + eps_z;
	double z2plus = z + 2.0*eps_z;
	double z2minus = z - 2.0*eps_z;
	double q, qplus, qminus, q2plus, q2minus;
	try {
		if (zminus < minimum_altitude()) {
			qplus = this->get( quantity, zplus );
			q2plus = this->get( quantity, z2plus );
			q = this->get( quantity, z );
			return (q2plus - 2.0*qplus + q) / std::pow(eps_z,2.0);
		} else if (zplus > maximum_altitude()) {
			q = this->get( quantity, z );
			qminus = this->get( quantity, zminus );
			return (q - 2.0*qminus + q2minus) / std::pow(eps_z,2.0);
		} else {
			qplus = this->get( quantity, zplus );
			qminus = this->get( quantity, zminus );
			q = this->get( quantity, z );
			return (qplus + qminus - 2.0*q) / std::pow(eps_z,2.0);
		}
	} catch (const std::out_of_bounds& oob) {
		throw std::out_of_bounds( "Requested altitude outside of valid bounds" );
	}
}


//    Unadjusted functions below here    //















// Calculates the magnitude of the wind speed vector.  Z units are in the
// specified output units and are passed on to the u() and v() functions.
// Wind speed is returned in the same units returned by u() and v().
double NCPA::AtmosphericProfile::wspeed( double z ) {
	return this->calculate_horizontal_wind_speed_( z );
	//return std::sqrt( std::pow( this->u( z ), 2.0 ) 
	//	+ std::pow( this->v( z ), 2.0 ) );
}

// Returns wind vector azimuth.  Z units are in the
// specified output units and are passed on the u() and v() functions.  Angle
// is returned in geophysical units (degrees clockwise from north)
double NCPA::AtmosphericProfile::wdirection( double z ) {
	return this->calculate_wind_direction_( z );
	//double angle = PI/2 - std::atan2( this->v(z), this->u(z) );
	//while (angle < 0) {
	//	angle += 2 * PI;
	//}
	//return NCPA::rad2deg(angle);
}

// Returns the component of the wind speed in the specified direction.  Direction
// is specified in degrees clockwise from north.  Z units are in the
// specified output units
double NCPA::AtmosphericProfile::wcomponent( double z, double phi ) {
	return this->calculate_wind_component_( z, phi );
	//return wspeed( z ) * std::cos( NCPA::deg2rad(phi) - NCPA::deg2rad(this->wdirection( z )) );
}

// Spatial derivatives of c0.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::dc0dz( double z ) {
	if (z-eps_z < z0()) {
		return (this->c0(z+eps_z) - this->c0(z)) / eps_z;
	} else {
		return (this->c0(z+eps_z) - this->c0(z-eps_z)) / (2.0*eps_z);
	}
}

// Spatial derivatives of effective sound speed.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::dceffdz( double z, double phi ) {
	if (z-eps_z < z0()) {
		return (this->ceff(z+eps_z,phi) - this->ceff(z,phi)) / eps_z;
	} else {
		return (this->ceff(z+eps_z,phi) - this->ceff(z-eps_z,phi)) / (2.0*eps_z);
	}
}

// Spatial second derivatives of sound speed.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::ddc0dzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->c0(z+2*eps_z) - 2*this->c0(z+eps_z) + this->c0(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->c0(z + eps_z) + this->c0(z - eps_z)- 2.0*this->c0(z))/(std::pow(eps_z,2.0));
	}
}


// Spatial second derivatives of effective sound speed.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::ddceffdzdz(double z, double phi) {
	if (z-eps_z < z0()) {
		return (this->ceff(z+2*eps_z,phi) - 2*this->ceff(z+eps_z,phi) + this->ceff(z,phi)) / (std::pow(eps_z,2.0));
	} else {
		return (this->ceff(z + eps_z, phi) + this->ceff(z - eps_z, phi)- 2.0*this->ceff(z,phi))/(std::pow(eps_z,2.0));
	}
}


// Spatial derivatives of temperature.  Z units are in the
// specified output units.  Returns the derivative in specified 
// temperature/distance units
double NCPA::AtmosphericProfile::dtdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->t(z+eps_z) - this->t(z)) / eps_z;
	} else {
		return (this->t(z+eps_z) - this->t(z-eps_z)) / (2.0*eps_z);
	}
}


// Spatial second derivatives of temperature.  Z units are in the
// specified output units.  Returns the derivative in specified 
// temperature/distance units
double NCPA::AtmosphericProfile::ddtdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->t(z+2*eps_z) - 2*this->t(z+eps_z) + this->t(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->t(z + eps_z) + this->t(z - eps_z)- 2.0*this->t(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of zonal winds.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::dudz( double z ) {
	if (z-eps_z < z0()) {
		return (this->u(z+eps_z) - this->u(z)) / eps_z;
	} else {
		return (this->u(z+eps_z) - this->u(z-eps_z)) / (2.0*eps_z);
	}
}

// Spatial second derivatives of zonal winds.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::ddudzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->u(z+2*eps_z) - 2*this->u(z+eps_z) + this->u(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->u(z + eps_z) + this->u(z - eps_z)- 2.0*this->u(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of meridional winds.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::dvdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->v(z+eps_z) - this->v(z)) / eps_z;
	} else {
		return (this->v(z+eps_z) - this->v(z-eps_z)) / (2.0*eps_z);
	}
}

// Spatial second derivatives of meridional winds.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::ddvdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->v(z+2*eps_z) - 2*this->v(z+eps_z) + this->v(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->v(z + eps_z) + this->v(z - eps_z)- 2.0*this->v(z))/(std::pow(eps_z,2.0));
	}
}


// Spatial derivatives of vertical winds.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::dwdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->w(z+eps_z) - this->w(z)) / eps_z;
	} else {
		return (this->w(z+eps_z) - this->w(z-eps_z)) / (2.0*eps_z);
	}
}

// Spatial second derivatives of vertical winds.  Z units are in the
// specified output units.  Returns the derivative in specified 
// speed/distance units
double NCPA::AtmosphericProfile::ddwdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->w(z+2*eps_z) - 2*this->w(z+eps_z) + this->w(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->w(z + eps_z) + this->w(z - eps_z) - 2.0*this->w(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of pressure.  Z units are in the
// specified output units.  Returns the derivative in specified 
// pressure/distance units
double NCPA::AtmosphericProfile::dpdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->p(z+eps_z) - this->p(z)) / eps_z;
	} else {
		return (this->p(z+eps_z) - this->p(z-eps_z)) / (2.0*eps_z);
	}
}


// Spatial second derivatives of pressure.  Z units are in the
// specified output units.  Returns the derivative in specified 
// pressure/distance units
double NCPA::AtmosphericProfile::ddpdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->p(z+2*eps_z) - 2*this->p(z+eps_z) + this->p(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->p(z + eps_z) + this->p(z - eps_z)- 2.0*this->p(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of density.  Z units are in the
// specified output units.  Returns the derivative in specified 
// density/distance units
double NCPA::AtmosphericProfile::drhodz( double z ) {
	if (z-eps_z < z0()) {
		return (this->rho(z+eps_z) - this->rho(z)) / eps_z;
	} else {
		return (this->rho(z+eps_z) - this->rho(z-eps_z)) / (2.0*eps_z);
	}
}


// Spatial second derivatives of density.  Z units are in the
// specified output units.  Returns the derivative in specified 
// density/distance units
double NCPA::AtmosphericProfile::ddrhodzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->rho(z+2*eps_z) - 2*this->rho(z+eps_z) + this->rho(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->rho(z + eps_z) + this->rho(z - eps_z)- 2.0*this->rho(z))/(std::pow(eps_z,2.0));
	}
}

double NCPA::AtmosphericProfile::z0() { return this->minimumAltitude(); }
double NCPA::AtmosphericProfile::t( double z ) { return this->get( NCPA::QUANTITY_TEMPERATURE, z ); }
double NCPA::AtmosphericProfile::u( double z ) { return this->get( NCPA::QUANTITY_WIND_SPEED_WEST_TO_EAST, z ); }
double NCPA::AtmosphericProfile::v( double z ) { return this->get( NCPA::QUANTITY_WIND_SPEED_SOUTH_TO_NORTH, z ); }
double NCPA::AtmosphericProfile::w( double z ) { return this->get( NCPA::QUANTITY_WIND_SPEED_VERTICAL, z ); }
double NCPA::AtmosphericProfile::p( double z ) { return this->get( NCPA::QUANTITY_AIR_PRESSURE, z ); }
double NCPA::AtmosphericProfile::rho( double z ) { return this->get( NCPA::QUANTITY_AIR_DENSITY, z ); }





/*
OLD, OBSOLETE CODE HERE, KEPT FOR COMPARISON AND TESTING PURPOSES
DELETE AFTER TESTING OF NEW CODE
*/