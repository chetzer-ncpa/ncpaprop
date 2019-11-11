#include <cmath>
#include "AtmosphericProfile.h"
#include "util.h"
#include "units.h"
#include <stdexcept>

#define GAM 1.4
#ifndef PI
#define PI 3.14159
#endif
#define R 287.0

NCPA::AtmosphericProfile::~AtmosphericProfile() { }

bool NCPA::AtmosphericProfile::good() {
	return good_;
}

bool NCPA::AtmosphericProfile::hasW() {
	return hasW_;
}

bool NCPA::AtmosphericProfile::hasP() {
	return hasP_;
}

bool NCPA::AtmosphericProfile::hasRho() {
	return hasRho_;
}

// Returns sound speed calculated from P and Rho if they are both present,
// from T otherwise.  Altitude z is assumed to be in the specified output units.
double NCPA::AtmosphericProfile::c0( double z ) {
	double zi = NCPA::convert_units( z, z_units_output_, z_units_internal_ );
        if (this->hasRho() && this->hasP()) {
        	return this->calculate_c0_using_p_( zi );
        } else {
        	return this->calculate_c0_using_t_( zi );
        }
}

// Calculate sound speed using sqrt( GAM * R * t ) as a function of z, 
// assuming z has already been converted to the internal altitude units
// Returns c in the specified output speed units
double NCPA::AtmosphericProfile::calculate_c0_using_t_( double z ) {
	double t = NCPA::convert_units( this->t( z ), 
		t_units_internal_, NCPA::UNITS_TEMPERATURE_K );
	return NCPA::convert_units( std::sqrt( GAM * R * t ), 
		NCPA::UNITS_SPEED_MPS, c_units_output_ );	
	//return 1.0e-3 * sqrt( GAM * R * t );
}

// Calculate sound speed using sqrt( GAM * P / Rho ) as a function of z, 
// assuming z has already been converted to the internal altitude units
// Returns c in the specified output speed units
//double NCPA::AtmosphericProfile::calculate_c0_using_p_( double p, double rho ) {
double NCPA::AtmosphericProfile::calculate_c0_using_p_( double z ) {
	double p = NCPA::convert_units( this->p( z ), 
		p_units_internal_, NCPA::UNITS_PRESSURE_PA );
	double rho = NCPA::convert_units( this->rho( z ), 
		rho_units_internal_, NCPA::UNITS_DENSITY_KGPM3 );
	// Add factor of 0.1 to correct for non-SI units
	//return 1.0e-3 * sqrt( 0.1 * GAM * p / rho );    
	return NCPA::convert_units( std::sqrt( GAM * p / rho ),
		NCPA::UNITS_SPEED_MPS, c_units_output_ );
}

// Calculate the effective sound speed approximation for the specified
// propagation direction.  Z units are in the specified output units
double NCPA::AtmosphericProfile::ceff( double z, double phi ) {
	return this->c0( z ) + this->wcomponent( z, phi );
}

// Calculates the magnitude of the wind speed vector.  Z units are in the
// specified output units and are passed on to the u() and v() functions.
// Wind speed is returned in the same units returned by u() and v().
double NCPA::AtmosphericProfile::wspeed( double z ) {
	return std::sqrt( std::pow( this->u( z ), 2.0 ) 
		+ std::pow( this->v( z ), 2.0 ) );
}

// Returns wind vector azimuth.  Z units are in the
// specified output units and are passed on the u() and v() functions.  Angle
// is returned in geophysical units (degrees clockwise from north)
double NCPA::AtmosphericProfile::wdirection( double z ) {
	double angle = PI/2 - std::atan2( this->v(z), this->u(z) );
	while (angle < 0) {
		angle += 2 * PI;
	}
	return NCPA::rad2deg(angle);
}

// Returns the component of the wind speed in the specified direction.  Direction
// is specified in degrees clockwise from north.  Z units are in the
// specified output units
double NCPA::AtmosphericProfile::wcomponent( double z, double phi ) {
	return wspeed( z ) * std::cos( NCPA::deg2rad(phi) - NCPA::deg2rad(this->wdirection( z )) );
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


void NCPA::setUnitsTemperature( NCPA::UNITS_TEMPERATURE u ) {
	t_units_output_ = u;
}
	
void NCPA::setUnitsAltitude( NCPA::UNITS_DISTANCE u ) {
	z_units_output_ = u;
}
	
void NCPA::setUnitsZonalWind( NCPA::UNITS_SPEED u ) {
	u_units_output_ = u;
}
	
void NCPA::setUnitsMeridionalWind( NCPA::UNITS_SPEED u ) {
	v_units_output_ = u;
}
	
void NCPA::setUnitsVerticalWind( NCPA::UNITS_SPEED u ) {
	w_units_output_ = u;
}
	
void NCPA::setUnitsDensity( NCPA::UNITS_DENSITY u ) {
	rho_units_output_ = u;
}
	
void NCPA::setUnitsPressure( NCPA::UNITS_PRESSURE u ) {
	p_units_output_ = u;
}
	
void NCPA::setUnitsSoundSpeed( NCPA::UNITS_SPEED u ) {
	c_units_output_ = u;
}