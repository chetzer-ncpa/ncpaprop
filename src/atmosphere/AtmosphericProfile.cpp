#include <cmath>
#include "AtmosphericProfile.h"
#include "util.h"
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

double NCPA::AtmosphericProfile::c0( double z ) {
        if (this->hasRho() && this->hasP()) {
        	return this->calculate_c0_using_p_( this->p(z), this->rho(z) );
        } else {
        	return this->calculate_c0_using_t_( this->t(z) );
        }
}

double NCPA::AtmosphericProfile::calculate_c0_using_t_( double t ) {
	return 1.0e-3 * sqrt( GAM * R * t );
}

double NCPA::AtmosphericProfile::calculate_c0_using_p_( double p, double rho ) {
	// Add factor of 0.1 to correct for non-SI units
	return 1.0e-3 * sqrt( 0.1 * GAM * p / rho );    
}

double NCPA::AtmosphericProfile::ceff( double z, double phi ) {
	return this->c0(z) + this->wcomponent(z,phi);
}


double NCPA::AtmosphericProfile::wspeed( double z ) {
	return std::sqrt( std::pow(this->u(z),2.0) 
	+ std::pow(this->v(z),2.0) );
}

// Returns azimuth in degrees clockwise from north
double NCPA::AtmosphericProfile::wdirection( double z ) {
	double angle = PI/2 - std::atan2( this->v(z), this->u(z) );
	while (angle < 0) {
		angle += 2 * PI;
	}
	return NCPA::rad2deg(angle);
}

double NCPA::AtmosphericProfile::wcomponent( double z, double phi ) {
	return wspeed( z ) * std::cos( NCPA::deg2rad(phi) - NCPA::deg2rad(this->wdirection( z )) );
}

// Spatial derivatives of c0
double NCPA::AtmosphericProfile::dc0dz( double z ) {
	if (z-eps_z < z0()) {
		return (this->c0(z+eps_z) - this->c0(z)) / eps_z;
	} else {
		return (this->c0(z+eps_z) - this->c0(z-eps_z)) / (2.0*eps_z);
	}
}

// Spatial derivatives of effective sound speed
double NCPA::AtmosphericProfile::dceffdz( double z, double phi ) {
	if (z-eps_z < z0()) {
		return (this->ceff(z+eps_z,phi) - this->ceff(z,phi)) / eps_z;
	} else {
		return (this->ceff(z+eps_z,phi) - this->ceff(z-eps_z,phi)) / (2.0*eps_z);
	}
}

// Spatial second derivatives of sound speed
double NCPA::AtmosphericProfile::ddc0dzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->c0(z+2*eps_z) - 2*this->c0(z+eps_z) + this->c0(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->c0(z + eps_z) + this->c0(z - eps_z)- 2.0*this->c0(z))/(std::pow(eps_z,2.0));
	}
}


// Spatial second derivatives of effective sound speed
double NCPA::AtmosphericProfile::ddceffdzdz(double z, double phi) {
	if (z-eps_z < z0()) {
		return (this->ceff(z+2*eps_z,phi) - 2*this->ceff(z+eps_z,phi) + this->ceff(z,phi)) / (std::pow(eps_z,2.0));
	} else {
		return (this->ceff(z + eps_z, phi) + this->ceff(z - eps_z, phi)- 2.0*this->ceff(z,phi))/(std::pow(eps_z,2.0));
	}
}


// Spatial derivatives of temperature
double NCPA::AtmosphericProfile::dtdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->t(z+eps_z) - this->t(z)) / eps_z;
	} else {
		return (this->t(z+eps_z) - this->t(z-eps_z)) / (2.0*eps_z);
	}
}


// Spatial second derivatives of temperature
double NCPA::AtmosphericProfile::ddtdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->t(z+2*eps_z) - 2*this->t(z+eps_z) + this->t(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->t(z + eps_z) + this->t(z - eps_z)- 2.0*this->t(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of zonal winds
double NCPA::AtmosphericProfile::dudz( double z ) {
	if (z-eps_z < z0()) {
		return (this->u(z+eps_z) - this->u(z)) / eps_z;
	} else {
		return (this->u(z+eps_z) - this->u(z-eps_z)) / (2.0*eps_z);
	}
}

double NCPA::AtmosphericProfile::ddudzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->u(z+2*eps_z) - 2*this->u(z+eps_z) + this->u(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->u(z + eps_z) + this->u(z - eps_z)- 2.0*this->u(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of meridional winds
double NCPA::AtmosphericProfile::dvdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->v(z+eps_z) - this->v(z)) / eps_z;
	} else {
		return (this->v(z+eps_z) - this->v(z-eps_z)) / (2.0*eps_z);
	}
}

double NCPA::AtmosphericProfile::ddvdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->v(z+2*eps_z) - 2*this->v(z+eps_z) + this->v(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->v(z + eps_z) + this->v(z - eps_z)- 2.0*this->v(z))/(std::pow(eps_z,2.0));
	}
}


// Spatial derivatives of vertical winds
double NCPA::AtmosphericProfile::dwdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->w(z+eps_z) - this->w(z)) / eps_z;
	} else {
		return (this->w(z+eps_z) - this->w(z-eps_z)) / (2.0*eps_z);
	}
}

double NCPA::AtmosphericProfile::ddwdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->w(z+2*eps_z) - 2*this->w(z+eps_z) + this->w(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->w(z + eps_z) + this->w(z - eps_z) - 2.0*this->w(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of pressure
double NCPA::AtmosphericProfile::dpdz( double z ) {
	if (z-eps_z < z0()) {
		return (this->p(z+eps_z) - this->p(z)) / eps_z;
	} else {
		return (this->p(z+eps_z) - this->p(z-eps_z)) / (2.0*eps_z);
	}
}


// Spatial second derivatives of pressure
double NCPA::AtmosphericProfile::ddpdzdz(double z) {
	if (z-eps_z < z0()) {
		return (this->p(z+2*eps_z) - 2*this->p(z+eps_z) + this->p(z)) / (std::pow(eps_z,2.0));
	} else {
		return (this->p(z + eps_z) + this->p(z - eps_z)- 2.0*this->p(z))/(std::pow(eps_z,2.0));
	}
}

// Spatial derivatives of density
double NCPA::AtmosphericProfile::drhodz( double z ) {
	if (z-eps_z < z0()) {
		return (this->rho(z+eps_z) - this->rho(z)) / eps_z;
	} else {
		return (this->rho(z+eps_z) - this->rho(z-eps_z)) / (2.0*eps_z);
	}
}


// Spatial second derivatives of density
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
