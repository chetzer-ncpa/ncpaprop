#include <cmath>
#include "AtmosphericSpecification.h"
#include "util.h"
//#include "ExceptionWithStack.h"
#include <stdexcept>


#define GAM 1.4
#define PI 3.14159
#define R 287.0


NCPA::AtmosphericSpecification::~AtmosphericSpecification() { 
	good_ = false;
	eps_x = 0.01;
	eps_z = 0.01;
}

bool NCPA::AtmosphericSpecification::good() {
	return good_;
}

bool NCPA::AtmosphericSpecification::hasW() {
	return hasW_;
}

bool NCPA::AtmosphericSpecification::hasP() {
	return hasP_;
}

bool NCPA::AtmosphericSpecification::hasRho() {
	return hasRho_;
}

double NCPA::AtmosphericSpecification::c0( double x, double y, double z ) {
    if (this->hasRho() && this->hasP()) {
    	return 1.0e-3 * sqrt( GAM * this->p(x,y,z) / this->rho(x,y,z) );
    } else {
    	return 1.0e-3 * sqrt(GAM*R*this->t(x,y,z));
    }
}

double NCPA::AtmosphericSpecification::ceff( double x, double y, double z, double phi ) {
    return this->c0(x,y,z) + this->wcomponent(x,y,z,phi);
}


double NCPA::AtmosphericSpecification::wspeed( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return std::sqrt( std::pow(this->u(x,y,z),2.0) 
                        + std::pow(this->v(x,y,z),2.0)
                        + std::pow(this->w(x,y,z),2.0) );
}

// Returns azimuth in degrees clockwise from north
double NCPA::AtmosphericSpecification::wdirection( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	double angle = PI/2 - std::atan2( this->v(x,y,z), this->u(x,y,z) );
	while (angle < 0) {
		angle += 2 * PI;
	}
	return NCPA::rad2deg(angle);
}

double NCPA::AtmosphericSpecification::wcomponent( double x, double y, double z, double phi ) {
//      if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return wspeed( x, y, z ) * std::cos( NCPA::deg2rad(phi) - NCPA::deg2rad(this->wdirection( x, y, z )) );
}

// Spatial derivatives of c0
double NCPA::AtmosphericSpecification::dc0dx( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->c0(x+eps_x,y,z) - this->c0(x-eps_x,y,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dc0dy( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->c0(x,y+eps_x,z) - this->c0(x,y-eps_x,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dc0dz( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->c0(x,y,z+eps_z) - this->c0(x,y,z)) / eps_z;
	} else {
		return (this->c0(x,y,z+eps_z) - this->c0(x,y,z-eps_z)) / (2.0*eps_z);
	}
}

// Spatial derivatives of effective sound speed
double NCPA::AtmosphericSpecification::dceffdx( double x, double y, double z, double phi ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->ceff(x+eps_x,y,z,phi) - this->ceff(x-eps_x,y,z,phi)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dceffdy( double x, double y, double z, double phi ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->ceff(x,y+eps_x,z,phi) - this->ceff(x,y-eps_x,z,phi)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dceffdz( double x, double y, double z, double phi ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->ceff(x,y,z+eps_z,phi) - this->ceff(x,y,z,phi)) / eps_z;
	} else {
		return (this->ceff(x,y,z+eps_z,phi) - this->ceff(x,y,z-eps_z,phi)) / (2.0*eps_z);
	}
}

// Spatial second derivatives of sound speed
double NCPA::AtmosphericSpecification::ddc0dxdx(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->c0(x + eps_x, y, z) + this->c0(x - eps_x, y, z)- 2.0*this->c0(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddc0dydy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->c0(x, y + eps_x, z) + this->c0(x, y - eps_x, z)- 2.0*this->c0(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddc0dzdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->c0(x,y,z+2*eps_z) - 2*this->c0(x,y,z+eps_z) + this->c0(x,y,z)) / (pow(eps_z,2.0));
	} else {
		return (this->c0(x,y,z + eps_z) + this->c0(x,y,z - eps_z)- 2.0*this->c0(x,y,z))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddc0dxdy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->c0(x + eps_x, y + eps_x, z) - this->c0(x + eps_x, y - eps_x, z) - this->c0(x - eps_x, y + eps_x, z) + this->c0(x - eps_x, y - eps_x, z))/(4.0*eps_x*eps_x);
}
double NCPA::AtmosphericSpecification::ddc0dxdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->c0(x + eps_x,y,z + eps_z) - this->c0(x - eps_x, y, z + eps_z) - this->c0(x + eps_x, y, z) + this->c0(x - eps_x, y, z) ) / (2.0 * eps_x * eps_z);
	} else {
	        return (this->c0(x + eps_x, y, z + eps_z) - this->c0(x + eps_x, y, z - eps_z) - this->c0(x - eps_x, y , z + eps_z) + this->c0(x - eps_x, y, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddc0dydz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->c0(x,y + eps_x,z + eps_z) - this->c0(x, y - eps_x, z + eps_z) - this->c0(x, y + eps_x, z) + this->c0(x, y - eps_x, z)) / (2.0 * eps_x * eps_z);
	} else {
	        return (this->c0(x, y + eps_x, z + eps_z) - this->c0(x, y - eps_x, z + eps_z) - this->c0(x, y + eps_x, z - eps_z) + this->c0(x, y - eps_x, z - eps_z))/(4.0*eps_x*eps_z);
	}
}


// Spatial second derivatives of effective sound speed
double NCPA::AtmosphericSpecification::ddceffdxdx(double x, double y, double z, double phi) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->ceff(x + eps_x, y, z, phi) + this->ceff(x - eps_x, y, z, phi)- 2.0*this->ceff(x, y, z, phi))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddceffdydy(double x, double y, double z, double phi) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->ceff(x, y + eps_x, z, phi) + this->ceff(x, y - eps_x, z, phi)- 2.0*this->ceff(x, y, z, phi))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddceffdzdz(double x, double y, double z, double phi) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->ceff(x,y,z+2*eps_z,phi) - 2*this->ceff(x,y,z+eps_z,phi) + this->ceff(x,y,z,phi)) / (pow(eps_z,2.0));
        } else {
	        return (this->ceff(x,y,z + eps_z, phi) + this->ceff(x,y,z - eps_z, phi)- 2.0*this->ceff(x,y,z,phi))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddceffdxdy(double x, double y, double z, double phi) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->ceff(x + eps_x, y + eps_x, z, phi) - this->ceff(x + eps_x, y - eps_x, z, phi) - this->ceff(x - eps_x, y + eps_x, z, phi) + this->ceff(x - eps_x, y - eps_x, z, phi))/(4.0*eps_x*eps_x);
}
double NCPA::AtmosphericSpecification::ddceffdxdz(double x, double y, double z, double phi) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->ceff(x + eps_x,y,z + eps_z,phi) - this->ceff(x - eps_x, y, z + eps_z,phi) - this->ceff(x + eps_x, y, z,phi) + this->ceff(x - eps_x, y, z,phi) ) / (2.0 * eps_x * eps_z);
        } else {
	        return (this->ceff(x + eps_x, y, z + eps_z,phi) - this->ceff(x + eps_x, y, z - eps_z,phi) - this->ceff(x - eps_x, y , z + eps_z,phi) + this->ceff(x - eps_x, y, z - eps_z,phi))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddceffdydz(double x, double y, double z, double phi) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->ceff(x,y + eps_x,z + eps_z,phi) - this->ceff(x, y - eps_x, z + eps_z,phi) - this->ceff(x, y + eps_x, z,phi) + this->ceff(x, y - eps_x, z,phi)) / (2.0 * eps_x * eps_z);
        } else {
	        return (this->ceff(x, y + eps_x, z + eps_z,phi) - this->ceff(x, y - eps_x, z + eps_z,phi) - this->ceff(x, y + eps_x, z - eps_z,phi) + this->ceff(x, y - eps_x, z - eps_z,phi))/(4.0*eps_x*eps_z);
	}
}


// Spatial derivatives of temperature
double NCPA::AtmosphericSpecification::dtdx( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->t(x+eps_x,y,z) - this->t(x-eps_x,y,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dtdy( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->t(x,y+eps_x,z) - this->t(x,y-eps_x,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dtdz( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->t(x,y,z+eps_z) - this->t(x,y,z)) / eps_z;
        } else {
		return (this->t(x,y,z+eps_z) - this->t(x,y,z-eps_z)) / (2.0*eps_z);
	}
}


// Spatial second derivatives of temperature
double NCPA::AtmosphericSpecification::ddtdxdx(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->t(x + eps_x, y, z) + this->t(x - eps_x, y, z)- 2.0*this->t(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddtdydy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->t(x, y + eps_x, z) + this->t(x, y - eps_x, z)- 2.0*this->t(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddtdzdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->t(x,y,z+2*eps_z) - 2*this->t(x,y,z+eps_z) + this->t(x,y,z)) / (pow(eps_z,2.0));
        } else {
		return (this->t(x, y, z + eps_z) + this->t(x, y, z - eps_z)- 2.0*this->t(x, y, z))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddtdydz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->t(x,y + eps_x,z + eps_z) - this->t(x, y - eps_x, z + eps_z) - this->t(x, y + eps_x, z) + this->t(x, y - eps_x, z)) / (2.0 * eps_x * eps_z);
        } else {
		return (this->t(x, y + eps_x, z + eps_z) - this->t(x, y - eps_x, z + eps_z) - this->t(x, y + eps_x, z - eps_z) + this->t(x, y - eps_x, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddtdxdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->t(x + eps_x,y,z + eps_z) - this->t(x - eps_x, y, z + eps_z) - this->t(x + eps_x, y, z) + this->t(x - eps_x, y, z) ) / (2.0 * eps_x * eps_z);
        } else {
		return (this->t(x + eps_x, y, z + eps_z) - this->t(x - eps_x, y, z + eps_z) - this->t(x + eps_x, y, z - eps_z) + this->t(x - eps_x, y, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddtdxdy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->t(x + eps_x, y + eps_x, z) - this->t(x - eps_x, y + eps_x, z) - this->t(x + eps_x, y - eps_x, z) + this->t(x - eps_x, y - eps_x, z))/(4.0*eps_x*eps_x);
}

// Spatial derivatives of zonal winds
double NCPA::AtmosphericSpecification::dudx( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->u(x+eps_x,y,z) - this->u(x-eps_x,y,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dudy( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->u(x,y+eps_x,z) - this->u(x,y-eps_x,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dudz( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->u(x,y,z+eps_z) - this->u(x,y,z)) / eps_z;
        } else {
		return (this->u(x,y,z+eps_z) - this->u(x,y,z-eps_z)) / (2.0*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddudxdx(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->u(x + eps_x, y, z) + this->u(x - eps_x, y, z)- 2.0*this->u(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddudydy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->u(x, y + eps_x, z) + this->u(x, y - eps_x, z)- 2.0*this->u(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddudzdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->u(x,y,z+2*eps_z) - 2*this->u(x,y,z+eps_z) + this->u(x,y,z)) / (pow(eps_z,2.0));
        } else {
		return (this->u(x, y, z + eps_z) + this->u(x, y, z - eps_z)- 2.0*this->u(x, y, z))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddudydz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->u(x,y + eps_x,z + eps_z) - this->u(x, y - eps_x, z + eps_z) - this->u(x, y + eps_x, z) + this->u(x, y - eps_x, z)) / (2.0 * eps_x * eps_z);
        } else {
		return (this->u(x, y + eps_x, z + eps_z) - this->u(x, y - eps_x, z + eps_z) - this->u(x, y + eps_x, z - eps_z) + this->u(x, y - eps_x, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddudxdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->u(x + eps_x,y,z + eps_z) - this->u(x - eps_x, y, z + eps_z) - this->u(x + eps_x, y, z) + this->u(x - eps_x, y, z) ) / (2.0 * eps_x * eps_z);
        } else {
		return (this->u(x + eps_x, y, z + eps_z) - this->u(x - eps_x, y, z + eps_z) - this->u(x + eps_x, y, z - eps_z) + this->u(x - eps_x, y, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddudxdy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->u(x + eps_x, y + eps_x, z) - this->u(x - eps_x, y + eps_x, z) - this->u(x + eps_x, y - eps_x, z) + this->u(x - eps_x, y - eps_x, z))/(4.0*eps_x*eps_x);
}


// Spatial derivatives of meridional winds
double NCPA::AtmosphericSpecification::dvdx( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->v(x+eps_x,y,z) - this->v(x-eps_x,y,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dvdy( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
    return (this->v(x,y+eps_x,z) - this->v(x,y-eps_x,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dvdz( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->v(x,y,z+eps_z) - this->v(x,y,z)) / eps_z;
        } else {
		return (this->v(x,y,z+eps_z) - this->v(x,y,z-eps_z)) / (2.0*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddvdxdx(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->v(x + eps_x, y, z) + this->v(x - eps_x, y, z)- 2.0*this->v(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddvdydy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->v(x, y + eps_x, z) + this->v(x, y - eps_x, z)- 2.0*this->v(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddvdzdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->v(x,y,z+2*eps_z) - 2*this->v(x,y,z+eps_z) + this->v(x,y,z)) / (pow(eps_z,2.0));
        } else {
		return (this->v(x, y, z + eps_z) + this->v(x, y, z - eps_z)- 2.0*this->v(x, y, z))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddvdydz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->v(x,y + eps_x,z + eps_z) - this->v(x, y - eps_x, z + eps_z) - this->v(x, y + eps_x, z) + this->v(x, y - eps_x, z)) / (2.0 * eps_x * eps_z);
        } else {
		return (this->v(x, y + eps_x, z + eps_z) - this->v(x, y - eps_x, z + eps_z) - this->v(x, y + eps_x, z - eps_z) + this->v(x, y - eps_x, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddvdxdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->v(x + eps_x,y,z + eps_z) - this->v(x - eps_x, y, z + eps_z) - this->v(x + eps_x, y, z) + this->v(x - eps_x, y, z) ) / (2.0 * eps_x * eps_z);
        } else {
		return (this->v(x + eps_x, y, z + eps_z) - this->v(x - eps_x, y, z + eps_z) - this->v(x + eps_x, y, z - eps_z) + this->v(x - eps_x, y, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddvdxdy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->v(x + eps_x, y + eps_x, z) - this->v(x - eps_x, y + eps_x, z) - this->v(x + eps_x, y - eps_x, z) + this->v(x - eps_x, y - eps_x, z))/(4.0*eps_x*eps_x);
}


// Spatial derivatives of vertical winds
double NCPA::AtmosphericSpecification::dwdx( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->w(x+eps_x,y,z) - this->w(x-eps_x,y,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dwdy( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->w(x,y+eps_x,z) - this->w(x,y-eps_x,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dwdz( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->w(x,y,z+eps_z) - this->w(x,y,z)) / eps_z;
        } else {
		return (this->w(x,y,z+eps_z) - this->w(x,y,z-eps_z)) / (2.0*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddwdxdx(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->w(x + eps_x, y, z) + this->w(x - eps_x, y, z)- 2.0*this->w(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddwdydy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->w(x, y + eps_x, z) + this->w(x, y - eps_x, z)- 2.0*this->w(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddwdzdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->w(x,y,z+2*eps_z) - 2*this->w(x,y,z+eps_z) + this->w(x,y,z)) / (pow(eps_z,2.0));
        } else {
		return (this->w(x, y, z + eps_z) + this->w(x, y, z - eps_z)- 2.0*this->w(x, y, z))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddwdydz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->w(x,y + eps_x,z + eps_z) - this->w(x, y - eps_x, z + eps_z) - this->w(x, y + eps_x, z) + this->w(x, y - eps_x, z)) / (2.0 * eps_x * eps_z);
        } else {
	        return (this->w(x, y + eps_x, z + eps_z) - this->w(x, y - eps_x, z + eps_z) - this->w(x, y + eps_x, z - eps_z) + this->w(x, y - eps_x, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddwdxdz(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->w(x + eps_x,y,z + eps_z) - this->w(x - eps_x, y, z + eps_z) - this->w(x + eps_x, y, z) + this->w(x - eps_x, y, z) ) / (2.0 * eps_x * eps_z);
        } else {
		return (this->w(x + eps_x, y, z + eps_z) - this->w(x - eps_x, y, z + eps_z) - this->w(x + eps_x, y, z - eps_z) + this->w(x - eps_x, y, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddwdxdy(double x, double y, double z) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
        return (this->w(x + eps_x, y + eps_x, z) - this->w(x - eps_x, y + eps_x, z) - this->w(x + eps_x, y - eps_x, z) + this->w(x - eps_x, y - eps_x, z))/(4.0*eps_x*eps_x);
}



// Spatial derivatives of pressure
double NCPA::AtmosphericSpecification::dpdx( double x, double y, double z ) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->p(x+eps_x,y,z) - this->p(x-eps_x,y,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dpdy( double x, double y, double z ) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->p(x,y+eps_x,z) - this->p(x,y-eps_x,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::dpdz( double x, double y, double z ) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->p(x,y,z+eps_z) - this->p(x,y,z)) / eps_z;
	} else {
		return (this->p(x,y,z+eps_z) - this->p(x,y,z-eps_z)) / (2.0*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddpdxdx(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->p(x + eps_x, y, z) + this->p(x - eps_x, y, z)- 2.0*this->p(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddpdydy(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->p(x, y + eps_x, z) + this->p(x, y - eps_x, z)- 2.0*this->p(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddpdzdz(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->p(x,y,z+2*eps_z) - 2*this->p(x,y,z+eps_z) + this->p(x,y,z)) / (pow(eps_z,2.0));
	} else {
		return (this->p(x, y, z + eps_z) + this->p(x, y, z - eps_z)- 2.0*this->p(x, y, z))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddpdydz(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->p(x,y + eps_x,z + eps_z) - this->p(x, y - eps_x, z + eps_z) - this->p(x, y + eps_x, z) + this->p(x, y - eps_x, z)) / (2.0 * eps_x * eps_z);
	} else {
		return (this->p(x, y + eps_x, z + eps_z) - this->p(x, y - eps_x, z + eps_z) - this->p(x, y + eps_x, z - eps_z) + this->p(x, y - eps_x, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddpdxdz(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->p(x + eps_x,y,z + eps_z) - this->p(x - eps_x, y, z + eps_z) - this->p(x + eps_x, y, z) + this->p(x - eps_x, y, z) ) / (2.0 * eps_x * eps_z);
	} else {
		return (this->p(x + eps_x, y, z + eps_z) - this->p(x - eps_x, y, z + eps_z) - this->p(x + eps_x, y, z - eps_z) + this->p(x - eps_x, y, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddpdxdy(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->p(x + eps_x, y + eps_x, z) - this->p(x - eps_x, y + eps_x, z) - this->p(x + eps_x, y - eps_x, z) + this->p(x - eps_x, y - eps_x, z))/(4.0*eps_x*eps_x);
}


// Spatial derivatives of vertical winds
double NCPA::AtmosphericSpecification::drhodx( double x, double y, double z ) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->rho(x+eps_x,y,z) - this->rho(x-eps_x,y,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::drhody( double x, double y, double z ) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->rho(x,y+eps_x,z) - this->rho(x,y-eps_x,z)) / (2.0*eps_x);
}
double NCPA::AtmosphericSpecification::drhodz( double x, double y, double z ) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->rho(x,y,z+eps_z) - this->rho(x,y,z)) / eps_z;
	} else {
		return (this->rho(x,y,z+eps_z) - this->rho(x,y,z-eps_z)) / (2.0*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddrhodxdx(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->rho(x + eps_x, y, z) + this->rho(x - eps_x, y, z)- 2.0*this->rho(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddrhodydy(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->rho(x, y + eps_x, z) + this->rho(x, y - eps_x, z)- 2.0*this->rho(x, y, z))/(pow(eps_x,2.0));
}
double NCPA::AtmosphericSpecification::ddrhodzdz(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->rho(x,y,z+2*eps_z) - 2*this->rho(x,y,z+eps_z) + this->rho(x,y,z)) / (pow(eps_z,2.0));
	} else {
		return (this->rho(x, y, z + eps_z) + this->rho(x, y, z - eps_z)- 2.0*this->rho(x, y, z))/(pow(eps_z,2.0));
	}
}
double NCPA::AtmosphericSpecification::ddrhodydz(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->rho(x,y + eps_x,z + eps_z) - this->rho(x, y - eps_x, z + eps_z) - this->rho(x, y + eps_x, z) + this->rho(x, y - eps_x, z)) / (2.0 * eps_x * eps_z);
	} else {
		return (this->rho(x, y + eps_x, z + eps_z) - this->rho(x, y - eps_x, z + eps_z) - this->rho(x, y + eps_x, z - eps_z) + this->rho(x, y - eps_x, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddrhodxdz(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	if (z-eps_z < z0(x,y)) {
		return (this->rho(x + eps_x,y,z + eps_z) - this->rho(x - eps_x, y, z + eps_z) - this->rho(x + eps_x, y, z) + this->rho(x - eps_x, y, z) ) / (2.0 * eps_x * eps_z);
	} else {
		return (this->rho(x + eps_x, y, z + eps_z) - this->rho(x - eps_x, y, z + eps_z) - this->rho(x + eps_x, y, z - eps_z) + this->rho(x - eps_x, y, z - eps_z))/(4.0*eps_x*eps_z);
	}
}
double NCPA::AtmosphericSpecification::ddrhodxdy(double x, double y, double z) {
	//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return (this->rho(x + eps_x, y + eps_x, z) - this->rho(x - eps_x, y + eps_x, z) - this->rho(x + eps_x, y - eps_x, z) + this->rho(x - eps_x, y - eps_x, z))/(4.0*eps_x*eps_x);
}


NCPA::AtmosphericPoint NCPA::AtmosphericSpecification::take_snapshot(double r1, double r2, double r3) {
	NCPA::AtmosphericPoint snap;
	snap.c0 = this->c0(r1,r2,r3);
	snap.c00 = this->c0(0,0,0);
	snap.u = this->u(r1,r2,r3);
	snap.v = this->v(r1,r2,r3);
	snap.w = this->w(r1,r2,r3);
	
	snap.dudx = this->dudz(r1,r2,r3);
	snap.dudy = this->dudy(r1,r2,r3);
	snap.dudz = this->dudz(r1,r2,r3);
	snap.dvdx = this->dvdx(r1,r2,r3);
	snap.dvdy = this->dvdy(r1,r2,r3);
	snap.dvdz = this->dvdz(r1,r2,r3);
	snap.dwdx = this->dwdx(r1,r2,r3);
	snap.dwdy = this->dwdy(r1,r2,r3);
	snap.dwdz = this->dwdz(r1,r2,r3);
	snap.dtdx = this->dwdx(r1,r2,r3);
	snap.dtdy = this->dwdy(r1,r2,r3);
	snap.dtdz = this->dwdz(r1,r2,r3);
	snap.dpdx = this->dwdx(r1,r2,r3);
	snap.dpdy = this->dwdy(r1,r2,r3);
	snap.dpdz = this->dwdz(r1,r2,r3);
	snap.drhodx = this->dwdx(r1,r2,r3);
	snap.drhody = this->dwdy(r1,r2,r3);
	snap.drhodz = this->dwdz(r1,r2,r3);
	snap.dc0dx = this->dc0dx(r1,r2,r3);
	snap.dc0dy = this->dc0dy(r1,r2,r3);
	snap.dc0dz = this->dc0dz(r1,r2,r3);
	
	snap.ddudxdx = this->ddudxdz(r1,r2,r3);
	snap.ddudxdy = this->ddudxdy(r1,r2,r3);
	snap.ddudxdz = this->ddudxdz(r1,r2,r3);
	snap.ddvdxdx = this->ddvdxdx(r1,r2,r3);
	snap.ddvdxdy = this->ddvdxdy(r1,r2,r3);
	snap.ddvdxdz = this->ddvdxdz(r1,r2,r3);
	snap.ddwdxdx = this->ddwdxdx(r1,r2,r3);
	snap.ddwdxdy = this->ddwdxdy(r1,r2,r3);
	snap.ddwdxdz = this->ddwdxdz(r1,r2,r3);
	snap.ddc0dxdx = this->ddc0dxdx(r1,r2,r3);
	snap.ddc0dxdy = this->ddc0dxdy(r1,r2,r3);
	snap.ddc0dxdz = this->ddc0dxdz(r1,r2,r3);
	
	snap.ddudydy = this->ddudydy(r1,r2,r3);
	snap.ddudydz = this->ddudydz(r1,r2,r3);
	snap.ddvdydy = this->ddvdydy(r1,r2,r3);
	snap.ddvdydz = this->ddvdydz(r1,r2,r3);
	snap.ddwdydy = this->ddwdydy(r1,r2,r3);
	snap.ddwdydz = this->ddwdydz(r1,r2,r3);
	snap.ddc0dydy = this->ddc0dydy(r1,r2,r3);
	snap.ddc0dydz = this->ddc0dydz(r1,r2,r3);
	
	snap.ddudzdz = this->ddudzdz(r1,r2,r3);
	snap.ddvdzdz = this->ddvdzdz(r1,r2,r3);
	snap.ddwdzdz = this->ddwdzdz(r1,r2,r3);
	snap.ddc0dzdz = this->ddc0dzdz(r1,r2,r3);
	
	snap.ddtdxdx = this->ddudxdz(r1,r2,r3);
	snap.ddtdxdy = this->ddudxdy(r1,r2,r3);
	snap.ddtdxdz = this->ddudxdz(r1,r2,r3);
	snap.ddtdydy = this->ddudydy(r1,r2,r3);
	snap.ddtdydz = this->ddudydz(r1,r2,r3);
	snap.ddtdzdz = this->ddudzdz(r1,r2,r3);
	
	snap.ddpdxdx = this->ddudxdz(r1,r2,r3);
	snap.ddpdxdy = this->ddudxdy(r1,r2,r3);
	snap.ddpdxdz = this->ddudxdz(r1,r2,r3);
	snap.ddpdydy = this->ddudydy(r1,r2,r3);
	snap.ddpdydz = this->ddudydz(r1,r2,r3);
	snap.ddpdzdz = this->ddudzdz(r1,r2,r3);
	
	snap.ddrhodxdx = this->ddudxdz(r1,r2,r3);
	snap.ddrhodxdy = this->ddudxdy(r1,r2,r3);
	snap.ddrhodxdz = this->ddudxdz(r1,r2,r3);
	snap.ddrhodydy = this->ddudydy(r1,r2,r3);
	snap.ddrhodydz = this->ddudydz(r1,r2,r3);
	snap.ddrhodzdz = this->ddudzdz(r1,r2,r3);
	
	return snap;
}

