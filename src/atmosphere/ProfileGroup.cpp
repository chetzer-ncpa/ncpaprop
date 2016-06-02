#include "ProfileGroup.h"
#include "AtmosphericProfile.h"
#include "geographic.h"
#include <stdexcept>

NCPA::ProfileGroup::~ProfileGroup() { }

NCPA::AtmosphericProfile* NCPA::ProfileGroup::getProfile( double lat, double lon, bool exact ) {

	// first check to see if the location has been checked before
	if ( !buffer_.empty() ) {
		for ( std::vector< NCPA::AtmosphericProfile * >::const_iterator it = buffer_.begin(); it != buffer_.end(); it++ ) {
			if (exact) {
				if (lat == ( *it )->lat() && lon == ( *it )->lon())
					return *it;
			} else {
				double range = NCPA::range ( lat, lon, ( *it )->lat(), ( *it )->lon() );
				if ( range < eps_x )
					return *it;
			}
		}
	}

	double minrange = 999999999999999;

	double currange;
	int index = -1;

	for ( unsigned int i = 0; i < profiles_.size(); i++ ) {
			currange = NCPA::range ( lat, lon, profiles_[i]->lat(), profiles_[i]->lon() );

			if ( currange < minrange ) {
					minrange = currange;
					index = i;
				}
		}

	buffer_.push_back ( profiles_[ index ] );

	return profiles_[ index ];
}

double NCPA::ProfileGroup::z0 ( double x, double y ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->z0();
}

double NCPA::ProfileGroup::t ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->t ( z );
}

double NCPA::ProfileGroup::u ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->u ( z );
}

double NCPA::ProfileGroup::v ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->v (z );
}

double NCPA::ProfileGroup::w ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->w ( z );

}

// Should default to US Std Atmosphere
double NCPA::ProfileGroup::rho ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->rho ( z );
}

double NCPA::ProfileGroup::p ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->p ( z );
}

void NCPA::ProfileGroup::setOrigin ( double lat0, double lon0 ) {
	lat0_ = lat0;
	lon0_ = lon0;
}

double NCPA::ProfileGroup::dtdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->dtdz ( z );
}

double NCPA::ProfileGroup::dudz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->dudz ( z );
}

double NCPA::ProfileGroup::dvdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->dvdz ( z );
}

double NCPA::ProfileGroup::dwdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->dwdz ( z );
}

double NCPA::ProfileGroup::dpdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->dpdz ( z );
}

double NCPA::ProfileGroup::drhodz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->drhodz ( z );
}

double NCPA::ProfileGroup::dceffdz( double x, double y, double z, double phi ) {
	NCPA::Location ll = NCPA::xy2latlon( x, y, lat0_, lon0_ );
	return getProfile( ll.lat(), ll.lon() )->dceffdz( z, phi );
}

double NCPA::ProfileGroup::dc0dz( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon( x, y, lat0_, lon0_ );
	return getProfile( ll.lat(), ll.lon() )->dc0dz( z );
}

double NCPA::ProfileGroup::ddceffdzdz ( double x, double y, double z, double phi ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddceffdzdz ( z, phi );
}

double NCPA::ProfileGroup::ddc0dzdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddc0dzdz ( z );
}

double NCPA::ProfileGroup::ddtdzdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddtdzdz ( z );
}

double NCPA::ProfileGroup::ddudzdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddudzdz ( z );
}

double NCPA::ProfileGroup::ddvdzdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddvdzdz ( z );
}

double NCPA::ProfileGroup::ddwdzdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddwdzdz ( z );
}

double NCPA::ProfileGroup::ddpdzdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddpdzdz ( z );
}

double NCPA::ProfileGroup::ddrhodzdz ( double x, double y, double z ) {
	NCPA::Location ll = NCPA::xy2latlon ( x, y, lat0_, lon0_ );
	return getProfile ( ll.lat(), ll.lon() )->ddrhodzdz ( z );
}

// Spatial derivatives of pressure
double NCPA::ProfileGroup::dpdx ( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return ( this->p ( x + eps_x, y, z ) - this->p ( x - eps_x, y, z ) ) / ( 2.0*eps_x );
}

double NCPA::ProfileGroup::dpdy ( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return ( this->p ( x, y + eps_x, z ) - this->p ( x, y - eps_x, z ) ) / ( 2.0*eps_x );
}

// Spatial derivatives of density
double NCPA::ProfileGroup::drhodx ( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return ( this->rho ( x + eps_x, y, z ) - this->rho ( x - eps_x, y, z ) ) / ( 2.0*eps_x );
}

double NCPA::ProfileGroup::drhody ( double x, double y, double z ) {
//        if (z < z0(x,y)) MY_THROW( "Invalid z value requested!" );
	return ( this->rho ( x, y + eps_x, z ) - this->rho ( x, y - eps_x, z ) ) / ( 2.0*eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddtdxdz ( double x, double y, double z ) {
	return ( this->dtdz ( x + eps_x, y, z ) - this->dtdz ( x - eps_x, y, z ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddtdydz ( double x, double y, double z ) {
	return ( this->dtdz ( x, y + eps_x, z ) - this->dtdz ( x, y - eps_x, z ) ) / ( 2.0 * eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddudxdz ( double x, double y, double z ) {
	return ( this->dudz ( x + eps_x, y, z ) - this->dudz ( x - eps_x, y, z ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddudydz ( double x, double y, double z ) {
	return ( this->dudz ( x, y + eps_x, z ) - this->dudz ( x, y - eps_x, z ) ) / ( 2.0 * eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddvdxdz ( double x, double y, double z ) {
	return ( this->dvdz ( x + eps_x, y, z ) - this->dvdz ( x - eps_x, y, z ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddvdydz ( double x, double y, double z ) {
	return ( this->dvdz ( x, y + eps_x, z ) - this->dvdz ( x, y - eps_x, z ) ) / ( 2.0 * eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddwdxdz ( double x, double y, double z ) {
	return ( this->dwdz ( x + eps_x, y, z ) - this->dwdz ( x - eps_x, y, z ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddwdydz ( double x, double y, double z ) {
	return ( this->dwdz ( x, y + eps_x, z ) - this->dwdz ( x, y - eps_x, z ) ) / ( 2.0 * eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddpdxdz ( double x, double y, double z ) {
	return ( this->dpdz ( x + eps_x, y, z ) - this->dpdz ( x - eps_x, y, z ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddpdydz ( double x, double y, double z ) {
	return ( this->dpdz ( x, y + eps_x, z ) - this->dpdz ( x, y - eps_x, z ) ) / ( 2.0 * eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddrhodxdz ( double x, double y, double z ) {
	return ( this->drhodz ( x + eps_x, y, z ) - this->drhodz ( x - eps_x, y, z ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddrhodydz ( double x, double y, double z ) {
	return ( this->drhodz ( x, y + eps_x, z ) - this->drhodz ( x, y - eps_x, z ) ) / ( 2.0 * eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddceffdxdz ( double x, double y, double z, double phi ) {
	return ( this->dceffdz ( x + eps_x, y, z, phi ) - this->dceffdz ( x - eps_x, y, z, phi ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddceffdydz ( double x, double y, double z, double phi ) {
	return ( this->dceffdz ( x, y + eps_x, z, phi ) - this->dceffdz ( x, y - eps_x, z, phi ) ) / ( 2.0 * eps_x );
}

// Horizontal mixed derivatives
double NCPA::ProfileGroup::ddc0dxdz ( double x, double y, double z ) {
	return ( this->dc0dz ( x + eps_x, y, z ) - this->dc0dz ( x - eps_x, y, z ) ) / ( 2.0 * eps_x );
}

double NCPA::ProfileGroup::ddc0dydz ( double x, double y, double z ) {
	return ( this->dc0dz ( x, y + eps_x, z ) - this->dc0dz ( x, y - eps_x, z ) ) / ( 2.0 * eps_x );
}
