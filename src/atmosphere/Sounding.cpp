#include <cmath>
#include "geographic.h"
#include "Sounding.h"
#include "SampledProfile.h"
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <fstream>

NCPA::Sounding::Sounding() {
	init_();
}

void NCPA::Sounding::init_() {
	good_ = false;
	profile_ = 0;
}

void NCPA::Sounding::clearOut() {
        delete profile_;
	good_ = false;
}

NCPA::Sounding::~Sounding() {
	clearOut();
}

NCPA::Sounding::Sounding( AtmosphericProfile *p ) {
	init_();
	profile_ = p;
	good_ = true;
}

NCPA::Sounding::Sounding( int nz, double *z, double *t, double *u, double *v, double *w, 
			double *rho, double *p, double z0 ) {
        init_();
	profile_ = new SampledProfile( nz, z, t, u, v, w, rho, p, z0 );
        good_ = profile_->good();
}

NCPA::Sounding::Sounding( double lat, double lon, int nz, double *z, double *t, double *u, double *v, double *w, 
			double *rho, double *p, double z0 ) {
        init_();
	profile_ = new SampledProfile( lat, lon, nz, z, t, u, v, w, rho, p, z0 );
        good_ = profile_->good();
}


bool NCPA::Sounding::stratified() const { 
	return true; 
}

NCPA::AtmosphericProfile *NCPA::Sounding::getProfile( double lat, double lon, bool exact ) {
	return profile_;
}

double NCPA::Sounding::lat() const { 
	return profile_->lat(); 
}
double NCPA::Sounding::lon() const { 
	return profile_->lon(); 
}
double NCPA::Sounding::z0( double x, double y ) { 
	return profile_->z0();
}

double NCPA::Sounding::t( double x, double y, double z ) {
	return profile_->t( z );
}

double NCPA::Sounding::u( double x, double y, double z ) {
	return profile_->u( z );
}

double NCPA::Sounding::v( double x, double y, double z ) {
	return profile_->v( z );
}

double NCPA::Sounding::w( double x, double y, double z ) {
	return profile_->w( z );
}

double NCPA::Sounding::p( double x, double y, double z ) {
	return profile_->p( z );
}

double NCPA::Sounding::rho( double x, double y, double z ) {
	return profile_->rho( z );
}

double NCPA::Sounding::c0( double x, double y, double z ) {
	return profile_->c0( z );
}

double NCPA::Sounding::dtdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dtdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddtdxdz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddtdydz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddtdxdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddtdydy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddtdxdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dceffdx( double x, double y, double z, double phi ) { return 0.0; }
double NCPA::Sounding::dceffdy( double x, double y, double z, double phi ) { return 0.0; }
double NCPA::Sounding::ddceffdxdz( double x, double y, double z, double phi ) { return 0.0; }
double NCPA::Sounding::ddceffdydz( double x, double y, double z, double phi ) { return 0.0; }
double NCPA::Sounding::ddceffdxdx( double x, double y, double z, double phi ) { return 0.0; }
double NCPA::Sounding::ddceffdydy( double x, double y, double z, double phi ) { return 0.0; }
double NCPA::Sounding::ddceffdxdy( double x, double y, double z, double phi ) { return 0.0; }
double NCPA::Sounding::dc0dx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dc0dy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddc0dydy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddc0dxdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddc0dxdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddc0dxdz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddc0dydz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dudx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dudy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddudxdz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddudydz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddudxdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddudydy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddudxdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dvdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dvdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddvdxdz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddvdydz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddvdxdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddvdydy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddvdxdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dwdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dwdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddwdxdz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddwdydz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddwdxdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddwdydy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddwdxdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dpdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::dpdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddpdxdz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddpdydz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddpdxdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddpdydy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddpdxdy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::drhodx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::drhody( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddrhodxdz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddrhodydz( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddrhodxdx( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddrhodydy( double x, double y, double z ) { return 0.0; }
double NCPA::Sounding::ddrhodxdy( double x, double y, double z ) { return 0.0; }


double NCPA::Sounding::dtdz( double x, double y, double z ) {
	return profile_->dtdz( z );
}

double NCPA::Sounding::dudz( double x, double y, double z ) {
	return profile_->dudz( z );
}
double NCPA::Sounding::dvdz( double x, double y, double z ) {
	return profile_->dvdz( z );
}
double NCPA::Sounding::dwdz( double x, double y, double z ) {
	return profile_->dwdz( z );
}
double NCPA::Sounding::dpdz( double x, double y, double z ) {
	return profile_->dpdz( z );
}
double NCPA::Sounding::drhodz( double x, double y, double z ) {
	return profile_->drhodz( z );
}
double NCPA::Sounding::dceffdz( double x, double y, double z, double phi ) {
	return profile_->dceffdz( z, phi );
}
double NCPA::Sounding::dc0dz( double x, double y, double z ) {
	return profile_->dc0dz( z );
}


double NCPA::Sounding::ddtdzdz( double x, double y, double z ) {
	return profile_->ddtdzdz( z );
}
double NCPA::Sounding::ddudzdz( double x, double y, double z ) {
	return profile_->ddudzdz( z );
}
double NCPA::Sounding::ddvdzdz( double x, double y, double z ) {
	return profile_->ddvdzdz( z );
}
double NCPA::Sounding::ddwdzdz( double x, double y, double z ) {
	return profile_->ddwdzdz( z );
}
double NCPA::Sounding::ddpdzdz( double x, double y, double z ) {
	return profile_->ddpdzdz( z );
}
double NCPA::Sounding::ddrhodzdz( double x, double y, double z ) {
	return profile_->ddrhodzdz( z );
}
double NCPA::Sounding::ddceffdzdz( double x, double y, double z, double phi ) {
	return profile_->ddceffdzdz( z, phi );
}
double NCPA::Sounding::ddc0dzdz( double x, double y, double z ) {
	return profile_->ddc0dzdz( z );
}

