#include "Slice.h"
#include "util.h"
#include "geographic.h"
#include "SampledProfile.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sstream>

void NCPA::Slice::init_() {
	good_ = false;
	eps_x = 0;
	eps_z = 0.01;
	
	profiles_.clear();
	pathAz_ = -1;
	ranges_.clear();
	strict_ = false;
	azTolerance_ = 5;
}

void NCPA::Slice::clearOut() {
	for (unsigned int i = 0; i < profiles_.size(); i++) {
		delete profiles_[ i ];
	}
	profiles_.clear();
	ranges_.clear();
}

NCPA::Slice::~Slice() {
	this->clearOut();
}

NCPA::Slice::Slice() {
	this->init_();
}

double NCPA::Slice::pathAzimuth() const {
	return pathAz_;
}

void NCPA::Slice::readSummaryFile( std::string summaryFile ) {
	this->clearOut();
	this->init_();
	
	// File should have the following format: 
	// lat <origin_lat>
	// lon <origin_lon>
	// azimuth <path_azimuth>
	// order <field_order>
	// header <skiplines>
	// <range>  <filename>
	// <range>  <filename>
	// ...
	std::ifstream summary( summaryFile.c_str(), std::ifstream::in );
	int headerLines = 0;
	double lat = 0, lon = 0;
	char fieldOrder[256];
	std::memset(fieldOrder,'\0',256);
	
	// read the first five lines
	bool hasLat = false, hasLon = false, hasAzimuth = false, hasOrder = false, hasHeader = false;
	std::string lineHolder;
	std::vector< std::string > tokens;
	while (!(hasLat && hasLon && hasAzimuth && hasOrder && hasHeader)) {
		NCPA::safe_getline( summary, lineHolder );
		lineHolder = NCPA::deblank( lineHolder );
		if (lineHolder[0] != '#' && lineHolder.length() > 0) {
			tokens = NCPA::split( lineHolder, " ,:\t" );
			if (tokens.size() != 2) {
				std::runtime_error e( "Malformed header line in Slice summary file." );
				throw e;
			}
			
			if (tokens[0].compare("lat") == 0) {
				lat = std::atof( tokens[1].c_str() );
				hasLat = true;
			} else if (tokens[0].compare("lon") == 0) {
				lon = std::atof( tokens[1].c_str() );
				hasLon = true;
			} else if (tokens[0].compare("azimuth") == 0) {
				pathAz_ = std::atof( tokens[1].c_str() );
				hasAzimuth = true;
			} else if (tokens[0].compare("order") == 0) {
				std::strcpy( fieldOrder, tokens[1].c_str() );
				hasOrder = true;
			} else if (tokens[0].compare("header") == 0) {
				headerLines = std::atoi( tokens[1].c_str() );
				hasHeader = true;
			} else {
				std::ostringstream oss("");
				oss << "Unknown header value " << tokens[0] << " specified in Slice summary file.";
				std::runtime_error e( oss.str() );
				throw e;
			}
		}
	}
	
	// Got all the header information, now we read in the lines
	NCPA::AtmosphericProfile *pro;
	while (NCPA::safe_getline( summary, lineHolder )) {
		lineHolder = deblank( lineHolder );
		if (lineHolder.length() > 0) {
			tokens = NCPA::split( lineHolder, " :,\t" );
			if (tokens.size() != 2) {
				std::runtime_error e( "Malformed data line in Slice summary file." );
				throw e;
			}
			double r = std::atof( tokens[ 0 ].c_str() );
			pro = new SampledProfile( tokens[ 1 ], fieldOrder, headerLines );
			pro->setOrigin(lat,lon);
			profiles_.push_back( pro );
			ranges_.push_back( r );
		}
	}
	good_ = true;
}



// Find the profiles on either side of the requested range
void NCPA::Slice::getBrackets_( double r, int &nearIndex, int &farIndex ) {
	
	if (r < 0) {
		std::runtime_error e("Range supplied to getBrackets_ cannot be negative!" );
		throw e;
	}
	
	nearIndex = -1;
	farIndex = -1;
	double negdiff = 1e15, posdiff = 1e15;
	for (unsigned int i = 0; i < profiles_.size(); i++) {
		double diff = ranges_[ i ] - r;
		
		// negative diff means the point is closer to the origin than r, so it goes in nearIndex
		if (diff < 0 && -diff < negdiff) {
			negdiff = -diff;
			nearIndex = i;
		// positive diff means the point is farther from the origin than r, so it goes in farIndex
		} else if (diff >= 0 && diff < posdiff) {
			posdiff = diff;
			farIndex = i;
		}
	}
	
}

double NCPA::Slice::t( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->t(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->t(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->t(z),
					ranges_[farIndex],profiles_[farIndex]->t(z),
					r );
	}
}

double NCPA::Slice::u( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->u(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->u(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->u(z),
					ranges_[farIndex],profiles_[farIndex]->u(z),
					r );
	}
}

double NCPA::Slice::v( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->v(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->v(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->v(z),
					ranges_[farIndex],profiles_[farIndex]->v(z),
					r );
	}
}

double NCPA::Slice::w( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->w(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->w(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->w(z),
					ranges_[farIndex],profiles_[farIndex]->w(z),
					r );
	}
}

double NCPA::Slice::p( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->p(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->p(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->p(z),
					ranges_[farIndex],profiles_[farIndex]->p(z),
					r );
	}
}

double NCPA::Slice::rho( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->rho(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->rho(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->rho(z),
					ranges_[farIndex],profiles_[farIndex]->rho(z),
					r );
	}
}

bool NCPA::Slice::stratified() const { return false; }

double NCPA::Slice::z0( double x, double y ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->z0();
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->z0();
	} else {
		return NCPA::max( profiles_[nearIndex]->z0(), profiles_[farIndex]->z0() );
	}
}

NCPA::AtmosphericProfile *NCPA::Slice::getProfile( double lat, double lon, bool exact ) {
	int ind = -1;
	int minRange = 1e5;
	for (unsigned int i = 0; i < profiles_.size(); i++) {
		double r = NCPA::range( lat, lon, profiles_[i]->lat(), profiles_[i]->lon() );
		if ( r < minRange ) {
			minRange = r;
			ind = i;
		}
	}
	return profiles_[ ind ];
}

void NCPA::Slice::verticalTolerance(double eps_z) {
	this->eps_z = eps_z;
}

void NCPA::Slice::strict( bool s ) {
	strict_ = s;
}

double NCPA::Slice::dudz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->dudz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->dudz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->dudz(z),
					ranges_[farIndex],profiles_[farIndex]->dudz(z),
					r );
	}
}

double NCPA::Slice::ddudzdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddudzdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddudzdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddudzdz(z),
					ranges_[farIndex],profiles_[farIndex]->ddudzdz(z),
					r );
	}
}

double NCPA::Slice::dtdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->dtdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->dtdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->dtdz(z),
					ranges_[farIndex],profiles_[farIndex]->dtdz(z),
					r );
	}
}

double NCPA::Slice::ddtdzdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddtdzdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddtdzdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddtdzdz(z),
					ranges_[farIndex],profiles_[farIndex]->ddtdzdz(z),
					r );
	}
}

double NCPA::Slice::dvdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->dvdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->dvdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->dvdz(z),
					ranges_[farIndex],profiles_[farIndex]->dvdz(z),
					r );
	}
}

double NCPA::Slice::ddvdzdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddvdzdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddvdzdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddvdzdz(z),
					ranges_[farIndex],profiles_[farIndex]->ddvdzdz(z),
					r );
	}
}

double NCPA::Slice::dwdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->dwdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->dwdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->dwdz(z),
					ranges_[farIndex],profiles_[farIndex]->dwdz(z),
					r );
	}
}

double NCPA::Slice::ddwdzdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddwdzdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddwdzdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddwdzdz(z),
					ranges_[farIndex],profiles_[farIndex]->ddwdzdz(z),
					r );
	}
}

double NCPA::Slice::dpdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->dpdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->dpdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->dpdz(z),
					ranges_[farIndex],profiles_[farIndex]->dpdz(z),
					r );
	}
}

double NCPA::Slice::ddpdzdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddpdzdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddpdzdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddpdzdz(z),
					ranges_[farIndex],profiles_[farIndex]->ddpdzdz(z),
					r );
	}
}

double NCPA::Slice::drhodz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->drhodz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->drhodz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->drhodz(z),
					ranges_[farIndex],profiles_[farIndex]->drhodz(z),
					r );
	}
}

double NCPA::Slice::ddrhodzdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddrhodzdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddrhodzdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddrhodzdz(z),
					ranges_[farIndex],profiles_[farIndex]->ddrhodzdz(z),
					r );
	}
}

double NCPA::Slice::dc0dz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->dc0dz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->dc0dz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->dc0dz(z),
					ranges_[farIndex],profiles_[farIndex]->dc0dz(z),
					r );
	}
}

double NCPA::Slice::ddc0dzdz( double x, double y, double z ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddc0dzdz(z);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddc0dzdz(z);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddc0dzdz(z),
					ranges_[farIndex],profiles_[farIndex]->ddc0dzdz(z),
					r );
	}
}

double NCPA::Slice::dceffdz( double x, double y, double z, double phi ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->dceffdz(z,phi);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->dceffdz(z,phi);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->dceffdz(z,phi),
					ranges_[farIndex],profiles_[farIndex]->dceffdz(z,phi),
					r );
	}
}

double NCPA::Slice::ddceffdzdz( double x, double y, double z, double phi ) {
	// Compute r from x and y.  Assume r is along the path of profiles
	double r = std::sqrt( x*x + y*y );
	if (strict_ && r > 0) {
		// check that r is along the path of profiles
		double requestedAz = 90 - NCPA::rad2deg(std::atan2( y, x ));
		if (!NCPA::checkAzimuthLimits(requestedAz, pathAz_, azTolerance_)) {
			std::cerr << "Warning: Requested point (" << x << "," << y << ") has azimuth " << requestedAz
				<< ", which differs from " << pathAz_ << " by more than " << azTolerance_ << " degrees." << std::endl;
 		}
	}
	
	// find the brackets for this r
	int nearIndex, farIndex;
	getBrackets_( r, nearIndex, farIndex );
	if (nearIndex < 0) {
		return profiles_[farIndex]->ddceffdzdz(z,phi);
	} else if (farIndex < 0) {
		return profiles_[nearIndex]->ddceffdzdz(z,phi);
	} else {
		return NCPA::linearInterp(	ranges_[nearIndex],profiles_[nearIndex]->ddceffdzdz(z,phi),
					ranges_[farIndex],profiles_[farIndex]->ddceffdzdz(z,phi),
					r );
	}
}

