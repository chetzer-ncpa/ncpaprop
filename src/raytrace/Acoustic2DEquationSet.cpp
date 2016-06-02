#include "Acoustic2DEquationSet.h"
#include "AcousticEquationSet.h"
#include "AtmosphericSpecification.h"
#include <cmath>
#include "util.h"

#ifndef Pi
#define Pi 3.14159
#endif

NCPA::Acoustic2DEquationSet::Acoustic2DEquationSet( NCPA::AtmosphericSpecification *p, 
		double takeoffangle, double propagation_azimuth,
		bool rangeDependent ) {
	profile = p;
	azimuth = propagation_azimuth;
	takeoff = takeoffangle;
	rangeDependent_ = rangeDependent;
}

int NCPA::Acoustic2DEquationSet::numberOfEquations() const { return 6; }

void NCPA::Acoustic2DEquationSet::setupInitialConditions( double *initialConditions, 
	double zmin, double c0 ) const {
	initialConditions[ 0 ] = 0.0;
	initialConditions[ 1 ] = zmin;
	initialConditions[ 2 ] = std::sin( takeoff ) / c0;
	initialConditions[ 3 ] = 0.0;
	initialConditions[ 4 ] = 0.0;
	initialConditions[ 5 ] = std::cos( takeoff ) / c0;
}

void NCPA::Acoustic2DEquationSet::setupInitialConditions( double *initialConditions, double zmin ) {
	double c0 = profile->ceff( 0, 0, zmin, azimuth );
	initialConditions[ 0 ] = 0.0;
	initialConditions[ 1 ] = zmin;
	initialConditions[ 2 ] = std::sin( takeoff ) / c0;
	initialConditions[ 3 ] = 0.0;
	initialConditions[ 4 ] = 0.0;
	initialConditions[ 5 ] = std::cos( takeoff ) / c0;
}

void NCPA::Acoustic2DEquationSet::results( double s, double *current_values, double *new_values ) {
	// Identify the current values by name, for ease of reading/debugging
	double r = current_values[ 0 ];
	double z = current_values[ 1 ];
	double Zeta = current_values[ 2 ];
	//double rt = current_values[ 3 ];
	double zt = current_values[ 4 ];
	double Zetat = current_values[ 5 ];
	
	// Calculate x and y values from r and azimuth
	double cart_angle = NCPA::deg2rad( 90 - this->azimuth );
	double x = 0.0;
	double y = 0.0;
	if (rangeDependent_) {
		x = r * std::cos( cart_angle );
		y = r * std::sin( cart_angle );
	}
	
	double c = this->profile->ceff( x, y, z, this->azimuth );
	double c0 = this->profile->ceff( x, y, profile->z0(x,y), this->azimuth );
	double dc = this->profile->dceffdz( x, y, z, this->azimuth );
	double ddc = this->profile->ddceffdzdz( x, y, z, this->azimuth );
	
	new_values[ 0 ] = std::cos( this->takeoff ) * c / c0;
	new_values[ 1 ] = c * Zeta;
	new_values[ 2 ] = -dc / std::pow( c, 2 );
	new_values[ 3 ] = (dc * zt * std::cos( this->takeoff ) / c0 )
				- (std::sin( this->takeoff ) * c / c0 );
	new_values[ 4 ] = dc * zt * Zeta + c * Zetat;
	new_values[ 5 ] = (2.0 / std::pow( c, 3.0 ) * std::pow( dc, 2.0 ) * zt)
				- (ddc / std::pow( c, 2.0 ) * zt);
	
}
	
	
/*
double NCPA::Acoustic2DEquationSet::result( double s, double *current_values, int which ) {

	// Identify the current values by name, for ease of reading/debugging
	double r = current_values[ 0 ];
	double z = current_values[ 1 ];
	double Zeta = current_values[ 2 ];
	//double rt = current_values[ 3 ];
	double zt = current_values[ 4 ];
	double Zetat = current_values[ 5 ];

	// Calculate x and y values from r and azimuth
	double cart_angle = NCPA::deg2rad( 90 - this->azimuth );
	double x = 0.0;
	double y = 0.0;
	if (rangeDependent_) {
		x = r * std::cos( cart_angle );
		y = r * std::sin( cart_angle );
	}

	double c = this->profile->ceff( x, y, z, this->azimuth );
	double c0 = this->profile->ceff( x, y, profile->z0(x,y), this->azimuth );
	double dc = this->profile->dceffdz( x, y, z, this->azimuth );
	double ddc = this->profile->ddceffdzdz( x, y, z, this->azimuth );

	switch (which) {
		case 0:
			return std::cos( this->takeoff ) * c / c0;
		case 1:
			return c * Zeta;
		case 2:
			return -dc / std::pow( c, 2 );
		case 3:
			return (dc * zt * std::cos( this->takeoff ) / c0 )
				- (std::sin( this->takeoff ) * c / c0 );
		case 4:
			return dc * zt * Zeta + c * Zetat;
		case 5:
			return (2.0 / std::pow( c, 3.0 ) * std::pow( dc, 2.0 ) * zt)
				- (ddc / std::pow( c, 2.0 ) * zt);
		default:
			throw new NCPA::RequestOutOfBoundsException( "Only 6 equations defined!" );
	}
}
*/

double NCPA::Acoustic2DEquationSet::calculateTravelTime( double **solution, int steps, double step_size, 
		double azimuth ) const {
	double tau = 0.0;
	double x, y;
	double phi = NCPA::deg2rad( 90 - azimuth );
	for (int i = 0; i < steps; i++) {
		x = solution[i][0] * std::cos( phi );
		y = solution[i][0] * std::sin( phi );
		tau += step_size / profile->ceff(x,y,solution[i][1],azimuth);
	}
	return tau;
}

double NCPA::Acoustic2DEquationSet::calculateAmplitude( double **solution, int steps ) const {

	double phi = NCPA::deg2rad( 90 - azimuth );
	double r = solution[steps][0];
	double x = r * std::cos( phi );
	double y = r * std::sin( phi );
	double z = solution[steps][1];
	double c = profile->ceff(x,y,z,azimuth);
	double c0 = profile->ceff(x,y,0,azimuth);

	double drds_dzdtheta = solution[steps][4] * std::cos( takeoff )	* c / c0;
	double dzds_drdtheta = c * solution[steps][2] * solution[steps][3];

	double D = r * (drds_dzdtheta - dzds_drdtheta);
	double amp = 1.0 / (4.0 * Pi) * std::sqrt( 
			std::fabs( profile->rho(x,y,z) * c * std::cos(takeoff) / (profile->rho(x,y,0.0) * c0 * D) ) 
						 );
	return amp;
}

double NCPA::Acoustic2DEquationSet::transmissionLoss( double **soln, int k, double &refDist ) const {
	
	// iterate through ranges until you get past refDist km.  Since we expect refDist to occur very early
	// in the solution, this is probably comparably efficient to a bisection search
	int ceiling = 0;
	double range = soln[ceiling][0];
	double bottomrange;
	while (range < refDist && ceiling <= k) {
		ceiling++;
		bottomrange = range;
		range = soln[ ceiling ][ 0 ];
	}
	int index;
	
	if (std::fabs( bottomrange - refDist ) < std::fabs( range - refDist )) {
		index = ceiling-1;
		refDist = bottomrange;
	} else {
		index = ceiling;
		refDist = range;
	}
	double A0 = this->calculateAmplitude(soln, index);
	double A = this->calculateAmplitude(soln, k);
	return 20*std::log10(A/A0);
}