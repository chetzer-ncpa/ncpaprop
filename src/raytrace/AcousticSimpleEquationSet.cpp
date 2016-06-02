#include "AcousticSimpleEquationSet.h"
#include "AtmosphericSpecification.h"
#include <cmath>

NCPA::AcousticSimpleEquationSet::AcousticSimpleEquationSet( NCPA::AtmosphericSpecification *p, 
		double takeoffangle, double propagation_azimuth ) {
	profile = p;
	azimuth = propagation_azimuth;
	takeoff = takeoffangle;
}

int NCPA::AcousticSimpleEquationSet::numberOfEquations() const { return 3; }

void NCPA::AcousticSimpleEquationSet::changeTakeoffAngle( double newAngle ) {
	takeoff = newAngle;
}

void NCPA::AcousticSimpleEquationSet::changeAzimuth( double newAngle ) {
	azimuth = newAngle;
}

double NCPA::AcousticSimpleEquationSet::result( double s, double *current_values, int which ) {

	// Identify the current values by name, for ease of reading/debugging
	//double r = current_values[ 0 ];
	double z = current_values[ 1 ];
	double Zeta = current_values[ 2 ];

	double c = this->profile->ceff( 0.0, 0.0, z, this->azimuth );
	double c0 = this->profile->ceff( 0.0, 0.0, 0.0, this->azimuth );
	double dc = this->profile->dceffdz( 0.0, 0.0, z, this->azimuth );

	switch (which) {
		case 0:
			return std::cos( this->takeoff ) * c / c0;
		case 1:
			return c * Zeta;
		case 2:
			return -dc / std::pow( c, 2 );
		default:
			throw new NCPA::RequestOutOfBoundsException( "Only 3 equations defined!" );
	}
}

