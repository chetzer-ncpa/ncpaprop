#include "AcousticEquationSet.h"

void NCPA::AcousticEquationSet::changeAzimuth( double newAngle ) {
	azimuth = newAngle;
}

void NCPA::AcousticEquationSet::changeTakeoffAngle( double newAngle ) {
	takeoff = newAngle;
}

NCPA::AcousticEquationSet::~AcousticEquationSet() { }

NCPA::AtmosphericSpecification *NCPA::AcousticEquationSet::getProfile() const { return profile; }