#include "AtmosphericBreakConditions.h"
#include "AtmosphericSpecification.h"
#include <stdexcept>

NCPA::HitGroundCondition::~HitGroundCondition() {  }

NCPA::HitGroundCondition::HitGroundCondition( NCPA::AtmosphericSpecification *profile, int xindex, int yindex, int zindex, std::string breakMessage ) {

	if (xindex < 0 && yindex < 0) {
		throw std::runtime_error( "Either an x or y index, or both, must be supplied!" );
	}

	profile_ = profile;
	zindex_ = zindex;
	xindex_ = xindex;
	yindex_ = yindex;
	message = breakMessage;
}

bool NCPA::HitGroundCondition::shouldBreak( double **solution, int currentIteration ) {

	if (xindex_ >= 0 && yindex_ >= 0) {
		return (solution[ currentIteration ][ zindex_ ] <= 
			profile_->z0( solution[ currentIteration ][ xindex_ ], 
				     solution[ currentIteration ][ yindex_ ] ) );
	} else if (xindex_ < 0 ) {
		return (solution[ currentIteration ][ zindex_ ] <= profile_->z0( 0, solution[ currentIteration ][ yindex_ ] ) );
	} else {
		return (solution[ currentIteration ][ zindex_ ] <= profile_->z0( solution[ currentIteration ][ xindex_ ], 0 ) );
	}
}


NCPA::MaximumRangeCondition::~MaximumRangeCondition() {}

NCPA::MaximumRangeCondition::MaximumRangeCondition( int xindex, int yindex, int zindex, double maxRange, std::string messageout ) {
	if (xindex < 0 && yindex & 0) {
                throw std::runtime_error( "Either an x or y index, or both, must be supplied!" );
        }

	zindex_ = zindex;
        xindex_ = xindex;
        yindex_ = yindex;
        message = messageout;
	maxRange_ = maxRange;
}

bool NCPA::MaximumRangeCondition::shouldBreak( double **solution, int currentIteration ) {
	double currentRange = 0.0;
	if (xindex_ >= 0 && yindex_ >= 0) {
		currentRange = std::sqrt( std::pow( solution[currentIteration][xindex_], 2.0 ) 
		+ std::pow( solution[currentIteration][yindex_], 2.0 ) );
	} else if (xindex_ < 0 ) {
		currentRange = solution[currentIteration][yindex_];
	} else {
		currentRange = solution[currentIteration][xindex_];
	}

	return currentRange >= maxRange_;
}

NCPA::UpwardRefractionCondition::~UpwardRefractionCondition() {}

NCPA::UpwardRefractionCondition::UpwardRefractionCondition( int zindex, std::string messageout, unsigned int maxt ) {

        zindex_ = zindex;
        message = messageout;
	turns = 0;
	maxTurns = maxt;
}

bool NCPA::UpwardRefractionCondition::shouldBreak( double **solution, int currentIteration ) {
	if (currentIteration < 2) return false;

	double diff1 = solution[currentIteration-1][zindex_] - solution[currentIteration-2][zindex_];
	if (diff1 > 0) return false;
	double diff2 = solution[currentIteration][zindex_] - solution[currentIteration-1][zindex_];
	if (diff2 >= 0) {
		if (maxTurns > 0 && ++turns == maxTurns) {
			turns = 0;
			return true;
		}
	}
	return false;
}
