#include "GeneralBreakConditions.h"
#include <iostream>


NCPA::MinimumBreakCondition::MinimumBreakCondition( int index, double flag ) {
	eqnInd = index;
	flagValue = flag;
	message = "";
}

NCPA::MinimumBreakCondition::MinimumBreakCondition( int index, double flag, std::string m ) {
        eqnInd = index;
        flagValue = flag;
	message = m;
}



bool NCPA::MinimumBreakCondition::shouldBreak( double **solution, int currentIteration ) {
	double testValue = solution[ currentIteration ][ eqnInd ];
	if (testValue < flagValue) { return true; } else { return false; }
}

NCPA::MaximumBreakCondition::MaximumBreakCondition( int index, double flag ) {
        eqnInd = index;
        flagValue = flag;
	message = "";
}

NCPA::MaximumBreakCondition::MaximumBreakCondition( int index, double flag, std::string m ) {
        eqnInd = index;
        flagValue = flag;
	message = m;
}


bool NCPA::MaximumBreakCondition::shouldBreak( double **solution, int currentIteration ) {
	double testValue = solution[ currentIteration ][ eqnInd ];
	if (testValue > flagValue) { return true; } else { return false; }
}


