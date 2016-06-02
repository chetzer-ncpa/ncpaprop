#ifndef __GENERALBREAKCONDITIONS_H__
#define __GENERALBREAKCONDITIONS_H__
#include <iostream>
#include "ODESystemBreakCondition.h"

namespace NCPA {

	/**
	  * Break condition for when a part of the solution decreases below a threshold.
	  * @version 2.0.0
	  */
        class MinimumBreakCondition : public ODESystemBreakCondition {

		protected:
			double flagValue;
			int eqnInd;

                public:
                        MinimumBreakCondition( int index, double flag, std::string message );
			MinimumBreakCondition( int index, double flag );
                        bool shouldBreak( double **solution, int currentIteration );
        };

	/**
          * Break condition for when a part of the solution increases above a threshold.
          * @version 2.0.0
          */
        class MaximumBreakCondition : public ODESystemBreakCondition {

		protected:
			double flagValue;
			int eqnInd;

                public:
                        MaximumBreakCondition( int index, double flag );
			MaximumBreakCondition( int index, double flag, std::string message );
                        bool shouldBreak( double **solution, int currentIteration );
        };
}


#endif  // #ifndef __GENERALBREAKCONDITIONS_H__
