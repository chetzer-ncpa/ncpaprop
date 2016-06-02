#ifndef __ATMOSPHERICBREAKCONDITIONS_H__
#define __ATMOSPHERICBREAKCONDITIONS_H__

#include <iostream>
#include <cmath>
#include "ODESystemBreakCondition.h"
#include "AtmosphericSpecification.h"

namespace NCPA {

	class HitGroundCondition : public ODESystemBreakCondition {

		protected:
			AtmosphericSpecification *profile_;
			int zindex_, xindex_, yindex_;

		public:
			HitGroundCondition( AtmosphericSpecification *profile, int xindex, int yindex, int zindex, std::string messageout = "" );
			~HitGroundCondition();

			bool shouldBreak( double **solution, int currentIteration );

	};

	class MaximumRangeCondition : public ODESystemBreakCondition {

		protected:
			int zindex_, xindex_, yindex_;
			double maxRange_;

		public:
			MaximumRangeCondition( int xindex, int yindex, int zindex, double maxRange, std::string messageout = "" );
			~MaximumRangeCondition();

			bool shouldBreak( double **solution, int currentIteration );
	};

	class UpwardRefractionCondition : public ODESystemBreakCondition {
		protected:
			int zindex_;
			unsigned int maxTurns;
			unsigned int turns;

		public:
			UpwardRefractionCondition( int zindex, std::string messageout = "", unsigned int turns = 0 );
			~UpwardRefractionCondition();

			bool shouldBreak( double **solution, int currentIteration );
	};
}
			
#endif
