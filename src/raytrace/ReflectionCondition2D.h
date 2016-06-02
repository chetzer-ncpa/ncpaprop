/**
 * Handles the 2D reflection condition
 * @version 1.0.2
 * @author Claus Hetzer, claus@olemiss.edu
 * @date 2012-04-02
 * 
 * Changelog:
 * 1.0.2: Added reset() function, triggered_ flag, triggered() function
 */

#ifndef __REFLECTIONCONDITION2D_H__
#define __REFLECTIONCONDITION2D_H__
#include <iostream>
#include "ODESystemBreakCondition.h"
#include "AtmosphericSpecification.h"

namespace NCPA {
	

        class ReflectionCondition2D : public ODESystemBreakCondition {

		protected:
			double propAzimuth;
			double launchAngle;
			AtmosphericSpecification *spec;
			unsigned int bounces, maxbounces;
			bool triggered_;

                public:
                        ReflectionCondition2D( AtmosphericSpecification *atmosphere, double propagationAzimuth, unsigned int maxbounces = 0 );
                        bool shouldBreak( double **solution, int currentIteration );
			void setLaunchAzimuthDegrees( double degrees );
			void setLaunchAngleRadians( double radians );
			unsigned int countBounces() const;
			bool triggered() const;
			void reset();
        };

	
}


#endif  // #ifndef __REFLECTIONCONDITION2D_H__
