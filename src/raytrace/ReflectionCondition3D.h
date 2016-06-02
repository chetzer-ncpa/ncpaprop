	
/**
  * Handles the 3D reflection condition
  * @author Claus Hetzer, claus@olemiss.edu
  * @date 2012-04-02
  * @version 1.0.0
  */

#ifndef __REFLECTIONCONDITION3D_H__
#define __REFLECTIONCONDITION3D_H__
#include <iostream>
#include "ODESystemBreakCondition.h"
#include "AtmosphericSpecification.h"

namespace NCPA {

        class ReflectionCondition3D : public ODESystemBreakCondition {

		protected:
			double propAzimuth;
			AtmosphericSpecification *spec;
			unsigned int bounces, maxbounces;
			bool triggered_;

                public:
                        ReflectionCondition3D( AtmosphericSpecification *atmosphere, double propagationAzimuth, unsigned int maxbounces = 0 );
                        bool shouldBreak( double **solution, int currentIteration );
			void setAzimuth( double az );
			bool triggered() const;
			void reset();
        };

	
}


#endif  // #ifndef __REFLECTIONCONDITION3D_H__
