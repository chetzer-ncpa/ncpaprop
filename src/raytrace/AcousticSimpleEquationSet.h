#ifndef __ACOUSTICSIMPLEEQUATIONSET_H__
#define __ACOUSTICSIMPLEEQUATIONSET_H__

#include "AtmosphericSpecification.h"
#include "EquationSet.h"

namespace NCPA {
	class AcousticSimpleEquationSet : public EquationSet {

		public:
			AcousticSimpleEquationSet( AtmosphericSpecification *profile, 
				double takeoffAngle, double propagation_azimuth );
			double result( double t, double *current_values, int which );
			void changeTakeoffAngle( double newAngle );
			void changeAzimuth( double newAzimuth );
			int numberOfEquations() const;

		protected:
			AtmosphericSpecification *profile;
			double azimuth, takeoff;
	};
}




#endif
