#ifndef __ACOUSTICEQUATIONSET_H__
#define __ACOUSTICEQUATIONSET_H__

#include "EquationSet.h"
#include "AtmosphericSpecification.h"


namespace NCPA {
	class AcousticEquationSet : public EquationSet {
		
		public:
			virtual ~AcousticEquationSet();
			virtual void changeTakeoffAngle( double newAngle );
			virtual void changeAzimuth( double newAzimuth );
			virtual double calculateTravelTime( double **solution, int steps,
				double step_size, double azimuth ) const = 0;
			virtual double calculateAmplitude( double **solution, int steps ) const = 0;
			virtual double transmissionLoss( double **solution, int steps, double &refDist ) const = 0;
			//virtual void setupInitialConditions( double *initialConditions, double zmin ) = 0;
			virtual AtmosphericSpecification *getProfile() const;
			
		protected:
			AtmosphericSpecification *profile;
			double azimuth, takeoff;
			bool rangeDependent_;

	};

}



#endif
