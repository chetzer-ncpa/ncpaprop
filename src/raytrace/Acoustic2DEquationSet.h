#ifndef __ACOUSTIC2DEQUATIONSET_H__
#define __ACOUSTIC2DEQUATIONSET_H__

#include "AtmosphericSpecification.h"
#include "AcousticEquationSet.h"

namespace NCPA {
	class Acoustic2DEquationSet : public AcousticEquationSet {

		public:
			Acoustic2DEquationSet( AtmosphericSpecification *profile, 
				double takeoffAngle, double propagation_azimuth,
				bool rangeDependent = true );
			//double result( double t, double *current_values, int which );
			void results( double t, double *current_values, double *new_values );
			//void changeTakeoffAngle( double newAngle );
			//void changeAzimuth( double newAzimuth );
			int numberOfEquations() const;
			double calculateTravelTime( double **solution, int steps, 
				double step_size, double azimuth ) const;
			double calculateAmplitude( double **solution, int steps ) const;
			double transmissionLoss( double **solution, int steps, double &refDist ) const;
			void setupInitialConditions( double *initialConditions, double zmin, double c0 ) const;
			void setupInitialConditions( double *initialConditions, double zmin );

	};
}




#endif
