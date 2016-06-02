#ifndef __ACOUSTIC3DEQUATIONSET_H__
#define __ACOUSTIC3DEQUATIONSET_H__

#include "AtmosphericSpecification.h"
#include "AcousticEquationSet.h"

namespace NCPA {
	class Acoustic3DEquationSet : public AcousticEquationSet {

		public:
			Acoustic3DEquationSet( AtmosphericSpecification *spec, 
				double takeoffAngle, double propagation_azimuth,
				bool rangeDependent = true );
			//double result( double t, double *current_values, int which );
			void results( double t, double *current_values, double *new_values );
			int numberOfEquations() const;

			double calculateTravelTime( double **solution, int steps, 
				double step_size, double azimuth = 0 ) const;
			double calculateAmplitude( double **solution, int steps ) const;
			double transmissionLoss( double **solution, int k, double &refRange ) const;
			double Jacobian( double **solution, int steps ) const;
			double calculateArrivalAzimuth( double **solution, int steps, int averages ) const;
			void setupInitialConditions( double *initialConditions, double zmin );


		protected:
			/*
			The following functions are neceessary to calculate ray path and Jacobian source functions in 
			3 dimensions with source functions varying in more than 1 dimension (non-stratified):
			        1.  Magnitude of nu vector
			        2.  Components and magnitude of group velocity vector
			        3.  Derivative of sound speed, winds, and group velocity with respect to angle
			*/
			double mag_nu( double r1, double r2, double r3,	double nu1, double nu2, double nu3 ) const;
			double cg1( double r1, double r2, double r3, double nu1, double nu2, double nu3 ) const;
			double cg2( double r1, double r2, double r3, double nu1, double nu2, double nu3 ) const;
			double cg3( double r1, double r2, double r3, double nu1, double nu2, double nu3 ) const;
			double cg_mag( double r1, double r2, double r3, double nu1, double nu2, double nu3 ) const;
			double dc_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle ) const;
			double du_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle ) const;
			double dv_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle ) const;
			double dw_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle ) const;
			double dCg_angle( double r1, double r2, double r3, double nu1, double nu2, double nu3, double dr1_angle, 
					double dr2_angle, double dr3_angle, double dnu1_angle, double dnu2_angle, double dnu3_angle ) const;
					
			/*
			 * This set of functions supplies the required atmospheric characteristics, instead of calculating them from
			 * the AtmosphericSpecification.
			 */
			double mag_nu(  double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double c0, double c0_0, double u, double v, double w ) const;
			double cg1( 	double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double c0, double c0_0, double u, double v, double w ) const;
			double cg2( 	double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double c0, double c0_0, double u, double v, double w ) const;
			double cg3( 	double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double c0, double c0_0, double u, double v, double w ) const;
			double cg_mag(  double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double c0, double c0_0, double u, double v, double w ) const;
			double dc_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					 double dc0dx, double dc0dy, double dc0dz ) const;
			double du_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					 double dudx, double dudy, double dudz ) const;
			double dv_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					 double dvdx, double dvdy, double dvdz ) const;
			double dw_angle( double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					 double dwdx, double dwdy, double dwdz ) const;
			double dCg_angle( double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double dr1_angle, double dr2_angle, double dr3_angle, 
					double dnu1_angle, double dnu2_angle, double dnu3_angle,
					double c0, double c0_0, double u, double v, double w,
					double dc0dx, double dc0dy, double dc0dz,
					double dudx, double dudy, double dudz,
					double dvdx, double dvdy, double dvdz,
					double dwdx, double dwdy, double dwdz ) const;

	};
}




#endif
