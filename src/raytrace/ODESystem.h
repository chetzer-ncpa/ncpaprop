#ifndef __ODESYSTEM_H__
#define __ODESYSTEM_H__

#include <iostream>
#include <vector>
#include "ODESystemBreakCondition.h"
#include "EquationSet.h"

namespace NCPA {

	/**
	A system of ordinary differential equations that is amenable to solution via
	the fourth-order Runge-Kutta method.
	@version 2.0
	@author Claus Hetzer
	@email claus@olemiss.edu
	@date May 11, 2011
	*/
	class ODESystem {

                public:

			// Base constructor
			ODESystem();
			ODESystem( EquationSet * );

			virtual int rk4( double **solution, int steps, double *initialconditions,
				double t0, double tend ) const;
			virtual int rk4( double **solution, int steps, double *initialconditions,
				double t0, double tend,
				std::vector< ODESystemBreakCondition * > conditions ) const;

			virtual int rk4( EquationSet *equations, double **solution, int steps,
				double *initialconditions, double t0, double tend ) const;
			virtual int rk4( EquationSet *equations, double **solution, int steps,
				double *initialconditions, double t0, double tend,
				std::vector< ODESystemBreakCondition * > conditions ) const;


		protected:
			EquationSet *equations_;

        };


}

#endif  // #define __ODESYSTEM_H__
