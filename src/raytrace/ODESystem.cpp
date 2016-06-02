#include <vector>
#include <stdexcept>
#include "ODESystem.h"

// Default constructor.  Sets up the object with a null equation set.
NCPA::ODESystem::ODESystem() {
	equations_ = 0;
}

// Basic constructor.  Sets up the pointer to the equation set
NCPA::ODESystem::ODESystem( NCPA::EquationSet *equationset ) {
	equations_ = equationset;
}

// Runs the 4th-order Runge-Kutta solver with no break conditions (i.e. each run is
// carried through to the full number of steps
int NCPA::ODESystem::rk4( double **solution, int steps, double *initialconditions, 
		double t0, double tend ) const {
	if (equations_ != 0)
		return this->rk4( equations_, solution, steps, initialconditions, t0, tend );
	else {
		std::runtime_error e( "No default equation set has been defined!" );
		throw e;
	}
}

// Runs the 4th-order Runge-Kutta solver with specified break conditions
int NCPA::ODESystem::rk4( double **solution, int steps, double *initialconditions, 
		double t0, double tend, 
		std::vector< NCPA::ODESystemBreakCondition * > conditions ) const {

	if (equations_ != 0)
		return this->rk4( equations_, solution, steps, initialconditions, t0, tend, conditions );
	else {
		std::runtime_error e( "No default equation set has been defined!" );
                throw e;
        }

}

int NCPA::ODESystem::rk4( NCPA::EquationSet *equations, double **solution, int steps,
		double *initialconditions, double t0, double tend ) const {
	std::vector< NCPA::ODESystemBreakCondition * > conditions;
	conditions.clear();
	return this->rk4( equations, solution, steps, initialconditions, t0, tend, conditions );
}



int NCPA::ODESystem::rk4( NCPA::EquationSet *equations, double **solution, int steps, 
		double *initialconditions, double t0, double tend, 
		std::vector< NCPA::ODESystemBreakCondition * > conditions ) const {

	// Set up the initial conditions in the solution matrix
	for (int i = 0; i < equations->numberOfEquations(); i++) {
		solution[0][i] = initialconditions[ i ];
	}

	// Get ready for the iteration procedure
	// Calculate the step size
	double h = (tend - t0) / steps;

	// Preallocate the iterator for the break conditions
	std::vector< NCPA::ODESystemBreakCondition * >::iterator bci;

	// The iteration variable is actually significant.  Since it tracks how many 
	// steps were taken, it allows us to reuse the same solutions matrix over
	// and over again.  So, declare it now and return it later.
	int k = 0;

	// Independent variable (we'll call it t for convenience although
	// it may or may not be time in any given solution)
	double t;
	
	// Allocate arrays to hold the various steps along the way
	double *temp0 = new double[ equations->numberOfEquations() ];
	double *temp1 = new double[ equations->numberOfEquations() ];
	double *temp2 = new double[ equations->numberOfEquations() ];
	double *temp3 = new double[ equations->numberOfEquations() ];
	double *temp4 = new double[ equations->numberOfEquations() ];
	double *partial1 = new double[ equations->numberOfEquations() ];
	double *partial2 = new double[ equations->numberOfEquations() ];
	double *partial3 = new double[ equations->numberOfEquations() ];


	for (k = 0; k < steps; k++) {
		t = t0 + k*h;  // Set up independent variable

		// Get the starting values
		for (int i = 0; i < equations->numberOfEquations(); i++) {
			temp0[ i ] = solution[ k ][ i ];
		}

		// Step 1 and prep for step 2
		equations->results( t, temp0, temp1 );
		for (int i = 0; i < equations->numberOfEquations(); i++) {
			//temp1[ i ] = h * equations->result( t, temp0, i );
			temp1[ i ] *= h;
			partial1[ i ] = solution[ k ][ i ] + temp1[ i ] / 2;
		}

		// Step 2 and prep for step 3
		equations->results( t + h/2, partial1, temp2 );
		for (int i = 0; i < equations->numberOfEquations(); i++) {
			//temp2[ i ] = h * equations->result( t + h/2, partial1, i );
			temp2[ i ] *= h;
			partial2[ i ] = solution[ k ][ i ] + temp2[ i ] / 2;
		}

		// Step 3 and prep for step 4
		equations->results( t + h/2, partial2, temp3 );
		for (int i = 0; i < equations->numberOfEquations(); i++) {
			//temp3[ i ] = h * equations->result( t + h/2, partial2, i );
			temp3[ i ] *= h;
			partial3[ i ] = solution[ k ][ i ] + temp3[ i ];
		}

		// Step 4 and calculate total of this iteration
		equations->results( t + h, partial3, temp4 );
		for (int i = 0; i < equations->numberOfEquations(); i++) {
			//temp4[ i ] = h * equations->result( t + h, partial3, i );
			temp4[ i ] *= h;
			solution[ k+1 ][ i ] = solution[ k ][ i ] 
				+ temp1[ i ] / 6 + temp2[ i ] / 3
				+ temp3[ i ] / 3 + temp4[ i ] / 6;
		}

		// Check to see if break conditions have been satisfied
		for (bci = conditions.begin(); bci != conditions.end(); bci++) {
                        if ((*bci)->shouldBreak( solution, k+1 )) {
                                (*bci)->printMessage();
                                return (k+1);
                        }
                }

	}

	// Deallocate memory to avoid leaks
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;
	delete [] temp4;
	delete [] partial1;
	delete [] partial2;
	delete [] partial3;

	return k;
}












/*
int NCPA::ODESystem::rk4( double idep_beg, double idep_end, int steps,
                        double *init, double **soln ) const {
	std::vector< NCPA::ODESystemBreakCondition * > conditions;
	conditions.clear();
	return this->rk4( idep_beg, idep_end, steps, init, soln, conditions );
}

int NCPA::ODESystem::rk4( double idep_beg, double idep_end, int steps,
                        double *init, double **soln,
                        std::vector< NCPA::ODESystemBreakCondition * > conditions ) const {

        double h = (idep_end - idep_beg) / steps;
        double *step1 = new double[ neqs ];
        double *step2 = new double[ neqs ];
        double *step3 = new double[ neqs ];
        double *step4 = new double[ neqs ];
        double *temp1 = new double[ neqs ];
        double *temp2 = new double[ neqs ];
        double *temp3 = new double[ neqs ];

        for (int eq = 0; eq < neqs; eq++) {
                soln[ 0 ][ eq ] = init[ eq ];
        }
	
	int i;
        for ( i = 0; i < steps; i++) {
                double idep_n = idep_beg + i*h;
                for (int k = 0; k < neqs; k++) {
                        step1[ k ] = h * sfn[ k ]( this, idep_n, soln[ i ] );
                        temp1[ k ] = soln[ i ][ k ] + step1[ k ] / 2;
                }
                for (int k = 0; k < neqs; k++) {
                        step2[ k ] = h * sfn[ k ]( this, idep_n + h/2, temp1 );
                        temp2[ k ] = soln[ i ][ k ] + step2[ k ] / 2;
                }
                for (int k = 0; k < neqs; k++) {
                        step3[ k ] = h * sfn[ k ]( this, idep_n + h/2, temp2 );
                        temp3[ k ] = soln[ i ][ k ] + step3[ k ];
                }
                for (int k = 0; k < neqs; k++) {
                        step4[ k ] = h * sfn[ k ]( this, idep_n + h, temp3 );
                        soln[ i+1 ][ k ] = soln[ i ][ k ]
                                                + step1[k]/6
                                                + step2[k]/3
                                                + step3[k]/3
                                                + step4[k]/6;
                }
        }

	// Clean up dynamically-allocated memory
	delete [] step1;
	delete [] step2;
	delete [] step3;
	delete [] step4;
	delete [] temp1;
	delete [] temp2;
	delete [] temp3;

        // check break conditions
        for (std::vector< ODESystemBreakCondition * >::iterator j = conditions.begin(); 
			j != conditions.end(); j++) {
		if ((*j)->shouldBreak( soln[i+1][(*j)->index()] )) {
			return (i+1);
		}
	}

	return i;
}
*/
