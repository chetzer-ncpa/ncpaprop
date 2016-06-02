#ifndef __EQUATIONSET_H__
#define __EQUATIONSET_H__

#include <exception>
#include <cstring>

namespace NCPA {
	class RequestOutOfBoundsException : public std::exception {

		public:
			RequestOutOfBoundsException( const char *m ) {
				strcpy( message, m );
			}

			virtual const char *what() const throw() {
				return message;
			}

		protected:
			char *message;
	};


	class EquationSet {

		public:
// Return the result of one of the coupled differential equations.  An instance
// of this class represents a set of coupled differential equations, 
//    d/dt(v_0) = ...
//    d/dt(v_1) = ...
// etc, where t is the independent variable.  This function accepts the
// independent variable and the current values of all of the dependent
// variables, and returns d/dt(v_n) where n is supplied as the third argument.
//
// @param indep The value of the independent variable
// @param current_values The current values of the dependent variables
// @param which The number of the equation for which we want the derivative
// @return The derivative of variable [which], as calculated from equation [which]
			//virtual double result( double indep, double *current_values, 
			//	int which ) = 0;
			virtual void results( double indep, double *current_values, double *new_values ) = 0;
			virtual int numberOfEquations() const = 0;
			

	};

}
#endif
