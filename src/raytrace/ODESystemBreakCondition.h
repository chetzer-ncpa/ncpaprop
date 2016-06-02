#ifndef __ODESYSTEMBREAKCONDITION_H__
#define __ODESYSTEMBREAKCONDITION_H__
#include <iostream>

namespace NCPA {

	 class ODESystemBreakCondition {

                protected:
			std::string message;

                public:
			virtual ~ODESystemBreakCondition();
                        virtual bool shouldBreak( double **solution, int currentIteration ) = 0;
			void printMessage() const;
        };
}


#endif  // #ifndef __ODESYSTEMBREAKCONDITION_H__
