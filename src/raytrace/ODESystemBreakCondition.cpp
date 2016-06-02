#include "ODESystemBreakCondition.h"
#include <iostream>

NCPA::ODESystemBreakCondition::~ODESystemBreakCondition() {}

void NCPA::ODESystemBreakCondition::printMessage() const {
	if (message.length() > 0) {
		std::cout << message << std::endl;
	}
}


