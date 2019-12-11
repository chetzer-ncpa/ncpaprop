#include <map>
#include "parameterset.h"



// constructor
NCPA::ParameterSet::ParameterSet() {
	type_map_.clear();
	string_map_.clear();
	int_map_.clear();
	double_map_.clear();
}

// destructor
NCPA::ParameterSet::~ParameterSet() { }

unsigned long NCPA::ParameterSet::removeKey_( std::string key ) {
	unsigned long total = type_map_.erase( key );
	total += string_map_.erase( key );
	total += int_map_.erase( key );
	total += double_map_.erase( key );
	return total;
}

void NCPA::ParameterSet::setStringValue( std::string key, std::string value ) {
	
	// see if there is already a value set
	type_map_.erase( key );
	string_map_.erase( key );
	
	
	
}