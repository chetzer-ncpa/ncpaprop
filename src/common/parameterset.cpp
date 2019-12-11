#include "parameterset.h"
#include <iostream>
#include <cstring>
#include <cstdlib>


/*
Code for ParameterSet class
*/
NCPA::ParameterSet::ParameterSet() {
	_params.clear();
	_unparsed.clear();
	_delims = ":";
	_strict = true;
}

NCPA::ParameterSet::~ParameterSet() {
	for ( std::vector< NCPA::GenericParameter * >::iterator it = _params.begin(); 
			it != _params.end(); ++it ) {
		delete *it;
	}
	_params.clear();
	_unparsed.clear();
}

void NCPA::ParameterSet::setDelimiters( std::string newdelim ) {
	_delims = newdelim;
}

void NCPA::ParameterSet::addParameter( GenericParameter *newParam ) {
	_params.push_back( newParam );
}

void NCPA::ParameterSet::setStrict( bool tf ) {
	_strict = tf;
}

bool NCPA::ParameterSet::validate() {
	bool allgood = true;
	
	for ( std::vector< NCPA::GenericParameter * >::iterator it = _params.begin(); 
			it != _params.end(); ++it ) {
		allgood = (it->validate()) && allgood;
	}
	
	return allgood;
}

unsigned int NCPA::ParameterSet::parseCommandLine( int argc, const char **argv ) {
	bool expectingArg = false;
	std::string lastarg;
	unsigned int nOptions = 0;
	
	for (unsigned int i = 0; i < argc; i++) {
		
		std::string currentarg = argv[ i ];
		if (_isLongOption( currentarg )) {
			i = _processLongOption( argc, argv, i );
			nOptions++;
		} else if (_isShortOption( currentarg )) {
			i = _processShortOption( argc, argv, i );
			nOptions++;
		} else {
			_unparsed.push_back( currentarg );
		}
	}
	return nOptions;
}

// returns true if the string is at least 3 characters and the first two are "--"
bool NCPA::ParameterSet::_isLongOption( std::string opt ) {
	if (opt.size() < 3) {
		return false;
	}
	if (opt.compare(0,2,"--") == 0) {
		return true;
	} else {
		return false;
	}
}

// returns true if the string is at least 2 characters, the first character is "-", and the
// second is NOT "-"
bool NCPA::ParameterSet::_isShortOption( std::string opt ) {
	if (opt.size() < 2) {
		return false;
	}
	if (opt.compare(0,1,"-") == 0 && ( ! opt.compare(1,1,"-") )) {
		return true;
	} else {
		return false;
	}
}

unsigned int NCPA::ParameterSet::_processLongOption( int argc, const char **argv, unsigned int i ) {
	
	std::string fullarg = argv[ i ];
	
	// first, strip off the leading hyphens
	std::string stripped = fullarg.substr( fullarg.find_first_not_of( "-" ) );
	
	// Two cases.  One argument with an '=' character, or two arguments with spaces
	std::size_t found = str.find_first_of("=");
	if ( found != std::string::npos ) {    // not found
		
	}
	
	
}