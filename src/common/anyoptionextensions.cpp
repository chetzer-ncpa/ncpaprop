#include "anyoptionextensions.h"
#include "anyoption.h"
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstring>
#include <algorithm>


/**********************************************************************
NCPA::AnyOptionValidator class methods here
**********************************************************************/


// Constructor.  Does nothing.
NCPA::AnyOptionValidator::AnyOptionValidator() {
	this->_criteria.clear();
	this->_failed.clear();
}

// Destructor.  Frees memory held by tests and clears out the criteria vector
NCPA::AnyOptionValidator::~AnyOptionValidator() {
	if (! _criteria.empty()) {
		std::vector< NCPA::OptionTest * >::iterator it;
		for ( it=_criteria.begin(); it != _criteria.end(); ++it ) {
			delete (*it);
		}
		_criteria.clear();
		_failed.clear();
	}
}

// Validates options, sets up the vector of failed options, and returns false if any options failed,
// true otherwise.  Empty criteria sets return true by default.
bool NCPA::AnyOptionValidator::validateOptions( AnyOption *opts ) {
	
	// Return true if there are no criteria to satisfy
	if ( _criteria.empty() ) {
		return true;
	}
	
	// Get ready
	_failed.clear();
	std::vector< NCPA::OptionTest * >::iterator it;
	for ( it=_criteria.begin(); it != _criteria.end(); ++it ) {
		std::cout << "Checking option: " << (*it)->description() << std::endl;
		try {
			if ( ! (*it)->validate( opts ) ) {
				_failed.push_back( (*it) );
				std::cout << "Failed!" << std::endl;
			} else {
				std::cout << "Passed!" << std::endl;
			}
		} catch (const std::invalid_argument& ia) {
			std::cerr << "Number formatting error: " << (*it)->optionName() << " = " 
				<< opts->getValue( (*it)->optionName().c_str() ) << " threw " 
				<< ia.what() << std::endl;
			_failed.push_back( (*it) );
	        } catch (const std::logic_error* le) {
	        	std::cerr << "Incomplete test setup: " << (*it)->description()
				<< " - parameters not fully defined." << std::endl;
			_failed.push_back( (*it) );
	        }
	}
	
	return _failed.empty();
}

// Returns any failure messages for processing by the user
std::vector< NCPA::OptionTest * > NCPA::AnyOptionValidator::getFailedTests() const {
	return _failed;
}

// Prints any failure messages to the provided output stream
void NCPA::AnyOptionValidator::printFailedTests( std::ostream *out ) const {
	if (_failed.empty()) {
		return;
	}
	
	std::vector< NCPA::OptionTest * >::const_iterator it;
	for ( it = _failed.begin(); it != _failed.end(); ++it ) {
		(*out) << (*it)->failureMessage() << std::endl;
	}
}

// Prints criteria descriptions to the provided output stream
void NCPA::AnyOptionValidator::printDescriptions( std::ostream *out ) const {
	std::vector< NCPA::OptionTest * >::const_iterator it;
	for ( it=_criteria.begin(); it != _criteria.end(); ++it ) {
		(*out) << (*it)->description() << std::endl;
	}
}


// addOption(): adds a test to the queue and returns a pointer to the test
NCPA::OptionTest * NCPA::AnyOptionValidator::addTest( const std::string option,
		OPTION_TEST_TYPE option_type ) {
			
	NCPA::OptionTest *crit;
	switch (option_type) {
		case OPTION_REQUIRED:
			crit = new NCPA::RequiredTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_RADIO_BUTTON:
			crit = new NCPA::RadioButtonTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_STRING_SET:
			crit = new NCPA::StringSetTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_POSITIVE:
			crit = new NCPA::IntegerGreaterThanTest( option );
			crit->addIntegerParameter( 0 );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_NEGATIVE:
			crit = new NCPA::IntegerLessThanTest( option );
			crit->addIntegerParameter( 0 );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_ZERO:
			crit = new NCPA::IntegerEqualToTest( option );
			crit->addIntegerParameter( 0 );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_NONZERO:
			crit = new NCPA::IntegerNotEqualToTest( option );
			crit->addIntegerParameter( 0 );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_POSITIVE:
			crit = new NCPA::FloatGreaterThanTest( option );
			crit->addFloatParameter( 0.0 );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_NEGATIVE:
			crit = new NCPA::FloatLessThanTest( option );
			crit->addFloatParameter( 0.0 );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_ZERO:
			crit = new NCPA::FloatEqualToTest( option );
			crit->addFloatParameter( 0.0 );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_NONZERO:
			crit = new NCPA::FloatNotEqualToTest( option );
			crit->addFloatParameter( 0.0 );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_GREATER_THAN:
			crit = new NCPA::IntegerGreaterThanTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_GREATER_THAN_OR_EQUAL:
			crit = new NCPA::IntegerGreaterThanOrEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_LESS_THAN:
			crit = new NCPA::IntegerLessThanTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_LESS_THAN_OR_EQUAL:
			crit = new NCPA::IntegerLessThanOrEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_EQUAL:
			crit = new NCPA::IntegerEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_NOT_EQUAL:
			crit = new NCPA::IntegerNotEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_GREATER_THAN:
			crit = new NCPA::FloatGreaterThanTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_GREATER_THAN_OR_EQUAL:
			crit = new NCPA::FloatGreaterThanOrEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_LESS_THAN:
			crit = new NCPA::FloatLessThanTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_LESS_THAN_OR_EQUAL:
			crit = new NCPA::FloatLessThanOrEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_EQUAL:
			crit = new NCPA::FloatEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_NOT_EQUAL:
			crit = new NCPA::FloatNotEqualToTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_STRING_MINIMUM_LENGTH:
			crit = new NCPA::StringMinimumLengthTest( option );
			_criteria.push_back( crit );
			break;
		case OPTION_STRING_MAXIMUM_LENGTH:
			crit = new NCPA::StringMaximumLengthTest( option );
			_criteria.push_back( crit );
			break;
		default:
			throw std::invalid_argument( "Undefined test requested" );
	}
	return crit;
}



/**********************************************************************
NCPA::OptionTest class methods here
**********************************************************************/

// Destructor for ABC, must not be pure virtual
NCPA::OptionTest::~OptionTest() { }

std::string NCPA::OptionTest::optionName() const {
	return _optName;
}

void NCPA::OptionTest::addIntegerParameter( int param ) { }
void NCPA::OptionTest::addFloatParameter( double param ) { }
void NCPA::OptionTest::addStringParameter( std::string param ) { }
bool NCPA::OptionTest::ready() const { return _ready; }



/**********************************************************************
NCPA::OptionTest derived class methods here
Each gets a constructor, a validate() method, a description, a failure message,
and optional add<Type>Parameter() method(s)
**********************************************************************/

NCPA::RequiredTest::RequiredTest( const std::string optionName ) {
	_optName = optionName;
	_ready = true;
}
std::string NCPA::RequiredTest::description() const {
	return _optName + " is present.";
}
std::string NCPA::RequiredTest::failureMessage() const {
	return _optName + " is not present.";
}
bool NCPA::RequiredTest::validate( AnyOption *opt )  {
	// Just check to see if it's been provided
	if (opt->getValue( _optName.c_str()) != NULL || opt->getFlag( _optName.c_str() ) ) {
		return true;
	} else {
		return false;
	}
}
std::string NCPA::RequiredTest::valueString() const { return ""; }






NCPA::RadioButtonTest::RadioButtonTest( const std::string optionName ) {
	_buttons.clear();
	_optName = optionName;
	_matched.clear();
}
std::string NCPA::RadioButtonTest::description() const {
	return _optName + ": One and only one of " + this->valueString() + " must be present.";
}
std::string NCPA::RadioButtonTest::failureMessage() const {
	ostringstream oss;
	oss << _optName << ": " << _matched.size() << " of " << this->valueString()
		<< " are present; must be one and only one.";
	return oss.str();
}
std::string NCPA::RadioButtonTest::valueString() const {
	ostringstream oss;
	oss << "{ ";
	for (std::vector<std::string>::const_iterator it = _buttons.begin();
			it != _buttons.end(); ++it) {
		if (it != _buttons.begin()) {
			oss << ", ";
		}
		oss << *it;
	}
	oss << " }";
	return oss.str();
}
bool NCPA::RadioButtonTest::validate( AnyOption *opt )  {
	_matched.clear();
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	for (std::vector<std::string>::const_iterator it = _buttons.begin();
			it != _buttons.end(); ++it) {
		if (opt->getValue( it->c_str()) != NULL || opt->getFlag( it->c_str() ) ) {
			_matched.push_back( *it );
		}
	}
	
	return (_matched.size() == 1);
}
void NCPA::RadioButtonTest::addStringParameter( const std::string newButton ) {
	std::string str = newButton;
	_buttons.push_back( str );
}
std::vector< std::string > NCPA::RadioButtonTest::lastMatched() const {
	std::vector< std::string > v( _matched );
	return v;
}
bool NCPA::RadioButtonTest::ready() const {
	return !( _buttons.empty() );
}





NCPA::IntegerGreaterThanTest::IntegerGreaterThanTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::IntegerGreaterThanTest::description() const {
	return _optName + " is greater than " + 
		this->valueString() + ".";
}
std::string NCPA::IntegerGreaterThanTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be greater than " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerGreaterThanTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val > _value );
}
std::string NCPA::IntegerGreaterThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerGreaterThanTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerGreaterThanOrEqualToTest::IntegerGreaterThanOrEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::IntegerGreaterThanOrEqualToTest::description() const {
	return _optName + " is greater than or equal to " + 
		this->valueString() + ".";
}
std::string NCPA::IntegerGreaterThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be greater than or equal to " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerGreaterThanOrEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val >= _value );
}
std::string NCPA::IntegerGreaterThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerGreaterThanOrEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerLessThanTest::IntegerLessThanTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::IntegerLessThanTest::description() const {
	return _optName + " is less than " + 
		this->valueString() + ".";
}
std::string NCPA::IntegerLessThanTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be less than " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerLessThanTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val < _value );
}
std::string NCPA::IntegerLessThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerLessThanTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerLessThanOrEqualToTest::IntegerLessThanOrEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::IntegerLessThanOrEqualToTest::description() const {
	return _optName + " is less than or equal to " + 
		this->valueString() + ".";
}
std::string NCPA::IntegerLessThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be less than or equal to " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerLessThanOrEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val <= _value );
}
std::string NCPA::IntegerLessThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerLessThanOrEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}





NCPA::IntegerEqualToTest::IntegerEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::IntegerEqualToTest::description() const {
	return _optName + " is equal to " + this->valueString() + ".";
}
std::string NCPA::IntegerEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be equal to " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val == _value );
}
std::string NCPA::IntegerEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}



NCPA::IntegerNotEqualToTest::IntegerNotEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::IntegerNotEqualToTest::description() const {
	return _optName + " is not equal to " + this->valueString() + ".";
}
std::string NCPA::IntegerNotEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must not be equal to " 
		+ this->valueString() + ".";
}
bool NCPA::IntegerNotEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val != _value );
}
std::string NCPA::IntegerNotEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::IntegerNotEqualToTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}






NCPA::FloatGreaterThanTest::FloatGreaterThanTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0.0;
	_ready = false;
}
std::string NCPA::FloatGreaterThanTest::description() const {
	return _optName + " is greater than " 
		 + this->valueString() + ".";
}
std::string NCPA::FloatGreaterThanTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be greater than " 
		+ this->valueString() + ".";
}
bool NCPA::FloatGreaterThanTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val > _value );
}
std::string NCPA::FloatGreaterThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatGreaterThanTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}



NCPA::FloatGreaterThanOrEqualToTest::FloatGreaterThanOrEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0.0;
	_ready = false;
}
std::string NCPA::FloatGreaterThanOrEqualToTest::description() const {
	return _optName + " is greater than or equal to " 
		 + this->valueString() + ".";
}
std::string NCPA::FloatGreaterThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be greater than or equal to " 
		+ this->valueString() + ".";
}
bool NCPA::FloatGreaterThanOrEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val >= _value );
}
std::string NCPA::FloatGreaterThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatGreaterThanOrEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}





NCPA::FloatLessThanTest::FloatLessThanTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0.0;
	_ready = false;
}
std::string NCPA::FloatLessThanTest::description() const {
	return _optName + " is less than " 
		 + this->valueString() + ".";
}
std::string NCPA::FloatLessThanTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be less than " 
		+ this->valueString() + ".";
}
bool NCPA::FloatLessThanTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val < _value );
}
std::string NCPA::FloatLessThanTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatLessThanTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}



NCPA::FloatLessThanOrEqualToTest::FloatLessThanOrEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0.0;
	_ready = false;
}
std::string NCPA::FloatLessThanOrEqualToTest::description() const {
	return _optName + " is less than " 
		 + this->valueString() + ".";
}
std::string NCPA::FloatLessThanOrEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be less than " 
		+ this->valueString() + ".";
}
bool NCPA::FloatLessThanOrEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val <= _value );
}
std::string NCPA::FloatLessThanOrEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatLessThanOrEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}




NCPA::FloatEqualToTest::FloatEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0.0;
	_ready = false;
}
std::string NCPA::FloatEqualToTest::description() const {
	return _optName + " is equal to " + this->valueString() + ".";
}
std::string NCPA::FloatEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be equal to " 
		+ this->valueString() + ".";
}
bool NCPA::FloatEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val == _value );
}
std::string NCPA::FloatEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}



NCPA::FloatNotEqualToTest::FloatNotEqualToTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0.0;
	_ready = false;
}
std::string NCPA::FloatNotEqualToTest::description() const {
	return _optName + " is not equal to " + this->valueString() + ".";
}
std::string NCPA::FloatNotEqualToTest::failureMessage() const {
	return _optName + " (" + _testedValue + ") must not be equal to " 
		+ this->valueString() + ".";
}
bool NCPA::FloatNotEqualToTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val != _value );
}
std::string NCPA::FloatNotEqualToTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::FloatNotEqualToTest::addFloatParameter( double param ) {
	_value = param;
	_ready = true;
}






NCPA::StringMinimumLengthTest::StringMinimumLengthTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::StringMinimumLengthTest::description() const {
	return _optName + " is at least " + this->valueString() + " characters.";
}
std::string NCPA::StringMinimumLengthTest::failureMessage() const {
	return _optName + " (\"" + _testedValue + "\") must be at least " 
		+ this->valueString() + " characters long.";
}
bool NCPA::StringMinimumLengthTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;

	return ( _testedValue.length() >= _value );
}
std::string NCPA::StringMinimumLengthTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::StringMinimumLengthTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}




NCPA::StringMaximumLengthTest::StringMaximumLengthTest( const std::string optionName ) {
	_optName = optionName;
	_testedValue.clear();
	_value = 0;
	_ready = false;
}
std::string NCPA::StringMaximumLengthTest::description() const {
	return _optName + " is at most " + this->valueString() + " characters.";
}
std::string NCPA::StringMaximumLengthTest::failureMessage() const {
	return _optName + " (\"" + _testedValue + "\") must be at most " 
		+ this->valueString() + " characters long.";
}
bool NCPA::StringMaximumLengthTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;

	return ( _testedValue.length() <= _value );
}
std::string NCPA::StringMaximumLengthTest::valueString() const {
	return this->ready() ? std::to_string(_value) : "";
}
void NCPA::StringMaximumLengthTest::addIntegerParameter( int param ) {
	_value = param;
	_ready = true;
}




NCPA::StringSetTest::StringSetTest( const std::string optionName ) {
	_choices.clear();
	_optName = optionName;
}
std::string NCPA::StringSetTest::description() const {
	return _optName + " must be in " + this->valueString() + ".";
}
std::string NCPA::StringSetTest::failureMessage() const {
	return _optName + ": " + '"' + _testedValue + '"' + " is not in " + this->valueString() + ".";
}
std::string NCPA::StringSetTest::valueString() const {
	ostringstream oss;
	oss << "{ ";
	for (std::vector<std::string>::const_iterator it = _choices.begin();
			it != _choices.end(); ++it) {
		if (it != _choices.begin()) {
			oss << ", ";
		}
		oss << '"' << *it << '"';
	}
	oss << " }";
	return oss.str();
}
bool NCPA::StringSetTest::validate( AnyOption *opt )  {
	if (! this->ready() ) {
		throw new std::logic_error( _optName + ": no options defined." );
	}
	
	_testedValue = opt->getValue( _optName.c_str() );
	std::vector< std::string >::const_iterator it = std::find( _choices.begin(), _choices.end(), _testedValue );
	return ( it != _choices.end() );
}
void NCPA::StringSetTest::addStringParameter( std::string newChoice ) {
	_choices.push_back( newChoice );
}
bool NCPA::StringSetTest::ready() const {
	return ! ( _choices.empty() );
}