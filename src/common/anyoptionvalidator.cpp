#include "anyoptionvalidator.h"
#include "anyoption.h"
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>

// Constructor.  Does nothing.
NCPA::AnyOptionValidator::AnyOptionValidator() {
	this->_criteria.clear();
	this->_failed.clear();
}

// Destructor.  Frees memory held by tests and clears out the criteria map
NCPA::AnyOptionValidator::~AnyOptionValidator() {
	if (! _criteria.empty()) {
		std::vector< NCPA::OptionValidationCriterion * >::iterator it;
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
	std::vector< NCPA::OptionValidationCriterion * >::iterator it;
	for ( it=_criteria.begin(); it != _criteria.end(); ++it ) {
		try {
			if ( ! (*it)->validate( opts ) ) {
				_failed.push_back( (*it) );
			}
		} catch (const std::invalid_argument& ia) {
			std::cerr << "Number formatting error: " << (*it)->optionName() << " = " 
				<< opts->getValue( (*it)->optionName().c_str() ) << " threw " 
				<< ia.what() << std::endl;
			_failed.push_back( (*it) );
	        }
	}
	
	return _failed.empty();
}

// Returns any failure messages for processing by the user
std::vector< NCPA::OptionValidationCriterion * > NCPA::AnyOptionValidator::getFailedChecks() const {
	return _failed;
}

// Prints any failure messages to the provided output stream
void NCPA::AnyOptionValidator::printFailedChecks( std::ostream *out ) const {
	if (_failed.empty()) {
		return;
	}
	
	std::vector< NCPA::OptionValidationCriterion * >::const_iterator it;
	for ( it = _failed.begin(); it != _failed.end(); ++it ) {
		(*out) << (*it)->failureMessage() << std::endl;
	}
}

// Prints criteria descriptions to the provided output stream
void NCPA::AnyOptionValidator::printDescriptions( std::ostream *out ) const {
	std::vector< NCPA::OptionValidationCriterion * >::const_iterator it;
	for ( it=_criteria.begin(); it != _criteria.end(); ++it ) {
		(*out) << (*it)->description() << std::endl;
	}
}


// addOption() functions for tests that do not have comparison criteria or for those whose
// criteria are implicit in the type of test (e.g. positive)
void NCPA::AnyOptionValidator::addOption( const std::string &option, OPTION_NOTYPE_TEST_TYPE option_type ) {
	NCPA::OptionValidationCriterion *crit;
	switch (option_type) {
		case OPTION_REQUIRED:
			crit = new NCPA::RequiredCriterion( option );
			_criteria.push_back( crit );
			break;
		default:
			throw std::invalid_argument( "Type error, this test is for required options regardless of type" );
	}
}


void NCPA::AnyOptionValidator::addOption( const std::string &option, OPTION_INTEGER_TEST_TYPE option_type ) {
	NCPA::OptionValidationCriterion *crit;
	switch (option_type) {
		case OPTION_INTEGER_POSITIVE:
			crit = new NCPA::IntegerGreaterThanCriterion( option, 0, false );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_NEGATIVE:
			crit = new NCPA::IntegerLessThanCriterion( option, 0, false );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_ZERO:
			crit = new NCPA::IntegerEqualsCriterion( option, 0 );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_NONZERO:
			crit = new NCPA::IntegerNotEqualsCriterion( option, 0 );
			_criteria.push_back( crit );
			break;
		default:
			throw std::invalid_argument( "Type error, this test is for integer tests with implicit conditions" );
	}
}

void NCPA::AnyOptionValidator::addOption( const std::string &option, OPTION_FLOAT_TEST_TYPE option_type ) {
	NCPA::OptionValidationCriterion *crit;
	switch (option_type) {
		case OPTION_FLOAT_POSITIVE:
			crit = new NCPA::FloatGreaterThanCriterion( option, 0.0, false );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_NEGATIVE:
			crit = new NCPA::FloatLessThanCriterion( option, 0.0, false );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_ZERO:
			crit = new NCPA::FloatEqualsCriterion( option, 0.0 );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_NONZERO:
			crit = new NCPA::FloatNotEqualsCriterion( option, 0.0 );
			_criteria.push_back( crit );
			break;
		default:
			throw std::invalid_argument( "Type error, this test is for floating-point tests with implicit conditions" );
	}
}







void NCPA::AnyOptionValidator::addOption( const std::string &option, OPTION_INTEGER_TEST_TYPE option_type, int boundary1 ) {
	NCPA::OptionValidationCriterion *crit;
	switch (option_type) {
		case OPTION_INTEGER_GREATER_THAN:
			crit = new NCPA::IntegerGreaterThanCriterion( option, boundary1, false );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_GREATER_THAN_EQUAL:
			crit = new NCPA::IntegerGreaterThanCriterion( option, boundary1, true );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_LESS_THAN:
			crit = new NCPA::IntegerLessThanCriterion( option, boundary1, false );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_LESS_THAN_EQUAL:
			crit = new NCPA::IntegerLessThanCriterion( option, boundary1, true );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_EQUAL:
			crit = new NCPA::IntegerEqualsCriterion( option, boundary1 );
			_criteria.push_back( crit );
			break;
		case OPTION_INTEGER_NOT_EQUAL:
			crit = new NCPA::IntegerNotEqualsCriterion( option, boundary1 );
			_criteria.push_back( crit );
			break;
		default:
			throw std::invalid_argument( "Type error, this test is for integer tests with one condition" );
	}
}

void NCPA::AnyOptionValidator::addOption( const std::string &option, OPTION_FLOAT_TEST_TYPE option_type, double boundary1 ) {
	NCPA::OptionValidationCriterion *crit;
	switch (option_type) {
		case OPTION_FLOAT_GREATER_THAN:
			crit = new NCPA::FloatGreaterThanCriterion( option, boundary1, false );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_GREATER_THAN_EQUAL:
			crit = new NCPA::FloatGreaterThanCriterion( option, boundary1, true );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_LESS_THAN:
			crit = new NCPA::FloatLessThanCriterion( option, boundary1, false );
			_criteria.push_back( crit );
			break;
		case OPTION_FLOAT_LESS_THAN_EQUAL:
			crit = new NCPA::FloatLessThanCriterion( option, boundary1, true );
			_criteria.push_back( crit );
			break;
		default:
			throw std::invalid_argument( "Type error, this test is for floating-point tests with one condition" );
	}
}

void NCPA::AnyOptionValidator::addOption( const std::string &option, OPTION_STRING_TEST_TYPE option_type, int boundary1 ) {
	NCPA::OptionValidationCriterion *crit;
	switch (option_type) {
		case OPTION_STRING_MINIMUM_LENGTH:
			crit = new NCPA::StringMinimumLengthCriterion( option, boundary1 );
			_criteria.push_back( crit );
			break;
		case OPTION_STRING_MAXIMUM_LENGTH:
			crit = new NCPA::StringMaximumLengthCriterion( option, boundary1 );
			_criteria.push_back( crit );
			break;
		default:
			throw std::invalid_argument( "Type error, this test is for string tests with one condition" );
	}
}


// Destructor for ABC, must not be pure virtual
NCPA::OptionValidationCriterion::~OptionValidationCriterion() { }

std::string NCPA::OptionValidationCriterion::optionName() const {
	return _optName;
}

/*
 * DEFINITIONS FOR INDIVIDUAL CRITERION CLASSES GO BELOW HERE.
 * Each gets a constructor, a validate() method, a description, and a failure message
 */
NCPA::RequiredCriterion::RequiredCriterion( const std::string optionName ) {
	_optName = optionName;
}
std::string NCPA::RequiredCriterion::description() const {
	return _optName + " is present.";
}
std::string NCPA::RequiredCriterion::failureMessage() const {
	return _optName + " is not present.";
}
bool NCPA::RequiredCriterion::validate( AnyOption *opt )  {
	// Just check to see if it's been provided
	if (opt->getValue( _optName.c_str()) != NULL || opt->getFlag( _optName.c_str() ) ) {
		return true;
	} else {
		return false;
	}
}





NCPA::IntegerGreaterThanCriterion::IntegerGreaterThanCriterion( const std::string optionName, int comparison, bool trueIfEquals = false ) {
	_optName = optionName;
	_trueIfEquals = trueIfEquals;
	_value = comparison;
	_testedValue.clear();
}
std::string NCPA::IntegerGreaterThanCriterion::description() const {
	return _optName + " is greater than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
std::string NCPA::IntegerGreaterThanCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be greater than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
bool NCPA::IntegerGreaterThanCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	if (_trueIfEquals) {
		return ( val >= _value );
	} else {
		return ( val > _value );
	}
}


NCPA::IntegerLessThanCriterion::IntegerLessThanCriterion( const std::string optionName, int comparison, bool trueIfEquals = false ) {
	_optName = optionName;
	_trueIfEquals = trueIfEquals;
	_value = comparison;
	_testedValue.clear();
	
}
std::string NCPA::IntegerLessThanCriterion::description() const {
	return _optName + " is less than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
std::string NCPA::IntegerLessThanCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be less than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
bool NCPA::IntegerLessThanCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	if (_trueIfEquals) {
		return ( val <= _value );
	} else {
		return ( val < _value );
	}
}



NCPA::IntegerEqualsCriterion::IntegerEqualsCriterion( const std::string optionName, int comparison ) {
	_optName = optionName;
	_value = comparison;
	_testedValue.clear();
}
std::string NCPA::IntegerEqualsCriterion::description() const {
	return _optName + " is equal to " + std::to_string(_value) + ".";
}
std::string NCPA::IntegerEqualsCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be equal to " + std::to_string(_value) + ".";
}
bool NCPA::IntegerEqualsCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val == _value );
}




NCPA::IntegerNotEqualsCriterion::IntegerNotEqualsCriterion( const std::string optionName, int comparison ) {
	_optName = optionName;
	_value = comparison;
	_testedValue.clear();
}
std::string NCPA::IntegerNotEqualsCriterion::description() const {
	return _optName + " is not equal to " + std::to_string(_value) + ".";
}
std::string NCPA::IntegerNotEqualsCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must not be equal to " + std::to_string(_value) + ".";
}
bool NCPA::IntegerNotEqualsCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	int val = std::stoi( valStr );
	return ( val != _value );
}







NCPA::FloatGreaterThanCriterion::FloatGreaterThanCriterion( const std::string optionName, double comparison, bool trueIfEquals = false ) {
	_optName = optionName;
	_trueIfEquals = trueIfEquals;
	_value = comparison;
	_testedValue.clear();
}
std::string NCPA::FloatGreaterThanCriterion::description() const {
	return _optName + " is greater than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
std::string NCPA::FloatGreaterThanCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be greater than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
bool NCPA::FloatGreaterThanCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	if (_trueIfEquals) {
		return ( val >= _value );
	} else {
		return ( val > _value );
	}
}


NCPA::FloatLessThanCriterion::FloatLessThanCriterion( const std::string optionName, double comparison, bool trueIfEquals = false ) {
	_optName = optionName;
	_trueIfEquals = trueIfEquals;
	_value = comparison;
	_testedValue.clear();
}
std::string NCPA::FloatLessThanCriterion::description() const {
	return _optName + " is less than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
std::string NCPA::FloatLessThanCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be less than " 
		+ (_trueIfEquals ? "or equal to " : "" ) + std::to_string(_value) + ".";
}
bool NCPA::FloatLessThanCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	if (_trueIfEquals) {
		return ( val <= _value );
	} else {
		return ( val < _value );
	}
}



NCPA::FloatEqualsCriterion::FloatEqualsCriterion( const std::string optionName, double comparison ) {
	_optName = optionName;
	_value = comparison;
	_testedValue.clear();
}
std::string NCPA::FloatEqualsCriterion::description() const {
	return _optName + " is equal to " + std::to_string(_value) + ".";
}
std::string NCPA::FloatEqualsCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must be equal to " + std::to_string(_value) + ".";
}
bool NCPA::FloatEqualsCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val == _value );
}




NCPA::FloatNotEqualsCriterion::FloatNotEqualsCriterion( const std::string optionName, double comparison ) {
	_optName = optionName;
	_value = comparison;
	_testedValue.clear();
}
std::string NCPA::FloatNotEqualsCriterion::description() const {
	return _optName + " is not equal to " + std::to_string(_value) + ".";
}
std::string NCPA::FloatNotEqualsCriterion::failureMessage() const {
	return _optName + " (" + _testedValue + ") must not be equal to " + std::to_string(_value) + ".";
}
bool NCPA::FloatNotEqualsCriterion::validate( AnyOption *opt ) {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		return true;
	}
	
	_testedValue = valStr;
	double val = std::stof( valStr );
	return ( val != _value );
}



NCPA::StringMinimumLengthCriterion::StringMinimumLengthCriterion( const std::string optionName, int minLength ) {
	_optName = optionName;
	_value = minLength;
	_testedValue.clear();
}
std::string NCPA::StringMinimumLengthCriterion::description() const {
	return _optName + " is at least " + std::to_string(_value) + " characters.";
}
std::string NCPA::StringMinimumLengthCriterion::failureMessage() const {
	return _optName + " (\"" + _testedValue + "\") must be at least " + std::to_string(_value) + " characters.";
}
bool NCPA::StringMinimumLengthCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		//std::cout << "Option " << _optName << " not specified, passing" << std::endl;
		return true;
	}
	
	_testedValue = valStr;
	//std::cout << "Length of '" << _testedValue << "' is " << _testedValue.length() << std::endl;
	return ( _testedValue.length() >= _value );
}

NCPA::StringMaximumLengthCriterion::StringMaximumLengthCriterion( const std::string optionName, int maxLength ) {
	_optName = optionName;
	_value = maxLength;
	_testedValue.clear();
}
std::string NCPA::StringMaximumLengthCriterion::description() const {
	return _optName + " is at most " + std::to_string(_value) + " characters.";
}
std::string NCPA::StringMaximumLengthCriterion::failureMessage() const {
	return _optName + " (\"" + _testedValue + "\") must be at most " + std::to_string(_value) + " characters.";
}
bool NCPA::StringMaximumLengthCriterion::validate( AnyOption *opt )  {
	
	char *valStr = opt->getValue( _optName.c_str() );
	if (valStr == NULL) {
		_testedValue.clear();
		//std::cout << "Option " << _optName << " not specified, passing" << std::endl;
		return true;
	}
	
	_testedValue = valStr;
	//std::cout << "Length of '" << _testedValue << "' is " << _testedValue.length() << std::endl;
	return ( _testedValue.length() <= _value );
}


