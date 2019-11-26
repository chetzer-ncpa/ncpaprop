/**
  * Extensions of and tools used with the AnyOption class presented in anyoption.h
  */


/*
Example:

	// Option parser and validator
	using namespace std;
	using namespace NCPA;
	AnyOption *opt = new AnyOption();
	AnyOptionValidator *validator = new AnyOptionValidator();

	// Required and must be positive
	opt->setOption( "positive_int" );  	
	validator->addTest( "positive_int", NCPA::OPTION_REQUIRED );
	validator->addTest( "positive_int", NCPA::OPTION_INTEGER_POSITIVE );

	// one and only one of these must exist
	opt->setFlag( "flagOne" );		
	opt->setFlag( "flagTwo" );		
	opt->setFlag( "flagThree" );
	OptionTest *radioTest = validator->addTest( "groupOne", NCPA::OPTION_RADIO_BUTTON );
	radioTest->addStringParameter( "flagOne" );
	radioTest->addStringParameter( "flagTwo" );
	radioTest->addStringParameter( "flagThree" );

	// Optional, but must be one of a set of values
	opt->setOption( "str_set" );
	OptionTest *setTest = validator->addTest( "str_set", NCPA::OPTION_STRING_SET );
	setCrit->addStringParameter( "this" );
	setCrit->addStringParameter( "is" );
	setCrit->addStringParameter( "a" );
	setCrit->addStringParameter( "test" );

	// get options from file and validate
	opt->processFile( "./tester.options" );
	if (validator->validateOptions( opt )) {
		cout << endl << "All options OK" << endl;
	} else {
		cout << endl << "Error.  Failed options:" << endl;
		validator->printFailedTests( &cout );
	}

*/

#ifndef __NCPA_ANYOPTION_EXTENSIONS__
#define __NCPA_ANYOPTION_EXTENSIONS__

#include "anyoption.h"
#include <string>
#include <vector>

namespace NCPA {
	
	/**
	  * These indicate tests that do not depend on the type of the tested
	  * value.
	  */
	enum OPTION_TEST_TYPE : unsigned int {
		
		/** This option/flag must be present. */
		OPTION_REQUIRED,
		
		/** 
		  * Designates a group of options, one and only one of which
		  * must be present.
		  */
		OPTION_RADIO_BUTTON,
		
		/** Integer that must be > 0 */
		OPTION_INTEGER_POSITIVE,
		
		/** Integer that must be < 0 */
		OPTION_INTEGER_NEGATIVE,
		
		/** Integer that must be > another integer */
		OPTION_INTEGER_GREATER_THAN,
		
		/** Integer that must be >= another integer */
		OPTION_INTEGER_GREATER_THAN_OR_EQUAL,
		
		/** Integer that must be < another integer */
		OPTION_INTEGER_LESS_THAN,
		
		/** Integer that must be <= another integer */
		OPTION_INTEGER_LESS_THAN_OR_EQUAL,
		
		/** Integer that must == 0 */
		OPTION_INTEGER_ZERO,
		
		/** Integer that must != 0 */
		OPTION_INTEGER_NONZERO,
		
		/** Integer that must == another integer */
		OPTION_INTEGER_EQUAL,
		
		/** Integer that must != another integer */
		OPTION_INTEGER_NOT_EQUAL,
		
		/** Double that must > 0.0 */
		OPTION_FLOAT_POSITIVE,
		
		/** Double that must < 0.0 */
		OPTION_FLOAT_NEGATIVE,
		
		/** Double that must > another double */
		OPTION_FLOAT_GREATER_THAN,
		
		/** Double that must >= another double */
		OPTION_FLOAT_GREATER_THAN_OR_EQUAL,
		
		/** Double that must < another double */
		OPTION_FLOAT_LESS_THAN,
		
		/** Double that must <= another double */
		OPTION_FLOAT_LESS_THAN_OR_EQUAL,
		
		/** Double that must == another double (standard floating point caveats apply) */
		OPTION_FLOAT_EQUAL,
		
		/** Double that must != another double (standard floating point caveats apply) */
		OPTION_FLOAT_NOT_EQUAL,
		
		/** Double that must == 0.0 (standard floating point caveats apply) */
		OPTION_FLOAT_ZERO,
		
		/** Double that must != 0.0 (standard floating point caveats apply) */
		OPTION_FLOAT_NONZERO,		
		
		/** string .size() must be >= an integer */
		OPTION_STRING_MINIMUM_LENGTH,
		
		/** string .size() must be <= an integer */
		OPTION_STRING_MAXIMUM_LENGTH,
		
		/** string must match one of a set of strings */
		OPTION_STRING_SET
	};

	/**
	  * Abstract base class for validation criteria.  See subclass descriptions
	  * for general usage, you won't instantiate this class directly.
	  */
	class OptionTest {
	public:
		/**
		  * Destructor.  Cleans up any dynamically allocated memory.
		  */
		virtual ~OptionTest();
		
		/**
		  * Run the validation checks specified for this option
		  * @param opts A pointer to the AnyOption object that has ingested
		  *             the command line and file options
		  * @return true if the test passes, false otherwise
		  */
		virtual bool validate( AnyOption *opts ) = 0;
		
		/**
		  * A text description of the test that is to be run.
		  * @return A std::string containing a description of the test
		  */
		virtual std::string description() const = 0;
		
		/**
		  * A text description of why the test failed, if applicable.
		  * @return A std::string containing a description of the failure
		  */
		virtual std::string failureMessage() const = 0;
		
		/**
		  * The name of the parameter to be checked
		  * @return A std::string containing the parameter name
		  */
		virtual std::string optionName() const;
		
		/**
		  * Add an integer parameter to the test, if applicable.
		  *
		  * Does nothing unless overridden.
		  * @param param	The parameter to add to the test
		  */
		virtual void addIntegerParameter( int param );
		
		/**
		  * Add an integer parameter to the test, if applicable.
		  *
		  * Does nothing unless overridden.
		  * @param param	The parameter to add to the test
		  */
		virtual void addFloatParameter( double param );
		
		/**
		  * Add an integer parameter to the test, if applicable.
		  *
		  * Does nothing unless overridden.
		  * @param param	The parameter to add to the test
		  */
		virtual void addStringParameter( const std::string param );
		
		/**
		  * Indicates if the test is ready to be run (i.e. any necessary 
		  * parameters have been supplied).
		  *
		  * @return true if the test can be run meaningfully, false otherwise
		  */
		virtual bool ready() const;
		
		/**
		  * Returns the value the test is checking for, as a string.
		  *
		  * @return the test value in string form, implementation-dependent
		  */
		virtual std::string valueString() const = 0;
		
	protected:
		
		/**
		  * The option name.
		  */
		std::string _optName;
		
		/**
		  * The value that was last checked against.
		  */
		std::string _testedValue;
		
		/**
		  * All information has been provided and the test can be run.
		  */
		bool _ready;
	};
	
	/**
	  * A class for simple validation of AnyOption parameters.
	  */
	class AnyOptionValidator {
	public:
		/**
		  * Default constructor.  Allocates internal objects.
		  */
		AnyOptionValidator();
		
		/**
		  * Destructor.  Clears internal objects and releases memory held by
		  * option tests.
		  */
		~AnyOptionValidator();
		
		/**
		  * Creates and stores a type-independent parameter test.
		  * @param option The name of the option or flag to be tested
		  * @param option_type The type of test to be run.
		  * @return A pointer to the test criterion, so that additional 
		  *         parameters may be specified.
		  */
		OptionTest * addTest( const std::string option, OPTION_TEST_TYPE option_type );
		
		/**
		  * Runs all indicated tests.
		  *
		  * Runs all tests in unspecified order.  Failed tests can be retrieved with the
		  * getFailedTests() or printed with printFailedTests().
		  *
		  * @param opts 	A pointer containing the parameters
		  * @return true if all tests pass, false otherwise
		  * @see getFailedTests()
		  * @see printFailedTests()
		  */
		bool validateOptions( AnyOption *opts );
		
		/**
		  * Retrieves the failed test objects.
		  *
		  * @return a vector of pointers to the test objects reporting failure
		  */
		std::vector< NCPA::OptionTest * > getFailedTests() const;
		
		/**
		  * Prints the failure messages of all tests reporting failure.
		  *
		  * @param out The string to which to print them.
		  */
		void printFailedTests( std::ostream *out ) const;
		
		/**
		  * Prints the description messages of all tests.
		  *
		  * @param out The string to which to print them.
		  */
		void printDescriptions( std::ostream *out ) const;
		
	private:
		std::vector< NCPA::OptionTest * > _criteria;
		std::vector< NCPA::OptionTest * > _failed;
	};
	
	
	/** Test whether an option is present */
	class RequiredTest : public OptionTest {
	public:
		RequiredTest( const std::string option_name );
		bool validate(  AnyOption *opts );
		std::string description() const;
		std::string failureMessage() const;
		std::string valueString() const;
	};
	
	/** Test whether one and only one of a set of options or flags is present */
	class RadioButtonTest : public OptionTest {
	public:
		RadioButtonTest( const std::string option_name );
		bool validate(  AnyOption *opts );
		std::string description() const;
		std::string failureMessage() const;
		void addStringParameter( const std::string param );
		std::string valueString() const;
		std::vector< std::string > lastMatched() const;
		bool ready() const;
	private:
		std::vector< std::string > _buttons;
		std::vector< std::string > _matched;
	};

	/** Test whether an integer option is greater than a parameter */
	class IntegerGreaterThanTest : public OptionTest {
	public:
		IntegerGreaterThanTest( const std::string option_name );
		bool validate(  AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value;
	};

	/** Test whether an integer option is greater than or equal to a parameter */
	class IntegerGreaterThanOrEqualToTest : public OptionTest {
	public:
		IntegerGreaterThanOrEqualToTest( const std::string option_name );
		bool validate(  AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value;
	};

	/** Test whether an integer option is less than a parameter */
	class IntegerLessThanTest : public OptionTest {
	public:
		IntegerLessThanTest( const std::string option_name );
		bool validate(  AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value;
	};
	
	// Test whether an integer option is less than
	class IntegerLessThanOrEqualToTest : public OptionTest {
	public:
		IntegerLessThanOrEqualToTest( const std::string option_name );
		bool validate(  AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		int _value;
	};
	
	// Test whether an integer option is equal to
	class IntegerEqualToTest : public OptionTest {
	public:
		IntegerEqualToTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		int _value;
	};
	
	// Test whether an integer option is not equal to
	class IntegerNotEqualToTest : public OptionTest {
	public:
		IntegerNotEqualToTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		int _value;
	};
	
	
	// Test whether a floating point option is greater than
	class FloatGreaterThanTest : public OptionTest {
	public:
		FloatGreaterThanTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		double _value;
	};
	
	// Test whether a floating point option is greater than or equal to
	class FloatGreaterThanOrEqualToTest : public OptionTest {
	public:
		FloatGreaterThanOrEqualToTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		double _value;
	};

	// Test whether a floating point option is less than
	class FloatLessThanTest : public OptionTest {
	public:
		FloatLessThanTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		double _value;
	};
	
	// Test whether a floating point option is less than or equal to
	class FloatLessThanOrEqualToTest : public OptionTest {
	public:
		FloatLessThanOrEqualToTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		bool _trueIfEquals;
		double _value;
	};
	
	// Test whether a floating point option is equal to (standard caveats apply)
	class FloatEqualToTest : public OptionTest {
	public:
		FloatEqualToTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value;
	};
	
	// Test whether a floating point option is not equal to (standard caveats apply)
	class FloatNotEqualToTest : public OptionTest {
	public:
		FloatNotEqualToTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addFloatParameter( double param );
		std::string valueString() const;
	private:
		double _value;
	};
	
	// Test whether the length of a string is at least N characters
	class StringMinimumLengthTest : public OptionTest {
	public:
		StringMinimumLengthTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		unsigned int _value;
	};
	
	// Test whether the length of a string is at most N characters
	class StringMaximumLengthTest : public OptionTest {
	public:
		StringMaximumLengthTest( const std::string option_name );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		void addIntegerParameter( int param );
		std::string valueString() const;
	private:
		unsigned int _value;
	};
	
	
	// test whether a string is a member of a set or not
	class StringSetTest : public OptionTest {
	public:
		StringSetTest( const std::string option_name );
		bool validate(  AnyOption *opts );
		std::string description() const;
		std::string failureMessage() const;
		void addStringParameter( const std::string choice_name );
		std::string valueString() const;
		bool ready() const;
	private:
		std::vector< std::string > _choices;
	};
	
}

#endif