#ifndef _NCPA_PARAMETERSET_H_
#define _NCPA_PARAMETERSET_H_

#include <vector>
#include <string>

namespace NCPA {
	
	/**
	  * These indicate tests that do not depend on the type of the tested
	  * value.
	  */
	enum OPTION_TEST_TYPE : unsigned int {
	
		/** This option/flag must be present. */
		OPTION_REQUIRED,
	
		/** 
		  * This option/flag must be present if at least one of a set of 
		  * other options is present. 
		*/
		OPTION_REQUIRED_IF,
	
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
	
	
	
	

	class GenericParameter {
	public:
		~GenericParameter();
		virtual bool validate() = 0;
		virtual OptionTest * addTest( OPTION_TEST_TYPE option_type );
		
	protected:
		std::string _key;
		std::vector< OptionTest * > _tests;
	};

	
	class IntegerParameter : public GenericParameter {
	public:
		IntegerParameter( std::string key );
		IntegerParameter( std::string key, int defaultValue );
		bool validate();
		void setValue( int newval );
	protected:
		int _value;
	};
	
	class FloatParameter : public GenericParameter {
	public:
		FloatParameter( std::string key );
		FloatParameter( std::string key, double defaultValue );
		bool validate();
		void setValue( double newval );
	protected:
		double _value;
	};
	
	class StringParameter : public GenericParameter {
	public:
		StringParameter( std::string key );
		StringParameter( std::string key, std::string defaultValue );
		bool validate();
		void setValue( std::string newval );
	protected:
		std::string _value;
	};
	
	class FlagParameter : public GenericParameter {
	public:
		FlagParameter( std::string key );
		FlagParameter( std::string key, bool defaultValue );
		bool validate();
		void setValue( bool newval );
	protected:
		bool _value;
	};


	class ParameterSet {
	public:
		ParameterSet();
		~ParameterSet();
		
		unsigned int parseCommandLine( int argc, const char **argv );
		void parseFile( std::string filename );
		void setDelimiters( std::string newdelims );
		void addParameter( GenericParameter *newParam );
		bool validate();
		std::vector< std::string > getUnparsedOptions();
		void setStrict( bool tf );
		
		int getIntegerValue( std::string key );
		double getFloatValue( std::string key );
		std::string getStringValue( std::string key );
		bool getBoolValue( std::string key );
		
	protected:
		std::vector< GenericParameter * > _params;
		std::string _delims;
		std::vector< std::string > _unparsed;
		bool _strict;
		
		bool _isLongOption( std::string );
		bool _isShortOption( std::string );
		unsigned int _processShortOption( int argc, const char **argv, unsigned int index );
		unsigned int _processLongOption( int argc, const char **argv, unsigned int index );
	};
}




#endif