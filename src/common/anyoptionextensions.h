/**
  * Extensions of and tools used with the AnyOption class presented in anyoption.h
  */

#ifndef __NCPA_ANYOPTION_EXTENSIONS__
#define __NCPA_ANYOPTION_EXTENSIONS__

#include "anyoption.h"
#include <map>
#include <string>
#include <vector>

namespace NCPA {
	
	/**
	  * These indicate tests that do not depend on the type of the tested
	  * value.
	  */
	enum OPTION_NOTYPE_TEST_TYPE : unsigned int {
		
		/** This option/flag must be present. */
		OPTION_NOTYPE_REQUIRED,
		
		/** 
		  * Designates a group of options, one and only one of which
		  * must be present.
		  */
		OPTION_NOTYPE_RADIO_BUTTON
	};
	
	/**
	  * These indicate tests that apply to integer values
	  */
	enum OPTION_INTEGER_TEST_TYPE : unsigned int {
		
		/** Integer that must be > 0 */
		OPTION_INTEGER_POSITIVE,
		
		/** Integer that must be < 0 */
		OPTION_INTEGER_NEGATIVE,
		
		/** Integer that must be > another integer */
		OPTION_INTEGER_GREATER_THAN,
		
		/** Integer that must be >= another integer */
		OPTION_INTEGER_GREATER_THAN_EQUAL,
		
		/** Integer that must be < another integer */
		OPTION_INTEGER_LESS_THAN,
		
		/** Integer that must be <= another integer */
		OPTION_INTEGER_LESS_THAN_EQUAL,
		
		/** Integer that must == 0 */
		OPTION_INTEGER_ZERO,
		
		/** Integer that must != 0 */
		OPTION_INTEGER_NONZERO,
		
		/** Integer that must == another integer */
		OPTION_INTEGER_EQUAL,
		
		/** Integer that must != another integer */
		OPTION_INTEGER_NOT_EQUAL
	};
	
	/**
	  * These indicate tests that apply to double values
	  */
	enum OPTION_FLOAT_TEST_TYPE : unsigned int {
		
		/** Double that must > 0.0 */
		OPTION_FLOAT_POSITIVE,
		
		/** Double that must < 0.0 */
		OPTION_FLOAT_NEGATIVE,
		
		/** Double that must > another double */
		OPTION_FLOAT_GREATER_THAN,
		
		/** Double that must >= another double */
		OPTION_FLOAT_GREATER_THAN_EQUAL,
		
		/** Double that must < another double */
		OPTION_FLOAT_LESS_THAN,
		
		/** Double that must <= another double */
		OPTION_FLOAT_LESS_THAN_EQUAL,
		
		/** Double that must == 0.0 (standard floating point caveats apply) */
		OPTION_FLOAT_ZERO,
		
		/** Double that must != 0.0 (standard floating point caveats apply) */
		OPTION_FLOAT_NONZERO			
	};
	
	/**
	  * These indicate tests that apply to string values
	  */
	enum OPTION_STRING_TEST_TYPE : unsigned int {
		
		/** string .size() must be >= an integer */
		OPTION_STRING_MINIMUM_LENGTH,
		
		/** string .size() must be <= an integer */
		OPTION_STRING_MAXIMUM_LENGTH
	};

	/**
	  * Abstract base class for validation criteria.  See subclass descriptions
	  * for general usage, you won't instantiate this class directly.
	  */
	class OptionValidationCriterion {
	public:
		/**
		  * Destructor.  Cleans up any dynamically allocated memory.
		  */
		virtual ~OptionValidationCriterion();
		
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
		
	protected:
		
		/**
		  * The option name.
		  */
		std::string _optName;
		
		/**
		  * The value that was last checked against.
		  */
		std::string _testedValue;
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
		OptionValidationCriterion * addOption( const std::string &option, OPTION_NOTYPE_TEST_TYPE option_type );
		OptionValidationCriterion * addOption( const std::string &option, OPTION_INTEGER_TEST_TYPE option_type );
		OptionValidationCriterion * addOption( const std::string &option, OPTION_FLOAT_TEST_TYPE option_type );
		
		// One argument required
		OptionValidationCriterion * addOption( const std::string &option, OPTION_INTEGER_TEST_TYPE option_type, int boundary1 );
		OptionValidationCriterion * addOption( const std::string &option, OPTION_FLOAT_TEST_TYPE option_type, double boundary1 );
		OptionValidationCriterion * addOption( const std::string &option, OPTION_STRING_TEST_TYPE option_type, int boundary1 );
		
		bool validateOptions( AnyOption *opts );
		std::vector< NCPA::OptionValidationCriterion * > getFailedChecks() const;
		void printFailedChecks( std::ostream *out ) const;
		void printDescriptions( std::ostream *out ) const;
	private:
		std::vector< NCPA::OptionValidationCriterion * > _criteria;
		std::vector< NCPA::OptionValidationCriterion * > _failed;
	};
	
	
	// Test whether an option is present
	class RequiredCriterion : public OptionValidationCriterion {
	public:
		RequiredCriterion( const std::string option_name );
		bool validate(  AnyOption *opts );
		std::string description() const;
		std::string failureMessage() const;
	};
	
	// Test whether one and only one of a set is 
	class RadioButtonCriterion : public OptionValidationCriterion {
	public:
		RadioButtonCriterion( const std::string option_name );
		bool validate(  AnyOption *opts );
		std::string description() const;
		std::string failureMessage() const;
		void addParameter( const std::string option_name );
		std::string joinedString() const;
		std::vector< std::string > lastMatched() const;
	private:
		std::vector< std::string > _buttons;
		std::vector< std::string > _matched;
	};

	// Test whether an integer option is greater than (optionally or equal to)
	class IntegerGreaterThanCriterion : public OptionValidationCriterion {
	public:
		IntegerGreaterThanCriterion( const std::string option_name, int comparison, bool trueIfEquals );
		bool validate(  AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		bool _trueIfEquals;
		int _value;
	};

	// Test whether an integer option is less than (optionally or equal to)
	class IntegerLessThanCriterion : public OptionValidationCriterion {
	public:
		IntegerLessThanCriterion( const std::string option_name, int comparison, bool trueIfEquals );
		bool validate(  AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		bool _trueIfEquals;
		int _value;
	};
	
	// Test whether an integer option is equal to
	class IntegerEqualsCriterion : public OptionValidationCriterion {
	public:
		IntegerEqualsCriterion( const std::string option_name, int comparison );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		int _value;
	};
	
	// Test whether an integer option is not equal to
	class IntegerNotEqualsCriterion : public OptionValidationCriterion {
	public:
		IntegerNotEqualsCriterion( const std::string option_name, int comparison );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		int _value;
	};
	
	// Test whether an integer option is greater than (optionally or equal to)
	class IntegerBetweenCriterion : public OptionValidationCriterion {
	public:
		IntegerBetweenCriterion( const std::string option_name, int minval, int maxval, bool trueIfEquals );
		bool validate(  AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		bool _trueIfEquals;
		int _minval, _maxval;
	};
	
	
	// Test whether a floating point option is greater than (optionally or equal to)
	class FloatGreaterThanCriterion : public OptionValidationCriterion {
	public:
		FloatGreaterThanCriterion( const std::string option_name, double comparison, bool trueIfEquals );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		bool _trueIfEquals;
		double _value;
	};

	// Test whether a floating point option is less than (optionally or equal to)
	class FloatLessThanCriterion : public OptionValidationCriterion {
	public:
		FloatLessThanCriterion( const std::string option_name, double comparison, bool trueIfEquals );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		bool _trueIfEquals;
		double _value;
	};
	
	// Test whether a floating point option is equal to (standard caveats apply)
	class FloatEqualsCriterion : public OptionValidationCriterion {
	public:
		FloatEqualsCriterion( const std::string option_name, double comparison );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		double _value;
	};
	
	// Test whether a floating point option is not equal to (standard caveats apply)
	class FloatNotEqualsCriterion : public OptionValidationCriterion {
	public:
		FloatNotEqualsCriterion( const std::string option_name, double comparison );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		double _value;
	};
	
	// Test whether the length of a string is at least N characters
	class StringMinimumLengthCriterion : public OptionValidationCriterion {
	public:
		StringMinimumLengthCriterion( const std::string option_name, int minLength );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		int _value;
	};
	
	// Test whether the length of a string is at most N characters
	class StringMaximumLengthCriterion : public OptionValidationCriterion {
	public:
		StringMaximumLengthCriterion( const std::string option_name, int maxLength );
		bool validate( AnyOption *opts ) ;
		std::string description() const;
		std::string failureMessage() const;
		
	private:
		int _value;
	};
	
}

#endif