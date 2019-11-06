#ifndef __NCPA_ANYOPTION_VALIDATOR__
#define __NCPA_ANYOPTION_VALIDATOR__

#include "anyoption.h"
#include <map>
#include <string>
#include <vector>

namespace NCPA {
	
	enum OPTION_NOTYPE_TEST_TYPE : unsigned int {
		OPTION_REQUIRED
	};
	
	enum OPTION_INTEGER_TEST_TYPE : unsigned int {
		// integer tests.  Options are assumed to be optional and will return true if not present.
		OPTION_INTEGER_POSITIVE,
		OPTION_INTEGER_NEGATIVE,
		OPTION_INTEGER_GREATER_THAN,
		OPTION_INTEGER_GREATER_THAN_EQUAL,
		OPTION_INTEGER_LESS_THAN,
		OPTION_INTEGER_LESS_THAN_EQUAL,
		OPTION_INTEGER_ZERO,
		OPTION_INTEGER_NONZERO,
		OPTION_INTEGER_EQUAL,
		OPTION_INTEGER_NOT_EQUAL
	};
	
	enum OPTION_FLOAT_TEST_TYPE : unsigned int {
		// floating point tests.  Options are assumed to be optional and will return true if not present.
		OPTION_FLOAT_POSITIVE,
		OPTION_FLOAT_NEGATIVE,
		OPTION_FLOAT_GREATER_THAN,
		OPTION_FLOAT_GREATER_THAN_EQUAL,
		OPTION_FLOAT_LESS_THAN,
		OPTION_FLOAT_LESS_THAN_EQUAL,
		OPTION_FLOAT_ZERO,
		OPTION_FLOAT_NONZERO			
	};
	
	enum OPTION_STRING_TEST_TYPE : unsigned int {
		// string tests.  Options are assumed to be optional and will return true if not present.
		OPTION_STRING_MINIMUM_LENGTH,
		OPTION_STRING_MAXIMUM_LENGTH
	};

	// Abstract base class for validation criteria
	class OptionValidationCriterion {
	public:
		virtual ~OptionValidationCriterion();
		virtual bool validate( AnyOption *opts ) = 0;
		virtual std::string description() const = 0;
		virtual std::string failureMessage() const = 0;
		virtual std::string optionName() const;
	protected:
		std::string _optName;
		std::string _testedValue;
	};
	
	
	class AnyOptionValidator {
	public:
		AnyOptionValidator();
		~AnyOptionValidator();
		
		// No argument required
		void addOption( const std::string &option, OPTION_NOTYPE_TEST_TYPE option_type );
		void addOption( const std::string &option, OPTION_INTEGER_TEST_TYPE option_type );
		void addOption( const std::string &option, OPTION_FLOAT_TEST_TYPE option_type );
		
		// One argument required
		void addOption( const std::string &option, OPTION_INTEGER_TEST_TYPE option_type, int boundary1 );
		void addOption( const std::string &option, OPTION_FLOAT_TEST_TYPE option_type, double boundary1 );
		void addOption( const std::string &option, OPTION_STRING_TEST_TYPE option_type, int boundary1 );
		
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