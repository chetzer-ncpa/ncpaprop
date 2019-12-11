#ifndef NCPA_PARAMETER_SET_H_
#define NCPA_PARAMETER_SET_H_

#include <map>

namespace NCPA {

	class ParameterSet {
	
	public:
		ParameterSet();
		~ParameterSet();
		
		void setStringValue( std::string key, std::string val );
		void setIntegerValue( std::string key, int val );
		void setDoubleValue( std::string key, double val );
		
		std::string getStringValue( const std::string &key ) const;
		int getIntegerValue( const std::string &key ) const;
		double getDoubleValue( const std::string &key ) const;
		
	protected:
		const unsigned int STRING_TYPE_ = 0, INT_TYPE_ = 1, DOUBLE_TYPE_ = 2;
		std::map< std::string, unsigned int > type_map_;
		std::map< std::string, int > int_value_map_;
		std::map< std::string, double > double_value_map_;
		std::map< std::string, std::string > string_value_map_;
		
	};

}









#endif