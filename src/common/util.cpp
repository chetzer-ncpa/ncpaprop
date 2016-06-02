#include <cstdio>
#include <cmath>
#include <ctime>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <complex>
#include <cctype>
#include <cstring>
#include <vector>
#include <sstream>
#include "util.h"

#define PI 3.14159

double NCPA::min( double a, double b ) {
	return a < b ? a : b;
}

double NCPA::max( double a, double b ) {
	return a > b ? a : b;
}

float NCPA::min( float a, float b ) {
	return a < b ? a : b;
}

float NCPA::max( float a, float b ) {
	return a > b ? a : b;
}

int NCPA::min( int a, int b ) {
	return a < b ? a : b;
}

int NCPA::max( int a, int b ) {
	return a > b ? a : b;
}


float NCPA::max( const float *vals, int size ) {

	float maxval = vals[ 0 ];
	for (int i = 1; i < size; i++) {
		maxval = vals[i] > maxval ? vals[ i ] : maxval;
	}
	return maxval;
}

float NCPA::min( const float *vals, int size ) {

	float minval = vals[ 0 ];
	for (int i = 1; i < size; i++) {
		minval = vals[i] < minval ? vals[ i ] : minval;
	}
	return minval;
}

double NCPA::max( const double *vals, int size ) {

	double maxval = vals[ 0 ];
	for (int i = 1; i < size; i++) {
		maxval = vals[i] > maxval ? vals[ i ] : maxval;
	}
	return maxval;
}

double NCPA::min( const double *vals, int size ) {

	double minval = vals[ 0 ];
	for (int i = 1; i < size; i++) {
		minval = vals[i] < minval ? vals[ i ] : minval;
	}
	return minval;
}

int NCPA::max( const int *vals, int size ) {

	int maxval = vals[ 0 ];
	for (int i = 1; i < size; i++) {
		maxval = vals[i] > maxval ? vals[ i ] : maxval;
	}
	return maxval;
}

int NCPA::min( const int *vals, int size ) {

	int minval = vals[ 0 ];
	for (int i = 1; i < size; i++) {
		minval = vals[i] < minval ? vals[ i ] : minval;
	}
	return minval;
}

std::string NCPA::timeAsString(double d) {
    time_t temptime = (time_t)d;
    tm* uttime = std::gmtime( &temptime );
    double ipart, fpart;
    fpart = std::modf( d, &ipart);

    char* holder = new char[ 28 ];
    std::sprintf(holder,"%4d-%02d-%02d %02d:%02d:%02d.%03d GMT",
            uttime->tm_year+1900, uttime->tm_mon+1, uttime->tm_mday,
            uttime->tm_hour, uttime->tm_min, uttime->tm_sec,
            (int)(round(fpart * 1000)) );
    std::string s = holder;
    return s;
}

bool NCPA::fexists( const char *filename ) {
	std::ifstream ifile( filename );
	bool tf = ifile.good();
	ifile.close();
	return tf;
}

std::string NCPA::deblank( const std::string orig ) {
        int index = orig.length();
        while (orig[(int)--index] == ' ') ;
        std::string temp = orig.substr(0,index+1);
        index = -1;
        while (temp[++index] == ' ') ;
        return temp.substr(index);
}

double NCPA::deg2rad( double d ) {
	return d * PI / 180.0;
}

double NCPA::rad2deg( double d ) {
	return d * 180.0 / PI;
}

// Acts like std::getline, but checks for all three permutations of EOL characters
std::istream &NCPA::safe_getline( std::istream &is, std::string &s ) {
	
	char ch;
	s.clear();
	
	// Keep going until you get to a linefeed character
	while (is.get(ch) && ch != '\n' && ch != '\r')
		s += ch;
	
	// DOS systems use 2 consecutive characters, so make sure there isn't another one lurking
	if (ch == '\r') {
		ch = is.peek();
		if (ch == '\n') {
			is.get(ch);
		}
	}
	return is;
}

// perform a simple linear interpolation
double NCPA::linearInterp( double x1, double y1, double x2, double y2, double x ) {
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

void NCPA::toLowerCase( char *in ) {
	unsigned int len = std::strlen(in);
	for (unsigned int i = 0; i < len; i++) {
		in[ i ] = (char)std::tolower( in[ i ] );
	}
}

void NCPA::toUpperCase( char *in ) {
	unsigned int len = std::strlen(in);
	for (unsigned int i = 0; i < len; i++) {
		in[ i ] = (char)std::toupper( in[ i ] );
	}
}

std::string NCPA::toLowerCase( const std::string in ) {
	std::ostringstream oss("");
	for (unsigned int i = 0; i < in.length(); i++) {
		oss << (char)std::tolower( in[ i ] );
	}
	return oss.str();
}

std::string NCPA::toUpperCase( const std::string in ) {
	std::ostringstream oss("");
	for (unsigned int i = 0; i < in.length(); i++) {
		oss << (char)std::toupper( in[ i ] );
	}
	return oss.str();
}

std::vector< std::string > NCPA::split( std::string input, std::string delimiters ) {
	// lean on strtok for this functionality
	char *holder = new char[ input.size() + 1 ];
	std::memset( holder, '\0', input.size() + 1 );
	std::strcpy( holder, input.c_str() );
	char *tmp = strtok( holder, delimiters.c_str() );
	std::vector< std::string > tokens;
	tokens.clear();
	while (tmp != NULL) {
		std::string s( tmp );
		tokens.push_back( s );
		tmp = strtok(NULL,delimiters.c_str());
	}
	delete [] holder;
	return tokens;
}

bool NCPA::checkAzimuthLimits( double toCheck, double target, double tolerance ) {
	
	return ( std::cos(target - toCheck) >= std::cos( tolerance ) );
}

double NCPA::normalizeAzimuth( double in ) {
	double out = in;
	while (out < 0) {
		out += 360;
	}
	while (out >= 360) {
		out -= 360;
	}
	return out;
}