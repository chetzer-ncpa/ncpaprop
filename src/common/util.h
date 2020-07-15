#ifndef _UTIL_H_
#define _UTIL_H_

#include <string>
#include <istream>
#include <vector>
#include <complex>

#ifndef PI
#define PI 3.141592653589793
#endif

namespace NCPA {
	float min( const float *, int );
	float max( const float *, int );
	double min( const double *, int );
	double max( const double *, int );
	int min( const int *, int );
	int max( const int *, int );
	double min( double, double );
	double max( double, double );
	float min( float, float );
	float max( float, float );
	int min( int, int );
	int max( int, int );
	std::string timeAsString( double d );
	bool fexists( const char *filename );
	std::string deblank( const std::string& orig );
	double deg2rad( double );
	double rad2deg( double );
	double mean( double*, size_t );
	
	std::istream &safe_getline( std::istream &is, std::string &s );
	
	// perform a simple linear interpolation
	double linearInterp( double x1, double y1, double x2, double y2, double x );
	void toUpperCase( char *in );
	void toLowerCase( char *in );
	std::string toUpperCase( const std::string in );
	std::string toLowerCase( const std::string in );
	
	// Split a string into more strings by tokenizing
	std::vector< std::string > split( std::string input, std::string delimiters = " \t\n" );
	std::string deblank( const std::string& str, const std::string& whitespace );

	bool checkAzimuthLimits( double toCheck, double target, double tolerance );
	double normalizeAzimuth( double in );
	
	// Utility functions
	double **dmatrix(long nr, long nc);
	int free_dmatrix(double** v, long nr, long nc);
	std::complex<double> **cmatrix(long nr, long nc);
	int free_cmatrix(std::complex<double> **v, long nr, long nc);
	std::complex<double> ***c3Darray(size_t xlen, size_t ylen, size_t zlen);
	void free_c3Darray(std::complex<double> ***data, size_t xlen, size_t ylen);

	//int count_rows_arbcol(const char *filename);
	int count_rows_arbcol( const std::string& filename );

	size_t find_closest_index( double *z, size_t NZ, double zs );
	
//	double phase( std::complex< double > );
	
	// Constants
	//double R = 287.0;
	//double GAM = 1.4;
	//double PI=3.14159;
	
}




#endif
