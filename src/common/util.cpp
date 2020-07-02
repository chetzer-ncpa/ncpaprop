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


#define PI 3.141592653589793

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

/*
std::string NCPA::deblank( const std::string orig ) {
        int index = orig.length();
        while (orig[(int)--index] == ' ') ;
        std::string temp = orig.substr(0,index+1);
        index = -1;
        while (temp[++index] == ' ') ;
        return temp.substr(index);
}
*/

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
	// // lean on strtok for this functionality
	// char *holder = new char[ input.size() + 1 ];
	// std::memset( holder, '\0', input.size() + 1 );
	// std::strcpy( holder, input.c_str() );
	// char *tmp = strtok( holder, delimiters.c_str() );
	// std::vector< std::string > tokens;
	// tokens.clear();
	// while (tmp != NULL) {
	// 	std::string s( tmp );
	// 	tokens.push_back( s );
	// 	tmp = strtok(NULL,delimiters.c_str());
	// }
	// delete [] holder;
	// return tokens;

	std::vector< std::string > tokens;
	tokens.clear();
	size_t firstind, lastind;
	firstind = input.find_first_not_of( delimiters );
	lastind = input.find_first_of( delimiters, firstind );
	while( firstind != std::string::npos ) {
		tokens.push_back( input.substr( firstind, lastind - firstind ) );
		firstind = input.find_first_not_of( delimiters, lastind );
		lastind = input.find_first_of( delimiters, firstind );
	}
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

// Following functions came from various <module>_lib files
double** NCPA::dmatrix(long nr, long nc) {
  // allocate a double matrix
  double **v;
  v = new double* [nr];
  for (long i=0; i<nr; i++) {
      v[i] = new double [nc];
  }
  return v;
}


int NCPA::free_dmatrix(double **v, long nr, long nc) {
  // free a double matrix
  for (long i=0; i<nr; i++) {
      delete v[i];
  }
  delete v;
  return 0;
}

std::complex<double> **NCPA::cmatrix(long nr, long nc) {
  // allocate a complex<double> matrix
  std::complex<double> **v;
  v = new std::complex<double>* [nr];
  for (long i=0; i<nr; i++) {
      v[i] = new std::complex<double> [nc];
  }
  return v;
}

int NCPA::free_cmatrix( std::complex<double> **v, long nr, long nc) {
  // free a complex<double> matrix
  for (long i=0; i<nr; i++) {
      delete v[i];
  }
  delete v;
  return 0;
}

std::complex<double> ***NCPA::c3Darray(size_t xlen, size_t ylen, size_t zlen)
{
  std::complex<double> ***p;
  size_t i, j;

  if ((p = (std::complex<double> ***) calloc(xlen, (sizeof(*p)))) == NULL) {
      perror("calloc 1");
      return NULL;
  }

  for (i=0; i < xlen; ++i)
      p[i] = NULL;

  for (i=0; i < xlen; ++i)
      if ((p[i] = (std::complex<double> **) calloc(ylen, (sizeof (*p[i])))) == NULL) {
          perror("calloc 2");
          NCPA::free_c3Darray(p, xlen, ylen);
          return NULL;
      }

  for (i=0; i < xlen; ++i)
      for (j=0; j < ylen; ++j)
          p[i][j] = NULL;

  for (i=0; i < xlen; ++i)
      for (j=0; j < ylen; ++j)
          if ((p[i][j] = (std::complex<double> *) calloc(zlen, (sizeof (*p[i][j])))) == NULL) {
              perror("calloc 3");
              NCPA::free_c3Darray(p, xlen, ylen);
              return NULL;
          }
  return p;
}


void NCPA::free_c3Darray(std::complex<double> ***data, size_t xlen, size_t ylen)
{
  size_t i, j;
  for (i=0; i < xlen; ++i) {
      if (data[i] != NULL) {
          for (j=0; j < ylen; ++j)
              free(data[i][j]);
          free(data[i]);
      }
  }
  free(data);
}

std::string NCPA::deblank( const std::string& str ) {
	return NCPA::deblank( str, " \t\n" );
}

std::string NCPA::deblank( const std::string& str, const std::string& whitespace ) {
	const size_t strBegin = str.find_first_not_of( whitespace );
	if (strBegin == std::string::npos) {
		return "";
	}

	const size_t strEnd = str.find_last_not_of( whitespace );
	const size_t strRange = strEnd - strBegin + 1;

	return str.substr( strBegin, strRange );
}

int NCPA::count_rows_arbcol(const std::string& filename) {
  int answer,c;
  FILE *f=fopen(filename.c_str(),"r");

  if(f==NULL) {
      std::ostringstream es;
      es << "file << " << filename << " could not be opened.\n";
      throw std::invalid_argument(es.str());
  }
  answer=0;
  //read_header(f);
  while((c=getc(f))!=EOF){
      if(c=='\n') answer=answer+1;
  }
  fclose(f);
  return answer;
}