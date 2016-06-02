#ifndef ATMLIB_H
#define ATMLIB_H

//#include <string>

#define PI acos(-1.0)
#define GAMMA 1.4
#define GASCONSTANT 287.058

#include <iostream>

typedef long long int64;

class AtmLibrary {

private:
// put here the local variables?
/*      static bool bigEndianNative;
        static bool endianChecked;
        float swap( float );
        int swap( int );
        double swap( double );
        short swap( short );
        int64 swap( int64 ); */

public:
        void readG2SBinary( std::ifstream*, double*, double*, double**, double**, double**, double**, double**);
        void readG2SAscii( std::ifstream*, double*, double*, double*, double*, double*, double*);
        void getNumberOfLines( std::ifstream*, int*);
        void getBinaryG2SDimensions( std::ifstream*, int*, int*);
        void readEOFs( std::ifstream*, int, int, double**);
        void makeProfile( int, int, double*, double**, double*);
        void makeZeroProfile( int, double*);
        void makeGaussianProfile( int, double, double, double, double*, double*);
        //void writeProfile(char*,int,double*,double*,double*,double*,double*,double*);
        void writeProfile(char*,int,double,double*,double*,double*,double*,double*,double*);
        //void getAbsorptionCoefficients(int,double,double*,double*,double*,double*,double*);
        void getAbsorptionCoefficients(int,double,double*,double*,double*,double*,std::string,double*);
        
        
        AtmLibrary();          // this is called 'the constructor'
};

#endif // ATMLIB_H
