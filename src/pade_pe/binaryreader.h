#ifndef BINARYREADER_H
#define BINARYREADER_H

#include <iostream>

typedef long long int64;

class BinaryReader {

private:
        static bool bigEndianNative;
        static bool endianChecked;
        float swap( float );
        int swap( int );
        double swap( double );
        short swap( short );
        int64 swap( int64 );

public:
        void readBigFloatArray( std::ifstream*, int, float* );
        void readLittleFloatArray( std::ifstream*, int, float* );
        void readFloatArray( std::ifstream*, int, float*, bool );
        void readBigIntArray( std::ifstream*, int, int* );
        void readLittleIntArray( std::ifstream*, int, int* );
        void readIntArray( std::ifstream*, int, int*, bool );
        void readBigDoubleArray( std::ifstream*, int, double* );
        void readLittleDoubleArray( std::ifstream*, int, double* );
        void readDoubleArray( std::ifstream*, int, double*, bool );
        void readBigShortArray( std::ifstream*, int, short* );
        void readLittleShortArray( std::ifstream*, int, short* );
        void readShortArray( std::ifstream*, int, short*, bool );
        void readBigLongArray( std::ifstream*, int, int64* );
        void readLittleLongArray( std::ifstream*, int, int64* );
        void readLongArray( std::ifstream*, int, int64*, bool );

        BinaryReader();

};

#endif // BINARYREADER_H
