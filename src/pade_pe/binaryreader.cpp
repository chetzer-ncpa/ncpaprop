#include "binaryreader.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
/*
    BinaryReader.cpp: Reads binary data from files and returns arrays of
    numbers.
*/

/*
    Typedef shortcut for unsigned chars
*/
typedef unsigned char byte;

bool BinaryReader::bigEndianNative = false;
bool BinaryReader::endianChecked = false;

BinaryReader::BinaryReader()
{
        if (!endianChecked) {
                // Cute test for the endian-ness of the system.  Big thanks to a program
                // by Promit Roy, who credits the source code for Quake 2.
                byte SwapTest[2] = { 1, 0 };

                if( *(short *) SwapTest == 1 )
                {
                        bigEndianNative = false;
                } else {
                        bigEndianNative = true;
                }
                endianChecked = true;
        }
        // InitEndian();
}

void BinaryReader::readFloatArray( ifstream* infile, int samples, float* buffer, bool bigEndian )
{
        if (bigEndian) {
                readBigFloatArray( infile, samples, buffer );
        } else {
                readLittleFloatArray( infile, samples, buffer );
        }
}


 void BinaryReader::readBigFloatArray( ifstream* infile, int samples, float* buffer ) {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void BinaryReader::readLittleFloatArray( ifstream* infile, int samples, float* buffer ) {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

float BinaryReader::swap( float f )
{
        union {
                byte b[4];
                float f;
        } u1,u2;

        u1.f = f;
        u2.b[0] = u1.b[3];
        u2.b[1] = u1.b[2];
        u2.b[2] = u1.b[1];
        u2.b[3] = u1.b[0];
        return u2.f;
}



  void BinaryReader::readIntArray( ifstream* infile, int samples, int* buffer, bool bigEndian )
{
        if (bigEndian) {
                readBigIntArray( infile, samples, buffer );
        } else {
                readLittleIntArray( infile, samples, buffer );
        }
}


 void BinaryReader::readBigIntArray( ifstream* infile, int samples, int* buffer ) {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void BinaryReader::readLittleIntArray( ifstream* infile, int samples, int* buffer ) {
        int nBytes = samples * 4;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

int BinaryReader::swap( int i )
{
        union {
                byte b[4];
                int i;
        } u1,u2;

        u1.i = i;
        u2.b[0] = u1.b[3];
        u2.b[1] = u1.b[2];
        u2.b[2] = u1.b[1];
        u2.b[3] = u1.b[0];

        return u2.i;
}

  void BinaryReader::readDoubleArray( ifstream* infile, int samples, double* buffer, bool bigEndian )
{
        if (bigEndian) {
                readBigDoubleArray( infile, samples, buffer );
        } else {
                readLittleDoubleArray( infile, samples, buffer );
        }
}


 void BinaryReader::readBigDoubleArray( ifstream* infile, int samples, double* buffer ) {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void BinaryReader::readLittleDoubleArray( ifstream* infile, int samples, double* buffer ) {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

double BinaryReader::swap( double d )
{
        union {
                byte b[8];
                double d;
        } u1,u2;
        u1.d = d;

        for (int k = 0; k < 8; k++) {
            u2.b[k] = u1.b[8-1-k];
        }

        return u2.d;
}

  void BinaryReader::readShortArray( ifstream* infile, int samples, short* buffer, bool bigEndian )
{
        if (bigEndian) {
                readBigShortArray( infile, samples, buffer );
        } else {
                readLittleShortArray( infile, samples, buffer );
        }
}


 void BinaryReader::readBigShortArray( ifstream* infile, int samples, short* buffer ) {
        int nBytes = samples * 2;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void BinaryReader::readLittleShortArray( ifstream* infile, int samples, short* buffer ) {
        int nBytes = samples * 2;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

short BinaryReader::swap( short s )
{
        union {
                byte b[2];
                short s;
        } u1,u2;
        u1.s = s;

        for (int k = 0; k < 2; k++) {
            u2.b[k] = u1.b[1-k];
        }

        return u2.s;
}

void BinaryReader::readLongArray( ifstream* infile, int samples, int64* buffer, bool bigEndian )
{
        if (bigEndian) {
                readBigLongArray( infile, samples, buffer );
        } else {
                readLittleLongArray( infile, samples, buffer );
        }
}


 void BinaryReader::readBigLongArray( ifstream* infile, int samples, int64* buffer ) {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (!bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

 void BinaryReader::readLittleLongArray( ifstream* infile, int samples, int64* buffer ) {
        int nBytes = samples * 8;
        infile->read((char*)buffer,nBytes);
        if (bigEndianNative) {
            for (int i = 0; i < samples; i++) {
                    buffer[i] = swap( buffer[i] );
                }
        }
}

int64 BinaryReader::swap( int64 s )
{
        union {
                byte b[8];
                int64 s;
        } u1,u2;
        u1.s = s;

        for (int k = 0; k < 8; k++) {
            u2.b[k] = u1.b[7-k];
        }

        return u2.s;
}
