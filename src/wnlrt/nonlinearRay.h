//
//  nonlinearRay.h
//  
//
//  Created by Joel Lonzaga on 5/27/12.
//  Copyright (c) 2012 University of Mississippi. All rights reserved.
//


#ifndef _nonlinearRay_h
#define _nonlinearRay_h
#include <fftw3.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace std;

struct nonray{
    //vector<double> xx, yy, zz, rh, ps, sb, tt, Ur, Ui, ff, cc, tr;
    vector<double> xx, yy, zz, rh, ps, sb, tt, Ur, Ui, ff, cc, tr, ss;  // DV
    vector< vector<double> > uu;
};

struct waveform{
    vector<double> tt, pp;
};

struct blastParam{
    double pks, ttd;
};

struct attRed{ // attenuation reduction by factor attfac above an altitude zzRed
    double zzRed, atfac;
};


void savewf( const char*, waveform );

waveform loadWaveform( const char*, double, double );

waveform loadBlastModel( const char*, int, double, double, double );

waveform generateWaveForm( int, double, double, double& );

waveform generateWaveForm( int, double, double& );

waveform blastmodel( int, double, double, double, linray );

waveform generateNWave( int, double, double, double );

blastParam pkOverpress( double, double, linray );

double posPhaseDur( double, double );

vector< vector<double> > dataReader( const char*, int );

double wfResize( int, waveform, double*, double* );

void causticParams( vector<double>, linray, double*, double*, double* );

double stepSizeCaustic( double, double, double, double, double, double& );

void atmAttn(int, double, double, double, double*, double* );

class spectralParams
{
    int nn;
    double *cf, *fltr, *attnEf1, *attnEf2, *usq;
    fftw_complex *Usq;
    
public:
    spectralParams( int );
    ~spectralParams();
    void frequency( double, double* );
    void initFFT( double, double*, double* );
    void integEulerMethod( double, double, double, double, double* );
    void integLeapfrog( double, double );
    void applyAttn( double, double, double, double* );
    void FFT( double*, double& );
    void hilbert();
    fftw_complex *U1, *U2, *U3;
};


class raypathParams
{
    double beta, NN, vr, om, ja;
    double *spl[9], ptNlp[9];
	  vector< vector<double> > nlp;
	  linray LR;
    
public:
	  raypathParams( linray );
	  raypathParams( linray, int );
    ~raypathParams();
	  void physicalParams( double );
	  void physicalParams( double, int );
    double pressCorrection( double );
    void normalizingFactor( double, double);
    void psDopScaBeta( double, double&, double&);
    double scaledBetaCaustic( double );
	  double xx, yy, zz, cc, ps, tr, rh;
	  gsl_interp_accel *acc_; // DV
    gsl_spline *splin[9];   // DV
};

nonray nonlinearRay( int, double, attRed, linray, waveform );


#endif
