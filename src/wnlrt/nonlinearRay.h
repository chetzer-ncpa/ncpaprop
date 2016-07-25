//
//  nonlinearRay.h
//  
//
//  Created by Joel Lonzaga on 5/27/12.
//  Copyright (c) 2012 University of Mississippi. All rights reserved.
//
// Most references are to the article: "Modelling waveforms of infrasound 
// arrivals from impulsive sources using weakly non-linear ray theory" 
// by Joel Lonzaga, Roger Waxler, Jelle Assink, Carrick Talmadge
// in Geophysical Journal International (2015), JGI 200, pg 1347-1361


#ifndef _nonlinearRay_h
#define _nonlinearRay_h
#include <fftw3.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace std;

struct nonray{
    //vector<double> xx, yy, zz, rh, ps, sb, tt, Ur, Ui, ff, cc, tr; // Joel
    vector<double> xx, yy, zz, rh, ps, sb, tt, Ur, Ui, ff, cc, tr, ss;  // DV
    // note: ps is the scaling factor for acoustic pressure: p = u/ps where p = acoustic pressure
    vector< vector<double> > uu; // scaled acoustic pressure uu = pp*ps
    // DV
    //vector< vector<double> > UUr; // real (spectrum)
    //vector< vector<double> > UUi; // imag (spectrum)
    
};

struct waveform{
    vector<double> tt, pp; // time, pressure
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


// The meaning of these variables:
// rho0 = ambient density
// beta = coefficient of non-linearity
// phi = eikonal
// grad(phi) = slowness
// 
// om = 1 - grad(phi)*v0 = the Doppler effect induced by the wind with velocity v0
// vr = v0 + c0*nv = magnitude of the ray velocity; 
// where nv = unit vector pointing along the ray slowness (grad(phi)
//
// NN = rho0_i*c0_i*om_i/(Pmax^2*abs(J_i)*Vr_i) // 'i' means initial or reference value


class raypathParams
{
    double beta, NN, vr, om, ja; // NN is the normalizing factor; set by method normalizingFactor()
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
