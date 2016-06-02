//
//  linearRay3DStrat.h
//  
//
//  Created by Joel Lonzaga on 6/23/12.
//  Copyright (c) 2012 University of Mississippi. All rights reserved.
//

#ifndef _linearRay3DStrat_h
#define _linearRay3DStrat_h

#include <fstream>
#include <iostream>
#include <cstdlib>
#include<vector>
#include <cmath>
#define PI 3.141592653589

using namespace std;

double* spline( vector<double>, vector<double>, double, double );

double splint( vector<double>, vector<double>, double*, double x);

void splint( vector<double>, vector< vector<double> >, double**, double, double*);

void splintDeriv( vector<double>, vector< vector<double> >, double**, double, double [][2] );

void deriv1EndPoints( double, vector<double>, double* );

vector<double> deriv1( double, vector<double> );

vector<double> deriv1( vector<double>, vector<double> );

vector<double> deriv2( double, vector<double> );

vector<double> deriv2( vector<double>, vector<double> );

vector<double> zeros(vector<double>, vector<double> );

struct profile 
{	vector<double> zz, cc, wx, wy, rh;
};

profile profileReader(const char *);

profile profileReader2(const char *);


struct linray{
	vector<double> xx, yy, zz, ss, cc, rh, vr, tr, om, ja, ps;
    double tv, th, ph, rsd;
    // tv =  1/sqrt( Sx*Sx + Sy*Sy );  i.e. 1/(initial slowness?) = trace velocity?
    // rsd = dist_sr; i.e. distance source receiver?
};

linray linearRay3DStrat( double, double, double, double, double, double, double, profile );

linray linearRay3DStrat( double, double, double, double, double, double, profile );


// DV
linray eigrayReader(const char *filename);


class odesys
{
    double cc, wx, wy, dcdz, dwxdz, dwydz, d2cdz, d2wxdz, d2wydz, om, QQ1, QQ1rt, QQ2, dom;
    double ff1, ff2, QQ4, QQ5, p11, p12,p21, p22, p31, p32, dxds, dyds, dzds, dSzds, dSzdSxds;
    double Sx, Sy, Sz, dSzdSyds, dxdSxds, dxdSyds, dydSxds, dydSyds, dzdSxds, dzdSyds;
    
public:
    void vars( double*, double [][2], double, double, double );
    void vars( double*, double, double, double );
    void rk4( double, double, double, double, double, double* );
};

linray eigenray( double, double, double, double, double, double, double, double, double, double, profile );

#endif
