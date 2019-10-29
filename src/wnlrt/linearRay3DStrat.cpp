//
//  linearRay3DStrat.cpp
//  
//
//  Created by Joel Lonzaga on 5/27/12.
//  Copyright (c) 2012 University of Mississippi. All rights reserved.
//


#include "linearRay3DStrat.h"

//----------------------------------------------------------------------------------------------
// LINEAR RAY 3D PROPAGATION IN A STRATIFIED MEDIUM
//----------------------------------------------------------------------------------------------
linray linearRay3DStrat( double ds, double theta, double phi, double zs, double rr, double zr, profile prf )
{	
    linray R;
	  double zzprev, range, dz, dxds, dyds, dzds, wx, wy, cef, cc, Sx, Sy, Sz, om;
    double dSzdSy, dzdSy, dxdSx, dxdSy, dydSx, dydSy, dSzdSx, dzdSx;
    double d1ccep[2], d1wxep[2], d1wyep[2], A[4][14], ptAtm[3], ptAtmDer[3][2];
    double *spl[3], *splrh;
    vector< vector<double> > atm;
    
    if( theta <= 0.0)
    {   cerr << "Error: The subroutine \"linearRay3DStrat\" cannot handle a ray with initial" << endl;
        cerr << "       inclination angle of less than or equal to 0 deg." << endl;
        exit(1);
    }
    
    dz = prf.zz[1]-prf.zz[0];
    deriv1EndPoints( dz, prf.cc, d1ccep );
    deriv1EndPoints( dz, prf.wx, d1wxep );
    deriv1EndPoints( dz, prf.wy, d1wyep );
    
    spl[0] = spline( prf.zz, prf.cc, d1ccep[0], d1ccep[1] );
    spl[1] = spline( prf.zz, prf.wx, d1wxep[0], d1wxep[1] );
    spl[2] = spline( prf.zz, prf.wy, d1wyep[0], d1wyep[1] );
    splrh = spline( prf.zz, prf.rh, 1.0e50, 1.0e50 );
    
    //Putting the vectors of sound speed, etc into a single matrix
	atm.push_back( prf.cc );	
	atm.push_back( prf.wx );
	atm.push_back( prf.wy ); 	

  //Values at source location (xs = ys = 0; zs)
  splint( prf.zz, atm, spl, zs, ptAtm );     //interpolating atmospheric profile
  
  cc = ptAtm[0];
  wx = ptAtm[1];
  wy = ptAtm[2];
	cef = cc + wx*cos(theta)*cos(phi) + wy*cos(theta)*sin(phi);
	Sx = cos(theta)*cos(phi)/cef;
	Sy = cos(theta)*sin(phi)/cef;
	Sz = sin(theta)/cef;
	om = 1-wx*Sx-wy*Sy;
  dSzdSx = -( wx*om/pow(cc,2) + Sx )/ Sz;
	dzdSx = 0;
	dSzdSy = -( wy*om/pow(cc,2) + Sy )/ Sz;
	dzdSy = 0;
	dxdSx = 0;
	dxdSy = 0;
	dydSx = 0;
	dydSy = 0;
	
  R.cc.push_back( cc );
  R.ss.push_back( 0 );
  R.xx.push_back( 0 );
  R.yy.push_back( 0 );
  R.zz.push_back( zs );
  R.rh.push_back( splint( prf.zz, prf.rh, splrh, zs ) );
  R.om.push_back( om );
  R.vr.push_back( sqrt(  ( pow(wx,2) + pow(wy,2) - pow(cc,2) )*pow(om,2) + 2*pow(cc,2)*om )/om );
  R.ja.push_back( 0 );
  R.tr.push_back( 0 );
  R.tv =  1/sqrt( Sx*Sx + Sy*Sy );
  R.th = theta*180/PI;
  R.ph = phi*180/PI;
    
    
    // Marching along the ray path and solving the system of ode using the 4th order Runge-Kutta routine
    
    splintDeriv( prf.zz, atm, spl, R.zz.back(), ptAtmDer );
    
    
    
    odesys ray;
    ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz );
    ray.rk4( ds, dzdSx, dzdSy, dSzdSx, dSzdSy, A[0] );

    zzprev = zs;
    range = 0;
    int k;
    
    while( !( R.zz.back() <= zr && zzprev > zr ) && range <= rr+10000  )
    {   
        for( k=1; k<4; k++)
        {   
            if(k==3) 
            {   zs = R.zz.back()+A[2][2];
                if( zs < prf.zz[0] || zs > prf.zz.back() ) break;
                splint( prf.zz, atm, spl, zs, ptAtm );
                splintDeriv( prf.zz, atm, spl, zs, ptAtmDer );
                ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz+A[2][3] );
                ray.rk4( ds, dzdSx+A[2][8], dzdSy+A[2][9], dSzdSx+A[2][10], dSzdSy+A[2][11], A[k] );
            }
            else
            {   zs = R.zz.back()+A[k-1][2]/2;
                if( zs < prf.zz[0] || zs > prf.zz.back() ) break;
                splint( prf.zz, atm, spl, zs, ptAtm );
                splintDeriv( prf.zz, atm, spl, zs, ptAtmDer );
                ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz+A[k-1][3]/2 );
                ray.rk4( ds, dzdSx+A[k-1][8]/2, dzdSy+A[k-1][9]/2, dSzdSx+A[k-1][10]/2, dSzdSy+A[k-1][11]/2, A[k] );
            }
        }
        
        if( zs < prf.zz[0] || zs > prf.zz.back() ) break;
        zzprev = R.zz.back();
        
        Sz =         Sz + 1.0/6*( A[0][3]  +  2*(A[1][3]+A[2][3])  + A[3][3]  );
        dxdSx =   dxdSx + 1.0/6*( A[0][4]  +  2*(A[1][4]+A[2][4])  + A[3][4]  );
        dxdSy =   dxdSy + 1.0/6*( A[0][5]  +  2*(A[1][5]+A[2][5])  + A[3][5]  );
        dydSx =   dydSx + 1.0/6*( A[0][6]  +  2*(A[1][6]+A[2][6])  + A[3][6]  );
        dydSy =   dydSy + 1.0/6*( A[0][7]  +  2*(A[1][7]+A[2][7])  + A[3][7]  );
        
        dzdSx =   dzdSx + 1.0/6*( A[0][8]  +  2*(A[1][8]+A[2][8])  + A[3][8]  );
        dzdSy =   dzdSy + 1.0/6*( A[0][9]  +  2*(A[1][9]+A[2][9])  + A[3][9]  );
        dSzdSx = dSzdSx + 1.0/6*( A[0][10] + 2*(A[1][10]+A[2][10]) + A[3][10]  );
        dSzdSy = dSzdSy + 1.0/6*( A[0][11] + 2*(A[1][11]+A[2][11]) + A[3][11]  );
        
        R.ss.push_back(  R.ss.back()+ds );
        R.xx.push_back(  R.xx.back() + 1.0/6*(A[0][0] + 2*(A[1][0]+A[2][0]) + A[3][0])   );
        R.yy.push_back(  R.yy.back() + 1.0/6*(A[0][1] + 2*(A[1][1]+A[2][1]) + A[3][1])   );
        R.zz.push_back(  R.zz.back() + 1.0/6*(A[0][2] + 2*(A[1][2]+A[2][2]) + A[3][2])   );
        
        splint( prf.zz, atm, spl, R.zz.back(), ptAtm );
        splintDeriv( prf.zz, atm, spl, R.zz.back(), ptAtmDer );
        ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz );
        ray.rk4( ds, dzdSx, dzdSy, dSzdSx, dSzdSy, A[0] );
        
        R.cc.push_back(  ptAtm[0] );
        R.om.push_back( A[0][12] );
        R.vr.push_back( A[0][13] );
        R.tr.push_back( R.tr.back() + ds/A[0][13] );
        R.rh.push_back( splint( prf.zz, prf.rh, splrh, R.zz.back() ) );
        
        dxds = A[0][0]/ds; 
        dyds = A[0][1]/ds;
        dzds = A[0][2]/ds;
        R.ja.push_back( dxds*(dydSx*dzdSy-dydSy*dzdSx) + dxdSx*(dydSy*dzds-dyds*dzdSy) + dxdSy*(dyds*dzdSx-dydSx*dzds) ); // Jacobian
        range = sqrt( pow( R.xx.back(),2 ) + pow( R.yy.back(), 2 ) ); 
    }

    for( int i=0; i<3; i++ )
        delete[] spl[i];
    delete[] splrh;

    return R;
}   

//----------------------------------------------------------------------------------------------
// LINEAR RAY 3D PROPAGATION IN A STRATIFIED MEDIUM
//----------------------------------------------------------------------------------------------
linray linearRay3DStrat( double ds, double theta, double phi, double zs, double xr, double yr, double zr, profile prf )
{	
    linray R;
	  double zzprev, range, rr, dz, dxds, dyds, dzds, wx, wy, cef, cc, Sx, Sy, Sz, om;
    double dSzdSy, dzdSy, dxdSx, dxdSy, dydSx, dydSy, dSzdSx, dzdSx;
    double d1ccep[2], d1wxep[2], d1wyep[2], A[4][14], ptAtm[3], ptAtmDer[3][2];
    double *spl[3], *splrh;
    vector< vector<double> > atm;
    
    if( theta <= 0.0)
    {   cerr << "Error: The subroutine \"linearRay3DStrat\" cannot handle a ray with initial" << endl;
        cerr << "       inclination angle of less than or equal to 0 deg." << endl;
        exit(1);
    }
    
    dz = prf.zz[1]-prf.zz[0];
    deriv1EndPoints( dz, prf.cc, d1ccep );
    deriv1EndPoints( dz, prf.wx, d1wxep );
    deriv1EndPoints( dz, prf.wy, d1wyep );
    
    spl[0] = spline( prf.zz, prf.cc, d1ccep[0], d1ccep[1] );
    spl[1] = spline( prf.zz, prf.wx, d1wxep[0], d1wxep[1] );
    spl[2] = spline( prf.zz, prf.wy, d1wyep[0], d1wyep[1] );
    splrh = spline( prf.zz, prf.rh, 1.0e50, 1.0e50 );
    
    //Putting the vectors of sound speed, etc into a single matrix
	atm.push_back( prf.cc );	
	atm.push_back( prf.wx );
	atm.push_back( prf.wy ); 	

  //Values at source location (xs = ys = 0; zs)
  splint( prf.zz, atm, spl, zs, ptAtm );     //interpolating atmospheric profile
  
  cc = ptAtm[0];
  wx = ptAtm[1];
  wy = ptAtm[2];
	cef = cc + wx*cos(theta)*cos(phi) + wy*cos(theta)*sin(phi);
	Sx = cos(theta)*cos(phi)/cef;
	Sy = cos(theta)*sin(phi)/cef;
	Sz = sin(theta)/cef;
	om = 1-wx*Sx-wy*Sy;
  dSzdSx = -( wx*om/pow(cc,2) + Sx )/ Sz;
	dzdSx = 0;
	dSzdSy = -( wy*om/pow(cc,2) + Sy )/ Sz;
	dzdSy = 0;
	dxdSx = 0;
	dxdSy = 0;
	dydSx = 0;
	dydSy = 0;
	
  R.cc.push_back( cc );
  R.ss.push_back( 0 );
  R.xx.push_back( 0 );
  R.yy.push_back( 0 );
  R.zz.push_back( zs );
  R.rh.push_back( splint( prf.zz, prf.rh, splrh, zs ) );
  R.om.push_back( om );
  R.vr.push_back( sqrt(  ( pow(wx,2) + pow(wy,2) - pow(cc,2) )*pow(om,2) + 2*pow(cc,2)*om )/om );
  R.ja.push_back( 0 );
  R.tr.push_back( 0 );
  R.tv =  1/sqrt( Sx*Sx + Sy*Sy );
  R.th = theta*180/PI;
  R.ph = phi*180/PI;
    
    
    // Marching along the ray path and solving the system of ode using the 4th order Runge-Kutta routine
    
    splintDeriv( prf.zz, atm, spl, R.zz.back(), ptAtmDer );
    
    
    
    odesys ray;
    ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz );
    ray.rk4( ds, dzdSx, dzdSy, dSzdSx, dSzdSy, A[0] );

    zzprev = zs;
    range = rr = sqrt( xr*xr + yr*yr );
    int k;
    
    while( !( R.zz.back() <= zr && zzprev > zr ) && range <= rr+10000  )
    {   
        for( k=1; k<4; k++)
        {   
            if(k==3) 
            {   zs = R.zz.back()+A[2][2];
                if( zs < prf.zz[0] || zs > prf.zz.back() ) break;
                splint( prf.zz, atm, spl, zs, ptAtm );
                splintDeriv( prf.zz, atm, spl, zs, ptAtmDer );
                ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz+A[2][3] );
                ray.rk4( ds, dzdSx+A[2][8], dzdSy+A[2][9], dSzdSx+A[2][10], dSzdSy+A[2][11], A[k] );
            }
            else
            {   zs = R.zz.back()+A[k-1][2]/2;
                if( zs < prf.zz[0] || zs > prf.zz.back() ) break;
                splint( prf.zz, atm, spl, zs, ptAtm );
                splintDeriv( prf.zz, atm, spl, zs, ptAtmDer );
                ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz+A[k-1][3]/2 );
                ray.rk4( ds, dzdSx+A[k-1][8]/2, dzdSy+A[k-1][9]/2, dSzdSx+A[k-1][10]/2, dSzdSy+A[k-1][11]/2, A[k] );
            }
        }
        
        if( zs < prf.zz[0] || zs > prf.zz.back() ) break;
        zzprev = R.zz.back();
        
        Sz =         Sz + 1.0/6*( A[0][3]  +  2*(A[1][3]+A[2][3])  + A[3][3]  );
        dxdSx =   dxdSx + 1.0/6*( A[0][4]  +  2*(A[1][4]+A[2][4])  + A[3][4]  );
        dxdSy =   dxdSy + 1.0/6*( A[0][5]  +  2*(A[1][5]+A[2][5])  + A[3][5]  );
        dydSx =   dydSx + 1.0/6*( A[0][6]  +  2*(A[1][6]+A[2][6])  + A[3][6]  );
        dydSy =   dydSy + 1.0/6*( A[0][7]  +  2*(A[1][7]+A[2][7])  + A[3][7]  );
        
        dzdSx =   dzdSx + 1.0/6*( A[0][8]  +  2*(A[1][8]+A[2][8])  + A[3][8]  );
        dzdSy =   dzdSy + 1.0/6*( A[0][9]  +  2*(A[1][9]+A[2][9])  + A[3][9]  );
        dSzdSx = dSzdSx + 1.0/6*( A[0][10] + 2*(A[1][10]+A[2][10]) + A[3][10]  );
        dSzdSy = dSzdSy + 1.0/6*( A[0][11] + 2*(A[1][11]+A[2][11]) + A[3][11]  );
        
        R.ss.push_back(  R.ss.back()+ds );
        R.xx.push_back(  R.xx.back() + 1.0/6*(A[0][0] + 2*(A[1][0]+A[2][0]) + A[3][0])   );
        R.yy.push_back(  R.yy.back() + 1.0/6*(A[0][1] + 2*(A[1][1]+A[2][1]) + A[3][1])   );
        R.zz.push_back(  R.zz.back() + 1.0/6*(A[0][2] + 2*(A[1][2]+A[2][2]) + A[3][2])   );
        
        splint( prf.zz, atm, spl, R.zz.back(), ptAtm );
        splintDeriv( prf.zz, atm, spl, R.zz.back(), ptAtmDer );
        ray.vars( ptAtm, ptAtmDer, Sx, Sy, Sz );
        ray.rk4( ds, dzdSx, dzdSy, dSzdSx, dSzdSy, A[0] );
        
        R.cc.push_back(  ptAtm[0] );
        R.om.push_back( A[0][12] );
        R.vr.push_back( A[0][13] );
        R.tr.push_back( R.tr.back() + ds/A[0][13] );
        R.rh.push_back( splint( prf.zz, prf.rh, splrh, R.zz.back() ) );
        
        dxds = A[0][0]/ds; 
        dyds = A[0][1]/ds;
        dzds = A[0][2]/ds;
        R.ja.push_back( dxds*(dydSx*dzdSy-dydSy*dzdSx) + dxdSx*(dydSy*dzds-dyds*dzdSy) + dxdSy*(dyds*dzdSx-dydSx*dzds) ); // Jacobian
        range = sqrt( pow( R.xx.back(),2 ) + pow( R.yy.back(), 2 ) ); 
    }

    for( int i=0; i<3; i++ )
        delete[] spl[i];
    delete[] splrh;

    return R;
}   


//----------------------------------------------------------------------------------------------
//System of ODE's for Runge-Kutta Approach
//----------------------------------------------------------------------------------------------    

void odesys::vars( double *atm, double atm1stDer[][2], double SSx, double SSy, double SSz )
{   
    cc = atm[0];
    wx = atm[1];
    wy = atm[2];
    dcdz = atm1stDer[0][0];
    dwxdz = atm1stDer[1][0];
    dwydz = atm1stDer[2][0];
    d2cdz = atm1stDer[0][1];
    d2wxdz = atm1stDer[1][1];
    d2wydz = atm1stDer[2][1];
    Sx = SSx;
    Sy = SSy;
    Sz = SSz;
    
    om = 1-wx*Sx-wy*Sy;    
    QQ1 = ( pow(wx,2)+pow(wy,2)-pow(cc,2) )* pow(om,2) + 2*pow(cc,2)*om;    
    QQ1rt = sqrt(QQ1); 
    QQ2 = ( pow(wx,2) + pow(wy,2) )*om + pow(cc,2)*(1-om);
    dom = Sx*dwxdz + Sy*dwydz;
    ff1 = pow(om,2)/cc*dcdz + om*dom;
    ff2 = 2*om/cc*dcdz + dom;
    QQ4 = (wx*dwxdz + wy*dwydz - cc*dcdz)*pow(om,2) + 2*cc*dcdz*om - QQ2*dom;
    QQ5 = pow(om,2)*( d2cdz/cc - pow(dcdz/cc,2) + 1/om*(Sx*d2wxdz+Sy*d2wydz) ) - ff2*dom;
    
    p11 = dwxdz*om + 2*cc*dcdz*Sx - wx*dom;
    p12 = (wx*om + pow(cc,2)*Sx)*QQ4;
    p21 = dwydz*om + 2*cc*dcdz*Sy - wy*dom;
    p22 = (wy*om + pow(cc,2)*Sy)*QQ4;
    p31 = 2*cc*dcdz*Sz;
    p32 = pow(cc,2)*Sz*QQ4;
}

void odesys::rk4( double ds, double dzdSx, double dzdSy, double dSzdSx, double dSzdSy, double *A )
{   
    dxds = (wx*om + pow(cc,2)*Sx)/QQ1rt;
    dyds = (wy*om + pow(cc,2)*Sy)/QQ1rt;
    dzds = pow(cc,2)*Sz/QQ1rt;
    
    dSzds = -om*( om*dcdz/cc + Sx*dwxdz + Sy*dwydz )/QQ1rt;
    dSzdSxds = -1/pow(QQ1,3.0/2) * ( (QQ5*QQ1-ff1*QQ4)*dzdSx + om*dwxdz*QQ1 + wx*( QQ2*ff1-QQ1*ff2) );
    dSzdSyds = -1/pow(QQ1,3.0/2) * ( (QQ5*QQ1-ff1*QQ4)*dzdSy + om*dwydz*QQ1 + wy*( QQ2*ff1-QQ1*ff2) );
    
    dxdSxds = ( (p11*QQ1-p12)*dzdSx +QQ1*( pow(cc,2) - pow(wx,2) ) + QQ2*wx* (wx*om+pow(cc,2)*Sx) )/pow(QQ1,3.0/2);
    dxdSyds = ( (p11*QQ1-p12)*dzdSy -QQ1*wx*wy + QQ2*wy*(wx*om+pow(cc,2)*Sx) )/pow(QQ1,3.0/2);
    dydSxds = ( (p21*QQ1-p22)*dzdSx -QQ1*wx*wy + QQ2*wx*(wy*om+pow(cc,2)*Sy) )/pow(QQ1,3.0/2);
    dydSyds = ( (p21*QQ1-p22)*dzdSy +QQ1*(pow(cc,2)-pow(wy,2) ) + QQ2*wy*(wy*om+pow(cc,2)*Sy) )/pow(QQ1,3.0/2);
    dzdSxds = ( (p31*QQ1-p32)*dzdSx + QQ1*pow(cc,2)*dSzdSx + QQ2*wx*pow(cc,2)*Sz )/pow(QQ1,3.0/2);
    dzdSyds = ( (p31*QQ1-p32)*dzdSy + QQ1*pow(cc,2)*dSzdSy + QQ2*wy*pow(cc,2)*Sz )/pow(QQ1,3.0/2);
    
    A[0] = ds*dxds;
    A[1] = ds*dyds;
    A[2] = ds*dzds;
    A[3] = ds*dSzds;
    A[4] = ds*dxdSxds;
    A[5] = ds*dxdSyds;
    A[6] = ds*dydSxds;
    A[7] = ds*dydSyds;
    A[8] = ds*dzdSxds;
    A[9] = ds*dzdSyds;
    A[10] = ds*dSzdSxds;
    A[11] = ds*dSzdSyds;
    A[12] = om;
    A[13] = QQ1rt/om;  // v_ray
}


//----------------------------------------------------------------------------------------------
// Reading the atmospheric specifications
//----------------------------------------------------------------------------------------------

// DV
linray eigrayReader(const char *filename)
{
  // read eigenray file outputted from Phil Blom code as adapted by Doru Velea
  linray LR;
  vector<double> Sx, Sy, Sz, wx, wy, wz;
  
  // struct linray  {
	//    vector<double> xx, yy, zz, ss, cc, rh, vr, tr, om, ja, ps;
  //    double tv, th, ph, rsd;
  // };

  
//%%   1       2       3          4                   5                    6             7       8       9           10     11           12           13        14       15   
//# x [km]	y [km]	z [km]	Geo. Atten. [dB]	Atmo. Atten. [dB]	Travel Time [s]	  rho  	   c [km/s]  	   u [km/s]  	   v [km/s]  	   w [km/s]  	Slowness_x [s/km]	Slowness_y [s/km]	Slowness_z [s/km]	Jacobian [km^2/rad^2]

	double dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11, dat12, dat13, dat14, dat15;
	char line[512]; // to store the headerline text
	ifstream indata;
    
	indata.open( filename ); 	// opens the file
	if(!indata)  				// file couldn't be opened
	{ 
		cerr << "Error: File " << filename << " could not be opened!" << endl;
        exit(1);
	}
	indata.getline(line, 512);
	indata.getline(line, 512);
	indata.getline(line, 512);
	indata.getline(line, 512);
	
	indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5 >> dat6 >> dat7 >> dat8 >> dat9 >> dat10 >> dat11 >> dat12 >> dat13 >> dat14 >> dat15;;
	
	//cout << dat1 << "  " << dat2 << endl;
	//printf("dat1=%g\n", dat1);
	//printf("dat2=%g\n", dat2);
	//printf("dat3=%g\n", dat3);		

	int it = 1;
	while ( !indata.eof() ) 		// keep reading until end-of-file
	{
	  // in SI units
	  LR.xx.push_back(dat1*1000.0);
	  LR.yy.push_back(dat2*1000.0);
	  LR.zz.push_back(dat3*1000.0);
	  LR.cc.push_back(dat8*1000.0);
	  LR.rh.push_back(dat7*1000.0);
	  LR.tr.push_back(dat6); // time [sec]
	  LR.ja.push_back(dat15*1.0e6); // Jacobian in m^2/rad^2
	  wx.push_back(dat9*1000.0); // wind - x-direction
	  wy.push_back(dat10*1000.0);
	  wz.push_back(dat11*1000.0);
	  Sx.push_back(dat12*1.0e-3); // slowness components in seconds/meter
	  Sy.push_back(dat13*1.0e-3);
	  Sz.push_back(dat14*1.0e-3);
		indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5 >> dat6 >> dat7 >> dat8 >> dat9 >> dat10 >> dat11 >> dat12 >> dat13 >> dat14 >> dat15;
		it = it + 1;
		//printf("it=%d\n", it);
	}
	indata.close();
	
	//printf("LR.zz[0]=%g\n", LR.zz[0]);
  //printf("LR.cc[0]=%g\n", LR.cc[0]);

	// populate LR.vr, LR.om and LR.ss
	int j, N = Sx.size();
	double Smag, ds, s = 0.0;
	LR.ss.push_back(0.0); // the raypath length vector initial value
	
	j = 0;
	Smag = sqrt(Sx[j]*Sx[j] + Sy[j]*Sy[j] + Sz[j]*Sz[j]);
	//LR.om.push_back(Smag*LR.cc[j]); // initial OMEGA
	LR.om.push_back(1.0 - wx[j]*Sx[j] - wy[j]*Sy[j] - wz[j]*Sz[j]); // initial OMEGA = 1-v0*Slowness
	LR.tv =  1.0/Smag; // initial tv? -  not clear what this is - DV
	
	//printf("Smag=%g\n", Smag);
	//printf("LR.om[0]=%g\n", LR.om[0]);
	//printf("wx[0]=%g\n", wx[0]);
	//printf("wy[0]=%g\n", wy[0]);
	//printf("Sx size N=%d\n", N);
	
	// initial v_ray (here in terms of OMEGA: sqrt(v^2-c^2 +2c^2/om) which should be the same as Joel's
	// R.vr.push_back( sqrt(  ( pow(wx,2) + pow(wy,2) - pow(cc,2) )*pow(om,2) + 2*pow(cc,2)*om )/om );
	//LR.vr.push_back( sqrt(  ( pow(wx[0],2) + pow(wy[0],2) - pow(LR.cc[0],2) )*pow(LR.om[0],2) + 2*pow(LR.cc[0],2)*LR.om[0] )/LR.om[0] ); // Joel's
	LR.vr.push_back( sqrt(wx[0]*wx[0] + wy[0]*wy[0] + wz[0]*wz[0] + LR.cc[0]*LR.cc[0]*(2.0/LR.om[0]-1.0)) );
	
//	vray(j) = sqrt(c(j)^2 + (norm(v0(j,:)))^2 + 2*c(j)*dot(v0(j,:),S(j,:)/Smag(j)));
//	vray(j) = wx[0]*wx[0] + wy[0]*wy[0] + LR.cc[0]*LR.cc[0] + 2*c(j)*dot(v0(j,:),S(j,:)/Smag(j)));
	
	for (int j=1; j<N; j++) {
	  ds = sqrt( (LR.xx[j]-LR.xx[j-1])*(LR.xx[j]-LR.xx[j-1]) + \
	            (LR.yy[j]-LR.yy[j-1])*(LR.yy[j]-LR.yy[j-1]) + \
	            (LR.zz[j]-LR.zz[j-1])*(LR.zz[j]-LR.zz[j-1]) );
	            
	  s = s + ds;
	  
	  LR.ss.push_back(s); // ray path length
	  
	  Smag = sqrt(Sx[j]*Sx[j] + Sx[j]*Sx[j] + Sx[j]*Sx[j]);
	  //LR.om.push_back(Smag*LR.cc[j]);
	  LR.om.push_back(1.0 - wx[j]*Sx[j] - wy[j]*Sy[j] - wz[j]*Sz[j]); // OMEGA
	  LR.vr.push_back( sqrt(wx[j]*wx[j] + wy[j]*wy[j] + wz[j]*wz[j] + LR.cc[j]*LR.cc[j]*(2.0/LR.om[j]-1.0)) );
	  // alternate vr would be to compute ds/dt; have to check if they agree !!!
	  // LR.vr.push_back(ds/(LR.tr[j]-LR.tr[j-1]));
	  
	}
	
	//printf("LR.ss[1]=%g\n", LR.ss[1]);
	//printf("LR.vr[0]=%g\n", LR.vr[0]);
  //printf("LR.vr[1]=%g\n", LR.vr[1]);
	
	LR.rsd = 0.0; // this is the distance between the ray last point and receiver 
	              // set it to zero since this is a true eigenray
	              
	              
	              
  // LR.ps is not set yet - not clear is needed	to be set here           
	              
	return LR;
};



 profile profileReader3(const char *filename)
{
	profile prf;
	double dat1, dat2, dat3, dat4, dat5;
	ifstream indata;
    
	indata.open( filename ); 	// opens the file
	if(!indata)  				// file couldn't be opened
	{ 
		cerr << "Error: File " << filename << " could not be opened!" << endl;
        exit(1);
	}
	indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5;
	while ( !indata.eof() ) 		// keep reading until end-of-file
	{
	  prf.zz.push_back(dat1);
		prf.cc.push_back(dat2);
		prf.wx.push_back(dat3);
		prf.wy.push_back(dat4);
		prf.rh.push_back(dat5);
		indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5;
	}
	indata.close();
	return prf;
};

profile profileReader2(const char *filename)
{
	profile prf;
	double dat1, dat2, dat3, dat4, dat5, dat6;
	ifstream indata;
    
	indata.open( filename ); 	// opens the file
	if(!indata)  				// file couldn't be opened
	{ 
		cerr << "Error: File " << filename << "could not be opened!" << endl;
        exit(1);
	}
    
    // G2S style
	indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5 >> dat6;
    dat5 = sqrt(0.14*dat5/dat6);
	while ( !indata.eof() ) 		// keep reading until end-of-file
	{
	  	prf.zz.push_back( 1e3*dat1);
		prf.cc.push_back( dat5 );
		prf.wx.push_back( 1e3*dat3 );
		prf.wy.push_back( 1e3*dat4 );
		prf.rh.push_back( 1e3*dat6 );
        indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5 >> dat6;
        dat5 = sqrt(0.14*dat5/dat6);
	}
	indata.close();
	return prf;
};


// read G2S style atmo file: zuvwtdp; z [km] | d [g/cm^3] | P [hectoPa]|
profile profileReader(const char *filename)
{
	profile prf;
	double dat1, dat2, dat3, dat4, dat5, dat6, dat7;
  double tmp1;
	ifstream indata;
    
	indata.open( filename ); 	// opens the file
	if(!indata)  				// file couldn't be opened
	{ 
		cerr << "Error: File " << filename << " could not be opened!" << endl;
        exit(1);
	}
    
    // G2S style
	indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5 >> dat6 >> dat7;
    tmp1 = sqrt(0.14*dat7/dat6);
    //tmp1 = sqrt(1.4*287.058*dat5);
    dat6 = dat6*1000;
    dat1 = dat1*1000;
	while ( !indata.eof() ) 		// keep reading until end-of-file
	{
	  prf.zz.push_back(dat1);
		prf.cc.push_back(tmp1);
		prf.wx.push_back(dat2);
		prf.wy.push_back(dat3);
		prf.rh.push_back(dat6);
    indata >> dat1 >> dat2 >> dat3 >> dat4 >> dat5 >> dat6 >> dat7;
    tmp1 = sqrt(0.14*dat7/dat6);
    //tmp1 = sqrt(1.4*287.058*dat5);
    dat6 = dat6*1000;
    dat1 = dat1*1000;
    //printf("dat1=%g\n", dat1);
	}
	indata.close();
	return prf;
};


//----------------------------------------------------------------------------------------------
// Numerical Derivative at the end points
//----------------------------------------------------------------------------------------------
void deriv1EndPoints( double dx, vector<double> yy, double *dydx )
{
    int N = yy.size();
    //dydx[0] = ( -11*yy[0] + 18*yy[1] - 9*yy[2] + 2*yy[3] )/(6*dx);
    dydx[0] = ( -yy[2] + 4*yy[1] - 3*yy[0] )/(2*dx);
    //dydx[1] = ( 11*yy[N-1] - 18*yy[N-2] + 9*yy[N-3] - 2*yy[N-4] )/(6*dx);
    dydx[1] = ( 3*yy[N-1] - 4*yy[N-2] + yy[N-3] )/(2*dx);
}


//----------------------------------------------------------------------------------------------
// Spline Interpolation
//----------------------------------------------------------------------------------------------

double* spline( vector<double> x, vector<double> y, double y1p, double ynp )
{
    int N = x.size()-1;
    double h;
    double alph[N+1], ell[N+1], zed[N+1], mu[N+1];
    double *coef = new double[N+1];
    
    h = x[1]-x[0];
    alph[0] = 3*( y[1]-y[0] )/h - 3*y1p;
    alph[N] = 3*ynp - 3*(y[N]-y[N-1])/h;
    
    if( y1p>0.99e50 )                               
    {   ell[0] = 1;
        mu[0] = zed[0] = 0;
    }
    else                                      
    {   ell[0] = 2*h;
        mu[0] = 0.5;
        zed[0] = alph[0]/ell[0];
    }
    
    for( int i=1; i<N; i++ )
    {
        alph[i] = 3/h*( y[i+1] - 2*y[i] + y[i-1] );
        ell[i] = 4*h - h*mu[i-1];
        mu[i] = h/ell[i];
        zed[i] = ( alph[i] - h*zed[i-1] )/ell[i];
    }
    
    if( y1p>0.99e50 )
    {   ell[N] = 1;
        zed[N] = 0;
    }
    else
    {   ell[N] = h*(2-mu[N-1]);
        zed[N] = ( alph[N] - h*zed[N-1] )/ell[N];
    }
    
    coef[N] = zed[N];
    for( int j = N-1; j>=0; j--)
        coef[j] = zed[j] - mu[j]*coef[j+1];
    
    return coef;
}

double splint( vector<double> xx, vector<double> yy, double *coef, double x)
{
    int k, klo, khi, n;
    double h, b, d;
    
    n = xx.size();
    klo=0;
    khi=n-1;
    
    if( x < xx.front() || x > xx.back() )
	{	cerr << "Error: (in splint(... double x) ... Interpolating point is outside of interval." << endl;
        printf("ss(i) is %.2f while ss(end) is %.2f\n", x, xx.back() );
        exit(1);
	}
    
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xx[k] > x) khi=k;
        else klo=k;
    }    
    h=xx[khi]-xx[klo];
    
    if (h == 0.0) 
    {   cerr << "Bad independent variable input to routine splint." << endl;
        exit(1);
    }
    
    b = ( yy[khi] - yy[klo] )/h - h*( coef[khi] +2*coef[klo] )/3;
    d = ( coef[khi] - coef[klo] )/( 3*h);
    return yy[klo] + b*(x-xx[klo]) + coef[klo]*pow(x-xx[klo],2) + d*pow(x-xx[klo],3);
}


void splint( vector<double> xx, vector< vector<double> > yy, double **coef, double x, double *spl)
{
    int k, klo, khi, n, M;
    double h, b, d;
    
    n = xx.size();
    M = yy.size();

    klo=0;
    khi=n-1;
    
    if( x < xx.front() || x > xx.back() )
	{	cerr << "Error: splint(... double *spl) ... Interpolating point is outside of interval." << endl;
        printf("ss(i) is %.2f while ss(end) is %.2f\n", x, xx.back() );
        exit(1);
	}

    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xx[k] > x) khi=k;
        else klo=k;
    }    
    h=xx[khi]-xx[klo];
    
    if (h == 0.0) 
    {   cerr << "Bad independent variable input to routine splint." << endl;
        exit(1);
    }
    
    for( int i=0; i<M; i++ )
    {   b = ( yy[i][khi] - yy[i][klo] )/h - h*( coef[i][khi] +2*coef[i][klo] )/3;
        d = ( coef[i][khi] - coef[i][klo] )/( 3*h);
        spl[i] = yy[i][klo] + b*(x-xx[klo]) + coef[i][klo]*pow(x-xx[klo],2) + d*pow(x-xx[klo],3);
    }
}



void splintDeriv( vector<double> xx, vector< vector<double> > yy, double **coef, double x, double deriv[][2] )
{
    int k, klo, khi, M;
    double h, b, d;
    
    khi = xx.size()-1;
    klo = 0;
    
    if( x < xx.front() || x > xx.back() )
	{	cerr << "Error: Interpolating point is outside of interval." << endl;
        printf("ss(i) is %.2f while ss(end) is %.2f\n", x, xx.back() );
        exit(1);
	}
    
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xx[k] > x) khi=k;
        else klo=k;
    }    
    h=xx[khi]-xx[klo];
    
    if (h == 0.0) 
    {   cerr << "Bad independent variable input to routine splintDeriv." << endl;
        exit(1);
    }

    M = yy.size();

    for( int i=0; i<M; i++ )
    {   b = ( yy[i][khi] - yy[i][klo] )/h - h*( coef[i][khi] +2*coef[i][klo] )/3;
        d = ( coef[i][khi] - coef[i][klo] )/( 3*h);
        deriv[i][0] = b + 2*coef[i][klo]*(x-xx[klo]) + 3*d*pow(x-xx[klo],2);    //first derivative
        deriv[i][1] = 2*coef[i][klo] + 6*d*(x-xx[klo]);                         //second derivative
    }
}

// searching for an eigenray
linray eigenray( double ds, double theta, double phi, double dth, double dph, double zs, double xr, double yr, double zr, double tol, profile prf )
{
    linray eiray,ray;
    double thMid, phMid, dist_sr, dist2, th[3], ph[3];
    int solfnd, repmx, rep;

    rep = 1;
    repmx = 20;
    solfnd = 0;
    dist_sr = pow(10.0,5);
    
    thMid = theta;
    phMid = phi;
    
    while( solfnd == 0 )
    {   for(int i=0; i<3; i++ )
        {   th[i] = thMid + (i-1)*dth;
            ph[i] = phMid + (i-1)*dph;
            if(th[i] <= 0) th[i] = 1e-10;
        }
    
        for(int i=0; i<3; i++)
        {   for( int j=0; j<3; j++ )
            {   
                //printf("%.3f %.3f\n", 180/PI*th[j], 180/PI*ph[i]) ;
                ray = linearRay3DStrat( ds, th[j], ph[i], zs, xr, yr, zr, prf ); 
                dist2 = sqrt( pow(ray.xx.back()-xr,2) + pow(ray.yy.back()-yr,2) + pow(ray.zz.back()-zr,2) );
                if (dist2 < dist_sr)
                {   dist_sr = dist2;
                    eiray = ray;
                    thMid = th[j];
                    phMid = ph[i];
                }
            }
        }
    
        dph = dph/2;
        dth = dth/2;

        if(dist_sr <= tol)
            solfnd = 1;
        rep++;
        if(rep == repmx)
        {   solfnd = 2;
            printf("Covergence to tol=%.0f failed. Closest distance is %.1f m \n", tol, dist_sr);
        }
    }
    eiray.rsd = dist_sr;
    return eiray;
}


// -------------------------------------------------------------------------------
// Finding zeros of a set of data points
// -------------------------------------------------------------------------------
template <typename T> int sgn(T val) 
{
    return (T(0) < val) - (val < T(0));
}


vector<double> zeros(vector<double> xx, vector<double> yy )
{   double tol, xxl, xxu, yyl, yyu, xxa, yytemp, *spl;
    vector<double> sc;
    
    tol = 1;
    spl = spline( xx, yy, 1e50, 1e50 );
    
    for(unsigned int i=1; i<xx.size()-1; i++)
    {   if( sgn(yy[i+1])-sgn(yy[i])!=0 )
    {   xxl = xx[i];
        xxu = xx[i+1];
        yyl = yy[i];
        yyu = yy[i+1];
        
        if( fabs(yyl)<tol )  sc.push_back(xxl);
        else if( fabs(yyu) < tol ) sc.push_back(xxu);
        else
        {   
            xxa = (xxl+xxu)/2;
            yytemp = splint( xx, yy, spl, xxa);
            while( tol < xxu-xxl )
            {   if( sgn(yytemp)-sgn(yyl)==0 )    
            {   xxl = xxa;
                yyl = yytemp;
            }
            else  
            {   xxu = xxa;
                yyu = yytemp;
            }
                xxa = (xxl+xxu)/2;
                yytemp = splint( xx, yy, spl, xxa);
            }
            sc.push_back( xxa );
        }
    }
    }
    delete [] spl;
    return sc;
}

