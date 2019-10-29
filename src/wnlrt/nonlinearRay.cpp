//
//  nonlinearRay.cpp
//  
//
//  Created by Joel Lonzaga on 5/27/12.
//  Copyright (c) 2012 University of Mississippi. All rights reserved.
//

// DV made some changes, formatting, some comments
// See the article: "Modelling waveforms of infrasound arrivals from impulsive 
// sources using weakly non-linear ray theory" 
// by Joel Lonzaga, Roger Waxler, Jelle Assink, Carrick Talmadge
// in Geophysical Journal International (2015), JGI 200, pg 1347-1361

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "linearRay3DStrat.h"
#include "nonlinearRay.h"

//#define PI 3.141592653589793;

const double GAMMA    = 1.40;
const double STEP_MAX = 100.0;   //maximum step size allowed
const int NUMPTS      = pow(2.0,12);  // number of poins in the time/frequency domains to save


//==============================================================================================
//Main Function
//==============================================================================================
nonray nonlinearRay( int nn, double ssInit, attRed attRed, linray LR, waveform wf )
{
  nonray NR;
  int    skip, nsc;
  double dop, scbt, ffmx, pmax, umax, ss, ds, ssMem, atfc, dzta;
  double scbc, tt[nn], pp[nn], uu[nn], attn[nn/2+1], ff[nn/2+1];
  vector<double> sc, uhold;

  //------------------------------------------------------------------------
  // Find caustics if they exist
  sc  = zeros( LR.ss, LR.ja );
  nsc = sc.size();
  double ssup[nsc], sslo[nsc], slope[nsc];
  causticParams( sc, LR, ssup, sslo, slope ); // find J'=slope - at points where J->0

  //------------------------------------------------------------------------
  //Raypath Parameters, Waveform, and Spectrum at 1 km from the source
  //pmax = wfResize( nn, wf, tt, pp ); // sets interpolated tt, pp from wf.tt and wf.pp
  //cout << "pmax = " << pmax << endl;
  
  if (0) {
      cout << "in nonlinearRay() after wfResize: saving pp\n";
      FILE *fid;
      fid = fopen( "pp.dat", "w");
      for (int j=0; j<nn; j++) {
        fprintf(fid, "%15.5E\n", pp[j]);
      }
      
      fclose(fid);
  }
  
  if (1) { // "in nonlinearRay() after undoing wfResize: saving pp\n";
      //cout << "in nonlinearRay() after undoing wfResize: saving pp\n";
      //FILE *fid;
      //fid = fopen( "pp.dat", "w");
      pmax = wf.pp[0];
      for (int j=0; j<nn; j++) {      
          pp[j] = wf.pp[j];
          tt[j] = wf.tt[j];
          if (pmax<pp[j]) {
            pmax=pp[j];
            //fprintf(fid, "%15.5E\n", pp[j]);
          }
      } 
      //fclose(fid);
  }

  ffmx = 1/( 2.0*(tt[1]-tt[0]) );		   // Nyquist frequency
  //printf("Nyquist: ffmx=%g\n", ffmx);
  //printf("number of caustics=%d\n", nsc);

  // construct an instance RP of class raypathParams
  raypathParams RP( LR ); // sets beta = 1+(GAMMA-1)/2, spline accelerators, copies LR.xx,yy, ... to RP
  RP.normalizingFactor( pmax, ssInit); // calls physicalParams(ssinit) and 
                                       // interpolates xx, om, vr... at ssInit; 
                                       // calculates initial normalizing factor NN
                                       // NN=rh0_i*om_i*c0_i^2/(pmax^2*J_i*vr_i)
                                       // 
                                       // sets ps in class RP; 
                                       // below in initFFT the scaled pressure is u = pp*ps
                                       
  //cout << "pmax 2= " << pmax << endl;
                     
  RP.psDopScaBeta( ssInit, dop, scbt ); // scbt= beta tilde = scaled beta at distance ssinit:
                                        // scbt = beta*sqrt( pow(om,3)/( NN*rh*cc*cc*fabs(ja)*pow(vr,3) ) );
                                        // see eq 4 in JGI (200) pg 1347
                                        // sets ps= sqrt( NN*fabs(ja)*vr/( rh*om*pow(cc,2) ) )
                                        // i.e. ps is a scaling factor for pressure in class raypathParams
                                        // computes dop  = cc*om/vr;

  spectralParams SP( nn );     // allocate spectral quantitities U, etc
  SP.initFFT( RP.ps, pp, uu ); // get initial scaled acoustic pressure uu=ps*pp,  
                               // and also (U1, Usq as part of class spectralParams)
  SP.frequency( ffmx, ff );    // since df = 2*ffmx/nn ff is the vector of 
                               // (nn/2+1) positive frequencies 0:df:Nyquist
                               // cf and fltr are set here as part of spectralParams class data
  
  // accumulate all inititial quantities in NR
  NR.xx.push_back( RP.xx );
  NR.yy.push_back( RP.yy );
  NR.zz.push_back( RP.zz );
  NR.cc.push_back( RP.cc );
  NR.tr.push_back( RP.tr );
  NR.ps.push_back( RP.ps );
  NR.sb.push_back( scbt*1 ); // initial scaled beta i.e. beta tilde in eq 4 
  NR.ss.push_back( ssInit ); // DV
  
  //cout << "RP.ps = " << RP.ps << endl;
  
  // store every 'skip' point for a total of NUMPTS
  skip = nn/NUMPTS;
  for(int i=0; i < nn; i+=skip) // store subsampled uu and tt in NR.uu, NR.tt
  {	
    uhold.push_back( uu[i] );
    NR.tt.push_back( tt[i] );
  }
  NR.uu.push_back( uhold ); 
  
//------------------------------------------------------------------------
  //Propagate by 1 step, integrate using the first-order Euler method
  ds = 1.0/( 2.0*PI*ffmx*scbt ); // initial ds
  if( ds > STEP_MAX )	ds = STEP_MAX;
  ssMem = ss = ssInit;
  ss    = ss + ds;
  RP.psDopScaBeta( ss, dop, scbt); // scaled beta at new ss

  if( RP.zz < attRed.zzRed ) atfc = 1; // setting the attenuation adjusting factor above zzRed height
  else atfc = attRed.atfac;	
  atmAttn( nn/2+1, RP.zz*1e-3, RP.rh, RP.cc, ff, attn );
  SP.integEulerMethod( ds, atfc, scbt, dop, attn ); // U2 = (U1 + ds*cf*Usq)*exp(-a*c0*OM/vr); U3=U2
  
  
  // save U2
  if (0) {
      cout << "in nonlinearRay() after integEulerMethod(): saving U1 U2\n";
      FILE *fid;
      fid = fopen( "U2.dat", "w");
      for (int j=0; j<nn/2+1; j++) {
        fprintf(fid, "%15.5E  %15.5E\n", SP.U2[j][0], SP.U2[j][1]);
      } 
      fclose(fid);
      
    //U2[i][1] = U1[i][1] + ds*cf[i]*scbt*Usq[i][0];  // imaginary part
		//attnEf2[i] = atfc * attn[i]*fltr[i]*dop;        // alfa*c0*Omega/v_ray; alfa is altered by a factor atfc
		//U2[i][0] = U2[i][0]*exp( - attnEf2[i]*ds  );    // apply attenuation; see paper page 1351
		//U2[i][1] = U2[i][1]*exp( - attnEf2[i]*ds  );
  }
  
  
  SP.FFT( uu, umax );  // obtain uu = ifft(U3) and Usq=fft(u^2); also umax is returned
  //------------------------------------------------------------------------
  //Marching along the ray path
  int j  = 0;
  int it = 1;
  
  
  FILE *fidU2;
  bool save_U2 = false;
  if (save_U2) {
    cout << "in nonlinearRay() after integEulerMethod(): saving  U2a\n";
    fidU2 = fopen( "U2a.dat", "w");
  }
   
  while( ss+ds <= LR.ss.back() )  {	
    // Saving/storing needed parameters every ~0.5 km of ray path		
    if( int(ss)%500 < int(ssMem)%500 )
    {	
        NR.xx.push_back( RP.xx );
        NR.yy.push_back( RP.yy );
        NR.zz.push_back( RP.zz );
        NR.cc.push_back( RP.cc );
        NR.tr.push_back( RP.tr );
        NR.ps.push_back( RP.ps );
        NR.sb.push_back( scbt*umax ); // note scbt*umax
        NR.ss.push_back( ss );        // DV
        uhold.clear();
        for(int i=0; i < nn; i+=skip) 
            uhold.push_back( uu[i] );
        NR.uu.push_back( uhold ); // store entire waveform at this distance
            
      if (save_U2) {
      cout << "in nonlinearRay() after integEulerMethod(): saving  U2\n";

      for (int j=0; j<nn/2+1; j++) {
        fprintf(fidU2, "%15.5E  ", SP.U1[j][0]);
      }
      fprintf(fidU2, "\n");
      } 
    }
            
    //Nonlinear Integration 
    if( ss < sslo[j] || ss > ssup[j] )      // away from caustic
    {	
        ds = 1.0/( 2.0*PI*ffmx * scbt * umax ); // variable ds size 
        if( ds > STEP_MAX )	ds = STEP_MAX;
        SP.integLeapfrog( ds, scbt );  // U3 = U1 + i*ds*w_n*scbt*Usq - see eq after eq. 12 in GJI article
    }
    else  // around caustic; we obtain a transformed Burger's equation free of singularity
    {			
        //printf("around caustic\n");
        if( ss >= sc[j] && ssMem<sc[j] )  // adding a phase shift of +/-PI/2 to pos/neg freq components
        {   
            SP.hilbert();
            SP.FFT( uu, umax );
        }
        scbc = RP.scaledBetaCaustic( slope[j] ); // beta_caustic = beta_scaled * sqrt(J/J'); here slope is |J'|
        ds = stepSizeCaustic( ffmx, scbc, umax, sc[j], ss, dzta );
        SP.integLeapfrog( dzta, scbc ); // leapfrog the transformed Burger's eq. 
        
        // DV
        if ( (j+1)<sc.size() ) { // increment the zero-based caustic counter j 
          j++;
        }
    }

    //Stepping
    ssMem = ss;
    ss    = ss+ds;
    if( ss >= LR.ss.back() ) break;
    
    RP.psDopScaBeta( ss, dop, scbt); // update RP.xx, yy... dop and scbt

    // Attenuation and Fourier transforming
    if( RP.zz < attRed.zzRed ) atfc = 1;
    else atfc = attRed.atfac;
    atmAttn( nn/2+1, RP.zz*1e-3, RP.rh, RP.cc, ff, attn );
    SP.applyAttn( ds, atfc, dop, attn ); // Update U1,U2,U3, and apply fltr
    SP.FFT( uu, umax );  // obtain uu = ifft(U3) and Usq=fft(u^2); also umax is returned

    // Displaying status
    if( it%1000 == 0) {
      printf("--> Currently at altitude = %.2f km; pathlength=%g km; (total is %g km)\n", \
                   RP.zz*1e-3, ss*1e-3, LR.ss.back()*1e-3);
      }
      
    it++;

  } 	//ending marching along ray
  
  
  
  if (save_U2) {
  cout << "closing fidU2\n";
      fclose(fidU2);
  }
    
    
  
  NR.xx.push_back( RP.xx );
  NR.yy.push_back( RP.yy );
  NR.zz.push_back( RP.zz );
  NR.cc.push_back( RP.cc );
  NR.tr.push_back( RP.tr );
  NR.ps.push_back( RP.ps );
  NR.sb.push_back( scbt*umax );
  NR.ss.push_back( ss ); // DV
  
  //printf("NR.xx.back = %g\n", NR.xx.back());
  //printf("NR.ss.back = %g\n", NR.ss.back());
  
  
  // store in NR the last waveform and its spectrum
  uhold.clear();
	for( int i=0; i < nn; i+=skip ) 
	{	
	  uhold.push_back( uu[i] );
  }
	NR.uu.push_back( uhold );  // store this last time step (uhold) into matrix NR.uu
	
  for( int i=0; i<nn/2+1; i++ )
  {   
      NR.Ur.push_back( SP.U2[i][0] );
      NR.Ui.push_back( SP.U2[i][1] );
      NR.ff.push_back( ff[i] );
  }
    
  cout << "--> Finishing the Nonlinear Ray Acoustic Calculation." << endl;
  
	return NR;

} // end of nonlinearRay




//---------------------------------------------------------
//Save waveform to disk
//---------------------------------------------------------
void savewf( const char* fname, waveform wf )
{
    FILE *wave;
    wave = fopen( fname, "w");
    int N = wf.tt.size();
    for( int i=0; i<N; i++)
        fprintf( wave, "%.3f %10.3f\n", wf.tt[i], wf.pp[i] );
    fclose( wave );
}


// ----------------------------------------------------------------------------
//Generate N-Waveform
// ----------------------------------------------------------------------------
waveform generateNWave( int nn, double period, double taud, double amp )
{   // T=period; sampling interval  = dt = T/nn; tt is from -T/2 to T/2 and
    // N wave extends from [-taud:taud]
    waveform wf;
    int j = -nn/2;
    for( int i=0; i<nn; i++ )
    {   
        wf.tt.push_back( period*j/nn );
        if( abs( wf.tt[i] ) <= taud ) 
            wf.pp.push_back( -amp*wf.tt[i]/taud ); // [-taud:taud]
        else 
            wf.pp.push_back( 0 );
        j++;
    }
    return wf;
}


// ----------------------------------------------------------------------------
//Generate Regular Smooth Waveform
// ----------------------------------------------------------------------------
waveform generateWaveForm( int nn, double period, double taud, double &amp )
{   
    waveform wf;
    double ppmx;  // CHH 191029: fac Unused
    int j = -nn/2;
    ppmx = 0;
    for( int i=0; i<nn; i++ )
    {   
        wf.tt.push_back( period*j/nn );
        wf.pp.push_back( -sin( wf.tt[i]) * exp( -1.50*pow(wf.tt[i],2) ) );
        
        if( wf.pp[i] > ppmx ) ppmx = wf.pp[i];
        j++;
    }
    
    for( int i=0; i<nn; i++ )
    {
        wf.pp[i] = wf.pp[i]*amp/ppmx;
        wf.tt[i] = taud * wf.tt[i]/0.5493;  // to "normalize" the positive phase duration
    }
    return wf;
}

waveform generateWaveForm( int nn, double period, double &amp )
{   
    waveform wf;
    double ppmx, fac;
    int j = -nn/2;
    ppmx = 0;
    fac = 1.00389070263;
    for( int i=0; i<nn; i++ )
    {   
        wf.tt.push_back( period*j/nn );
        wf.pp.push_back( -3 * amp *fac * sin( wf.tt[i]) * exp( -1.50*pow(wf.tt[i],2) ) );
        
        if( wf.pp[i] > ppmx ) ppmx = wf.pp[i];
        j++;
    }
    amp = ppmx;
    return wf;
}


// -------------------------------------------------------------------------------
//Generate Blast Waveform
// -------------------------------------------------------------------------------
waveform blastmodel( int nn, double period, double posPhaseDur, double pkOverpress, linray LR )
{
    waveform wf;
    //double Po = LR.rh[0] * pow( LR.cc[0],2 )/1.4;  // CHH 191029: Unused
    int j = -nn/2;
    for( int i=0; i<nn; i++ )
    {
        wf.tt.push_back( period*j/nn );
        if( wf.tt[i] >= 0 ) 
            wf.pp.push_back( pkOverpress* ( 1-wf.tt[i]/posPhaseDur )*exp( -0.9*wf.tt[i]/posPhaseDur ) );
        else 
            wf.pp.push_back( 0 );
        wf.tt[i] = wf.tt[i]-posPhaseDur;
        j++;
    }
    return wf;
}


//----------------------------------------------------------------------------------------------
//Loading Waveform
//----------------------------------------------------------------------------------------------
waveform loadWaveform( const char *wfFile, double period, double amp )
{
    double dat1, dat2, pmx, pmn, tmx, tmn, tlf, trg, tav;
    vector<double> tt, pp;
    waveform wf;
    ifstream indata;
    indata.open( wfFile ); 	
    if( !indata )  				
    { 
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    trg = period/2;
    tlf = -trg;
    pmx = 0;
    pmn = 0;
    indata >> dat1 >> dat2;
    while ( !indata.eof() ) 		
    {  	tt.push_back( dat1 );
        pp.push_back( dat2 );
        if( pmx<dat2 )
        {   pmx = dat2;
            tmx = dat1;
        }
        if( pmn>dat2 )
        {   pmn = dat2;
            tmn = dat1;
        }
        indata >> dat1 >> dat2;
    }	
    tav = ( tmn+tmx )/2;
    if( abs(pmn) > abs(pmx) ) pmx = pmn;
    for( unsigned int i=0; i<tt.size(); i++ )
    {   dat1 = tt[i]-tav;
        if( dat1 >= tlf ) 
        {   wf.tt.push_back( dat1 );
            wf.pp.push_back( amp/pmx*pp[i] );
        }
        if( dat1 > trg ) break;
    }
    return wf;
}


//----------------------------------------------------------------------------------------------
//Loading BlastModel
//----------------------------------------------------------------------------------------------
waveform loadBlastModel( const char *wfFile, int nn, double period, double posPhasDur, double amp )
{
    double dat1, dat2;
    vector<double> tt, pp;
    waveform wf;
    ifstream indata;
    indata.open( wfFile ); 	
    if( !indata )  				
    { 
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    
    indata >> dat1 >> dat2;
    while ( !indata.eof() ) 		
    {  	tt.push_back( posPhasDur*dat1 );
        pp.push_back( amp*dat2 );
        indata >> dat1 >> dat2;
    }	
    double *spl = spline( tt, pp, 1e50, 1e50 );
    int j = -nn/2;
    for( int i=0; i<nn; i++ )
    {
        wf.tt.push_back( period*j/nn );
        wf.pp.push_back( splint( tt, pp, spl, wf.tt[i] ) );
        j++;
    }
    delete [] spl;
    return wf;
}


// -------------------------------------------------------------------------------
//Calculate peak overpressure and positive phase duration using the Kinsler and Kinney blast models
// -------------------------------------------------------------------------------
// Kinsler Model
blastParam pkOverpress( double distance, double yield, linray LR )
{
    blastParam bp;
    double R, fact1, fact2, fact3, fact4;    
    R = distance/ pow( 1.8*yield,1.0/3 );
    
    fact1 = pow( R/0.915, 25.0/4 );
    fact2 = pow( R/8.54, 21.0/8 );
    fact3 = sqrt( 1 + fact1 );
    fact4 = pow( 1+fact2, 1.0/3 );
    bp.ttd = 0.85*pow( R,4 )/( fact3*fact4 );
    bp.ttd = 1e-3 * bp.ttd * pow( 1.8*yield, 1.0/3 );
    
    fact1 = pow( R/7.2, 3 );
    fact2 = pow( 1 + fact1, 5.0/12 );
    fact3 = pow( R, 9.0/4 );
    bp.pks = 1.0/1.4*LR.rh[0]* pow( LR.cc[0], 2 ) * 9.7*fact2/fact3;
    
    //printf("R=%g\n", R);
    //printf("LR.rh[0]=%g\n", LR.rh[0]);
    //printf("LR.cc[0]=%g\n", LR.cc[0]);
    //printf("bp.pks=%g\n", bp.pks);

    return bp;
}

// Kinney-Graham Model
/*
 blastParam pkOverpress( double distance, double yield, linray LR )
{
    blastParam bp;
    double R, fact1, fact2, fact3, fact4;    
    R = distance/ pow( 1.8*yield,1.0/3 );
    
    fact1 = 808*( 1 + pow( R/4.5, 2 ) );
    fact2 = sqrt( 1 + pow( R/0.048, 2 ) );
    fact3 = sqrt( 1 + pow( R/0.32, 2 ) );
    fact4 = sqrt( 1 + pow( R/1.35, 2 ) );
    bp.pks = 1.0/1.4*LR.rh[0]* pow( LR.cc[0], 2 )* fact1/( fact2*fact3*fact4 );
    
    fact1 = 980*( 1 + pow( R/0.54, 10) );
    fact2 = 1 + pow( R/0.02, 3 );
    fact3 = 1 + pow( R/0.74, 6 );
    fact4 = sqrt( 1 + pow( R/6.9, 2 ) );
    bp.ttd = 1.0/1000*pow( 2*yield, 1.0/3 ) * fact1/ ( fact2*fact3*fact4 );
    
    return bp;
}
*/

// -------------------------------------------------------------------------------
//Calculate positive-phase duration
// -------------------------------------------------------------------------------
double posPhaseDur( double distance, double yield )
{
    double fact = 2;    // for hemispherical blast
    double R = distance/ pow( fact*yield, 1.0/3.0 );
    double taud = pow( 10.0, -2.75 + 0.27*log10(R) );
    return pow(fact*yield, 1.0/3 ) * taud;
}


// -------------------------------------------------------------------------------
// Loading a file with data to read with any number of columns represented by N
// -------------------------------------------------------------------------------
vector< vector<double> > dataReader( const char *fn, int N )
{
    // fn is the filename to read. N is the number of columns of data to read.
    double tmp1[N];
    vector<double> tmp2;
    vector< vector<double> > dat;
    ifstream indata;
    
    indata.open( fn );
    if(!indata)  				// file couldn't be opened
	{ 
		cerr << "Error: File" << fn << " could not be opened!" << endl;
        exit(1);
	}
    
    for( int j=0; j<N; j++ )
        indata >> tmp1[j];
    
    while( !indata.eof() )
    {
        for( int j=0; j<N; j++ )
        {
            tmp2.push_back( tmp1[j] );
            indata >> tmp1[j];
        }
        dat.push_back( tmp2 );
        tmp2.clear();
    }
    return dat;
}


// -------------------------------------------------------------------------------
//Waveform Resizing
// -------------------------------------------------------------------------------
double wfResize( int nn, waveform wf, double *tt, double *pp )
{  
    double pmax, dt;
    double *spl = spline( wf.tt, wf.pp, 1e50, 1e50 );
    
    
  if (1) { //DV
      cout << "in wfResize() saving pp\n";
      FILE *fid;
      fid = fopen( "pp1.dat", "w");
      for (int j=0; j<nn; j++) {
        fprintf(fid, "%15.5E\n", wf.pp[j]);
      }
      
      fclose(fid);
  }
    
    dt    = ( wf.tt.back() - wf.tt.front() )/nn;    
    tt[0] = wf.tt.front();
    pp[0] = wf.pp.front();	
    pmax  = wf.pp.front();
    
    for(int i=0; i<nn-1; i++)
    {   
        tt[i+1] = tt[i] + dt;
        pp[i+1] = splint( wf.tt, wf.pp, spl, tt[i+1] );
        if( pp[i+1] > pmax ) pmax = pp[i+1];
    }
    delete [] spl;
    return pmax;
}


// -------------------------------------------------------------------------------
//Finding Caustics Parameters
// -------------------------------------------------------------------------------
void causticParams( vector<double> sc, linray LR, double *ssup, double *sslo, double *slope )
{   
    double zc, xc, yc, ss_add, nsc, locSlo, actSlo;
    nsc = sc.size();
    if( nsc!=0 )
    {	
        double *spl[4];
        spl[0] = spline( LR.ss, LR.ja, 1e50, 1e50 );
        spl[1] = spline( LR.ss, LR.zz, 1e50, 1e50 );
        spl[2] = spline( LR.ss, LR.xx, 1e50, 1e50 );
        spl[3] = spline( LR.ss, LR.yy, 1e50, 1e50 );
        
        for(int j=0; j < nsc; j++)
        {	
            ss_add = STEP_MAX;
            while( sc[j]+ss_add > LR.ss.back() ) ss_add = ss_add/2;
            
            slope[j] = ( splint(LR.ss,LR.ja, spl[0], sc[j] + ss_add) - splint(LR.ss,LR.ja, spl[0], sc[j] - ss_add) )/( 2*ss_add );
            bool err = true;
            int k = 0;
            while( err )
            {	ss_add = (k+1)*STEP_MAX;
                locSlo = -slope[j]*( ss_add );
                actSlo = splint( LR.ss, LR.ja, spl[0], sc[j]-ss_add );
                err = 0.01 > fabs( ( locSlo - actSlo )/ actSlo );
                k = k+1;
            } 
            ssup[j] = sc[j]+ss_add;
            sslo[j] = sc[j]-ss_add;
            
            zc = 1e-3*splint(LR.ss, LR.zz, spl[1],  sc[j]);
            xc = 1e-3*splint(LR.ss, LR.xx, spl[2], sc[j]);
            yc = 1e-3*splint(LR.ss, LR.yy, spl[3], sc[j]);
            printf("--> Found caustic near \n\taltitude = %.1f km, \n\t%.1f km due East, and \n\t%.1f km due North.\n",zc,xc,yc);
        }
        for( int i=0; i<4; i++ )
            delete [] spl[i];
    }
    else
    {	sslo[0] = LR.ss.back()+1000;
        ssup[0]  = sslo[0];
        slope[0] = 0;
        printf("No caustic detected along the path.\n");
    }
}


// -------------------------------------------------------------------------------
// Step size around a caustic
// -------------------------------------------------------------------------------
double stepSizeCaustic( double ffmx, double scbt, double umax, double sc, double ss, double &dzta )
{   
    // note that dzta is also updated in here
    double ds, stcn;
    stcn = 1/(2*PI*ffmx*scbt*umax);
    dzta = stcn;
    ds = sc - ss - pow( sqrt(fabs(sc-ss)) - dzta/2, 2 ) ;
    if( ds <= 0 )   
    {	
      dzta = -stcn;
      ds = sc - ss + pow( sqrt(fabs(sc-ss)) - dzta/2, 2 ) ;
    }
    if( ds > STEP_MAX ) ds = STEP_MAX;
    dzta = abs(dzta);
    return ds;
}



// -------------------------------------------------------------------------------
// Declaration of the contructor and destructor for the class spectralParams
// -------------------------------------------------------------------------------

spectralParams::spectralParams( int N )
{   
  nn = N;
  usq     = new double[nn];
  cf      = new double[nn/2+1];
  fltr    = new double[nn/2+1];
  attnEf1 = new double[nn/2+1];
  attnEf2 = new double[nn/2+1];
    
	U1  = new fftw_complex[nn/2+1];
	U2  = new fftw_complex[nn/2+1];
	U3  = new fftw_complex[nn/2+1];
	Usq = new fftw_complex[nn/2+1];
}

spectralParams::~spectralParams()
{   delete [] usq;
    delete [] cf;
    delete [] fltr;
    delete [] attnEf1;
    delete [] attnEf2;
    delete [] U1;
    delete [] U2;
    delete [] U3;
    delete [] Usq;
}


// -------------------------------------------------------------------------------
//Frequency and Such
// -------------------------------------------------------------------------------

/*
void spectralParams::frequency( double ffmx, double *ff )	
{
  // since df = 2*ffmx/nn ff is the vector of positive frequencies 0:df:Nyquist
  // cf and fltr are set here as part of spectralParams class data
	for( int i=0; i< nn/2+1; i++ )
	{	
		ff[i] = ffmx*2*i/nn;    // ffmx=f_Nyquist; positive frequencies 0:df:Nyquist
		cf[i] = 0.5*2*PI*ff[i];	// sets cf = (ang freq)/2
		if( i< nn/4 ) {
		  fltr[i] = 1;
		}
		else {
		  fltr[i] = pow( ff[i]/ff[nn/4], 14 ); // 10 was an initial value
		  //if (i < nn/3) {
		  //    printf("fltr[%d]=%g\n", i, fltr[i]);
		  //}
		}
	}
	
	if (1) {
      cout << "in nonlinearRay() in spectralParams::frequency: saving fltr in fltr.dat\n";
      FILE *fid; 
      fid = fopen( "fltr.dat", "w");
      for (int j=0; j<nn/2+1; j++) {
        fprintf(fid, "%15.5E  %15.5E\n", ff[j], fltr[j]);
      }
      
      fclose(fid);
  }
	
}
*/


// DV - 
void spectralParams::frequency( double ffmx, double *ff )	
{
  // since df = 2*ffmx/nn ff is the vector of positive frequencies 0:df:Nyquist
  // cf and fltr are set here as part of spectralParams class data
	for( int i=0; i< nn/2+1; i++ )
	{	
		ff[i] = ffmx*2*i/nn;    // ffmx=f_Nyquist; positive frequencies 0:df:Nyquist
		cf[i] = 0.5*2*PI*ff[i];	// sets cf = (ang freq)/2
		if( i< round(nn/4) ) {
		  fltr[i] = 1;
		}
		else {
		  fltr[i] = pow( ff[i]/ff[(int) round(nn/4)], 20 ); // 100 was Joe's original value
		  //if (i < nn/3) {
		  //    printf("fltr[%d]=%g\n", i, fltr[i]);
		  //}
		}
	}
	
	if (0) {
      cout << "in nonlinearRay() in spectralParams::frequency: saving fltr in fltr.dat\n";
      FILE *fid; 
      fid = fopen( "fltr.dat", "w");
      for (int j=0; j<nn/2+1; j++) {
        fprintf(fid, "%15.5E  %15.5E\n", ff[j], fltr[j]);
      }
      
      fclose(fid);
  }
	
}


// -------------------------------------------------------------------------------
//Initial Fourier Fourier Transform
// -------------------------------------------------------------------------------
void spectralParams::initFFT( double ps, double *pp, double *u )
{
  fftw_plan planFor, psq;
	for(int i=0; i < nn; i++) 
	{ 	
	    u[i]   = pp[i]*ps;    // obtain initial u - see eq. 2 in GJI paper
      usq[i] = pow(u[i],2);	// u^2
	}
	planFor = fftw_plan_dft_r2c_1d(nn,u,U1,FFTW_ESTIMATE); // obtain U1; it's data in class spectralParams
	fftw_execute(planFor);
	psq = fftw_plan_dft_r2c_1d(nn,usq,Usq,FFTW_ESTIMATE);  // obtain Usq
	fftw_execute(psq);
  fftw_destroy_plan( planFor );
  fftw_destroy_plan( psq );
}


/*
// -------------------------------------------------------------------------------
//Initial integration using the Euler method
// -------------------------------------------------------------------------------
void spectralParams::integEulerMethod( double ds, double atfc, double scbt, double dop, double *attn )
{
  // See eq. 12 in Non-lin paper GJI 200 pg 1347 (2015)
  // Note: cf is defined as omega_n/2 (i.e. half the angular frequency - nn/2+1 components
  // U2 = (U1 + ds*cf*U^2)*exp(-a*c0*OM/vr)
  // U3 = U2;

	for( int i=0; i< nn/2+1; i++ )
	{	U2[i][0] = U1[i][0] - ds*cf[i]*scbt*Usq[i][1]; 	// real part of the Fourier component  
		U2[i][1] = U1[i][1] + ds*cf[i]*scbt*Usq[i][0];  // imaginary part
		attnEf2[i] = atfc * attn[i]*fltr[i]*dop;        // alfa*c0*Omega/v_ray; alfa is altered by a factor atfc
		U2[i][0] = U2[i][0]*exp( - attnEf2[i]*ds  );    // apply attenuation; see paper page 1351
		U2[i][1] = U2[i][1]*exp( - attnEf2[i]*ds  );
    U3[i][0] = U2[i][0];
		U3[i][1] = U2[i][1];
	}
}
*/


// DV
// -------------------------------------------------------------------------------
//Initial integration using the Euler method
// -------------------------------------------------------------------------------
void spectralParams::integEulerMethod( double ds, double atfc, double scbt, double dop, double *attn )
{
  // See eq. 12 in Non-lin paper GJI 200 pg 1347 (2015)
  // Note: cf is defined as omega_n/2 (i.e. half the angular frequency - nn/2+1 components
  // U2 = (U1 + ds*cf*U^2)*exp(-a*c0*OM/vr)
  // U3 = U2;

	for( int i=0; i< nn/2+1; i++ )
	{	U2[i][0] = U1[i][0] - ds*cf[i]*scbt*Usq[i][1]; 	// real part of the Fourier component  
		U2[i][1] = U1[i][1] + ds*cf[i]*scbt*Usq[i][0];  // imaginary part
		attnEf2[i] = atfc * attn[i]*fltr[i]*dop;        // alfa*c0*Omega/v_ray; alfa is altered by a factor atfc
		U2[i][0] = U2[i][0]*exp( - attnEf2[i]*ds  );    // apply attenuation; see paper page 1351
		U2[i][1] = U2[i][1]*exp( - attnEf2[i]*ds  );
    U3[i][0] = U2[i][0];
		U3[i][1] = U2[i][1];
	}
}




// -------------------------------------------------------------------------------
// Integration Using the Leapfrog method
// -------------------------------------------------------------------------------
void spectralParams::integLeapfrog( double ds, double scbt )
{
    // U3 = U1 + i*ds*w_n*scbt*Usq - see eq after eq. 12 in GJI article
    for( int i=0; i< nn/2+1; i++ )
    {	
      U3[i][0] = U1[i][0] - 2*ds*cf[i]*scbt*Usq[i][1]; 	  
      U3[i][1] = U1[i][1] + 2*ds*cf[i]*scbt*Usq[i][0];
    }            
}

/* // original
// -------------------------------------------------------------------------------
//Application of Attenuation
// -------------------------------------------------------------------------------
void spectralParams::applyAttn( double ds, double atfc, double dop, double *attn )
{   
    for( int i=0; i< nn/2+1; i++ )
    {   attnEf1[i] = attnEf2[i];
        attnEf2[i] = atfc * attn[i]*fltr[i]*dop;
        U3[i][0] = U3[i][0]*exp(- (attnEf1[i] + attnEf2[i])*ds/2  );
        U3[i][1] = U3[i][1]*exp(- (attnEf1[i] + attnEf2[i])*ds/2  );
        U1[i][0] = U2[i][0];
        U1[i][1] = U2[i][1];
        U2[i][0] = U3[i][0];				
        U2[i][1] = U3[i][1];
    }  
}
*/

// DV - add a filter here -such as a sine filter, etc
// -------------------------------------------------------------------------------
//Application of Attenuation
// -------------------------------------------------------------------------------
void spectralParams::applyAttn( double ds, double atfc, double dop, double *attn )
{   
    for( int i=0; i< nn/2+1; i++ )
    {   attnEf1[i] = attnEf2[i];
        attnEf2[i] = atfc * attn[i]*fltr[i]*dop;
        //U3[i][0] = U3[i][0]*exp(- (attnEf1[i] + attnEf2[i])*ds/2  );
        //U3[i][1] = U3[i][1]*exp(- (attnEf1[i] + attnEf2[i])*ds/2  );
        U3[i][0] = U3[i][0]*exp(- (attnEf1[i] + attnEf2[i])*ds/2  )*sin(PI*(nn/2-i)/(nn)); // DV
        U3[i][1] = U3[i][1]*exp(- (attnEf1[i] + attnEf2[i])*ds/2  )*sin(PI*(nn/2-i)/(nn)); // DV       
        U1[i][0] = U2[i][0];
        U1[i][1] = U2[i][1];
        U2[i][0] = U3[i][0];				
        U2[i][1] = U3[i][1];
    }  
}


// -------------------------------------------------------------------------------
// Fourier Transform
// -------------------------------------------------------------------------------
void spectralParams::FFT( double *u, double &umax )
{   
    // obtain u from U3; set u^2 and U^2; update umax (passed by reference here)
    fftw_plan planBac, planFor;
    planBac = fftw_plan_dft_c2r_1d(nn,U3,u,FFTW_ESTIMATE); // obtain u from U3
    fftw_execute( planBac );	
	
    umax = 0;
    for(int i=0; i < nn; i++) 
    { 	
        u[i] = u[i]/nn;               // proper normalization for u from fft algorithm
        usq[i] = pow(u[i],2);         // u^2
		    if( u[i] > umax ) umax = u[i];
        if( u[i] >= 5.0 )
        {	cout << "\nError! Numerical instability is detected." << endl;
            cout << "Please increase the frequency range.\n" << endl;
            exit(1);
        }
    }
    
    planFor = fftw_plan_dft_r2c_1d(nn,usq,Usq,FFTW_ESTIMATE); // obtain Usq = fft(u^2)
    fftw_execute( planFor );
    fftw_destroy_plan( planBac );
    fftw_destroy_plan( planFor );
}



// -------------------------------------------------------------------------------
//Hilbert transform
// -------------------------------------------------------------------------------
void spectralParams::hilbert()
{
    // U3 is used to store the Hilbert transform of U1 first then of U2
    // U1, U2 end up storing their own transform 
    cout << "--> Passing through a caustic." << endl;
    for(int i=0; i<nn/2+1; ++i)
    {   U3[i][0] = -U1[i][1];
        U3[i][1] = U1[i][0];
        U1[i][0] = U3[i][0];
        U1[i][1] = U3[i][1];
        U3[i][0] = -U2[i][1];
        U3[i][1] = U2[i][0];
        U2[i][0] = U3[i][0];
        U2[i][1] = U3[i][1];
    }
}


//DV
// -----------------------------------------------------------------------------
//Definition of the class raypathParameters and its members
// -----------------------------------------------------------------------------
raypathParams::raypathParams( linray LR )
{	

  //printf("In raypathParams( linray LRinput )\n");

  beta = 1+(GAMMA-1)/2; // set beta

  // do spline interpolation 
  //gsl_interp_accel *acc_;
  //gsl_spline *splin[9];
  
  acc_ = gsl_interp_accel_alloc();

  for (int j=0; j<9; j++) {
    splin[j] = gsl_spline_alloc(gsl_interp_cspline, LR.ss.size());
    //printf("splin[%d]=%p\n", j, splin[j]);
  }
  
  gsl_spline_init(splin[0], LR.ss.data(), LR.xx.data(), LR.ss.size());
  gsl_spline_init(splin[1], LR.ss.data(), LR.yy.data(), LR.ss.size());
  gsl_spline_init(splin[2], LR.ss.data(), LR.zz.data(), LR.ss.size());
  gsl_spline_init(splin[3], LR.ss.data(), LR.cc.data(), LR.ss.size());
  gsl_spline_init(splin[4], LR.ss.data(), LR.vr.data(), LR.ss.size());
  gsl_spline_init(splin[5], LR.ss.data(), LR.om.data(), LR.ss.size());
  gsl_spline_init(splin[6], LR.ss.data(), LR.rh.data(), LR.ss.size());
  gsl_spline_init(splin[7], LR.ss.data(), LR.ja.data(), LR.ss.size());
  gsl_spline_init(splin[8], LR.ss.data(), LR.tr.data(), LR.ss.size());

  //for (int j=1; j<9; j++) {
  //  gsl_spline_free(splin);
  //}
  //gsl_interp_accel_free(acc_);

	nlp.push_back( LR.xx );
	nlp.push_back( LR.yy );
	nlp.push_back( LR.zz );
	nlp.push_back( LR.cc );
	nlp.push_back( LR.vr );
	nlp.push_back( LR.om );
	nlp.push_back( LR.rh );
	nlp.push_back( LR.ja );
  nlp.push_back( LR.tr );
  
  //printf("LR.xx.size=%d   last elem=%g\n", LR.xx.size(), LR.xx.back());
  //exit(1);
  //printf("At end of raypathParams( linray LRinput )\n");  
}


// DV
void raypathParams::physicalParams( double ss )
{

  //printf("In physicalParams( double ss ): ss=%g\n", ss);
  //printf("beta=%g\n", beta);


  // do spline interpolation 
  //gsl_interp_accel *acc_;
  //acc_ = gsl_interp_accel_alloc();
  
  //printf("acc_=%p\n", acc_);
  //int j=0;
  //printf("splin[%d]=%p\n", j, splin[j]);
  //j=1;
  //printf("splin[%d]=%p\n", j, splin[j]);
  //j=2;
  //printf("splin[%d]=%p\n", j, splin[j]);
  //printf("In physicalParams( double ss ) - here 0\n");
  
  xx = gsl_spline_eval(splin[0], ss, acc_ ); // ray's x coord
  yy = gsl_spline_eval(splin[1], ss, acc_ ); // y
  zz = gsl_spline_eval(splin[2], ss, acc_ ); // z
  cc = gsl_spline_eval(splin[3], ss, acc_ ); // c - sound speed
  vr = gsl_spline_eval(splin[4], ss, acc_ ); // ray velocity 
  om = gsl_spline_eval(splin[5], ss, acc_ ); // OMEGA = 1-grad(phi)*v_0 - Doppler effect from wind v_0
  rh = gsl_spline_eval(splin[6], ss, acc_ ); // density
  ja = gsl_spline_eval(splin[7], ss, acc_ ); // Jacobian
  tr = gsl_spline_eval(splin[8], ss, acc_ ); // reduced time

  //gsl_interp_accel_free(acc_);
       
  //printf("vr=%g\n", vr);
  //int j=6;
  //printf("LR.ss[0]=%g;   ss=%g\n", LR.ss[0], ss);
  //printf("j=%d;  LR.ss[j]=%g; LR.xx[j]=%g; LR.zz[j]=%g \n", j, LR.ss[j], LR.xx[j], LR.zz[j]);
  //printf("xx=%g;  yy=%g;  zz=%g\n", xx, yy, zz);
  
  
  //double splint( vector<double> xx, vector<double> yy, double *coef, double x)
  //double x = splint( LR.ss, LR.xx, spl[0], ss);
  //printf("x=%g\n", x);
  
  //exit(1);

  //printf("zz=%g\n", zz);
  //printf("cc=%g\n", cc);
  //printf("LR.ss[1]=%g\n", LR.ss[1]);
  //printf("LR.vr[0]=%g\n", LR.vr[0]);
  //printf("LR.vr[1]=%g\n", LR.vr[1]);
}


// This is the old Joel constructor - DV added int arg to distinguish it from the newer constructor
// -------------------------------------------------------------------------------
//Definition of the class raypathParameters and its members
// -------------------------------------------------------------------------------
raypathParams::raypathParams( linray LRinput, int j )
{	
// routine that is called once to set the coefficients spl used subsequently by splint; 
  double der1 = 1e50;
  //double der1 = 100.0; // DV
  beta = 1+(GAMMA-1)/2;
	LR = LRinput;
    
  spl[0] = spline( LR.ss, LR.xx, der1, der1 );
  spl[1] = spline( LR.ss, LR.yy, der1, der1 );
  spl[2] = spline( LR.ss, LR.zz, der1, der1 );
  spl[3] = spline( LR.ss, LR.cc, der1, der1 );
  spl[4] = spline( LR.ss, LR.vr, der1, der1 );
  spl[5] = spline( LR.ss, LR.om, der1, der1 );
  spl[6] = spline( LR.ss, LR.rh, der1, der1 );
  spl[7] = spline( LR.ss, LR.ja, der1, der1 );
  spl[8] = spline( LR.ss, LR.tr, der1, der1 );
    
	nlp.push_back( LR.xx );
	nlp.push_back( LR.yy );
	nlp.push_back( LR.zz );
	nlp.push_back( LR.cc );
	nlp.push_back( LR.vr );
	nlp.push_back( LR.om );
	nlp.push_back( LR.rh );
	nlp.push_back( LR.ja );
  nlp.push_back( LR.tr );
  
  printf("LR.xx.size=%lu   last elem=%g\n", LR.xx.size(), LR.xx.back());
  

  //exit(1);
}




// DV
raypathParams::~raypathParams()
{
  //printf("At start of raypathParams destructor\n");
  //    for( int i=0; i<9; i++ )
  //        delete [] spl[i];

  //printf("In raypathParams destructor 1 \n");       
  for (int j=0; j<9; j++) {
    gsl_spline_free(splin[j]);
  }   
  //printf("In raypathParams destructor 2\n");
  gsl_interp_accel_free(acc_);     
  //printf("At end of raypathParams destructor\n");
}  
  
// Joel's old destructor
//raypathParams::~raypathParams()
//{
//    for( int i=0; i<9; i++ )
//        delete [] spl[i];
//}
    
void raypathParams::physicalParams( double ss, int j )
{
  splint( LR.ss, nlp, spl, ss, ptNlp );
	xx = ptNlp[0];
	yy = ptNlp[1];
	zz = ptNlp[2];
	cc = ptNlp[3];
	vr = ptNlp[4];
	om = ptNlp[5];
	rh = ptNlp[6];
	ja = ptNlp[7];
  tr = ptNlp[8];
  
  //printf("vr=%g\n", vr);
  //int j=6;
  //printf("LR.ss[0]=%g;   ss=%g\n", LR.ss[0], ss);
  //printf("j=%d;  LR.ss[j]=%g; LR.xx[j]=%g; LR.zz[j]=%g \n", j, LR.ss[j], LR.xx[j], LR.zz[j]);
  //printf("xx=%g;  yy=%g;  zz=%g\n", xx, yy, zz);
  
  
  //double splint( vector<double> xx, vector<double> yy, double *coef, double x)
  //double x = splint( LR.ss, LR.xx, spl[0], ss);
  //printf("x=%g\n", x);
  
  //exit(1);

  //printf("zz=%g\n", zz);
  //printf("cc=%g\n", cc);
  //printf("LR.ss[1]=%g\n", LR.ss[1]);
  //printf("LR.vr[0]=%g\n", LR.vr[0]);
  //printf("LR.vr[1]=%g\n", LR.vr[1]);
}

double raypathParams::pressCorrection( double ssInit )
{
    double p1k, pch, j1k;
    physicalParams( 1000 ); // this distance is in meters; 
                            // ensure the eigenray is longer than 1000 m

    j1k = ja;
    p1k = rh*om*pow(cc,2)/( vr*pow( ssInit/1000,2 ) );
    
    physicalParams( ssInit ); // set initial ray physical parameters, x,y,z, om, etc
    pch = rh*om*pow(cc,2)/( vr*ja/j1k ); 
    
    //printf("pch=%g\n", pch);
    //printf("p1k=%g\n", p1k);
    //printf("j1k=%g\n", j1k);
    //printf("vr=%g\n", vr);
    //exit(1);
    
    return sqrt( pch/p1k );
}

// set normalizing factor NN - 
void raypathParams::normalizingFactor( double pmax, double ssInit)
{
    physicalParams( ssInit ); // initial physical params: interpolates xx,om, vr... at ssInit
	  NN = rh*om*pow(cc,2)/( pow(pmax,2)*fabs(ja)*vr ); // calc with initial ray param values
	  //printf("in normalizingFactor; NN=%g ja=%g; vr=%g om=%g rho=%g, cc=%g\n", NN, ja, vr, om, rh, cc);
}

void raypathParams::psDopScaBeta( double ss, double &dop, double &scbt)
{
    physicalParams( ss ); // interpolates xx, yy, ...om, vr Jac at ss
	  ps   = sqrt( NN*fabs(ja)*vr/( rh*om*pow(cc,2) ) ); // sets ps (p_scaled) in class raypathParams
	  dop  = cc*om/vr;
    scbt = beta*sqrt( pow(om,3)/( NN*rh*cc*cc*fabs(ja)*pow(vr,3) ) ); // eq. (4) in JGI (200) pg 1347
    //printf("in psDopScaBeta; scbt=%g; ja=%g; vr=%g om=%g ss=%g, rho=%g, cc=%g\n", scbt, ja, vr, om, ss, rh, cc);
}

double raypathParams::scaledBetaCaustic( double slope )
{ // beta_caustic = beta_scaled * sqrt(J/J'); here slope is |J'|
	return beta*sqrt( pow(om,3)/( NN*rh*cc*cc*fabs(slope)*pow(vr,3) ) );
}



// -------------------------------------------------------------------------------
// Atmospheric Attenuation
// -------------------------------------------------------------------------------

//-------- Gas fraction polynomial fits -----------------------------------
void fracMole( double z, double *X )
{
    double exp0;
    
    // O2 profile
    if (z < 90.)                                         
        exp0 = -0.67887;
    else  
        exp0 = 49.296 - 1.5524*z + 1.8714e-2*pow(z,2) - 1.1069e-4*pow(z,3) + 3.199e-7*pow(z,4) - 3.6211e-10*pow(z,5);
    X[0] = pow(10.0, exp0 );
    
    // N2 profile
    if (z < 76.)                                         
        exp0 = -0.10744;
    else  
        exp0 = 1.3972E-1 - 5.6269E-3*z + 3.9407E-5*pow(z,2) - 1.0737E-7*pow(z,3);
    X[1] = pow(10.0, exp0);
    
    // CO2 profile
    X[2]  = pow(10.0,-3.3979);                             
    
    // O3 profile
    if (z < 80.) 
        exp0 = -19.027 + 1.3093*z - 4.6496E-2*pow(z,2) + 7.8543E-4*pow(z,3) - 6.5169E-6*pow(z,4) + 2.1343E-8*pow(z,5);
    else
        exp0 = -4.234 - 3.0975E-2*z;              
    X[3] = pow(10.0, exp0);
    
    // O profile
    if (z < 95. )                        
        exp0 = -11.195 + 1.5408E-1*z - 1.4348E-3*pow(z,2) + 1.0166E-5*pow(z,3);
    else
        exp0 = -3.2456 + 4.6642E-2*z - 2.6894E-4*pow(z,2) + 5.264E-7*pow(z,3);
    X[4] = pow(10.0, exp0);
    
    // N profile 
    exp0 = -53.746 + 1.5439*z - 1.8824E-2*pow(z,2) + 1.1587E-4*pow(z,3) - 3.5399E-7*pow(z,4) + 4.2609E-10*pow(z,5);
    X[5]  = pow(10.0, exp0);
    
    // H2O profile
    if (z < 30. )                                        
        exp0 = -1.7491 + 4.4986E-2*z - 6.8549E-2*pow(z,2) + 5.4639E-3*pow(z,3) - 1.5539E-4*pow(z,4) + 1.5063E-06*pow(z,5);
    else if(z < 100)
        exp0 = -4.2563 + 7.6245E-2*z - 2.1824E-3*pow(z,2) - 2.3010E-6*pow(z,3) + 2.4265E-7*pow(z,4) - 1.2500E-09*pow(z,5);
    else
        exp0 = -0.62534 - 0.083665*z;
    X[6] = pow(10.0, exp0);
    X[7] = .012*X[1];
}


//---------Vibrational relaxation-------------------------------------------
void vibRelaxnViscosity( double TT, double PP, double &mu, double *X, double *fvib, double *A_max )
{
    double C_R, pressMuCf;
    double theta[4];
    theta[0]= 2239.1;                                 // Charact. temperature (O2)
    theta[1]= 3352;                                   // Charact. temperature (N2)
    theta[2]= 915;                                    // Charact. temperature (CO2)
    theta[3]= 1037;                                   // Charact. temperature (O3)
    double Cp_R[4], Cv_R[4];
    Cv_R[0] = 5.0/2; //3.0/2;                                    // Heat capacity|volume (O2)
    Cv_R[1] = 5.0/2; //3.0/2;                                    // Heat capacity|volume (N2)
    Cv_R[2] = 3.0;      //3.0/2;                                      // Heat capacity|volume (CO2)
    Cv_R[3] = 3.0;    //3.0/2;                                      // Heat capacity|volume (O3)
    Cp_R[0] = 7.0/2;  //5.0/2;                                    // Heat capacity|pressure (O2)
    Cp_R[1] = 7.0/2;   //5.0/2;                                    // Heat capacity|pressure (N2)
    Cp_R[2] = 4.0;        //5.0/2;                                      // Heat capacity|pressure (CO2)
    Cp_R[3] = 4.0;      //5.0/2;                                      // Heat capacity|pressure (O3)
    
    double mu_o  = 18.192e-6;       // Reference viscosity [kg/(m*s)]
    double P_o = 101325;
    double T_o = 293.15;
    double S  = 117; 
    double Tr = pow(TT/T_o,-1./3.)-1;
    double A1 = (X[0]+X[1])*24*exp(-9.16*Tr);
    double A2 = (X[4]+X[5])*2400;
    double B  = 40400*exp(10*Tr);
    double C  = 0.02*exp(-11.2*Tr);
    double D  = 0.391*exp(8.41*Tr);
    double E  = 9*exp(-19.9*Tr);
    double F  = 60000;
    double G  = 28000*exp(-4.17*Tr);
    double H  = 22000*exp(-7.68*Tr);
    double I  = 15100*exp(-10.4*Tr);
    double J  = 11500*exp(-9.17*Tr);
    double K  = (8.48E08)*exp(9.17*Tr);
    double L  = exp(-7.72*Tr);
    double ZZ = H*X[2]+I*(X[0]+0.5*X[4])+J*(X[1]+0.5*X[5])+K*(X[6]+X[3]);
    double hu = 100*(X[3]+X[6]);
    
    mu = mu_o*sqrt(TT/T_o) * ( 1+S/T_o )/( 1+S/TT ); // Viscosity [kg/(m*s)]
    pressMuCf = (PP/P_o)*mu_o/mu;
    fvib[0] = pressMuCf*(A1+A2+B*hu*(C+hu)*(D+hu));    //O2
    fvib[1] = pressMuCf*(E+F*X[3]+G*X[6]);             //N2
    fvib[2] = pressMuCf*ZZ;                            //C02
    fvib[3] = pressMuCf*(1.2E5)*L;                     //O3
    
    
    for (int m=0; m<4; m++)
    {   C_R          = ((pow(theta[m]/TT,2))*exp(-theta[m]/TT))/(pow(1.0-exp(-theta[m]/TT),2));
        A_max[m]    = (X[m]*(PI/2)*C_R)/( Cp_R[m]*( Cv_R[m]+C_R ) );
    }
}


// -------------------------------------------------------------------------------
// Molar weigth of air as function of altitude
void molGammaTemp( double cc, double *X, double &M, double &gamma, double &TT )
{
    double R = 8.31448;  // universal gas constant (J.mol-1.K)
    double Agamma[]= {1.371, 2.46e-4, -6.436e-7, 5.2e-10, -1.796e-13, 2.182e-17};
    double MX, XX;
    double Mi[8];
    Mi[0] = 0.0319988;          //O2
    Mi[1] = 0.0280134;          //N2
    Mi[2] = 0.0440095;          //CO2
    Mi[3] = 0.0479982;          //O3
    Mi[4] = 0.0159994;          //O
    Mi[5] = 0.0140067;          //N
    Mi[6] = 0.0180153;          //H2O
    Mi[7] = 0.039948;           //Ar
    
    // Mole fraction
    MX = 0;
    XX = 0;
    for( int i=0; i<8; i++ )
    {   MX += X[i]*Mi[i];
        XX += X[i];
    }
    M = MX/XX;
    
    // Gamma and Temperature
    gamma = GAMMA;
    for( int i=0; i<2; i++ )
    {   
        TT = pow(cc,2)*M/(gamma*R);
        gamma = 0;
        for( int j=0; j<6; j++ )
            gamma += Agamma[j]*pow(TT,j);
    }
}


// -------------------------------------------------------------------------------
// Rotational Collision Number
double rotCollisionNumber( double TT, double *X )
{   double AA, Zrot1, Zrot2, Zrot;
    AA = pow(TT,-1./3.);
    Zrot1 = 54.1*exp( -17.3* AA );
    Zrot2 = 63.3*exp( -16.7* AA );
    Zrot   = 1./( ( X[1]/Zrot2 ) + ( X[0]/Zrot1 ) );
    return (4./5.)*sqrt(3./7.)*Zrot;
}


// -------------------------------------------------------------------------------
// Atmospheric Attenuation
void atmAttn(int N, double z, double rh, double cc, double *freq, double *alpha )
{   
    double X[8], fvib[4], A_max[4], M, gamma, TT, PP, mu, nn, X_ON, sig, nu, chi, cchi;
    double pownu2, bt0, bt1bt0, bt2bt0, alph1, alph2, attn_cl, attn_vib, attn_rot;
    
    fracMole( z, X);
    molGammaTemp( cc, X, M, gamma, TT );        // input X, output: M, TT, gamma
    PP = rh*pow(cc,2)/gamma;
    vibRelaxnViscosity( TT, PP, mu, X, fvib, A_max );
    nn = rotCollisionNumber( TT, X );
    X_ON = (X[0] + X[1])/0.9903;
    sig = 5./sqrt(21.);
    
    for(int i=0; i<N; ++i)
    {   
        nu = 8*PI*freq[i]*mu/(3*PP);                   // Nondimensional frequency
        chi = 3.*nn*nu/4.;
        cchi = 2.36*chi;
        pownu2 = pow(nu,2);
        //---------Classical + rotational loss/dispersion--------------------------
        bt0  = 2*PI*freq[i]/cc; 
        bt1bt0  = sqrt( 0.5*( sqrt( 1+ pownu2 )+1 )/( 1+ pownu2 ) );
        bt2bt0  = sqrt( ( 1+pow(cchi,2) )/( 1+pow(sig*cchi,2) ) ); 
        alph1 = bt0 * sqrt( 0.5*( sqrt( 1+pownu2 ) - 1 )/( 1+pownu2 ) );
        alph2 = bt0 * ( sig/2-1/(2*sig) )*chi / sqrt( ( 1+pow(cchi,2) )*( 1+pow(sig*cchi,2) ) );
        attn_cl  = alph1*bt2bt0;
        attn_rot  = alph2 * bt1bt0 * X_ON;
        
        //---------Vibrational absorption-------------------------------------------
        attn_vib = 0.;
        for (int m=0; m<4; m++)
            attn_vib  += A_max[m]/cc * 2*pow(freq[i],2)/fvib[m] /( 1 + pow(freq[i]/fvib[m],2) );
        
        // Total absorption
        alpha[i] = attn_cl*1.003 + attn_rot + attn_vib;
    }
}

