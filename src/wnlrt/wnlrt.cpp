//
// wnlrt.cpp
// Much of this code is Joel Lonzaga's code adapted by Doru Velea
//  
//  To compile:
//       g++ -o wnlrt wnlrt.cpp nonlinearRay.cpp linearRay3DStrat.cpp -lfftw3 -lgsl -lgslcblas

//  To run:
// ./wnlrt --eigenrayfile ToyAtmo_Eigenray-0.dat 

// Another option to run is the follwoing but it has not been decided if it will make 
// it in the final cut of the code
// ./wnlrt --inclin 10 --azimuth 50 --eigenrayfile ToyAtmo_Eigenray-0.dat --shootray --atmosfile NCPA_canonical_profile_zuvwtdp.dat --src_z 0 --range 1000 --rcv_z 1

#include <cstring>
#include <stdexcept>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "linearRay3DStrat.h"
#include "nonlinearRay.h"
#include "anyoption.h"
#include "ProcessOptionsNRT.h"

using namespace NCPA;
using namespace std;

void save2disk( const char*, double, double, linray, nonray );
void save2disk( const char *fn, double atfac, linray ER, nonray NR );

void saveRay(linray);

// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

// convert to low case
void lowchar(std::string &s);

int main(int argc,char **argv) {

    // parse options from the command line as well as an options file
    AnyOption *opt = parseInputOptions( argc, argv ); 

    // object to process the options
    ProcessOptionsNRT *oNRT = new ProcessOptionsNRT(opt);
    //cout << "theta=" << oNRT->getInclination() << endl;
    //cout << "azi =" << oNRT->getAzimuth() << endl;
    
    double theta = oNRT->getInclination();
    double phi   = oNRT->getAzimuth();
    double zs    = oNRT->getSrcz();
    double xr    = oNRT->getRcvx();
    double yr    = oNRT->getRcvy();
    double zr    = oNRT->getRcvz();    
    double dth   = oNRT->getDth();
    double dph   = oNRT->getAzimuth();;
    double tol   = oNRT->getTol();
      
    double ds    = 100; // meters
    string fn;
    //const char *wfFile   = argv[8]; 
    
    // convert to radians   
    theta        = theta*PI/180;
    phi          = phi*PI/180;
    dth          = dth*PI/180;
    dph          = dph*PI/180;    

    clock_t start, end; // set, start timer
    start = clock();
   
    //----------------------------------------------------------
	  // Parameters
    vector<double> yield;
    yield.push_back( 100e3 );
    
    double ssInit = 1000;  // reference (start) point on ray [meters]
    int nn = pow(2.0,14);  // number of points in the waveform
    double period = 80.0;  // reduced time ranges over this period from -T/2 to T/2 
    attRed attRed;         // attenuation reduction factor structure
    attRed.zzRed = 90000;  // aplly attn. red. factor above this height
    attRed.atfac = 1 ;

    profile prf;
    linray ER;
    //-----------------------------------------------------------

    // Ray Path and Linear Ray Amplitude Calculation
    
    // Option 1: provide the eigenray file
    if (opt->getValue( "eigenrayfile" ) != NULL)  {
      fn = oNRT->getEigenfile();
      printf("--> Reading the eigenray from file: \"%s\".\n", fn.c_str());
      ER = eigrayReader(fn.c_str());
      if( ER.rsd > tol ) return 0;
    }
    else {  // Options 2,3 - will see if they survive in the final version of the package
	    //Reading the atmospheric specifications
	    fn = oNRT->getAtmosfile();
	    printf("--> Loading profile %s\n", fn.c_str());
      prf = profileReader( fn.c_str() );
	    if ( opt->getFlag( "findeigenray") ) {
        //-----------------------------------------------------------
        // Ray Path and Linear Ray Amplitude Calculation
        printf("--> Determining an eigenray for the profile \"%s\".\n", fn.c_str());
        ER = eigenray( ds, theta, phi, dth, dph, zs, xr, yr, zr, tol, prf );
        if( ER.rsd > tol ) return 0;
        
      } else if (opt->getFlag( "shootray")) {
        // determine the ray 
        ER = linearRay3DStrat( ds, theta, phi, zs, oNRT->getRange(), zr, prf ); 
      } else {
         cout << "--> Please specify one of 3 options: eigenrayfile, findeigenray or shootray"  << endl;
         return(0);
      }
    }
    
    //----------------------------------------------------------
	  // Determining pressure correction when initial distance is 
	  // other than 1 km from source. Otherwise, it is just 1.
    raypathParams RP( ER );
    //double presCor = RP.pressCorrection( ssInit );  // CHH 191029: Unused
    
    //----------------------------------------------------------
    //Generating model waveform
    waveform wf; // waveform with fields reduced time (tt) and pressure (pp)
    
    // wf struct has fields tt, pp;
    // T=period; sampling interval = dt = T/nn; 
    // tt is from -T/2 to T/2
    // N wave extends from [-ttd:ttd]
                                                         
    //blastParam bp = pkOverpress( ssInit, yield[0], ER );
    //double pks = presCor*bp.pks;
    //double ttd = bp.ttd;
    
    //waveform wf = loadBlastModel( wfFile, nn, period, ttd, pks );
    //waveform wf = blastmodel( nn, period, ttd, pks, ER );
    
    //double amp = 1000;
    //waveform wf = generateWaveForm(nn, period, amp ); // older version
    
    //double pks = 5000; // Pa
    //double ttd = 0.1;  // seconds
    //printf("ttd=%g \t pks=%g\n", ttd, pks);
    //cout << "oNRT->getWftype()=" << oNRT->getWftype() << endl;
    //cout << "(oNRT->getWftype()).compare() = " << (oNRT->getWftype()).compare("Nwave") << endl;
    
    // generate waveform with fields reduced time (tt) and pressure (pp)
    if ( (oNRT->getWftype()).compare("Nwave")==0) {
      //cout << "waveform type: " << oNRT->getWftype() <<  endl;
      wf = generateNWave( nn, period, oNRT->getWfduration(), oNRT->getWfampl() ); 
    } 
    else if ( (oNRT->getWftype()).compare("pulse")==0 ) {
      double ampl = oNRT->getWfampl();
      cout << "--> Initial amplitude = " << ampl << " Pa" << endl;
      wf = generateWaveForm(nn, period, ampl );
      //cout << "after generateWaveForm: ampl = " << ampl << endl;
    }                                      

    string init_wavef = "initwf.dat"; 
    savewf( init_wavef.c_str(), wf );
    cout << "--> Initial waveform saved to file: " << init_wavef << " with 2 columns" << endl;
    cout << "    [reduced_time, waveform] " << endl;       
    
    //----------------------------------------------------------
    // Nonlinear ray calculation
    //
    // printf("Starting the nonlinear ray acoustic calculation for yield=%.0f tonnes, tv=%.0f m/s, and attnReduxn=%.1f.\n", 1e-3*yield[0], ER.tv, attRed.atfac);
    printf("--> Starting the nonlinear ray acoustic calculation\n");
        
    // main work is done here
    nonray NR = nonlinearRay( nn, ssInit, attRed, ER, wf );
    
    //---------- end of Nonlinear ray calculation ---------------
    
    // Saving needed parameters
    //save2disk( fn.c_str(), 1e-3*yield[0], attRed.atfac, ER, NR );
    
    printf("--> Saving to disk ...\n");    
    save2disk( fn.c_str(), attRed.atfac, ER, NR );
    
    if (0) {
      saveRay(ER); // other output; not for final release
    }    
    
    end = clock();
    cout << "--> ...done. (Run time = " << (end-start)/CLOCKS_PER_SEC << " s.)" << endl << endl;
    
    delete opt;
    delete oNRT;
    return 0;
    
};  // end of main()


// Some additional functions
//----------------------------------------------------------
// Script to save to disk needed data: several functions with different parameter lists
//----------------------------------------------------------

void save2disk( const char *fn, double atfac, linray ER, nonray NR )
{
    char filename1[256], filename2[256], filename3[256];
       
    sprintf(filename1, "pressure_wf_evolution.dat");
    sprintf(filename2, "ray_params.dat");
    sprintf(filename3, "final_waveform_spectrum.dat");
    //sprintf(filename4, "%s_Cef_%03.0f_%.3f_%.1fEx.txt", tmp, yield, ER.tv, atfac );    
    
    FILE *acoPress, *acoParam, *acoSpect;
    acoPress = fopen( filename1, "w");
    acoParam = fopen( filename2, "w");
    acoSpect = fopen( filename3, "w");
    //acoCeffe = fopen( filename4, "w");
    
    //double cef;            // CHH 191029: Unused
    //double theta = ER.th;  // CHH 191029: Unused
    //double phi   = ER.ph;  // CHH 191029: Unused
    
    // save effective sound speed
    //for(unsigned int i = 0; i<prf.zz.size(); i++)
    //{   cef = prf.cc[i] + prf.wx[i]*cos(theta)*cos(phi) + prf.wy[i]*cos(theta)*sin(phi);
    //    fprintf( acoCeffe, "%.5E %15.5E %15.5E %15.5E %15.5E\n", prf.zz[i], prf.cc[i], prf.wx[i], prf.wy[i], cef );
    //}
    
    int vecsiz = NR.tt.size();
    int stpsiz = NR.uu.size();
    
    //printf("stpsiz=%d\n", stpsiz);
    //printf("vecsiz=%d\n", vecsiz);
   
   // save ray path and pressure scaling factor ps?
   printf("--> Saving ray info and the pressure scaling factor to file %s\n", filename2);
   printf("    with columns: x, y, z, raypath_length, travel_time, pressure scaling factor\n");
    for(int i = 0; i < stpsiz; i++)  
       fprintf (acoParam, "%15.5E %15.5E %15.5E %15.5E %15.5E %15.5E\n", NR.xx[i], NR.yy[i], NR.zz[i], NR.ss[i], NR.tr[i], NR.ps[i]); 
    
    //---------------------------------------------------------------------
    // Saving waveform along the ray path
    
    printf("--> Saving waveform evolution along the ray path to %s\n", filename1);
    printf("    with %d columns:\n", stpsiz+1);
    printf("    [ reduced_time, pressure waveform_at_step1, waveform_at_step2, 3, 4, ...etc.]");   
    printf("    See the time steps in column 5 in the accompanying file %s.\n", filename2);
    
     for(int i = 0; i < vecsiz; i++ )
     {   
       fprintf( acoPress, "%.5E", NR.tt[i] ); // save reduced time in the first column
       for( int j = 0; j < stpsiz; j++ ) {
         //fprintf( acoPress, "%15.5E", NR.uu[j][i] );
         fprintf( acoPress, "%15.5E", NR.uu[j][i]/NR.ps[j] );
       }
       fprintf( acoPress, "\n");
     }
    /*
    //---------------------------------------------------------------------
    // Saving waveform and spectrum at the receiver
    for(int i = 0; i < vecsiz; i++ )
        fprintf( acoPress, "%.5E %15.5E\n", NR.tt[i], NR.uu[stpsiz-1][i] );
    */
    vecsiz = NR.ff.size();
    printf("--> Saving waveform spectrum at receiver to %s\n", filename3);
    printf("    with 3 columns: [frequency, real part, imag part].\n");
    for(int i = 0; i<vecsiz; i++) {
        //fprintf( acoSpect, "%.5E %15.5E %15.5E\n", NR.ff[i], NR.Ur[i], NR.Ui[i] );
        fprintf( acoSpect, "%.5E %15.5E %15.5E\n", NR.ff[i], NR.Ur[i]/NR.ps[stpsiz-1], NR.Ui[i]/NR.ps[stpsiz-1] );
    }
    
    fclose( acoParam );
    fclose( acoPress );
    fclose( acoSpect );
    //fclose( acoCeffe );
}



// save contents of a ray
void saveRay(linray R) {

FILE *fp;
string rayf = "currentRay.dat";
fp = fopen(rayf.c_str(), "w");

printf("Saving ray parameters in file %s\n", rayf.c_str());
printf("with columns: [ x, y, z, raypath_length, travel_time, OMEGA, Jacobian ]\n");

//printf("R.xx.size()=%d\n",R.xx.size());
for(int i = 0; i < R.xx.size(); i++)  {
  //printf("i=%d\n",i);
  fprintf (fp, "%15.5E %15.5E %15.5E %15.5E %15.5E %15.5E %15.5E\n", R.xx[i], R.yy[i], R.zz[i], R.ss[i], R.tr[i], R.om[i], R.ja[i]); 
}
printf("... done saving ray.\n");

}

// convert to lower case 
void lowchar(std::string &s) {
  int i=0; 
  char c;
  while (s[i]) {
    c = s[i];
    //putchar(tolower(c));
    s[i] = tolower(c);
    i++;
  }
}


AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();
  
  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "|                 NCPA Infrasound Propagation Package                      |" );
  opt->addUsage( "|            wnlrt - Weakly non-linear ray tracing module                  |" );
  opt->addUsage( "|               Authors: Joel Lonzaga (main); Doru Velea                   |" );
  opt->addUsage( "|                       Coordinator: Roger Waxler                          |" );    
  opt->addUsage( "|                           September 2015                                 |" );  
  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "This code uses a weakly non-linear ray tracing algorithm to propagate a signal" );  
  opt->addUsage( "from an impulsive source to a receiver through a stratified atmosphere with " );
  opt->addUsage( "winds. For the theory and algorithm refer to: " );  
  opt->addUsage( "Joel B. Lonzaga, Roger M. Waxler, Jelle D. Assink and Carrick L. Talmadge," ); 
  opt->addUsage( "\"Modeling waveforms of infrasound arrivals from impulsive sources using " );  
  opt->addUsage( "weakly non-linear ray theory\", Geophys. J. Int. (2015) 200, 1337-1361." );  
  opt->addUsage( "" );  
  opt->addUsage( "Usage: " );
  opt->addUsage( "./wnlrt [--option1 val1] [--option2 val2] [--flag1] [...]" );  
  opt->addUsage( "" );  
  opt->addUsage( "The options below can be specified at the command line or in a colon-separated" );
  opt->addUsage( "file \"wnlrt.options\". Command-line options override file options." ); 
  opt->addUsage( "Be sure to precede all options with two minuses (--). The program options can" );
  opt->addUsage( "be of two kinds: pairs of [--option_name value] or flags. The values can be " );
  opt->addUsage( "numbers or strings (arrays of characters). The flags do not take a value on the" );  
  opt->addUsage( "command line; they are boolean switches signaling the program for a certain " );  
  opt->addUsage( "action to occur." );  
  opt->addUsage( "" );  
  opt->addUsage( " --help -h            Print this message and exit" );
  opt->addUsage( "" );
  opt->addUsage( "PREREQUISITE" );  
  opt->addUsage( "This code requires the user to provide an eigenray file obtained from an external" );  
  opt->addUsage( "ray-tracing package GeoAc1.1.1 authored by Phil Blom, currently at Sandia" );
  opt->addUsage( "National Lab, and whose work began at NCPA under the direction of Roger Waxler." );  
  opt->addUsage( "(As an example see the included ToyAtmo_Eigenray-0.dat file.)" );  
  opt->addUsage( "The 15-column eigenray file has a 4-line header providing information about " );
  opt->addUsage( "the source and receiver locations, the ray launch angles and the column names" );
  opt->addUsage( "as in the following example:" ); 
  opt->addUsage( "  Source Location (kilometers) : (0, 0, 0)." );  
  opt->addUsage( "  Receiver Location (kilometers) : (350, 0, 0)." );  
  opt->addUsage( "  theta = 18.0303, phi = 90,  c0(zground) = 0.340322 km/s" );  
  opt->addUsage( "  # x [km]  y [km]  z [km]  Geo. Atten. [dB]  Atmo. Atten. [dB] Travel Time [s]" );  
  opt->addUsage( "    rho [gm/cm^3]  c [km/s]  u [km/s]  v [km/s]  w [km/s]" );  
  opt->addUsage( "    Slowness_x [s/km]  Slowness_y [s/km]  Slowness_z [s/km]  Jacobian [km^2/rad^2]" );  
  opt->addUsage( "" );
  opt->addUsage( "REQUIRED options:" );
  opt->addUsage( " --eigenrayfile        Provide name of previously obtained eigenray file." );
  opt->addUsage( "                       Currently this eigenray file is obtained by running" );
  opt->addUsage( "                       the modified version of GeoAc1.1.1 that outputs " );
  opt->addUsage( "                       relevant parameters along the ray path." );    	
  opt->addUsage( "" ); 
  	
  opt->addUsage( "OPTIONAL options [defaults]:" );
  opt->addUsage( " --waveform            Provide the type of waveform at the source. " );
  opt->addUsage( "                       Can be ""Nwave"" or ""pulse"" [Nwave]. " );
  opt->addUsage( " --ampl                Provide the initial pressure waveform amplitude. [1 Pa]" ); 
  opt->addUsage( " --duration            Provide the initial waveform duration [0.5 secs]" );
  opt->addUsage( "" );
  
  
  opt->addUsage( "OUTPUT text files:" );
  opt->addUsage( " ray_params.dat" );
  opt->addUsage( "                       Contains ray info and acoustic pressure [Pa]; 6 columns " );
  opt->addUsage( "                       [ x, y, z, raypath_length, travel_time, acoustic pressure ]" );
  opt->addUsage( "" );  
  opt->addUsage( " pressure_wf_evolution.dat" );
  opt->addUsage( "                       Stores N pressure waveforms at N time steps along the ray." );
  opt->addUsage( "                       The pressure units are Pascals." );
  opt->addUsage( "                       The file contains (N+1) columns in the following order: " );  
  opt->addUsage( "                       [ reduced_time, pressure waveform_at_step_1 ..." );
  opt->addUsage( "                       [ pressure waveform_at_step_2 ... waveform_at_step_N ]" );  
  opt->addUsage( "" );
  opt->addUsage( " final_waveform_spectrum.dat" );
  opt->addUsage( "                       Stores waveform spectrum at the ray's end point " );
  opt->addUsage( "                       Contains 3 columns: [ frequency, real part, imag part ]" );
  opt->addUsage( "" );   

  opt->addUsage( "QUICK-START EXAMPLES:" );  
  opt->addUsage( "./wnlrt --eigenrayfile ToyAtmo_Eigenray-0.dat" );  
  opt->addUsage( "./wnlrt --eigenrayfile ToyAtmo_Eigenray-0.dat --waveform Nwave --ampl 500 --duration 0.5" ); 
  opt->addUsage( "" ); 

  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  
  //opt->setFlag( "use_eigenrayfile" );
  opt->setFlag( "findeigenray" );
  opt->setFlag( "shootray" );

  opt->setOption( "eigenrayfile" );
  opt->setOption( "atmosfile" );
  opt->setOption( "atmosfileorder" );
  
  opt->setOption( "inclin" );
  opt->setOption( "azimuth" );
  opt->setOption( "src_z" );
  opt->setOption( "rcv_x" );
  opt->setOption( "rcv_y" );
  opt->setOption( "rcv_z" );
  opt->setOption( "range" );
  opt->setOption( "dth" );
  opt->setOption( "daz" );
  opt->setOption( "tol" );
  
  opt->setOption( "waveform" );
  opt->setOption( "ampl" );
  opt->setOption( "duration" );

  // Process the command-line arguments
  opt->processFile( "../samples/wnlrt.options" );
  opt->processCommandArgs( argc, argv );

  if( ! opt->hasOptions()) { // print usage if no options
		  opt->printUsage();
		  delete opt;
		  exit( 1 );
  }

  // Check to see if help text was requested
  if ( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) {
	  opt->printUsage();
	  exit( 1 );
  }

  return opt;
}

