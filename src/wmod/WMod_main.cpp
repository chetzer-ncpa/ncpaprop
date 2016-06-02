#include "Atmosphere.h"
#include "anyoption.h"
#include "ProcessOptionsNB.h"
#include "SolveWMod.h"
#include "WMod_lib.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

/*
 * This is Jelle Assink's Normal Modes (wide angle)- moved to C++
 * Contains some of Roger Waxler's original design also.
 * @version 0
 * @date 2012-09
 * @authors Jelle Assink; Roger Waxler; Claus Hetzer; Doru Velea;
 * 
 * Changelog:
 * 20130326: DV added turnoff_WKB flag
 * 201305  : DV added use_attn_file option (to allow atten. coeff loaded from a text file)
 * 201306  : DV modified the get() functions that return strings
 */

// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

//
// main
//
int main( int argc, char **argv ) {

  // Physical values are usually in SI (System International) units 
  // unless mentioned otherwise.
  
  //parse options from the command line as well as an options file
  AnyOption *opt = parseInputOptions( argc, argv ); 

  // object to process the options
  ProcessOptionsNB *oNB = new ProcessOptionsNB(opt);

  string atmosfile      = "";          // stores the atmospheric profile name
  string atmosfileorder = "";          // column order e.g. 'zuvwtdp'
  string wind_units     = "mpersec";   // m/s
  int    skiplines      = 0;           // skiplines in "atmosfile"
  
  // set up to measure the duration of this run
  time_t tm1 = time(NULL);

  // obtain the parameter values from the user's options
  atmosfile      = oNB->getAtmosfile();
  atmosfileorder = oNB->getAtmosfileorder();
  wind_units     = oNB->getWindUnits();
  //gnd_imp_model  = oNB->getGnd_imp_model(); 
  skiplines      = oNB->getSkiplines();  
  
  // get atmospheric profile object; the azimuth is set inside SolveWMod
  SampledProfile *atm_profile = new SampledProfile( atmosfile, atmosfileorder.c_str(), skiplines );

  // get solver object 
  SolveWMod *a = new SolveWMod(oNB, atm_profile);
  a->printParams();                    
	 
	//   					 
  // compute modes - main action happens here					 
  //
  a->computeModes();
	
  // plot with gnuplot - probably not part of the final release; 
  // good for a quick check of results				 	
  //if (oNB->getPlot_flg()) {
  //    //plotwGNUplot_from_script();
  //    //plotwGNUplot(freq, write_2D_TLoss, 0, write_phase_speeds);
  //    plotwGNUplot(oNB);
  //}
  
  a->printParams();
  
  // save atm. profile if requested
  if (oNB->getWriteAtmProfile()) {
      //saveAtm_profile(atm_profile);
      saveAtm_profile(atm_profile, wind_units);
      //printf(" write_atm_profile flag : %d\n", write_atm_profile);
  }	
  else { printf(" write_atm_profile flag : %d\n", oNB->getWriteAtmProfile()); }
  	  
  delete a;		  	 
  delete atm_profile;
  delete opt;
  delete oNB;
	cout << "\n ... main() is done." << endl;
	time_t tm2 = time(NULL);
  cout << "Run duration: " << difftime(tm2,tm1) << " seconds." << endl;	

	return 0;
} // end of main();


AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();

  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "|                             NCPA Infrasound                              |" );  
  opt->addUsage( "|                Wide Angle High-Mach Number Normal Modes                  |" );
  opt->addUsage( "|                            Single Frequency                              |" );
  opt->addUsage( "----------------------------------------------------------------------------" );	
  opt->addUsage( "Usage: " );
  opt->addUsage("By default the program computes the 1D transmission loss (TL)" );
  opt->addUsage("at the ground and saves the data to 2 files:" );
  opt->addUsage( "   file wtloss_1d.nm - considering attenuation in the atmosphere" );
  opt->addUsage( "   file wtloss_1d.lossless.nm  - no attenuation" );
  opt->addUsage( "Additionally, if the flag --write_2D_TLoss is present on the command line the 2D TL is saved to file wtloss2d.nm" );
  opt->addUsage( "The user can also choose to propagate in N different directions" );
  opt->addUsage( "i.e. (N by 2D mode) by using the option --Nby2Dprop ." );
  opt->addUsage( "" );  
  opt->addUsage( "The options below can be specified in a colon-separated file \"WMod.options\" or at the command line.  Command-line options override file options." );
  opt->addUsage( " --help -h                Print these instructions and exit" );
  opt->addUsage( "" );
  opt->addUsage( "To use an arbitrary 1-D atmospheric profile in ASCII format" );
  opt->addUsage( "(space or comma-separated) the following options apply:" );	
  opt->addUsage( "" );	
  opt->addUsage( "REQUIRED (no default values):" );
  opt->addUsage( " --atmosfile  <filename>  Uses an ASCII atmosphere file" );
  opt->addUsage( "                          referenced to Mean Sea Level (MSL)." );  
  opt->addUsage( " --atmosfileorder         The order of the (z,t,u,v,w,p,d) fields" );
  opt->addUsage( "                          in the ASCII file (Ex: 'zuvwtdp')" );
  opt->addUsage( "                          The units assumed in the ASCII file are" );
  opt->addUsage( "                          z[km], t [kelvin], d [g/cm^3], p [hectoPa]" );
  opt->addUsage( "                          The wind speeds are in m/s by default;" );
  opt->addUsage( "                          however if the winds are given in km/s then use" ); 
  opt->addUsage( "                          option --wind_units kmpersec" );
  opt->addUsage( " --skiplines              Lines at the beginning of the ASCII file to skip [0]" );
  opt->addUsage( " --azimuth                Degrees in range [0,360), clockwise from North" );
  opt->addUsage( " --freq                   Frequency [Hz]" );	
  opt->addUsage( "" );
  opt->addUsage( "REQUIRED for propagation in one direction (no default values):" );
  opt->addUsage( " --azimuth                Degrees in range [0,360], clockwise from North" );
  opt->addUsage( "" );
  opt->addUsage( "REQUIRED for propagation in N directions i.e. (N by 2D) (no default values):" );
  opt->addUsage( " --azimuth_start          Start azimuth ([0,360] degrees, clockwise from North)" );
  opt->addUsage( " --azimuth_end            End azimuth ([0,360] degrees, clockwise from North)" );
  opt->addUsage( " --azimuth_step           Step by which the azimuth is changed (in degrees)" );    	
  opt->addUsage( "" ); 	
  opt->addUsage( "OPTIONAL [defaults]:" );  
  opt->addUsage( " --maxheight_km           Calculation grid height in km above MSL [150 km]" );
  opt->addUsage( " --zground_km             Height of the ground level above MSL [0 km]" );  
  opt->addUsage( " --Nz_grid                Number of points on the z-grid from ground to maxheight [20000]" );  
  opt->addUsage( " --sourceheight_km        Source height in km Above Ground Level (AGL) [0]" );
  opt->addUsage( " --receiverheight_km      Receiver height in km AGL [0]" );
  opt->addUsage( " --maxrange_km            Maximum horizontal propagation distance from origin [1000 km]" );  
  opt->addUsage( " --Nrng_steps             Number of range steps to propagate [1000]" );	
  opt->addUsage( " --ground_impedance_model Name of the ground impedance models to be employed:" );
  opt->addUsage( "                          [rigid], others TBD" );
  opt->addUsage( " --Lamb_wave_BC           If ==1 it sets admittance = -1/2*dln(rho)/dz; [ 0 ]" );
 	opt->addUsage( " --wind_units             Use it to specify 'kmpersec' if the winds are given in km/s [mpersec]" );
	opt->addUsage( " --use_attn_file          Use it to specify a file name containing user-provided" );
	opt->addUsage( "                          attenuation coefficients to be loaded instead of " );
	opt->addUsage( "                          the default Sutherland-Bass attenuation. " ); 
	opt->addUsage( "                          The text file should contain two columns: " );
	opt->addUsage( "                              height (km AGL) and " );
	opt->addUsage( "                              attenuation coefficients in np/m." ); 	
  opt->addUsage( "" );	
  opt->addUsage( "FLAGS (no value required):" );
  opt->addUsage( " --write_2D_TLoss         Outputs the 2D transmission loss to" );
  opt->addUsage( "                          default file: wtloss2D.nm" );	
  opt->addUsage( " --write_phase_speeds     Output the phase speeds to" );
  opt->addUsage( "                          default file: wphasespeeds.nm" );
  opt->addUsage( " --write_modes            Output the modes to default files:" );
  opt->addUsage( "                          wmode_<mode_count>.nm" );
  opt->addUsage( " --write_dispersion       Output the mode dispersion to" );
  opt->addUsage( "                          default file: wdispersion_<freq>.nm" );
  opt->addUsage( " --Nby2Dprop              Flag to perform (N by 2D) propagation i.e." );
  opt->addUsage( "                          propagation in N directions specified by" ); 
  opt->addUsage( "                          options: azimuth_start, azimuth_end, azimuth_step " );
  opt->addUsage( "                          The output is saved into default files: " );  
  opt->addUsage( "                              Nby2D_wtloss_1d.lossless.nm" );
  opt->addUsage( "                              Nby2D_wtloss_1d.nm" );   
  opt->addUsage( " --write_atm_profile      Save the interpolated atm. profile to" );
  opt->addUsage( "                          default file: atm_profile.nm" );
  opt->addUsage( " --turnoff_WKB            Turn off the WKB least phase speed estimation" );
  opt->addUsage( "                          an approx. that speeds-up ground-to-ground propag." ); 
  opt->addUsage( "                          It has the value 1 (true) if any of the flags" );
  opt->addUsage( "                          write_2D_TLoss, write_phase_speeds, write_modes" );
  opt->addUsage( "                          or write_dispersion are true." );    
  opt->addUsage( "" ); 
  opt->addUsage( "" );
  opt->addUsage( " OUTPUT Files:  Format description (column order):" );
  opt->addUsage( "  wtloss_1d.nm:           r, 4*PI*Re(P), 4*PI*Im(P), (incoherent TL)" );
  opt->addUsage( "  wtloss_1d.lossless.nm:" );
  opt->addUsage( "  wtloss_2d.nm:           r, z, 4*PI*Re(P), 4*PI*Im(P)" );
  opt->addUsage( "  Nby2D_wtloss_1d.nm:     r, theta, 4*PI*Re(P), 4*PI*Im(P), (incoherent TL)" );
  opt->addUsage( "  Nby2D_wtloss_1d.lossless.nm:" );
  opt->addUsage( "" );    
  opt->addUsage( "  wphasespeeds.nm:        Mode#, phase speed [m/s], imag(k)" );
  opt->addUsage( "  wmode_<mode_count>.nm   z, (Mode amplitude)" );
  opt->addUsage( "  wdispersion_<freq>.nm   Contains one line with entries:" );
  opt->addUsage( "                          freq, (# of modes), rho(z_src)," );
  opt->addUsage( "                          followed for each mode 'i' by quadruples: " );
  opt->addUsage( "                          real(k(i)), imag(k(i)), Mode(i)(z_src), Mode(i)(z_rcv)" );  
  opt->addUsage( "" );   
  opt->addUsage( "  atm_profile.nm           z,u,v,w,t,d,p,c,c_eff" );  
  opt->addUsage( "" ); 
  
  opt->addUsage( " Examples (run from 'samples' directory):" );
  opt->addUsage( "    ../bin/WMod --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1" );
  opt->addUsage( "" );  
  opt->addUsage( "    ../bin/WMod --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1 --write_2D_TLoss --sourceheight_km 60 --receiverheight_km 60" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/WMod --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --freq 0.1 --Nby2Dprop --azimuth_start 0 --azimuth_end 360 --azimuth_step 1" );
  opt->addUsage( "" );  
  
  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  opt->setFlag( "write_2D_TLoss" );
  opt->setFlag( "write_phase_speeds" ); 
  opt->setFlag( "write_modes" );
  opt->setFlag( "write_dispersion");
  opt->setFlag( "Nby2Dprop");  
  opt->setFlag( "write_atm_profile");
  opt->setFlag( "turnoff_WKB");
  opt->setFlag( "plot" );

  opt->setOption( "atmosfile" );
  opt->setOption( "atmosfileorder" );
  opt->setOption( "wind_units" );
  opt->setOption( "skiplines" );		
  opt->setOption( "azimuth" );
  opt->setOption( "azimuth_start" );
  opt->setOption( "azimuth_end" );
  opt->setOption( "azimuth_step" );  
  opt->setOption( "freq" );
  opt->setOption( "maxrange_km" );
  opt->setOption( "sourceheight_km" );
  opt->setOption( "receiverheight_km" );
  opt->setOption( "maxheight_km" );
  opt->setOption( "zground_km" );  
  opt->setOption( "stepsize" );
  opt->setOption( "Nz_grid" );
  opt->setOption( "Nrng_steps" );
  opt->setOption( "write_2D_TLoss" );
  opt->setOption( "ground_impedance_model" );
  opt->setOption( "Lamb_wave_BC" );
  opt->setOption( "use_attn_file" );

  // Process the command-line arguments
  opt->processFile( "./WMod.options" );
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


 
