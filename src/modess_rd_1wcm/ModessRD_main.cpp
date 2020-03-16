#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <complex>

#include "Atmosphere.h"
#include "anyoption.h"
#include "ProcessOptionsNBRD.h"
#include "SolveModNB.h"
#include "ModessRD_lib.h"

#include "slepceps.h"
#include "slepcst.h"

#include "binaryreader.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

/*
 * Range-Dependent Normal Modes One-Way Coupled Modes - Effective Sound Speed.
 * 
 * @version 0
 * @date 2012-09
 * @authors Doru Velea; Jelle Assink; Roger Waxler; Claus Hetzer; 
 * Follows Doru Velea's notes: Normal Modes for Range-Dependent Environments: 
 * One-way Coupled Modes i.e. NMRD-OWCM
 * See also the theory in "Computational Ocean Acoustics" Section 5.9, page 315, 1994 ed.
 * 
 * Changelog:
 * 20130326: DV added turnoff_WKB flag
 * 201305  : DV added use_attn_file option (to allow atten. coeff loaded from a text file)
 * 201306  : DV modified the get() functions that return strings
 * 20131111: DV modified the code to handle regions having different number of modes 
 *           i.e. matrices R and C are now allowed to be rectangular. Makes for 
 *           slightly more optimized code.
 */

// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

//
// main
//
int main( int argc, char **argv ) {

  // Physical values are usually in SI (System International) units 
  // unless mentioned otherwise.
  
  // set up timer to measure the duration of this run
  time_t tm1 = time(NULL);  
  
  // parse options from the command line as well as an options file
  AnyOption *opt = parseInputOptions( argc, argv ); 
  
  // object to process the options
  ProcessOptionsNB *oNB = new ProcessOptionsNB(opt);
  
  // these strings may be populated by the run-time user options
  string prf_ranges_km    = "none";          // string specifying profile ranges
  string atmosfile        = "";          // stores the atmospheric profile name
  string atmosfileorder   = "";          // column order e.g. 'zuvwtdp'
  string atm_profile_dir  = "no_dir";    // the directory where atm profiles resides
  string wind_units       = "mpersec";   // m/s

  // obtain parameter values from the user's options
  atmosfile       = oNB->getAtmosfile();
  atmosfileorder  = oNB->getAtmosfileorder();
  wind_units      = oNB->getWindUnits();
  //gnd_imp_model  = oNB->getGnd_imp_model(); 
  atm_profile_dir = oNB->getAtm_profile_dir();
  bool inMPS = 0;
  if ( strcmp( wind_units.c_str(), "mpersec" ) == 0) {
    inMPS = 1;
  }
    
  
  int    Nz_grid        = oNB->getNz_grid();        // number of points on the z-grid [20000]
  int    Nrng_steps     = oNB->getNrng_steps();     // number of range steps [1000]
  int    filetype       = oNB->getFiletype();       // type: atmosfile, g2senvfile, etc.
  int    skiplines      = oNB->getSkiplines();      // skiplines in "atmosfile" [0]
  bool   write_2D_TLoss = oNB->getWrite_2D_TLoss(); // flag to output 2D transmission loss [0]
  
  //double freq           = oNB->getFreq();           // Hz [0.5] 
  double z_min          = oNB->getZ_min();          // meters MSL [0]
  double maxheight      = oNB->getMaxheight();      // meters MSL [150000]
  double sourceheight   = oNB->getSourceheight();   // meters AGL [0]
  double receiverheight = oNB->getReceiverheight(); // meters AGL [0]
  
  int    i, j, m, Nm, Nm_prev, Nmax; 
  int    Nprofiles, n_zsrc, n_zrcv, stepj;
  double rng, rng_step, rng0, RR;
  double dz, dz_km;
  complex<double> I (0.0, 1.0);
  complex<double> eIpir;
  SolveModNB *a;
  SampledProfile  *atm_profile;
  PetscErrorCode ierr;
  PetscMPIInt    rank, size;
  
  // populate vector Rv holding the ranges at which the atm profiles are requested
  // note that by convention R[0]=R[1]
  // also note that the ranges at which the profiles are available may not 
  // coincide with the ones requested in Rv. This implementation chooses 
  // the left-closest profile to the Rv elements. See diagram below
  /*
  // The 2D domain is partitioned into regions. Region jth ends at range Rj.
  // Each region has its atm profile defined at left side of the region 
  // e.g. the atm profile in Region 2 is the one defined at r = R1 or 
  // the left-closest to R1.  
  
  |   Region     |   Region     |      Region      |         |  Region  |
  |              |              |                  |         |          |
  |      1       |      2       |        3         |         |    N+1   |
  |              |              |                  |   ...   |          |
  |              |              |                  |         |          |
  |              |              |                  |         |          |
  0----x----x----R1-----x-------R2-----------------R3- ... --RN--------max
       r1   r2          rn                                             range
       
   R0 is a special case and R0=R1 therefore vector Rv contains: R1, R1, R2, R3, ..., RN+1
   We compute the pressure at ranges r1, r2, r3 , ... rM<= max range
  */

  vector<double> Rv(10,0.0);
  Nprofiles = 0;
  prf_ranges_km = oNB->getProfileRanges();
                      
  //cout << "oNB->getProfile_ranges_given_flag() = " << oNB->getProfile_ranges_given_flag() << endl;
  //cout << "prf_ranges_km: " << prf_ranges_km << endl;
  //cout << "oNB->getMaxrange() = " << oNB->getMaxrange() << endl;
  //cout << "oNB->getReq_profile_step() = " << oNB->getReq_profile_step() << endl;
                        
  getRegionBoundaries(oNB->getProfile_ranges_given_flag(), oNB->getMaxrange(), \
                      oNB->getReq_profile_step(), prf_ranges_km, Nprofiles, &Rv);                      
  
  //cout << "In ModessRD_main.cpp: " << endl;
  //cout << "Nprofiles = " << Nprofiles << endl;
  //cout << "Rv size   = " << Rv.size() << endl;
  //for (unsigned int i=0; i<Rv.size(); i++) {
  //  printf("Rv[%i] = %g\n", i, Rv[i]);
  //}

  // Initialize Slepc
  SlepcInitialize(PETSC_NULL,PETSC_NULL,(char*)0,PETSC_NULL);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  //
  // Process Region 1 - r in [0, R1]
  // 
  if (!atm_profile_dir.empty()) { // get ascii 1D profiles from files in directory 'atm_profile_dir'
      atm_profile = get_RngDepndProfiles_ascii(1, atmosfileorder, skiplines, atm_profile_dir, "profile", inMPS );
  }
  else if (filetype==1){ // get profiles from the .env file
      atm_profile = get_RngDepnd_profile(atmosfile, 0.0); // the very first profile
  }
  else if (filetype==0) { // get a single ascii file - this will just force a range-independent run
    atm_profile = new SampledProfile( atmosfile, atmosfileorder.c_str(), skiplines, inMPS );
  }
  else {
      throw invalid_argument( "Unrecognized atmospheric file type: try --g2senvfile, --use_1D_profiles_from_dir or --atmosfile" );
  }

  printf("\nRegion 1 (%g to %g km)\n", 0.0, Rv[1]/1000.0);
  
  //
  // !!! ensure maxheight is less than the max height covered by the provided atm profile
  // may avoid some errors asociated with the code thinking that it goes above 
  // the max height when in fact the height may only differ from max height by
  // a rounding error. This should be revisited.
  //
  if (maxheight/1000.0 >= atm_profile->z(atm_profile->nz()-1)) {
      maxheight = (atm_profile->z(atm_profile->nz()-1) - 1e-9)*1000.0; // slightly less
		  cout << "\nmaxheight adjusted to: " << maxheight 
		       << " meters (max. available in the profile file)" << endl;
  }
  
  rng_step = oNB->getMaxrange()/Nrng_steps;         // range step [meters]
  dz = (maxheight - z_min)/Nz_grid;	      // the z-grid spacing
  n_zsrc = (int) floor(sourceheight/dz);  // index of source on z-grid
  n_zrcv = (int) receiverheight/dz;       // index of receiver on z-grid
  
  // get the atmospheric profile in Region 1
  a = new SolveModNB(oNB, atm_profile);
  
  // print the parameters so far
  //oNB->printParams();
  a->printParams();
  if (filetype==2) {
      printf("atmospheric profile dir : %s\n", atm_profile_dir.c_str());
  } else {
      printf("    atmospheric profile : %s\n", atmosfile.c_str());
  }                 				       					 
  // compute modes				 
  a->computeModes();

  // number of modes; remember Nm may not stay constant over all profiles
  Nm      = a->getNumberOfModes();
  Nm_prev = Nm;
  
  // set maximum number of modes; say, some factor more than what the first profile gives us.
  //Nmax = Nm + (int) round(0.5*Nm); 
  Nmax = min((Nm + 3*Nm), MAX_MODES); 
    
  //
  // Copy the modes and wavenumbers
  //
  complex<double> PP, PP_ll, *A_curr, *A_prev, *k_prev, *kp;
  complex<double> *A_prev_ll, *A_curr_ll;
  double **v_prev, **vp;
  double *rho_prev;
  double sqrtrho_s, sqrtrho_r;
   
  rho_prev  = new double [Nz_grid];
  A_curr    = new complex<double> [Nmax];
  A_prev    = new complex<double> [Nmax];
  A_prev_ll = new complex<double> [Nmax];
  A_curr_ll = new complex<double> [Nmax];
  k_prev    = new complex<double> [Nmax];
  v_prev    = dmatrix(Nz_grid, Nmax);
  
  // zero entries in case the compiler doesn't do it
    for (m=0; m<Nmax; m++) {
      A_curr[m]    = 0.0; A_curr_ll[m] = 0.0;
      A_prev[m]    = 0.0; A_prev_ll[m] = 0.0;
      k_prev[m]    = 0.0;
  }
  
  for (j=0; j<Nz_grid; j++) {   // copy wavevectors
      for (m=0; m<Nmax; m++) { 
          v_prev[j][m] = 0.0;
      }
  }
  
  vp = a->getWavevectors();     // just copy pointers
  kp = a->getWavenumbers();
  
  //
  // Copy the "previous" state: wavevectors + wavenumbers
  //
  for (j=0; j<Nz_grid; j++) {   // copy wavevectors
      rho_prev[j] = atm_profile->rho(j*dz/1000.0)*1000.0; // store density in rho_prev
      for (m=0; m<Nm; m++) { 
          v_prev[j][m] = vp[j][m];
      }
  }
  
  for (m=0; m<Nm; m++) {        // copy wavenumbers
      k_prev[m] = kp[m];
  }

  //
  // compute A_m coefficients: eq. 8 in NMRD-OWCM
  //
  //printf("atm_profile->rho(sourceheight/1000.0)=%g\n", atm_profile->rho(sourceheight/1000.0));
  //double rho_s     = atm_profile->rho(sourceheight/1000.0)*1000.0;             // in kg/m^3
  sqrtrho_s = sqrt(atm_profile->rho(sourceheight/1000.0)*1000.0);   // in kg/m^3
  sqrtrho_r = sqrt(atm_profile->rho(receiverheight/1000.0)*1000.0); // in kg/m^3
  for (m=0; m<Nm; m++) {
      //A_prev[m] = exp(I*Pi/4.0)/sqrt(8*Pi*Rv[0])*exp(I*k_prev[m]*Rv[0])/sqrt(k_prev[m])*v_prev[n_zsrc][m]/sqrtrho_s; // see eq. 8 in DV notes NMRD-OWCM
      //A_prev_ll[m] = exp(I*Pi/4.0)/sqrt(8*Pi*Rv[0])*exp(I*real(k_prev[m])*Rv[0])/sqrt(real(k_prev[m]))*v_prev[n_zsrc][m]/sqrtrho_s; // lossless case
      
      A_prev[m] = exp(I*Pi/4.0)/sqrt(8*Pi*Rv[0])*exp(I*k_prev[m]*Rv[0])/sqrt(k_prev[m])*v_prev[n_zsrc][m]; // see eq. 8 in DV notes NMRD-OWCM
      A_prev_ll[m] = exp(I*Pi/4.0)/sqrt(8*Pi*Rv[0])*exp(I*real(k_prev[m])*Rv[0])/sqrt(real(k_prev[m]))*v_prev[n_zsrc][m]; // lossless case
  }

  //
  // compute/save the 1D pressure field (or TL) at points r in region 1 (from 0 to R[1])
  // note that rng is updated inside saveTLoss1D()
  //
  rng = rng_step; // initialize current (marching) range
  cout << "Writing 1D TL to file: tloss_rd_1d.lossless.nm" << endl;
  cout << "Writing 1D TL to file: tloss_rd_1d.nm" << endl;
  saveTLoss1D ( Nm, n_zsrc, sqrtrho_s, n_zrcv, sqrtrho_r, 1, Rv, rng_step, v_prev, k_prev, A_prev, A_prev_ll, "w", &rng);
  
  //
  // Compute 2D Transmission Loss if requested
  //
  dz_km = dz/1000.0;
  stepj = Nz_grid/200;  // controls the vertical sampling of 2D data saved
  if (write_2D_TLoss) {                   
      rng   = rng_step; // re-initialize current (marching) range
      if (stepj==0) {
          stepj = 1;	  // ensure it's never less than 1; necessary for the loop inside next function
      }
      cout << "Writing 2D transmission loss to file: tloss_rd_2d.nm." << endl;  
      saveTLoss2D( Nm, Nz_grid, n_zsrc, dz_km, stepj,sqrtrho_s, 1, Rv, \
                    rng_step, rho_prev, v_prev, k_prev, A_prev, A_prev_ll, "w", &rng);          
  }  

  delete atm_profile;
  delete a;
  
  // ------ End processing Region 1 --------------------------------------------

  //
  // BIG LOOP - for each atm profile in regions 2, 3...N
  //
  for (i=2; i<=Nprofiles; i++) {

      RR = Rv[i-1];
      printf("\nRegion %d (%g to %g km)\n", i, Rv[i-1]/1000.0, Rv[i]/1000.0);
      
      if (!atm_profile_dir.empty()) { // get ascii 1D profiles from files in directory 'atm_profile_dir'
          atm_profile = get_RngDepndProfiles_ascii(i, atmosfileorder, skiplines, atm_profile_dir, "profile", inMPS );
      }
      else if (filetype==1) { // get profiles from the .env file
          atm_profile = get_RngDepnd_profile(atmosfile, RR); //get left-closest profile in the .env file
      }
      else if (filetype==0) { // get a single ascii file - this will just force a range-independent run
          atm_profile = new SampledProfile( atmosfile, atmosfileorder.c_str(), skiplines, inMPS );
      }
      else {
          throw invalid_argument( "Unrecognized atmospheric file type: try --g2senvfile, --use_1D_profiles_from_dir or --atmosfile" );
      }

      // get normal mode solver object valid from Ri to R_i+1 
                          
      a = new SolveModNB(oNB, atm_profile);
                          				 				       					 
      // compute modes					 
      a->computeModes();
      Nm = a->getNumberOfModes();

      // Evaluate A_curr = R1*A_prev = 1/2(C_tilde + C_hat)*H1*A_prev    
      getAcurr( Nz_grid, Nm_prev, Nm, dz, A_prev, v_prev, k_prev, rho_prev, Rv[i-1], Rv[i-2], a->getWavevectors(), a->getWavenumbers(), atm_profile, A_curr); 
      
      getAcurr_ll( Nz_grid, Nm_prev, Nm, dz, A_prev_ll, v_prev, k_prev, rho_prev, Rv[i-1], Rv[i-2], a->getWavevectors(), a->getWavenumbers(), atm_profile, A_curr_ll);      
      
      // we have A_curr; now store the current wavenumbers and wavevectors
      // into the "previous" state - and get ready for the next iteration
      vp = a->getWavevectors();     // just copy pointers
      kp = a->getWavenumbers();
      
      //printf("copying wavevectors ...\n");      
             
      // copy current wavevectors to "previous" state
      for (j=0; j<Nz_grid; j++) {   
          rho_prev[j] = atm_profile->rho(j*dz/1000.0)*1000.0;
          for (m=0; m<Nm; m++) { 
              v_prev[j][m] = vp[j][m];
          }
      }
      
      //printf("copying wavenumbers ...\n"); 
      
      // copy wavenumbers
      for (m=0; m<Nm; m++) {        
          k_prev[m]    = kp[m];
          A_prev[m]    = A_curr[m];  A_prev_ll[m] = A_curr_ll[m];
      }
      
      // update Nm_prev
      //if (Nm>Nm_prev) { Nm_prev = Nm; }
      Nm_prev = Nm;     
      
      
      /*
      printf("zero the remaining elements (from Nm to Nmax) of k_prev, A_prev, A_curr, v_prev...\n"); 
      // zero the remaining elements (from Nm to Nmax) of k_prev, A_prev, A_curr, v_prev
      for (m=Nm; m<Nmax; m++) {
          k_prev[m]    = 0.0;      
          A_curr[m]    = 0.0; A_prev[m]    = 0.0;
          A_curr_ll[m] = 0.0; A_prev_ll[m] = 0.0;
      }
      
  
      for (j=0; j<Nz_grid; j++) {
          for (m=Nm; m<Nmax; m++) { 
              v_prev[j][m] = 0.0;
          }
      }
      */
 
      //
      // Compute/save the 1D pressure field (or TL); note that rng is updated inside saveTLoss1D()
      //
      rng0 = rng;  // remember the starting range
      //printf("going to saveTLoss1D\n");
      saveTLoss1D ( Nm, n_zsrc, sqrtrho_s, n_zrcv, sqrtrho_r, i, Rv, rng_step, v_prev, k_prev, \
                    A_prev, A_prev_ll, "a", &rng);

      //
      // Compute 2D Transmission Loss if requested
      //
      if (write_2D_TLoss) { 
          rng = rng0; // reset rng                      
          saveTLoss2D( Nm, Nz_grid, n_zsrc, dz_km, stepj,sqrtrho_s, i, Rv, \
                       rng_step, rho_prev, v_prev, k_prev, A_prev, A_prev_ll, "a", &rng);      
      }                

      delete atm_profile;
      if (i==Nprofiles) { // print parameters after the last iteration
          a->printParams();
          if (filetype==2) {
              printf("atmospheric profile dir : %s\n", atm_profile_dir.c_str());
          } else {
              printf("    atmospheric profile : %s\n", atmosfile.c_str());
          } 
      }
      delete a;
  } // End of BIG LOOP
  
  cout << "File tloss_rd_1d.lossless.nm written." << endl;
  cout << "File tloss_rd_1d.nm written." << endl;
  if (write_2D_TLoss) { 
      cout << "File tloss_rd_2d.nm written." << endl;
  }

  ierr = SlepcFinalize();CHKERRQ(ierr);	// Finalize Slepc.
  
  //// plot with gnuplot - probably not part of the final release; 
  //// good for a quick check of results				 
  //if (oNB->getPlot_flg()) {
	//    plotwGNUplot(freq, write_2D_TLoss);
  //}		
   
  //oNB->printParams();
  
  delete [] A_curr;    delete [] A_prev;
  delete [] A_curr_ll; delete [] A_prev_ll;
  delete [] k_prev;    delete [] rho_prev;
  free_dmatrix(v_prev, Nz_grid, Nmax);
  delete opt;
  delete oNB;

  cout << "Maximum number of modes used: " << Nm_prev << endl;
  cout << "Number of atm. profiles used in this run: " << Nprofiles << endl;
  cout << "\n ... main() is done." << endl;
  time_t tm2 = time(NULL);
  cout << "Run duration: " << difftime(tm2,tm1) << " seconds." << endl;	

  return (0);
} // -------------- end of main(); -----------------------------------------


//
// Function to parse the input options (both command lines and in the options file ModESSRD1WCM.options)
//
AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();

  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "|                             NCPA Infrasound                              |" );  
  opt->addUsage( "|               Normal Modes for Range-Dependent Environments              |" );
  opt->addUsage( "|                          One-Way Coupled Modes                           |" );  
  opt->addUsage( "|           Single Frequency; Effective Sound Speed Approximation          |" );
  opt->addUsage( "----------------------------------------------------------------------------" );	
  opt->addUsage( "Usage: " );
  opt->addUsage( "By default the program computes the 1D transmission loss (TL)" );
  opt->addUsage( "at the ground or the specified receiver height and saves the data to 2 files:" );
  opt->addUsage( "   file tloss_rd_1d.nm - considering attenuation in the atmosphere" );
  opt->addUsage( "   file tloss_rd_1d.lossless.nm  - no attenuation" );
	opt->addUsage( "Additionally, if the flag --write_2D_TLoss is present on the command line the 2D TL is saved to file tloss_rd_2d.nm" );
  opt->addUsage( "The options below can be specified in a colon-separated file \"ModESSRD1WCM.options\" or at the command line.  Command-line options override file options." );
  opt->addUsage( " --help -h                Print this message and exit" );
  opt->addUsage( "" );
  opt->addUsage(  " The atmosphere can be specified from one of 3 different sources:");
  opt->addUsage( "    1. An .env file containing the atmospheric specifications at certain ranges:" );
  opt->addUsage( "       use option --g2senvfile <filename>" );
  opt->addUsage( "    2. Several ASCII files stored in a given directory:" );
  opt->addUsage( "       use option --use_1D_profiles_from_dir <mydirname>" );
  opt->addUsage( "    3. A single ASCII file. This will just force a range-independent run." );
  opt->addUsage( "       use option --atmosfile <filename>" );  
  opt->addUsage( "" );  
  opt->addUsage( "The options available are:" );	
  opt->addUsage( "" );	
  opt->addUsage( "REQUIRED (no default values):" );
  opt->addUsage( " --atmosfileorder         The order of the (z,t,u,v,w,p,d) fields in the file");
  opt->addUsage( "                          (Ex: 'ztuvpd')" );
  opt->addUsage( " --skiplines              Lines at the beginning of the ASCII file to skip" );
  opt->addUsage( " --azimuth                Degrees in range [0,360), clockwise from North" );
  opt->addUsage( " --freq                   Frequency [Hz]" );

  opt->addUsage( " --g2senvfile <filename>  Uses an .env binary file (for range-dependent code)" );   
  opt->addUsage( " --use_1D_profiles_from_dir" );
  opt->addUsage( "                          e.g. --use_1D_profiles_from_dir <myprofiles>" );
  opt->addUsage( "                          This option allows to use the ascii profiles stored in" );
  opt->addUsage( "                          the specified directory. The profiles must have names" );
  opt->addUsage( "                          'profiles0001', 'profiles0002', etc. and will be" );
  opt->addUsage( "                          used in alphabetical order at the provided ranges" );
  opt->addUsage( "                          e.g. in conjunction with either" );
  opt->addUsage( "                          option  '--use_profile_ranges_km' " );
  opt->addUsage( "                          or option '--use_profiles_at_steps_km'" );
  opt->addUsage( "                          If there are more requested ranges than existing" );
  opt->addUsage( "                          profiles then the last profile is used repeatedly" );
  opt->addUsage( "                          as necessary." );  
  opt->addUsage( "    Example: >> ../bin/ModESSRD1WCM --atmosfileorder zuvwtdp --skiplines 1" );
  opt->addUsage( "                --azimuth 90 --freq 0.1 --use_1D_profiles_from_dir myprofiles" );
  opt->addUsage( "                --use_profile_ranges_km 100_300_500_600_700" );
  opt->addUsage( "" );  
  opt->addUsage( " --atmosfile  <filename>  Uses an ASCII atmosphere file." ); 
  opt->addUsage( "                          In this case the run will just be range-independent" );   	
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
  
  opt->addUsage( " --use_profile_ranges_km" );
  opt->addUsage( "                          e.g. --use_profile_ranges_km  20_50_80.5_300     " );   
  opt->addUsage( "                          The profiles at certain ranges specified by numbers" );
  opt->addUsage( "                          (in km) in a string such as 20_50_80.5_300 are");
  opt->addUsage( "                          requested in the propagation. Note that underscores" );
  opt->addUsage( "                          are necessary to separate the numbers." );
  opt->addUsage( "                          Note also that these are requested ranges;" );
  opt->addUsage( "                          however the left-closest profile available" );
  opt->addUsage( "                          in the .env file will actually be used; " );
  opt->addUsage( "                          for instance we request the profile at 300 km " );
  opt->addUsage( "                          but in the .env file the left-closest profile" );
  opt->addUsage( "                          may be available at 290 km and it is the one used." );  
  opt->addUsage( "" ); 
  opt->addUsage( " --use_profiles_at_steps_km" );
  opt->addUsage( "                          e.g. --use_profiles_at_steps_km 100" );
  opt->addUsage( "                          The profiles are requested at equidistant intervals " );
  opt->addUsage( "                          specified by this option [1000]" );
  opt->addUsage( "" );                        
  opt->addUsage( "" );  
  opt->addUsage( "FLAGS (no value required):" );
  opt->addUsage( " --write_2D_TLoss         Outputs the 2D transmission loss to" );
  opt->addUsage( "                          default file: tloss_rd_2D.nm" );
  opt->addUsage( " --turnoff_WKB            Turn off the WKB least phase speed estimation" );
  opt->addUsage( "                          an approx. that speeds-up ground-to-ground propag." ); 
  opt->addUsage( "                          It has the value 1 (true) if any of the flags" );
  opt->addUsage( "                          write_2D_TLoss, write_phase_speeds, write_modes" );
  opt->addUsage( "                          or write_dispersion are true." );  
  opt->addUsage( "" );
  opt->addUsage( "" );
  opt->addUsage( " The format of the output files are as follows (column order):" );
  opt->addUsage( "  tloss_rd_1d.nm:           r, 4*PI*Re(P), 4*PI*Im(P)" );
  opt->addUsage( "  tloss_rd_1d.lossless.nm:" );
  opt->addUsage( "  tloss_rd_2d.nm:           r, z, 4*PI*Re(P), 4*PI*Im(P)" );   
  opt->addUsage( "" );
  opt->addUsage( "  Examples to run (from 'samples' directory):" );
  opt->addUsage( "      ../bin/ModessRD1WCM --use_1D_profiles_from_dir profiles --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --use_profile_ranges_km 100_200" );   
  opt->addUsage( "" );  
  opt->addUsage( "      ../bin/ModessRD1WCM --use_1D_profiles_from_dir profiles --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --use_profiles_at_steps_km 200" );
  opt->addUsage( "" );
  opt->addUsage( "      ../bin/ModessRD1WCM --g2senvfile g2sgcp2011012606L.jordan.env --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --use_profile_ranges_km 50_100_150_200_250_300" );
  opt->addUsage( "" );
  opt->addUsage( "      ../bin/ModessRD1WCM  --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1" );
  opt->addUsage( "" );
  opt->addUsage( "  This last example is in fact not a range-dependent case since it just loads a single 1D atmospheric profile." ); 
  opt->addUsage( "" );
  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  opt->setFlag( "write_2D_TLoss" );
  opt->setFlag( "turnoff_WKB");
  opt->setFlag( "plot" );

  opt->setOption( "atmosfile" );
  opt->setOption( "atmosfileorder" );
  opt->setOption( "wind_units" );
  opt->setOption( "use_1D_profiles_from_dir" );
  opt->setOption( "slicefile" );  
  opt->setOption( "g2senvfile" ); 
  opt->setOption( "skiplines" );		
  opt->setOption( "azimuth" );
  opt->setOption( "freq" );
  opt->setOption( "maxrange_km" );
  opt->setOption( "sourceheight_km" );
  opt->setOption( "receiverheight_km" );
  opt->setOption( "maxheight_km" );
  opt->setOption( "zground_km" );  
  opt->setOption( "stepsize" );
  opt->setOption( "Nz_grid" );
  opt->setOption( "Nrng_steps" );
  opt->setOption( "ground_impedance_model" );
  opt->setOption( "Lamb_wave_BC" );
  opt->setOption( "use_profile_ranges_km" ); 
  opt->setOption( "use_profiles_at_steps_km" );
  opt->setOption( "use_attn_file" ); 

  // Process the command-line arguments
  opt->processFile( "./ModessRD1WCM.options" );
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


 
