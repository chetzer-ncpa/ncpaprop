#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <complex>

#include <sys/stat.h>
#include <ctime>

#include "Atmosphere.h"
#include "anyoption.h"
#include "ProcessOptionsNBRDCM.h"
#include "SolveModNBRDCM.h"
#include "ModessRDCM_lib.h"

// On Doru's computer only:
//Note 1!!!
// use these includes (with <slepceps.h>, <slepcst.h>)
//#include <slepceps.h>
//#include <slepcst.h>
// with petsc and slepc versions 3.3

//Note 2!!!
// use these includes (with "slepceps.h", "slepcst.h")
#include "slepceps.h"
#include "slepcst.h"
// with petsc and slepc versions 3.2
// /home/doru/lib/petsc-3.2-p5/arch-linux2-c-debug/lib
// /home/doru/lib/slepc-3.2-p3/arch-linux2-c-debug/lib

#include "binaryreader.h"
//#include <petscksp.h>

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

/*
 * Range-Dependent Normal Modes - (Two-Way) Coupled Modes - Effective Sound Speed.
 * 
 * @version 0
 * @date 2012-12
 * @authors Doru Velea; Jelle Assink; Roger Waxler; Claus Hetzer; 
 * Follows Doru Velea's notes: Normal Modes for Range-Dependent Environments: 
 * Two-Way Coupled Modes i.e. NMRD-CM
 * See also the theory in "Computational Ocean Acoustics" Section 5.9, page 315, 1994 ed.
 * 
 * Changelog:
 * 201306  : DV modified the get() functions that return strings
 */

// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

//
// main
//
int main( int argc, char **argv ) {

  // parse options from the command line as well as an options file
  AnyOption *opt = parseInputOptions( argc, argv ); 
  
  // object to process the options
  ProcessOptionsNB *oNB;
  oNB = new ProcessOptionsNB(opt);
  
  //
  // Physical values are usually in SI (System International) units 
  // unless mentioned otherwise.
  // The defaults below are typical. They are reset by the options provided at run time.
  //
  double freq             = 0.5;         // Hz
  double azi              = 0.0;         // degrees
  double maxrange         = 1.0E6;       // meters
  double maxheight        = 150000.0;    // meters
  double z_min            = 0.0;         // meters
  double sourceheight     = 0.0;         // meters
  double receiverheight   = 0.0;         // meters
  double tol              = 1.0E-8;      // tolerance for Slepc calculations	
  int   Nz_grid           = 20000;       // number of points on the z-grid
  int   Nrng_steps        = 1000;        // number of range steps	 
  int   skiplines         = 0;           // skiplines in "atmosfile"
  int   Lamb_wave_BC      = 0;           // 1 to enforce the Lamb wave BC
  string atmosfile        = "";          // stores the atmospheric profile name
  string atm_profile_dir  = "no_dir";    // the directory where atm profiles resides
  string wind_units       = "mpersec";   // m/s
  string gnd_imp_model    = "rigid";     // rigid ground
  
  string prf_ranges_km;                  // string specifying profile ranges
  double req_profile_step = maxrange;    // distance wanted between new atm profiles
  
	
  // The default output is the 1D Transmission Loss at the ground
  // set defaults for some flags determining the various other outputs
  bool write_2D_TLoss     = 0;           // flag to output 2D transmission loss
  bool write_phase_speeds = 0;           // write phase speeds flag
  bool write_modes        = 0;           // write modes flag
  bool write_dispersion   = 0;           // write wavenumbers flag

  SampledProfile  *atm_profile;
  string           atmosfileorder;
  
  // set up to measure the duration of this run
  std::time_t tm1 = std::time(NULL);

  // obtain the parameter values from the user's options
  atmosfile       = oNB->getAtmosfile();
  atmosfileorder  = oNB->getAtmosfileorder();
  gnd_imp_model   = oNB->getGnd_imp_model(); 
  atm_profile_dir = oNB->getAtm_profile_dir(); 
  wind_units      = oNB->getWindUnits();

  skiplines          = oNB->getSkiplines();
  z_min              = oNB->getZ_min();
  freq               = oNB->getFreq();  
  azi                = oNB->getAzimuth();
  maxrange           = oNB->getMaxrange();
  maxheight          = oNB->getMaxheight();
  sourceheight       = oNB->getSourceheight();
  receiverheight     = oNB->getReceiverheight();
  tol				   			 = oNB->getSlepcTolerance();    
  Nz_grid            = oNB->getNz_grid();
  Nrng_steps         = oNB->getNrng_steps();
  Lamb_wave_BC       = oNB->getLamb_wave_BC();
  req_profile_step   = oNB->getReq_profile_step(); //default is maxrange

  write_2D_TLoss     = oNB->getWrite_2D_TLoss();
  write_phase_speeds = oNB->getWrite_phase_speeds();
  write_modes        = oNB->getWrite_modes();
  write_dispersion   = oNB->getWrite_dispersion();
  
  
  int    j, m, n, n_z, Nmin, Nmax, Noptmax; 
  int    Nprofiles, n_zsrc;
  double rng_curr, rng_step, rng0;
  double dz, dz_km, sqrtrho_z;
  complex<double> I (0.0, 1.0);
  complex<double> eIpir;
  SolveModNBRDCM *a;
  PetscErrorCode ierr;
  PetscMPIInt    rank, size;
  
  
  rng_step = maxrange/Nrng_steps;           // range step [meters]
  dz = 0; // just an initialization //dz       = (maxheight - z_min)/Nz_grid;	  // the z-grid spacing
  n_zsrc   = (int) floor(sourceheight/dz);  // index of source on z-grid
  
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
  
  oNB->printParams();
  
  vector<double> Rv(10,0.0);
  Nprofiles = 0;
  prf_ranges_km = oNB->getProfileRanges();
  getRegionBoundaries(oNB->getProfile_ranges_given_flag(), maxrange, \
                      req_profile_step, prf_ranges_km, &Nprofiles, &Rv);
  
  //cout << "Nprofiles = " << Nprofiles << endl;
  //cout << "Rv size   = " << Rv.size() << endl;
  //for (i=0; i<Rv.size(); i++) {
  //  printf("Rv[%i] = %g\n", i, Rv[i]);
  //} 
 
  if (Rv.size()==2) {
      cout << "\nOnly one region is defined - this is a range-independent case ...\n";
  }

  string subdir, filename;
  stringstream oss;
  string filen_stub = "eigvalvecs"; // will save to eigvalvecs1,2,3...
  //
  // data will be saved in a new directory using date_time format
  // YYYY_MM_DD_HH_MM_SS
  //
  makeYYYY_MM_DD_subdir(&subdir);

  // Initialize Slepc
  SlepcInitialize(PETSC_NULL,PETSC_NULL,(char*)0,PETSC_NULL);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  Nmax    = 0;
  Noptmax = 0;
  Nmin    = MAX_MODES;
  

  // first pass to obtain global ceffmin, max
  
  double ceffMin, ceffMax, k_min, k_max;
  getGlobalceffMinMax(Nprofiles, Rv, azi, atm_profile_dir, atmosfile, atmosfileorder, skiplines, &ceffMin, &ceffMax);
  k_min = 2.0*Pi*freq/ceffMax;
  k_max = 2.0*Pi*freq/ceffMin;

  // a vector to store the optimal number of modes in each region
  vector<int> Nmopt (Nprofiles, 0); // initialize to zero
  // a vector to store the computed number of modes in each region
  vector<int> Nj (Nprofiles, 0); // initialize to zero
  bool pflg = 1;

  //
  // loop over all regions to save the wavenumbers and modes
  //
  for (j=1; j<=Nprofiles; j++) {
      if (j==1) {
          printf("\nRegion %d (%g to %g km)\n", j, 0.0, Rv[j]/1000.0);
      }
      else {
          printf("\nRegion %d (%g to %g km)\n", j, Rv[j-1]/1000.0, Rv[j]/1000.0);
      }
      if (!atm_profile_dir.empty()) { // get ascii 1D profiles from files in directory 'atm_profile_dir'
          atm_profile = get_RngDepndProfiles_ascii(j, atmosfileorder, skiplines, atm_profile_dir, "profile", pflg);
          pflg = 0;
      }
      else { // get profiles from the .env file
          if (j==1) {
              atm_profile = get_RngDepnd_profile(atmosfile, 0.0); // if the very first profile is required in the first region
          }
          else {
              atm_profile = get_RngDepnd_profile(atmosfile, Rv[j-1]);
          }
      }
                          
      a = new SolveModNBRDCM( oNB, atm_profile);                         

      a->printParams();
      if (!atm_profile_dir.empty()) {
          printf("atmospheric profile dir : %s\n", atm_profile_dir.c_str());
      } else {
          printf("    atmospheric profile : %s\n", atmosfile.c_str());
      } 
      // compute modes; we request all modes with speeds between omega/k_max and omega/k_max m/s
      // then we will use only the necessary number of modes for all profiles				 
      a->computeModes(k_min, k_max); 

      Nj[j-1] = a->getNumberOfModes();

      if (Nmin>Nj[j-1]) {
          Nmin = Nj[j-1]; // track min of computed modes over all regions
      } 
      
      if (Nmax<Nj[j-1]) {
          Nmax = Nj[j-1]; // track max of computed modes over all regions
      }            

      Nmopt[j-1] = a->getOptimalNumberOfModes();
      if (Noptmax<Nmopt[j-1]) {
          Noptmax = Nmopt[j-1]; // track max of optimal number of modes over all regions
      }	

      //     
      // save modes to files: eigvalvecs_j.dat
      //      
      oss << subdir << "/" << filen_stub << "_" << j << ".dat"; // full name for the file string eigs
      a->writeEigenValVecs(oss.str(), a->getNumberOfModes());
      //cout << "modes saved to files: " << oss.str() << endl;  
      oss.str(""); oss.clear();  // flush/prepare oss to be rewritten
      
      delete a;

  } // end loop computing/saving the eigenval/vecs

  // Nmin is the number of modes to be used = maximum of the optimal number of modes
  if (Noptmax <=Nmin) { Nmin = Noptmax; }
  //printf("Nmin = %d\n", Nmin);
  //printf("Noptmax = %d\n", Noptmax);

  complex<double> *k_curr, *k_next;
  double *rho_curr, *rho_next;
  double **v_curr, **v_next;
  int Nm_curr = 0; 
  int Nm_next = 0;
  
  // R matrix + lossless version
  complex<double> **Rmx, **Rmx_ll;
  Rmx    = cmatrix(2*Nmin, 2*Nmin);
  Rmx_ll = cmatrix(2*Nmin, 2*Nmin);
  
  // assemble big matrix S
  complex<double> **S_curr, **S_next, **S_curr_ll, **S_next_ll;
  S_curr    = cmatrix(2*Nmin, 2*Nmin);
  S_next    = cmatrix(2*Nmin, 2*Nmin);
  S_curr_ll = cmatrix(2*Nmin, 2*Nmin);
  S_next_ll = cmatrix(2*Nmin, 2*Nmin);
  
  // fill with zeros in case the compiler doesn't do it
  for (j=0; j<2*Nmin; j++) {   
      for (m=0; m<2*Nmin; m++) {
          S_curr[j][m]    = 0.0 + 0.0*I;
          S_curr_ll[j][m] = 0.0 + 0.0*I;
      }
  }
  
  // set big matrix S_curr to the identity matrix
  for (j=0; j<2*Nmin; j++) {   
      S_curr[j][j]    = 1.0;
      S_curr_ll[j][j] = 1.0;
  }

  // 
  // read the data in the Region 1
  //
  oss << subdir << "/" << filen_stub << "_" << 1 << ".dat"; // full name for the file string eigs
  getNmodesNz(oss.str(), &Nm_curr, &Nz_grid, &dz_km);       // read from eigenvalvec_1.dat file
  dz = dz_km*1000;

  // allocate containers 
  rho_curr  = new double [Nz_grid];
  k_curr    = new complex<double> [Nmin];
  v_curr    = dmatrix(Nz_grid, Nmin);
  
  rho_next  = new double [Nz_grid]; 
  k_next    = new complex<double> [Nmin];
  v_next    = dmatrix(Nz_grid, Nmin);
  
  readEigenValVecs(oss.str(), k_curr, rho_curr, v_curr, Nmin); 
  
  // Compute the diagonal matrix D and column vector ss
  //double sqrtrho_s = sqrt(rho_curr[n_zsrc]);
  complex<double> *DD, *ss;
  complex<double> *DD_ll, *ss_ll;  // lossless case
  DD    = new complex<double> [Nmin];
  ss    = new complex<double> [Nmin];
  DD_ll = new complex<double> [Nmin];
  ss_ll = new complex<double> [Nmin];  
  
  for (m=0; m<Nmin; m++) {
      DD[m]    = exp(2.0*I*(k_curr[m]*Rv[0] - Pi/4.0)); // e.g. eq 5.171 in Oc. Ac.
      DD_ll[m] = exp(2.0*I*(real(k_curr[m])*Rv[0] - Pi/4.0));
      ss[m]    = I*exp(-I*Pi/4.0)/(sqrt(8*Pi*k_curr[m]*Rv[0]))*v_curr[n_zsrc][m]*exp(I*k_curr[m]*Rv[0]); // note: the modes in the Oc. Acoust. book = our modes times sqrt(rho)
      ss_ll[m] = I*exp(-I*Pi/4.0)/(sqrt(8*Pi*real(k_curr[m])*Rv[0]))*v_curr[n_zsrc][m]*exp(I*real(k_curr[m])*Rv[0]); // note: the modes in the Oc. Acoust. book = our modes times sqrt(rho) 
      
      //// use the next 2 lines only if the modes need scaling by sqrt(rho)
      //ss[m]    = I*exp(-I*Pi/4.0)/(rho_curr[n_zsrc]*sqrt(8*Pi*k_curr[m]*Rv[0]))*sqrtrho_s*v_curr[n_zsrc][m]*exp(I*k_curr[m]*Rv[0]); // note: the modes in the Oc. Acoust. book = our modes times sqrt(rho)
      //ss_ll[m] = I*exp(-I*Pi/4.0)/(rho_curr[n_zsrc]*sqrt(8*Pi*real(k_curr[m])*Rv[0]))*sqrtrho_s*v_curr[n_zsrc][m]*exp(I*real(k_curr[m])*Rv[0]); // note: the modes in the Oc. Acoust. book = our modes times sqrt(rho)      
  }

  //
  // loop over regions 2-end to save R matrices and update S
  //
  for (n=2; n<=Nprofiles; n++) {
      oss.str(""); oss.clear();
      oss << subdir << "/" << filen_stub << "_" << n << ".dat"; // full name for the file string eigs
   
      getNmodesNz(oss.str(), &Nm_next, &Nz_grid, &dz_km);
      dz = dz_km*1000;

      readEigenValVecs(oss.str(), k_next, rho_next, v_next, Nmin);
      
      //
      // evaluate matrix Rmx and save it in a file  Rmat_%d.mat
      //     
      oss.str(""); oss.clear();
      oss << subdir << "/Rmat_" << n-1 << ".dat";     
      getRRmats(oss.str(), Rv[n-2], Rv[n-1], Nmin, Nz_grid, dz,  \
            k_curr, rho_curr, v_curr, k_next, rho_next, v_next, \
            Rmx);
 
      // lossless version         
      oss.str(""); oss.clear();
      oss << subdir << "/Rmat_ll_" << n-1 << ".dat";
      getRRmats_ll(oss.str(), Rv[n-2], Rv[n-1], Nmin, Nz_grid, dz,  \
            k_curr, rho_curr, v_curr, k_next, rho_next, v_next, \
            Rmx_ll);                          
             
      // update S_next = Rmx*S_curr
      MatCMultiply(Rmx, S_curr, S_next, 2*Nmin, 2*Nmin, 2*Nmin, 2*Nmin);  
      MatCMultiply(Rmx_ll, S_curr_ll, S_next_ll, 2*Nmin, 2*Nmin, 2*Nmin, 2*Nmin);  
      
      // copy S_next into S_curr to prepare for the next iteration
      for (j=0; j<2*Nmin; j++) {   
          for (m=0; m<2*Nmin; m++) {
              S_curr[j][m] = S_next[j][m];
              S_curr_ll[j][m] = S_next_ll[j][m];
          }
      }
      
      // update k_curr, rho_curr, v_curr      
      for (m=0; m<Nmin; m++) {
          k_curr[m]   = k_next[m];  
      } 
      
      for (j=0; j<Nz_grid; j++) {
          rho_curr[j] = rho_next[j];
          for (m=0; m<Nmin; m++) {
              v_curr[j][m] = v_next[j][m];
          }
      }                                   
  } // end big loop computing matrices R and S

  //
  // solve (S4+S3*D)b1 = -S3*ss - eq 3 in Evans JASA 74 (1) 1983 pg 191
  //
  // form the matrix A = S4+S3*D
  complex<double> **S4plusS3D, *S3ss, *ab_curr, *ab_next, *bj;
  S3ss      = new complex<double> [Nmin];
  bj        = new complex<double> [Nmin];
  ab_curr   = new complex<double> [2*Nmin];
  ab_next   = new complex<double> [2*Nmin];
  S4plusS3D = cmatrix(Nmin, Nmin);
  
  complex<double> **S4plusS3D_ll, *S3ss_ll, *ab_curr_ll, *ab_next_ll, *bj_ll;
  S3ss_ll      = new complex<double> [Nmin];
  bj_ll        = new complex<double> [Nmin];
  ab_curr_ll   = new complex<double> [2*Nmin];
  ab_next_ll   = new complex<double> [2*Nmin];
  S4plusS3D_ll = cmatrix(Nmin, Nmin);

  for (j=0; j<Nmin; j++) {
  	S3ss[j]    = 0.0;
  	S3ss_ll[j] = 0.0;
    for (m=0; m<Nmin; m++) {
        S4plusS3D[j][m]    = S_curr[j+Nmin][m+Nmin] + S_curr[j+Nmin][m]*DD[m]; // = S4+S3*D
        S3ss[j]            = S3ss[j] - S_curr[j+Nmin][m]*ss[m]; // -S3*ss
        S4plusS3D_ll[j][m] = S_curr_ll[j+Nmin][m+Nmin] + S_curr_ll[j+Nmin][m]*DD_ll[m]; // = S4+S3*D
        S3ss_ll[j]         = S3ss_ll[j] - S_curr_ll[j+Nmin][m]*ss_ll[m]; // -S3*ss        
    }
  }
  
  //
  // solve for b1
  //
  PetscInitialize(NULL,NULL,(char *)0,(char *)0);
  SolveLinSys2(S4plusS3D,S3ss,bj, Nmin);
  SolveLinSys2(S4plusS3D_ll,S3ss_ll,bj_ll, Nmin);
  
  // ab_curr is the column vector ( aj )
  //                              ( bj )
  //
  for (m=0; m<Nmin; m++) {
      ab_curr[m] = ss[m] + DD[m]*bj[m]; //a1 = ss + D*b1
      ab_curr[m+Nmin] = bj[m];
      
      ab_curr_ll[m] = ss_ll[m] + DD_ll[m]*bj_ll[m]; //a1 = ss + D*b1
      ab_curr_ll[m+Nmin] = bj_ll[m];       
  }

  //
  // compute pressure and TL in Region 1
  //
  printf("Region 1 (0 to %g km): computing 1D pressure and TL\n", Rv[1]/1000);
  rng_curr = rng_step; // start at the first rng_step
  // form full file name storing eigenvals and vecs 
  oss.str(""); oss.clear(); // flush/prepare oss to be rewritten
  oss << subdir << "/" << filen_stub.c_str() << "_1.dat";

  // read from eigenvalvec_j.dat file
  readEigenValVecs(oss.str(), k_curr, rho_curr, v_curr, Nmin);

  n_z = (int) receiverheight/dz;
  sqrtrho_z = sqrt(rho_curr[n_z]);
  //printf("n_z=%d; n_zsrc=%d\n", n_z, n_zsrc);

  rng0 = rng_curr;
  computePressure_1D( 1, Nmin, n_z, Rv, rng_step, sqrtrho_z, k_curr, v_curr, ab_curr, ab_curr_ll, "w", &rng_curr);
  
  int stepn = Nz_grid/500;  // controls the vertical sampling of 2D data saved
  if (stepn==0) { stepn = 10; }  // ensure it's never 0; necessary for the loop inside next function
  if (write_2D_TLoss) {
      printf("Region 1 (0 to %g km): computing 2D pressure and TL\n", Rv[1]/1000);
      rng_curr = rng0; // reset rng_curr 
      computePressure_2D( 1, Nmin, Nz_grid, stepn, Rv, rng_step, dz, rho_curr, k_curr,  v_curr, ab_curr, ab_curr_ll, "w", &rng_curr);
  }

  printf("...done\n");
  // ------------ end computing in region 1 -----------------------------------

  //cout << "from ajbj find aj+1 bj+1 recursively" << endl;
  for (n=2; n<=Nprofiles; n++) {     
      oss.str(""); oss.clear(); // flush/prepare oss to be rewritten
      oss << subdir << "/Rmat_" << n-1 << ".dat";

      readRRmat(oss.str(), Rmx, 2*Nmin);
      MatVecCMultiply(Rmx, ab_curr, ab_next, 2*Nmin, 2*Nmin, 2*Nmin);

      //printf("Region %d - a%d,b%d computed from R%d*(a%d,b%d)\n", n, n,n,n-1,n-1,n-1);
      // update ab_curr
      for (m=0; m<2*Nmin; m++) {
          ab_curr[m] = ab_next[m];
      }
      
      oss.str(""); oss.clear();
      oss << subdir << "/Rmat_ll_" << n-1 << ".dat";      
      readRRmat(oss.str(), Rmx_ll, 2*Nmin);
      MatVecCMultiply(Rmx_ll, ab_curr_ll, ab_next_ll, 2*Nmin, 2*Nmin, 2*Nmin);

      //printf("Region %d - a%d,b%d computed from R%d*(a%d,b%d)\n", n, n,n,n-1,n-1,n-1);
      // update ab_curr_ll
      for (m=0; m<2*Nmin; m++) {
          ab_curr_ll[m] = ab_next_ll[m];
      }      

      //
      // Compute/save the 1D pressure field (or TL); note that rng is updated inside saveTLoss1D()
      //
      oss.str(""); oss.clear();
      oss << subdir << "/" << filen_stub.c_str() << "_" << n << ".dat";
      // read from eigenvalvec_j.dat file
      readEigenValVecs(oss.str(), k_curr, rho_curr, v_curr, Nmin);
  
      // remember the starting range;
      rng0 = rng_curr;  // needed if we compute 2D TL where rng is reset to rng0
      
      printf("Region %d (%g to %g km): computing 1D pressure and TL\n",n, Rv[n-1]/1000, Rv[n]/1000);
      computePressure_1D( n, Nmin, n_z, Rv, rng_step, sqrtrho_z, \
                          k_curr, v_curr, ab_curr, ab_curr_ll, "a", &rng_curr);
                          
      if (write_2D_TLoss) {
          printf("Region %d (%g to %g km): computing 2D pressure and TL\n",n, Rv[n-1]/1000, Rv[n]/1000);
          rng_curr = rng0; // reset rng_curr 
          computePressure_2D( n, Nmin, Nz_grid, stepn, Rv, rng_step, dz, rho_curr, \
                              k_curr,  v_curr, ab_curr, ab_curr_ll, "a", &rng_curr);
      }                          
      printf("...done\n");
  }
  
  //// plot?
  //if (oNB->getPlot_flg()) {
  //    plotwGNUplot(freq, write_2D_TLoss);
  //}
  //oNB->printParams();
  
  ierr = PetscFinalize();CHKERRQ(ierr);

  free_cmatrix(Rmx,    2*Nmin, 2*Nmin);
  free_cmatrix(S_curr, 2*Nmin, 2*Nmin);
  free_cmatrix(S_next, 2*Nmin, 2*Nmin);

  free_cmatrix(Rmx_ll,    2*Nmin, 2*Nmin);
  free_cmatrix(S_curr_ll, 2*Nmin, 2*Nmin);
  free_cmatrix(S_next_ll, 2*Nmin, 2*Nmin);
    
  free_cmatrix(S4plusS3D, Nmin, Nmin);
  free_cmatrix(S4plusS3D_ll, Nmin, Nmin);
  
  delete [] DD;
  delete [] ss;
  delete [] S3ss;
  delete [] bj;
  delete [] ab_curr;
  delete [] ab_next;
  
  delete [] DD_ll;
  delete [] ss_ll;  
  delete [] S3ss_ll;
  delete [] bj_ll;
  delete [] ab_curr_ll;
  delete [] ab_next_ll;  
  
  delete [] k_curr;
  delete [] rho_curr;
  free_dmatrix(v_curr, Nmin, Nmin);
   
  delete [] k_next;
  delete [] rho_next;
  free_dmatrix(v_next, Nmin, Nmin);
  delete opt;
  delete oNB;

  cout << "Intermediary files are saved in subdirectory: " << subdir << "." << endl;
  cout << "File tloss_rd2wcm_1d.lossless.nm written." << endl;
  cout << "File tloss_rd2wcm_1d.nm written." << endl;
  if (write_2D_TLoss) {
      cout << "File tloss_rd2wcm_2d.nm written." << endl;
  }
  cout << "\n ... main() is done." << endl;
  time_t tm2 = time(NULL);
  cout << "Run duration: " << difftime(tm2,tm1) << " seconds." << endl;	

  return (0);
} // -------------- end of main(); -----------------------------------------


//
// Function to parse the input options (both command lines and in the options file ModessRD.options)
//
AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();

  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "|                             NCPA Infrasound                              |" );  
  opt->addUsage( "|               Normal Modes for Range-Dependent Environments              |" );
  opt->addUsage( "|                      Two-Way Coupled Mode Solution                       |" );  
  opt->addUsage( "|           Single Frequency: Effective Sound Speed Approximation          |" );
  opt->addUsage( "----------------------------------------------------------------------------" );	
  opt->addUsage( "Usage: " );
  opt->addUsage( "By default the program computes the 1D transmission loss (TL)" );
  opt->addUsage( "at the ground or the specified receiver height and saves the data to 2 files:" );
  opt->addUsage( "   file tloss_rd2wcm_1d.nm - considering attenuation in the atmosphere" );
  opt->addUsage( "   file tloss_rd2wcm_1d.lossless.nm  - no attenuation" );
	opt->addUsage( "Additionally, if the flag --write_2D_TLoss is present on the command line the 2D TL is saved to file tloss_rd_2d.nm" );
  opt->addUsage( "The options below can be specified in a colon-separated file \"Modess.options\" or at the command line.  Command-line options override file options." );
  opt->addUsage( " --help -h                Print this message and exit" );
  opt->addUsage( "" );
  opt->addUsage(  " The atmosphere can be specified from one of 2 different sources:");
  opt->addUsage( "    1. An .env file containing the atmospheric specifications at certain ranges:" );
  opt->addUsage( "       use option --g2senvfile <filename>" );
  opt->addUsage( "    2. Several ASCII files stored in a given directory:" );
  opt->addUsage( "       use option --use_1D_profiles_from_dir <mydirname>" );
  //opt->addUsage( "The program requires an .env file containing the atmospheric specifications at certain ranges" );
  opt->addUsage( "The following options apply:" );	
  opt->addUsage( "" );	
  opt->addUsage( "REQUIRED (no default values):" );
  //opt->addUsage( " --atmosfile  <filename>  Uses an ASCII atmosphere file" );
  opt->addUsage( " --g2senvfile <filename>  Uses an .env binary file containing multiple 1D profiles" );
  opt->addUsage( " --atmosfileorder         The order of the (z,t,u,v,w,p,d) fields in the ASCII file (Ex: 'ztuvpd')" );
  opt->addUsage( " --skiplines              Lines at the beginning of the ASCII file to skip" );
  opt->addUsage( " --azimuth                Degrees in range [0,360), clockwise from North" );
  opt->addUsage( " --freq                   Frequency [Hz]" );	
  opt->addUsage( "" );	
  opt->addUsage( "OPTIONAL [defaults]:" ); 
  opt->addUsage( " --maxheight_km           Calculation grid height in km above MSL [150 km]" );
  opt->addUsage( " --zground_km             Height of the ground level above MSL [0 km]" );  
  opt->addUsage( " --Nz_grid                Number of points on the z-grid from ground to maxheight [20000]" );  
  opt->addUsage( " --sourceheight_km        Source height in km Above Ground Level (AGL) [0]" );
  opt->addUsage( " --receiverheight_km      Receiver height in km AGL [0]" );
  opt->addUsage( " --maxrange_km            Maximum horizontal distance from origin to propagate [1000 km]" );
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
  opt->addUsage( "                          e.g. --use_profile_ranges_km  0_20_50_80.5_300     " );   
  opt->addUsage( "                          The profiles at certain ranges specified by numbers" );
  opt->addUsage( "                          (in km) in a string such as 0_20_50_80.5_300 are");
  opt->addUsage( "                          requested in the propagation. Note that underscores" );
  opt->addUsage( "                          are necessary to separate the numbers." );
  opt->addUsage( "                          Note also that these are requested ranges;" );
  opt->addUsage( "                          however the left-closest profile available" );
  opt->addUsage( "                          in the .env file will actually be used; " );
  opt->addUsage( "                          for instance we request the profile at 300 km " );
  opt->addUsage( "                          but in the .env file the left-closest profile" );
  opt->addUsage( "                          may be available at 290 km and it is the one used." );
  opt->addUsage( "    Example: >>  ../bin/ModessRD2WCM --atmosfile g2sgcp2011012606L.jordan.env ");
  opt->addUsage( "                 --atmosfileorder zuvwtdp --azimuth 90 --freq 0.01 ");
  opt->addUsage( "                 --use_profiles_ranges_km 100_200_250 --maxrange_km 500 ");
  opt->addUsage( "" ); 
  opt->addUsage( " --use_profiles_at_steps_km" );
  opt->addUsage( "                          e.g. --use_profiles_at_steps_km 100" );
  opt->addUsage( "                          The profiles are requested at equidistant intervals " );
  opt->addUsage( "                          specified by this option [1000]" );
  opt->addUsage( "" );
  opt->addUsage( " --use_1D_profiles_from_dir" );
  opt->addUsage( "                          e.g. --use_1D_profiles_from_dir myprofiles" );
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
  opt->addUsage( "    Example: >> ../bin/ModessRD2WCM --atmosfileorder zuvwtdp --skiplines 1" );
  opt->addUsage( "                --azimuth 90 --freq 0.1 --use_1D_profiles_from_dir myprofiles" );
  opt->addUsage( "                --use_profile_ranges_km 0_100_300_500 " );                        
  opt->addUsage( "" );  
  opt->addUsage( "FLAGS (no value required):" );
  opt->addUsage( " --write_2D_TLoss         Outputs the 2D transmission loss to" );
  opt->addUsage( "                          default file: tloss_rd2wcm_2d.nm" );	
  opt->addUsage( "" );
  opt->addUsage( "" );
  opt->addUsage( " The format of the output files are as follows (column order):" );
  opt->addUsage( "  tloss_rd2wcm_1d.nm:           r, 4*PI*Re(P), 4*PI*Im(P), (incoherent TL)" );
  opt->addUsage( "  tloss_rd2wcm_1d.lossless.nm:" );
  opt->addUsage( "  tloss_rd2wcm_2d.nm:           r, z, 4*PI*Re(P), 4*PI*Im(P)" );   
  opt->addUsage( "" );
  opt->addUsage( "  Examples to run (in the 'samples' directory):" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/ModessRD2WCM --use_1D_profiles_from_dir profiles --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --use_profile_ranges_km 0_100_200_300 --maxrange_km 500" );
  opt->addUsage( "" );  
  opt->addUsage( "    ../bin/ModessRD2WCM --g2senvfile g2sgcp2011012606L.jordan.env --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --use_profile_ranges_km 0_100_200_250 --maxrange_km 500" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/ModessRD2WCM --use_1D_profiles_from_dir profiles --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --maxrange_km 500" );
  opt->addUsage( "" ); 
  opt->addUsage( "    Note: if options --use_profile_ranges_km/--use_profiles_at_steps_km are not used then we fall back on the range-independent case using the first available atm. profile." );
  opt->addUsage( "" );

  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  opt->setFlag( "write_2D_TLoss" );
  opt->setFlag( "plot" );

  opt->setOption( "atmosfile" );
  opt->setOption( "atmosfileorder" );
  opt->setOption( "g2senvfile" );
  opt->setOption( "wind_units" );  
  opt->setOption( "use_1D_profiles_from_dir" );
  opt->setOption( "slicefile" );  
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
  opt->processFile( "./ModessRD2WCM.options" );
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


 
