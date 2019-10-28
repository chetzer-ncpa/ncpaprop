#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <complex>
#include <sys/stat.h>
#include <ctime>
//#include <omp.h>

#include "Atmosphere.h"
#include "anyoption.h"
#include "SolveModBB.h"
#include "slepceps.h"
#include "slepcst.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;


// constructor 1
NCPA::SolveModBB::SolveModBB( \
          string filename, string atmosfile, string wind_units, NCPA::SampledProfile *atm_profile, \
          int Nfreq, double f_min, double f_step, double f_max, int Nz_grid, double azi, \
          double z_min, double maxheight, double sourceheight, double receiverheight, \
          string gnd_imp_model, int Lamb_wave_BC, \
          bool out_dispersion, bool out_disp_src2rcv, bool usemodess_flg, \
          bool turnoff_WKB)
{
  setParams(  filename, atmosfile, wind_units, atm_profile, Nfreq, f_min, f_step, f_max, \
              Nz_grid, azi, z_min, maxheight, sourceheight, receiverheight, \
              gnd_imp_model, Lamb_wave_BC, out_dispersion, out_disp_src2rcv, \
              usemodess_flg, turnoff_WKB);
}

//
// constructor 2
//
NCPA::SolveModBB::SolveModBB(ProcessOptionsBB *oBB, SampledProfile *atm_profile)
{
  setParams( oBB, atm_profile );                  
}


//destructor
NCPA::SolveModBB::~SolveModBB()
{
  delete [] Hgt;
  delete [] zw;
  delete [] mw;
  delete [] T;
  delete [] rho;
  delete [] Pr;
  //printf("SolveModNB destructor done.\n");
}

// setParams
void NCPA::SolveModBB::setParams( \
          string filename1, string atmosfile1, string wind_units1, \
          NCPA::SampledProfile *atm_profile1, \
          int Nfreq1, double f_min1, double f_step1, double f_max1, int Nz_grid1, double azi1, \
          double z_min1, double maxheight1, double sourceheight1, double receiverheight1, \
          string gnd_imp_model1, int Lamb_wave_BC1, \
          bool out_dispersion1, bool out_disp_src2rcv1, bool usemodess_flg1, \
          bool turnoff_WKB1)						
{			
  tol               = 1.0E-8; // tolerance for Slepc computations
  filename          = filename1;
  atm_profile       = atm_profile1;
  atmosfile         = atmosfile1;
  wind_units        = wind_units1;
  Nfreq             = Nfreq1;
  f_min             = f_min1;
  f_step            = f_step1;
  f_max             = f_max1;
  Nz_grid           = Nz_grid1;
  z_min             = z_min1;
  azi               = azi1;
  maxheight         = maxheight1;
  sourceheight      = sourceheight1;
  receiverheight    = receiverheight1;
  gnd_imp_model     = gnd_imp_model1;
  Lamb_wave_BC      = Lamb_wave_BC1;
  out_dispersion    = out_dispersion1;
  out_disp_src2rcv  = out_disp_src2rcv1;
  usemodess_flg     = usemodess_flg1;
  turnoff_WKB       = turnoff_WKB1;
  
  
  // get Hgt, zw, mw, T, rho, Pr in SI units; deleted in destructor
  Hgt = new double [Nz_grid];
  zw  = new double [Nz_grid];
  mw  = new double [Nz_grid];
  T   = new double [Nz_grid];
  rho = new double [Nz_grid];
  Pr  = new double [Nz_grid];
  
  //
  // !!! ensure maxheight is less than the max height covered by the 
  // provided atm profile may avoid some errors asociated with 
  // the code thinking that it goes above the max height when in fact 
  // the height may only differ from max height by a rounding error. 
  //
  if (maxheight/1000.0 >= atm_profile->z(atm_profile->nz()-1)) {
      maxheight = (atm_profile->z(atm_profile->nz()-1) - 1e-9)*1000.0; // slightly less
      //cout << "maxheight adjusted to: " << maxheight << " meters" << endl;
  }  
  
  // fill and convert to SI units
  double dz       = (maxheight - z_min)/Nz_grid;	// the z-grid spacing
  double z_min_km = z_min/1000.0;
  double dz_km    = dz/1000.0;
  double kmps2mps = 1.0;
  if (!wind_units.compare("kmpersec")) {
      kmps2mps = 1000.0;
  }

  // Note: the rho, Pr, T, zw, mw are computed wrt ground level i.e.
  // the first value is at the ground level e.g. rho[0] = rho(z_min)  
  for (int i=0; i<Nz_grid; i++) {
      Hgt[i] = (z_min_km + i*dz_km)*1000.0; // Hgt[0] = zground MSL
      rho[i] = atm_profile->rho(z_min_km + i*dz_km)*1000.0;
      Pr[i]  = atm_profile->p(z_min_km + i*dz_km)*100.0;
      T[i]   = atm_profile->t(z_min_km + i*dz_km);
      zw[i]  = atm_profile->u(z_min_km + i*dz_km)*kmps2mps;
      mw[i]  = atm_profile->v(z_min_km + i*dz_km)*kmps2mps;
  }

  std::time_t tm1 = std::time(NULL);
  if (0) {
      // data will be saved in a new directory using date_time format
      // YYYY_MM_DD_HH_MM_SS

      // form a string containing time in the desired format
      //std::time_t tm1 = std::time(NULL);
      char sdir[25];
      if (std::strftime(sdir, 25, "%F_%H_%M_%S", std::localtime(&tm1))) {
	        //cout << sdir << endl;
      }
      else {
          throw invalid_argument("time string YYYY_MM_DD_HH_MM_SS needed for dir name could not be formed.");
      }
      subdir  = string(sdir);             // subdirectory name
      disp_fn = subdir + "/" + filename;  // full name for the dispersion file
  }
  else { // no special directory for output
      subdir  = string("");   // blanc
      disp_fn = filename;     // full name for the dispersion file  
  }	
}


// setParams() prototype 2
void NCPA::SolveModBB::setParams(ProcessOptionsBB *oBB, SampledProfile *atm_profile1)						
{			
/*
  tol               = 1.0E-8; // tolerance for Slepc computations
  filename          = filename1;
  atm_profile       = atm_profile1;
  atmosfile         = atmosfile1;
  wind_units        = wind_units1;
  Nfreq             = Nfreq1;
  f_min             = f_min1;
  f_step            = f_step1;
  f_max             = f_max1;
  Nz_grid           = Nz_grid1;
  z_min             = z_min1;
  azi               = azi1;
  maxheight         = maxheight1;
  sourceheight      = sourceheight1;
  receiverheight    = receiverheight1;
  gnd_imp_model     = gnd_imp_model1;
  Lamb_wave_BC      = Lamb_wave_BC1;
  out_dispersion    = out_dispersion1;
  out_disp_src2rcv  = out_disp_src2rcv1;
  usemodess_flg     = usemodess_flg1;
  turnoff_WKB       = turnoff_WKB1;
  
  
            string filename1, string atmosfile1, string wind_units1, \
          NCPA::SampledProfile *atm_profile1, \
          int Nfreq1, double f_min1, double f_step1, double f_max1, int Nz_grid1, double azi1, \
          double z_min1, double maxheight1, double sourceheight1, double receiverheight1, \
          string gnd_imp_model1, int Lamb_wave_BC1, \
          bool out_dispersion1, bool out_disp_src2rcv1, bool usemodess_flg1, \
          bool turnoff_WKB1
          
          
  */
  
  
  
  
  
    // obtain the parameter values from the user's options
    
      tol            = 1.0E-8; // tolerance for Slepc computations 
     //tol				   			 = oNB->getSlepcTolerance(); 
      atmosfile      = oBB->getAtmosfile();
      atm_profile    = atm_profile1;
      //atmosfileorder = oBB->getAtmosfileorder();
      wind_units     = oBB->getWindUnits();
      gnd_imp_model  = oBB->getGnd_imp_model();      
      filename       = oBB->getOut_disp_file();
      usrattfile     = oBB->getUsrAttFile();
      
      usemodess_flg  = oBB->getUsemodess_flg();
          
      //skiplines      = oBB->getSkiplines();
      z_min          = oBB->getZ_min();
      azi            = oBB->getAzimuth();
      //maxrange       = oBB->getMaxrange();
      maxheight      = oBB->getMaxheight();
      sourceheight   = oBB->getSourceheight();
      receiverheight = oBB->getReceiverheight();   
      Nz_grid        = oBB->getNz_grid();
      Nrng_steps     = oBB->getNrng_steps();
      Nfreq          = oBB->getNfreq();
      f_min          = oBB->getF_min();
      f_step         = oBB->getF_step();			
      f_max          = oBB->getF_max();
      Lamb_wave_BC   = oBB->getLamb_wave_BC();
      out_dispersion = oBB->getW_disp_flg();
      out_disp_src2rcv = oBB->getW_disp_src2rcv_flg();
  
      turnoff_WKB    = oBB->getTurnoff_WKB();
      //cout << "turnoff_WKB = " << turnoff_WKB << endl;

      // default values for c_min, c_max and wvnum_filter_flg
      c_min = 0.0;
      c_max = 0.0;
      wvnum_filter_flg = 0;

      // set c_min, c_max if wavenumber filtering is on
      wvnum_filter_flg = oBB->getWvnum_filter_flg();     
      if (wvnum_filter_flg==1) {
          c_min = oBB->getC_min();
          c_max = oBB->getC_max();
      };


  // get Hgt, zw, mw, T, rho, Pr in SI units; deleted in destructor
  Hgt = new double [Nz_grid];
  zw  = new double [Nz_grid];
  mw  = new double [Nz_grid];
  T   = new double [Nz_grid];
  rho = new double [Nz_grid];
  Pr  = new double [Nz_grid];
  
  //
  // !!! ensure maxheight is less than the max height covered by the 
  // provided atm profile may avoid some errors asociated with 
  // the code thinking that it goes above the max height when in fact 
  // the height may only differ from max height by a rounding error. 
  //
  if (maxheight/1000.0 >= atm_profile->z(atm_profile->nz()-1)) {
      maxheight = (atm_profile->z(atm_profile->nz()-1) - 1e-9)*1000.0; // slightly less
      //cout << "maxheight adjusted to: " << maxheight << " meters" << endl;
  }  
  
  // fill and convert to SI units
  double dz       = (maxheight - z_min)/Nz_grid;	// the z-grid spacing
  double z_min_km = z_min/1000.0;
  double dz_km    = dz/1000.0;
  double kmps2mps = 1.0;
  if (!wind_units.compare("kmpersec")) {
      kmps2mps = 1000.0;
  }

  // Note: the rho, Pr, T, zw, mw are computed wrt ground level i.e.
  // the first value is at the ground level e.g. rho[0] = rho(z_min)  
  for (int i=0; i<Nz_grid; i++) {
      Hgt[i] = (z_min_km + i*dz_km)*1000.0; // Hgt[0] = zground MSL
      rho[i] = atm_profile->rho(z_min_km + i*dz_km)*1000.0;
      Pr[i]  = atm_profile->p(z_min_km + i*dz_km)*100.0;
      T[i]   = atm_profile->t(z_min_km + i*dz_km);
      zw[i]  = atm_profile->u(z_min_km + i*dz_km)*kmps2mps;
      mw[i]  = atm_profile->v(z_min_km + i*dz_km)*kmps2mps;
  }

  std::time_t tm1 = std::time(NULL);
  if (0) {
      // data will be saved in a new directory using date_time format
      // YYYY_MM_DD_HH_MM_SS

      // form a string containing time in the desired format
      //std::time_t tm1 = std::time(NULL);
      char sdir[25];
      if (std::strftime(sdir, 25, "%F_%H_%M_%S", std::localtime(&tm1))) {
	        //cout << sdir << endl;
      }
      else {
          throw invalid_argument("time string YYYY_MM_DD_HH_MM_SS needed for dir name could not be formed.");
      }
      subdir  = string(sdir);             // subdirectory name
      disp_fn = subdir + "/" + filename;  // full name for the dispersion file
  }
  else { // no special directory for output
      subdir  = string("");   // blanc
      disp_fn = filename;     // full name for the dispersion file  
  }
}


// utility to print the parameters to the screen
void NCPA::SolveModBB::printParams() {
  printf(" Normal Modes - broadband - run info:\n");
  printf("                  Nfreq : %d\n", Nfreq);  
  printf("               freq_min : %g\n", f_min);
  printf("              freq_step : %g\n", f_step);  
  printf("               freq_max : %g\n", f_max);  
  printf("                azimuth : %g\n", azi);
  printf("                Nz_grid : %d\n", Nz_grid); 
  printf("      z_min (meters MSL): %g\n", z_min);
  printf("      maxheight_km (MSL): %g\n", maxheight/1000.0);
  printf("   sourceheight_km (AGL): %g\n", sourceheight/1000.0);
  printf(" receiverheight_km (AGL): %g\n", receiverheight/1000.0);  
  printf("          gnd_imp_model : %s\n", gnd_imp_model.c_str());
  printf("Lamb wave boundary cond : %d\n", Lamb_wave_BC);
  printf("  SLEPc tolerance param : %g\n", tol);
  printf("    atmospheric profile : %s\n", atmosfile.c_str());
  printf("       turnoff_WKB flag : %d\n", turnoff_WKB);
  if (!usrattfile.empty()) {
  printf("  User attenuation file : %s\n", usrattfile.c_str());
  }
  if (out_disp_src2rcv) {
  printf("          data saved in : %s\n", disp_fn.c_str());
  } else {
  printf("   data saved in subdir : %s\n", subdir.c_str());
  }
}


int NCPA::SolveModBB::computeModes(bool usemodess_flg) {
    if (usemodess_flg) {
    computeModESS();
    }
    else {
    computeWmodes();
    }
    return 0;
}


int NCPA::SolveModBB::computeModESS() {	
  //
  // Declarations related to Slepc computations
  //
  Mat             A;           /* problem matrix */
  EPS             eps;         /* eigenproblem solver context */
  ST              stx;
  KSP             kspx;
  PC              pcx;
  EPSType         type;        // CHH 191029: Removed const 
  PetscReal       re, im;
  PetscScalar     kr, ki, *xr_;
  Vec             xr, xi;
  PetscInt        Istart, Iend, col[3], its, maxit, nconv;
  PetscBool       FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
  PetscScalar     value[3];	
  PetscErrorCode  ierr;
  PetscMPIInt     rank, size;

  int    i, j, select_modes, nev;
  int    NN, Nz_subgrid;
  double freq, admittance, z_min_km, delZ, h2;
  //double dz_km;
  double k_min, k_max;			
  double *kreal, *kim;
  double *alpha, *diag, *k2, *k_s, **v, **v_s;	
  complex<double> *k_pert;

  alpha  = new double [Nz_grid];
  diag   = new double [Nz_grid];
  k2     = new double [MAX_MODES];
  k_s    = new double [MAX_MODES];
  kreal  = new double [MAX_MODES];
  kim    = new double [MAX_MODES];
  k_pert = new complex<double> [MAX_MODES];
  v      = dmatrix(Nz_grid,MAX_MODES);
  v_s    = dmatrix(Nz_grid,MAX_MODES); 	

  k_min = 0; 
  k_max = 0;
  nev   = 0;

  Nfreq = (int) round((f_max-f_min)/f_step) + 1;
  cout << "Nfreq = " << Nfreq << endl;	

  dz = (maxheight - z_min)/Nz_grid;	// the z-grid spacing
  h2 = dz*dz;
  //dz_km    = dz/1000.0;
  z_min_km = z_min/1000.0;
  
  //// get the density array
  //for (i=0; i<Nz_grid; i++) {
  //    rho[i] = atm_profile->rho(z_min_km + i*dz_km)*1000.0;
  //}

  // will save data on a z-subgrid with spacing delZ = round(c/(4*f_max)/dz)*dz;
  // this z-subgrid is a subset of the altitudes at which the modes are computed
  // the z-subgrid has to adequately sample the pulse spatially 
  NN         = (int) round(340.0/(4*f_max)/dz); 
  delZ       = (double) NN*dz;		// delZ = NN*dz
  Nz_subgrid = (int) floor(Nz_grid*dz/delZ);    // total number of points on the z-subgrid 

  if (0) { // if making a special dir for output is requested
  // make a dir
  // if the wrong directory permissions are set for this dir look up
  // man 2 umask and man 2 stat
  if (mkdir(subdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO)==-1) {
      std::ostringstream es;
      es << "attempt to make dir << " << subdir << " failed.";
      throw invalid_argument(es.str());
  }
  }

  // Initialize Slepc
  SlepcInitialize(PETSC_NULL,PETSC_NULL,(char*)0,PETSC_NULL);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
  
  FILE *fp;
  if (out_disp_src2rcv) {
    // open dispersion file for writing
    fp = fopen(disp_fn.c_str(),"w");
  }

  //
  // big loop over frequencies; no parallelization here yet
  //
  for (int ii = 0; ii<Nfreq; ii++) {

      freq = ii*f_step + f_min;
      cout << "Now processing frequency = " << freq << " Hz" << endl;

      ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
      ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Nz_grid,Nz_grid);CHKERRQ(ierr);
      ierr = MatSetFromOptions(A);CHKERRQ(ierr);
      
      // the following Preallocation call is needed in PETSc version 3.3
      ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL);CHKERRQ(ierr);

      // compute absorption
      //getAbsorption(Nz_grid, dz, atm_profile, freq, alpha);
      getAbsorption(Nz_grid, dz, atm_profile, freq, usrattfile, alpha);

      //
      // ground impedance model (to be implemented)
      //
      // at the ground the BC is: Phi' = (a - 1/2*dln(rho)/dz)*Phi
      // for a rigid ground a=0; and the BC is the Lamb Wave BC:
      // admittance = -1/2*dln(rho)/dz
      //  
      if ((gnd_imp_model.compare("rigid")==0) && Lamb_wave_BC) {
          admittance = -atm_profile->drhodz(z_min/1000.0)/1000.0/atm_profile->rho(z_min_km)/2.0; // SI units
      }
      else if (gnd_imp_model.compare("rigid")==0) {
          admittance = 0.0; // no Lamb_wave_BC
      }
      else {
          std::ostringstream es;
          es << "This ground impedance model is not implemented yet: " \
             << gnd_imp_model;
          throw invalid_argument(es.str()); 
      }	   

      // Get the main diagonal
      cout << "2 atm_profile->getPropagationAzimuth() = " << atm_profile->getPropagationAzimuth() << endl;
      i = getModalTraceModESS(Nz_grid, z_min, sourceheight, receiverheight, dz, \
                              atm_profile, admittance, freq, azi, diag, &k_min, &k_max, turnoff_WKB);
  
      // if wavenumber filtering is on, redefine k_min, k_max
      if (wvnum_filter_flg) {
          k_min = 2*Pi*freq/c_max;
          k_max = 2*Pi*freq/c_min;
      }

      i = getNumberOfModes(Nz_grid,dz,diag,k_min,k_max,&nev);

      printf("Number of modes found: nev = %d\n", nev);

      // skip this iteration if no modes are found
      if ( nev==0 ) {
          printf (" -> No modes found for frequency %g Hz and sound speed range from %6.2f m/s to %6.2f m/s\n", \
                  freq, 2*Pi*freq/k_max, 2*Pi*freq/k_min);
          printf(" -> Check your input values.\n");
          if (wvnum_filter_flg) {
              printf(" -> Check the frequency and/or perhaps choose different c_min and c_max values.\n");
          }
          continue;
          //printf("Aborted\n");
          //exit(1);
      }  

      printf ("______________________________________________________________________\n\n");
      printf (" -> Normal mode solution at %6.3f Hz and %6.2f deg (%d modes)...\n", freq, azi, nev);
      printf (" -> Discrete spectrum: %6.2f m/s to %6.2f m/s\n", 2*Pi*freq/k_max, 2*Pi*freq/k_min);

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         Compute the operator matrix that defines the eigensystem, Ax=kx
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

      // Make the A matrix 
      ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
      if (Istart==0) FirstBlock=PETSC_TRUE;
      if (Iend==Nz_grid) LastBlock=PETSC_TRUE;
      value[0]=1.0/h2; value[2]=1.0/h2;
      for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {
          value[1] = -2.0/h2 + diag[i];
          col[0]=i-1; col[1]=i; col[2]=i+1;
          ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
      }
      if (LastBlock) {
          i=Nz_grid-1; col[0]=Nz_grid-2; col[1]=Nz_grid-1;
          ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
      }
      if (FirstBlock) {
          i=0; col[0]=0; col[1]=1; value[0]=-2.0/h2 + diag[0]; value[1]=1.0/h2;
          ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
      }

      ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

      // CHH 191028: MatGetVecs deprecated, changed to MatCreateVecs
      ierr = MatCreateVecs(A,PETSC_NULL,&xr);CHKERRQ(ierr);
      ierr = MatCreateVecs(A,PETSC_NULL,&xi);CHKERRQ(ierr);

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    Create the eigensolver and set various options
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      /* 
         Create eigensolver context
      */
      ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

      /* 
         Set operators. In this case, it is a standard eigenvalue problem
      */
      ierr = EPSSetOperators(eps,A,PETSC_NULL);CHKERRQ(ierr);
      ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

      /*
         Set solver parameters at runtime
      */
      ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
      ierr = EPSSetType(eps,"krylovschur"); CHKERRQ(ierr);
      ierr = EPSSetDimensions(eps,10,PETSC_DECIDE,PETSC_DECIDE); CHKERRQ(ierr);
      //ierr = EPSSetTarget(eps,sigma); CHKERRQ(ierr);
      ierr = EPSSetTolerances(eps,tol,PETSC_DECIDE); CHKERRQ(ierr);

      ierr = EPSGetST(eps,&stx); CHKERRQ(ierr);
      ierr = STGetKSP(stx,&kspx); CHKERRQ(ierr);
      ierr = KSPGetPC(kspx,&pcx); CHKERRQ(ierr);
      ierr = STSetType(stx,"sinvert"); CHKERRQ(ierr);
      ierr = KSPSetType(kspx,"preonly");
      ierr = PCSetType(pcx,"cholesky");
      //ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);
      ierr = EPSSetWhichEigenpairs(eps,EPS_ALL); CHKERRQ(ierr);
      ierr = EPSSetInterval(eps,pow(k_min,2),pow(k_max,2)); CHKERRQ(ierr);

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                          Solve the eigensystem
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

      ierr = EPSSolve(eps);CHKERRQ(ierr);
      /*
         Optional: Get some information from the solver and display it
      */
      ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
      ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
      ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev);CHKERRQ(ierr);
      ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);

      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        Display solution and clean up
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      /* 
         Get number of converged approximate eigenpairs
      */
      ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

      if (nconv>0) {
          for (i=0;i<nconv;i++) {
              ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
              //ierr = EPSComputeRelativeError(eps,i,&error);CHKERRQ(ierr);

        #if defined(PETSC_USE_COMPLEX)
              re = PetscRealPart(kr);
              im = PetscImaginaryPart(kr);
        #else
              re = kr;
              im = ki;
        #endif 
              k2[i] = re;
              ierr = VecGetArray(xr,&xr_);CHKERRQ(ierr);
              for (j = 0; j < Nz_grid; j++) {
               		v[j][i] = xr_[j]/sqrt(dz);
              }
              ierr = VecRestoreArray(xr,&xr_);CHKERRQ(ierr);
          }
      }		

      // select and perturb the modes
      doSelectModESS(Nz_grid,nconv,k_min,k_max,k2,v,k_s,v_s,&select_modes); 
      doPerturb2(Nz_grid, z_min, dz, select_modes, freq, atm_profile, k_s, v_s, alpha, k_pert, kreal, kim);
        
      //
      // Output data: saving to file(s)
      //
      if (out_dispersion) {
          // printf ("Writing to file: %s at freq = %8.3f Hz...\n", filen.c_str(), freq);
          writeDispersion_bb_bin3( disp_fn, freq, Nfreq, f_step, Nz_grid, z_min, \
                                   Nz_subgrid, delZ, NN, select_modes, dz, \
                                   sourceheight, rho, kreal, kim, v_s);															
      }
      else if (out_disp_src2rcv) {			
          // source-to-receiver will be written to one ascii file
          //writeDispersion_bb_ascii( disp_fn, select_modes, dz, sourceheight, receiverheight, freq, rho, k_pert, v_s);
          writeDispersion_bb_ascii( fp, select_modes, dz, \
                                    sourceheight, receiverheight, freq, rho, k_pert, v_s);
      }			
      cout << "frequency = " << freq << " Hz processed" << endl;			

      // Free work space
      ierr = EPSDestroy(&eps);CHKERRQ(ierr);
      ierr = MatDestroy(&A);  CHKERRQ(ierr);
      ierr = VecDestroy(&xr); CHKERRQ(ierr);
      ierr = VecDestroy(&xi); CHKERRQ(ierr); 	
  }	// end BIG LOOP over frequencies
  
  if (out_disp_src2rcv) {
    fclose(fp); // close the dispersion file
  }

  // Clean up
  delete[] alpha;
  delete[] diag;
  delete[] k2;
  delete[] k_s;  
  delete[] k_pert;
  delete[] kreal;
  delete[] kim;			
  free_dmatrix(v,Nz_grid,MAX_MODES);
  free_dmatrix(v_s,Nz_grid,MAX_MODES);

  // SlepcFinalize
  ierr = SlepcFinalize();CHKERRQ(ierr);	
	
  return 0;		
}  // end of SolveModBB::computeModess




int NCPA::SolveModBB::computeWmodes() {	
  //
  // Declarations related to Slepc computations
  //  
  Mat            A, B;       
  EPS            eps;  		// eigenproblem solver context      
  ST             stx;
  EPSType  type;                
  PetscReal      re, im;
  PetscScalar    kr, ki, *xr_;
  Vec            xr, xi;
  PetscInt       Istart, Iend, col[3], its, maxit, nconv;
  PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
  PetscScalar    value[3];
  PetscErrorCode ierr;
  PetscMPIInt    rank, size;
  Vec            V_SEQ;
  VecScatter     ctx;

  int    i, j, select_modes, nev;
  int    NN, Nz_subgrid;
  double freq, dz, z_min_km, admittance, h2, delZ;
  //double dz_km;
  double k_min, k_max, sigma;			
  double *alpha, *diag, *kd, *md, *cd, *kH, *k_s, *kreal, *kim, **v, **v_s;	
  complex<double> *k_pert;

  alpha  = new double [Nz_grid];
  diag   = new double [Nz_grid];
  kd     = new double [Nz_grid];
  md     = new double [Nz_grid];
  cd     = new double [Nz_grid];
  kH     = new double [MAX_MODES];
  k_s    = new double [MAX_MODES];
  kreal  = new double [MAX_MODES];
  kim    = new double [MAX_MODES];
  k_pert = new complex<double> [MAX_MODES];
  v      = dmatrix(Nz_grid,MAX_MODES);
  v_s    = dmatrix(Nz_grid,MAX_MODES); 	

  nev   = 0;
  k_min = 0; 
  k_max = 0;
  
  Nfreq = (int) round((f_max-f_min)/f_step) + 1;
  cout << "Nfreq = " << Nfreq << endl;	

  dz = (maxheight - z_min)/Nz_grid;	// the z-grid spacing
  h2 = dz*dz;
  //dz_km    = dz/1000.0;
  z_min_km = z_min/1000.0;
  
  //// get the density array
  //for (i=0; i<Nz_grid; i++) {
  //    rho[i] = atm_profile->rho(z_min_km + i*dz_km)*1000.0;
  //}

  // will save data on a z-subgrid with spacing delZ = round(c/(4*f_max)/dz)*dz;
  // this z-subgrid is a subset of the altitudes at which the modes are computed
  // the z-subgrid has to adequately sample the pulse spatially 
  NN         = (int) round(340.0/(4*f_max)/dz); 
  delZ       = (double) NN*dz;		// delZ = NN*dz
  Nz_subgrid = (int) floor(Nz_grid*dz/delZ);    // total number of points on the z-subgrid 

  if (0) { // if making a special dir for output is requested
  // make a dir
  // if the wrong directory permissions are set for this dir look up
  // man 2 umask and man 2 stat
  if (mkdir(subdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO)==-1) {
      std::ostringstream es;
      es << "attempt to make dir << " << subdir << " failed.";
      throw invalid_argument(es.str());
  }
  }

  // Initialize Slepc
  SlepcInitialize(PETSC_NULL,PETSC_NULL,(char*)0,PETSC_NULL);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  FILE *fp;
  if (out_disp_src2rcv) {
    // open dispersion file for writing
    fp = fopen(disp_fn.c_str(),"w");
  }
  
  //
  // big loop over frequencies; no parallelization here yet
  //
  for (int ii = 0; ii<Nfreq; ii++) {

      freq = ii*f_step + f_min;
      cout << "Now processing frequency = " << freq << " Hz" << endl;

      ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
      ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Nz_grid,Nz_grid);CHKERRQ(ierr);
      ierr = MatSetFromOptions(A);CHKERRQ(ierr);
      
      // the following Preallocation call is needed in PETSc version 3.3
      ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL);CHKERRQ(ierr);

      // compute absorption
      //getAbsorption(Nz_grid, dz, atm_profile, freq, alpha);
      getAbsorption(Nz_grid, dz, atm_profile, freq, usrattfile, alpha);
     
      
      //
      // ground impedance model
      //
      // at the ground the BC is: Phi' = (a - 1/2*dln(rho)/dz)*Phi
      // for a rigid ground a=0; and the BC is the Lamb Wave BC:
      // admittance = -1/2*dln(rho)/dz
      //  
      admittance = 0.0; // default
      if ((gnd_imp_model.compare("rigid")==0) && Lamb_wave_BC) { //Lamb_wave_BC
          admittance = -atm_profile->drhodz(z_min/1000.0)/1000.0/atm_profile->rho(z_min_km)/2.0; // SI units
      }
      else if (gnd_imp_model.compare("rigid")==0) {
          admittance = 0.0; // no Lamb_wave_BC
      }
      else {
          std::ostringstream es;
          es << "This ground impedance model is not implemented yet: " << gnd_imp_model;
          throw invalid_argument(es.str());
      }      

      
    //
    // Get the main diagonal and the number of modes (WMod)
    //		
    i = getModalTraceWMod(Nz_grid, z_min, sourceheight, receiverheight, dz, atm_profile, admittance, \
                          freq, diag, kd, md, cd, &k_min, &k_max, turnoff_WKB);

    // if wavenumber filtering is on, redefine k_min, k_max
    if (wvnum_filter_flg) {
        k_min = 2*Pi*freq/c_max;
        k_max = 2*Pi*freq/c_min;
    }

    i = getNumberOfModes(Nz_grid,dz,diag,k_min,k_max,&nev);

    // skip this iteration if no modes are found
    if ( nev==0 ) {
        printf (" -> No modes found for frequency %g Hz and sound speed range from %6.2f m/s to %6.2f m/s\n", \
                  freq, 2*Pi*freq/k_max, 2*Pi*freq/k_min);
        printf(" -> Check your input values.\n");
        if (wvnum_filter_flg) {
            printf(" -> Check and perhaps choose different c_min and c_max values.\n");
        }
        continue;
        //printf("Aborted\n");
        //exit(1);
    }      
              
    //
    // set up parameters for eigenvalue estimation; double dimension of problem for linearization
    //

    //sigma = pow(0.5*(k_min+k_max),2);
    //sigma = (k_min*k_min + k_max*k_max)/2; //dv: !!! is this sigma better?  
    sigma     = 0.5*(k_min+k_max);
    
    int nev_2 = nev*2;
    int n_2   = Nz_grid*2;	

    if (rank == 0) {
        printf ("______________________________________________________________________\n\n");
        printf (" -> Solving wide-angle problem at %6.3f Hz and %6.2f deg (%d modes)...\n", freq, azi, nev_2);
        printf (" -> Discrete spectrum: %6.2f m/s to %6.2f m/s\n", 2*Pi*freq/k_max, 2*Pi*freq/k_min);
        printf (" -> Quadratic eigenvalue problem  - double dimensionality.\n");
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    //   Compute the operator matrices that define the eigensystem, A*x = k.B*x
    //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n_2,n_2);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
    ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n_2,n_2);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);CHKERRQ(ierr);
    
    // the following Preallocation calls are needed in PETSc version 3.3
    ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(B, 2, PETSC_NULL);CHKERRQ(ierr);

    /*
    We solve the quadratic eigenvalue problem (Mk^2 + Ck +D)v = 0 where k are the 
    eigenvalues and v are the eigenvectors and M , C, A are N by N matrices.
    We linearize this by denoting u = kv and arrive at the following generalized
    eigenvalue problem:

     / -D  0 \  (v )      / C  M \  (v ) 
    |         | (  ) = k |        | (  )
     \ 0   M /  (kv)      \ M  0 /  (kv)
     
     ie. of the form
     A*x = k.B*x
     so now we double the dimensions of the matrices involved: 2N by 2N.
     
     Matrix D is tridiagonal and has the form
     Main diagonal     : -2/h^2 + omega^2/c(1:N)^2 + F(1:N)
     Upper diagonal    :  1/h^2
     Lower diagonal    :  1/h^2
     Boundary condition:  A(1,1) =  (1/(1+h*beta) - 2)/h^2 + omega^2/c(1)^2 + F(1)
     
     where 
     F = 1/2*rho_0"/rho_0 - 3/4*(rho_0')^2/rho_0^2 where rho_0 is the ambient
     stratified air density; the prime and double prime means first and second 
     derivative with respect to z.
     beta  = alpha - 1/2*rho_0'/rho_0
     alpha is given in
     Psi' = alpha*Psi |at z=0 (i.e. at the ground). Psi is the normal mode.
     If the ground is rigid then alpha is zero. 
     
     Matrix M is diagonal:
     Main diagonal  : u0(1:N)^2/c(1:N)^2 - 1
      u0 = scalar product of the wind velocity and the horizontal wavenumber k_H
      u0 = v0.k_H
      
     Matrix C is diagonal:
     Main diagonal: 2*omega*u0(1:N)/c(1:N)^2 
     
    */

    // Assemble the A matrix (2N by 2N)
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
    if (Istart==0) FirstBlock=PETSC_TRUE;
    if (Iend==n_2) LastBlock =PETSC_TRUE;

    // matrix -D is placed in the first N by N block
    // kd[i]=(omega/c_T)^2
    for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? (Iend/2)-1: Iend/2); i++ ) {
        value[0]=-1.0/h2; value[1] = 2.0/h2 - kd[i]; value[2]=-1.0/h2;
        col[0]=i-1; col[1]=i; col[2]=i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (LastBlock) {
        i=(n_2/2)-1; col[0]=(n_2/2)-2; col[1]=(n_2/2)-1; value[0]=-1.0/h2; value[1]=2.0/h2 - kd[(n_2/2)-1];
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    // boundary condition
    if (FirstBlock) {
        i=0; col[0]=0; col[1]=1; value[0]=2.0/h2 - kd[0]; value[1]=-1.0/h2;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    // Insert matrix M into the lower N by N block of A
    // md is u0^2/c^2-1
    for ( i=(n_2/2); i<n_2; i++ ) {
        ierr = MatSetValue(A,i,i,md[i-(n_2/2)],INSERT_VALUES); CHKERRQ(ierr); 
    }

    // Assemble the B matrix
    for ( i=0; i<(n_2/2); i++ ) {
        col[0]=i; col[1]=i+(n_2/2); value[0]=cd[i]; value[1]=md[i];
        ierr = MatSetValues(B,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    for ( i=(n_2/2); i<n_2; i++ ) {
        ierr = MatSetValue(B,i,i-(n_2/2),md[i-(n_2/2)],INSERT_VALUES); CHKERRQ(ierr); 
    }

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    // CHH 191028: MatGetVecs deprecated, changed to MatCreateVecs
    ierr = MatCreateVecs(A,PETSC_NULL,&xr);CHKERRQ(ierr);
    ierr = MatCreateVecs(A,PETSC_NULL,&xi);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  Create the eigensolver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Create eigensolver context
    */
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

    /* 
       Set operators. In this case, it is a quadratic eigenvalue problem
    */
    ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);

    /*
       Set solver parameters at runtime
    */
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
    ierr = EPSSetType(eps,"krylovschur"); CHKERRQ(ierr);
    ierr = EPSSetDimensions(eps,nev_2,PETSC_DECIDE,PETSC_DECIDE); CHKERRQ(ierr);
    ierr = EPSSetTarget(eps,sigma); CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps,tol,PETSC_DECIDE); CHKERRQ(ierr);

    ierr = EPSGetST(eps,&stx); CHKERRQ(ierr);
    ierr = STSetType(stx,"sinvert"); CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        Solve the eigensystem
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = EPSSolve(eps);CHKERRQ(ierr);
    /*
       Optional: Get some information from the solver and display it
    */
    ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
    ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps,&nev_2,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev_2);CHKERRQ(ierr);
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Get number of converged approximate eigenpairs
    */
    ierr = EPSGetConverged(eps,&nconv); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %d\n\n",nconv);CHKERRQ(ierr);

    if (nconv>0) {
        for( i=0; i<nconv; i++ ) {
            ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
            #ifdef PETSC_USE_COMPLEX
                re = PetscRealPart(kr);
                im = PetscImaginaryPart(kr);
            #else
                re = kr;
                im = ki;
            #endif 
            kH[i] = re;
            ierr = VecScatterCreateToAll(xr,&ctx,&V_SEQ);CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx,xr,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx,xr,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
            if (rank == 0) {
                ierr = VecGetArray(V_SEQ,&xr_);CHKERRQ(ierr);
                for (j = 0; j < Nz_grid; j++) {
                    v[j][i] = xr_[j]/sqrt(dz);
                }
                ierr = VecRestoreArray(V_SEQ,&xr_);CHKERRQ(ierr);
            }
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            ierr = VecDestroy(&V_SEQ);CHKERRQ(ierr);
        }
    }

      // select and perturb wavenumbers to add attenuation as imaginary part
      doSelectWMod(Nz_grid, nconv, k_min, k_max, kH, v, k_s, v_s, &select_modes); 
      doPerturb2(Nz_grid, z_min, dz, select_modes, freq, atm_profile, k_s, v_s, alpha, k_pert, kreal, kim);
      
        
      //
      // Output data: saving to file(s)
      //
      if (out_dispersion) {
          // printf ("Writing to file: %s at freq = %8.3f Hz...\n", filen.c_str(), freq);
          writeDispersion_bb_bin3( disp_fn, freq, Nfreq, f_step, Nz_grid, z_min, \
                                  Nz_subgrid, delZ, NN, select_modes, dz, \
                                  sourceheight, rho, kreal, kim, v_s);															
      }
      else if (out_disp_src2rcv) {			
          // source-to-receiver will be written to one ascii file
          //writeDispersion_bb_ascii( disp_fn, select_modes, dz, sourceheight, receiverheight, freq, rho, k_pert, v_s);
          writeDispersion_bb_ascii( fp, select_modes, dz, \
                                    sourceheight, receiverheight, freq, rho, k_pert, v_s);
      }			
      cout << "frequency = " << freq << " Hz processed" << endl;			

      // Free work space
      ierr = EPSDestroy(&eps);CHKERRQ(ierr);
      ierr = MatDestroy(&A);  CHKERRQ(ierr);
      ierr = MatDestroy(&B);  CHKERRQ(ierr);
      ierr = VecDestroy(&xr); CHKERRQ(ierr);
      ierr = VecDestroy(&xi); CHKERRQ(ierr);       
      
      	
  }	// end BIG LOOP over frequencies
  
  if (out_disp_src2rcv) {
    fclose(fp); // close the dispersion file
  }
  // free the rest of locally dynamically allocated space
  delete[] alpha;
  delete[] diag;
  delete[] kd;
  delete[] md;
  delete[] cd;
  delete[] kH;
  delete[] k_s;
  delete[] k_pert;
  delete[] kreal;
  delete[] kim;	  
  free_dmatrix(v, Nz_grid, MAX_MODES);
  free_dmatrix(v_s, Nz_grid, MAX_MODES);  

  // SlepcFinalize
  ierr = SlepcFinalize();CHKERRQ(ierr);	
	
  return 0;		
}  // end of SolveModesBB::computeWmodes



/*
// updated getAbsorption function: bug fixed by Joel and Jelle - Jun 2012
int NCPA::SolveModBB::getAbsorption(int n, double dz, SampledProfile *p, double freq, double *alpha)
{
  // Expressions based on Bass and Sutherland, JASA 2004
  // Computes alpha(freq) for given G2S output
  // Subroutine can easily be modified to include dispersion effects
  // In that case, make alpha a complex array and
  // use the real part as absorption and the imaginary part for dispersion
  int    m, ii;
  double T_o, P_o, S, z;
  double X[7], X_ON, Z_rot[2], Z_rot_;
  double sigma, nn, chi, cchi, mu, nu, mu_o;
  double beta_0, beta_1, beta_2, alpha_1, alpha_2;
  double a_cl, a_rot, a_diff, a_vib;
  double T_z, P_z, c_snd_z, gamma;
  double A1, A2, B, C, D, E, F, G, H, I, J, K, L, ZZ, hu;
  double f_vib[4], a_vib_c[4], Cp_R[4], Cv_R[4], theta[4], C_R, A_max, Tr;

  // Atmospheric composition constants
  mu_o  = 18.192E-6;             // Reference viscosity [kg/(m*s)]
  T_o   = T[0];         // Reference temperature [K]
  P_o   = Pr[0];   // Reference pressure [Pa]
  S     = 117;   	   			       // Sutherland constant [K]       

  Cv_R[0] = 5/2;                                    // Heat capacity|volume (O2)
  Cv_R[1] = 5/2;                                    // Heat capacity|volume (N2)
  Cv_R[2] = 3;                                      // Heat capacity|volume (CO2)
  Cv_R[3] = 3;                                      // Heat capacity|volume (O3)
  Cp_R[0] = 7/2;                                    // Heat capacity|pressure (O2)
  Cp_R[1] = 7/2;                                    // Heat capacity|pressure (N2)
  Cp_R[2] = 4;                                      // Heat capacity|pressure (CO2)
  Cp_R[3] = 4;                                      // Heat capacity|pressure (O3)
  theta[0]= 2239.1;                                 // Charact. temperature (O2)
  theta[1]= 3352;                                   // Charact. temperature (N2)
  theta[2]= 915;                                    // Charact. temperature (CO2)
  theta[3]= 1037;                                   // Charact. temperature (O3)

	//gamma   = 1.371 + 2.46E-04*T_z - 6.436E-07*pow(T_z,2) + 5.2E-10*pow(T_z,3) - 1.796E-13*pow(T_z,4) + 2.182E-17*pow(T_z,5);
	gamma   = 1.4;
			 
  for (ii=0; ii<n; ii++) {
			z       = ii*dz/1000.0;		// km	
			T_z     = T[ii];	        // K
			P_z     = Pr[ii];	        // Pa;
			c_snd_z = sqrt(gamma*P_z/rho[ii]);  // in m/s					 
			mu      = mu_o*sqrt(T_z/T_o)*((1+S/T_o)/(1+S/T_z)); // Viscosity [kg/(m*s)]
			nu      = (8*Pi*freq*mu)/(3*P_z);                   // Nondimensional frequency
			 
			//-------- Gas fraction polynomial fits -----------------------------------
			if (z > 90.)                                         // O2 profile
				X[0] = pow(10,49.296-(1.5524*z)+(1.8714E-2*pow(z,2))-(1.1069E-4*pow(z,3))+(3.199E-7*pow(z,4))-(3.6211E-10*pow(z,5)));
			else
				X[0] = pow(10,-0.67887);

			if (z > 76.)                                         // N2 profile
				X[1] = pow(10,(1.3972E-1)-(5.6269E-3*z)+(3.9407E-5*pow(z,2))-(1.0737E-7*pow(z,3)));
			else
				X[1] = pow(10,-0.10744);

			X[2] = pow(10,-3.3979);                              // CO2 profile

			if (z > 80. )                                        // O3 profile
				X[3] = pow(10,-4.234-(3.0975E-2*z));
			else
				X[3] = pow(10,-19.027+(1.3093*z)-(4.6496E-2*pow(z,2))+(7.8543E-4*pow(z,3))-(6.5169E-6*pow(z,4))+(2.1343E-8*pow(z,5)));

			if (z > 95. )                                        // O profile
				X[4] = pow(10,-3.2456+(4.6642E-2*z)-(2.6894E-4*pow(z,2))+(5.264E-7*pow(z,3)));
			else
				X[4] = pow(10,-11.195+(1.5408E-1*z)-(1.4348E-3*pow(z,2))+(1.0166E-5*pow(z,3)));

							                                             // N profile 
			X[5]  = pow(10,-53.746+(1.5439*z)-(1.8824E-2*pow(z,2))+(1.1587E-4*pow(z,3))-(3.5399E-7*pow(z,4))+(4.2609E-10*pow(z,5)));

			if (z > 30. )                                         // H2O profile
				X[6] = pow(10,-4.2563+(7.6245E-2*z)-(2.1824E-3*pow(z,2))-(2.3010E-6*pow(z,3))+(2.4265E-7*pow(z,4))-(1.2500E-09*pow(z,5)));
			else 
			{
				if (z > 100.)
				 X[6] = pow(10,-0.62534-(8.3665E-2*z));
				else
				 X[6] = pow(10,-1.7491+(4.4986E-2*z)-(6.8549E-2*pow(z,2))+(5.4639E-3*pow(z,3))-(1.5539E-4*pow(z,4))+(1.5063E-06*pow(z,5)));
			}
			X_ON = (X[0] + X[1])/0.9903;

			//-------- Rotational collision number-------------------------------------
			Z_rot[0] = 54.1*exp(-17.3*(pow(T_z,-1./3.)));   // O2
			Z_rot[1] = 63.3*exp(-16.7*(pow(T_z,-1./3.)));   // N2
			Z_rot_   = 1./((X[1]/Z_rot[1])+(X[0]/Z_rot[0]));

			//-------- Nondimensional atmospheric quantities---------------------------
			sigma = 5./sqrt(21.);
			nn    = (4./5.)*sqrt(3./7.)*Z_rot_;
			chi   = 3.*nn*nu/4.;
			cchi  = 2.36*chi;

			//---------Classical + rotational loss/dispersion--------------------------
			beta_0  = 2*Pi*freq/c_snd_z; 
			beta_1  = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/(1+pow(nu,2)));
			beta_2  = beta_0*sqrt((1+pow(chi,2))/(1+pow((sigma*chi),2))); 
			alpha_1 = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))-1)/(1+pow(nu,2)));
			alpha_2 = beta_0*(((sigma/2-1/(2*sigma))*chi)/(sqrt((1+pow(chi,2))*(1+pow(sigma*chi,2)))));
			//a_cl    = alpha_1*(beta_2/beta_0);
			//a_rot = alpha_2*(beta_1/beta_0)*X_ON;

			a_cl    = (2*Pi*freq/c_snd_z)*sqrt(0.5*(sqrt(1+pow(nu,2))-1)*(1+pow(cchi,2))/((1+pow(nu,2))*(1+pow(sigma*cchi,2))));
			a_rot   = (2*Pi*freq/c_snd_z)*X_ON*((pow(sigma,2)-1)*chi/(2*sigma))*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/((1+pow(nu,2))*(1+pow(cchi,2))));
			a_diff  = 0.003*a_cl;

			//---------Vibrational relaxation-------------------------------------------
			Tr = pow(T_z/T_o,-1./3.)-1;
			A1 = (X[0]+X[1])*24*exp(-9.16*Tr);
			A2 = (X[4]+X[5])*2400;
			B  = 40400*exp(10*Tr);
			C  = 0.02*exp(-11.2*Tr);
			D  = 0.391*exp(8.41*Tr);
			E  = 9*exp(-19.9*Tr);
			F  = 60000;
			G  = 28000*exp(-4.17*Tr);
			H  = 22000*exp(-7.68*Tr);
			I  = 15100*exp(-10.4*Tr);
			J  = 11500*exp(-9.17*Tr);
			K  = (8.48E08)*exp(9.17*Tr);
			L  = exp(-7.72*Tr);
			ZZ = H*X[2]+I*(X[0]+0.5*X[4])+J*(X[1]+0.5*X[5])+K*(X[6]+X[3]);
			hu = 100*(X[3]+X[6]);
			f_vib[0] = (P_z/P_o)*(mu_o/mu)*(A1+A2+B*hu*(C+hu)*(D+hu));
			f_vib[1] = (P_z/P_o)*(mu_o/mu)*(E+F*X[3]+G*X[6]);
			f_vib[2] = (P_z/P_o)*(mu_o/mu)*ZZ;
			f_vib[3] = (P_z/P_o)*(mu_o/mu)*(1.2E5)*L;

			a_vib = 0.;
			for (m=0; m<4; m++)
			{
				C_R        = ((pow(theta[m]/T_z,2))*exp(-theta[m]/T_z))/(pow(1-exp(-theta[m]/T_z),2));
				A_max      = (X[m]*(Pi/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
				a_vib_c[m] = (A_max/c_snd_z)*((2*(pow(freq,2))/f_vib[m])/(1+pow(freq/f_vib[m],2)));
				a_vib      = a_vib + a_vib_c[m];
			}

			alpha[ii] = a_cl + a_rot + a_diff + a_vib;
  }
  return 0;
}
*/





// updated getAbsorption function: bug fixed by Joel and Jelle - Jun 2012
// updated 20160330: will accept attenuation coeff. loaded from a file
int NCPA::SolveModBB::getAbsorption(int n, double dz, SampledProfile *p, double freq, string usrattfile, double *alpha)
{

  if (usrattfile.empty()) {
    // Expressions based on Bass and Sutherland, JASA 2004
    // Computes alpha(freq) for given G2S output
    // Subroutine can easily be modified to include dispersion effects
    // In that case, make alpha a complex array and
    // use the real part as absorption and the imaginary part for dispersion
    int    m, ii;
    double T_o, P_o, S, z;
    double X[7], X_ON, Z_rot[2], Z_rot_;
    double sigma, nn, chi, cchi, mu, nu, mu_o;
    //double beta_0, beta_1, beta_2, alpha_1, alpha_2;
    double a_cl, a_rot, a_diff, a_vib;
    double T_z, P_z, c_snd_z, gamma;
    double A1, A2, B, C, D, E, F, G, H, I, J, K, L, ZZ, hu;
    double f_vib[4], a_vib_c[4], Cp_R[4], Cv_R[4], theta[4], C_R, A_max, Tr;
   
    //tweak absorption  - DV 20130905 
    double tweak_abs = 1; //1; //0.3; // tweak absorption alpha by this factor
    printf("Using Sutherland-Bass absorption tweaked by a factor of %g\n", tweak_abs);    

    // Atmospheric composition constants
    mu_o  = 18.192E-6;    // Reference viscosity [kg/(m*s)]
    T_o   = T[0];         // Reference temperature [K]
    P_o   = Pr[0];        // Reference pressure [Pa]
    S     = 117;          // Sutherland constant [K]       

    Cv_R[0] = 5/2;                                    // Heat capacity|volume (O2)
    Cv_R[1] = 5/2;                                    // Heat capacity|volume (N2)
    Cv_R[2] = 3;                                      // Heat capacity|volume (CO2)
    Cv_R[3] = 3;                                      // Heat capacity|volume (O3)
    Cp_R[0] = 7/2;                                    // Heat capacity|pressure (O2)
    Cp_R[1] = 7/2;                                    // Heat capacity|pressure (N2)
    Cp_R[2] = 4;                                      // Heat capacity|pressure (CO2)
    Cp_R[3] = 4;                                      // Heat capacity|pressure (O3)
    theta[0]= 2239.1;                                 // Charact. temperature (O2)
    theta[1]= 3352;                                   // Charact. temperature (N2)
    theta[2]= 915;                                    // Charact. temperature (CO2)
    theta[3]= 1037;                                   // Charact. temperature (O3)

	  //gamma   = 1.371 + 2.46E-04*T_z - 6.436E-07*pow(T_z,2) + 5.2E-10*pow(T_z,3) - 1.796E-13*pow(T_z,4) + 2.182E-17*pow(T_z,5);
	  gamma   = 1.4;
			   
    for (ii=0; ii<n; ii++) {
			  z       = ii*dz/1000.0;	// km	AGL		
			  T_z     = T[ii];	      // K
			  P_z     = Pr[ii];	      // Pa;
			  c_snd_z = sqrt(gamma*P_z/rho[ii]);  // in m/s					 
			  mu      = mu_o*sqrt(T_z/T_o)*((1+S/T_o)/(1+S/T_z)); // Viscosity [kg/(m*s)]
			  nu      = (8*Pi*freq*mu)/(3*P_z);                   // Nondimensional frequency
			   
			  //-------- Gas fraction polynomial fits -----------------------------------
			  if (z > 90.)                                         // O2 profile
				  X[0] = pow(10,49.296-(1.5524*z)+(1.8714E-2*pow(z,2))-(1.1069E-4*pow(z,3))+(3.199E-7*pow(z,4))-(3.6211E-10*pow(z,5)));
			  else
				  X[0] = pow(10,-0.67887);

			  if (z > 76.)                                         // N2 profile
				  X[1] = pow(10,(1.3972E-1)-(5.6269E-3*z)+(3.9407E-5*pow(z,2))-(1.0737E-7*pow(z,3)));
			  else
				  X[1] = pow(10,-0.10744);

			  X[2] = pow(10,-3.3979);                              // CO2 profile

			  if (z > 80. )                                        // O3 profile
				  X[3] = pow(10,-4.234-(3.0975E-2*z));
			  else
				  X[3] = pow(10,-19.027+(1.3093*z)-(4.6496E-2*pow(z,2))+(7.8543E-4*pow(z,3))-(6.5169E-6*pow(z,4))+(2.1343E-8*pow(z,5)));

			  if (z > 95. )                                        // O profile
				  X[4] = pow(10,-3.2456+(4.6642E-2*z)-(2.6894E-4*pow(z,2))+(5.264E-7*pow(z,3)));
			  else
				  X[4] = pow(10,-11.195+(1.5408E-1*z)-(1.4348E-3*pow(z,2))+(1.0166E-5*pow(z,3)));

							                                               // N profile 
			  X[5]  = pow(10,-53.746+(1.5439*z)-(1.8824E-2*pow(z,2))+(1.1587E-4*pow(z,3))-(3.5399E-7*pow(z,4))+(4.2609E-10*pow(z,5)));

			  if (z > 30. )                                         // H2O profile
				  X[6] = pow(10,-4.2563+(7.6245E-2*z)-(2.1824E-3*pow(z,2))-(2.3010E-6*pow(z,3))+(2.4265E-7*pow(z,4))-(1.2500E-09*pow(z,5)));
			  else 
			  {
				  if (z > 100.)
				   X[6] = pow(10,-0.62534-(8.3665E-2*z));
				  else
				   X[6] = pow(10,-1.7491+(4.4986E-2*z)-(6.8549E-2*pow(z,2))+(5.4639E-3*pow(z,3))-(1.5539E-4*pow(z,4))+(1.5063E-06*pow(z,5)));
			  }
			  X_ON = (X[0] + X[1])/0.9903;

			  //-------- Rotational collision number-------------------------------------
			  Z_rot[0] = 54.1*exp(-17.3*(pow(T_z,-1./3.)));   // O2
			  Z_rot[1] = 63.3*exp(-16.7*(pow(T_z,-1./3.)));   // N2
			  Z_rot_   = 1./((X[1]/Z_rot[1])+(X[0]/Z_rot[0]));

			  //-------- Nondimensional atmospheric quantities---------------------------
			  sigma = 5./sqrt(21.);
			  nn    = (4./5.)*sqrt(3./7.)*Z_rot_;
			  chi   = 3.*nn*nu/4.;
			  cchi  = 2.36*chi;

			  //---------Classical + rotational loss/dispersion--------------------------
			  //beta_0  = 2*Pi*freq/c_snd_z; 
			  //beta_1  = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/(1+pow(nu,2)));
			  //beta_2  = beta_0*sqrt((1+pow(chi,2))/(1+pow((sigma*chi),2))); 
			  //alpha_1 = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))-1)/(1+pow(nu,2)));
			  //alpha_2 = beta_0*(((sigma/2-1/(2*sigma))*chi)/(sqrt((1+pow(chi,2))*(1+pow(sigma*chi,2)))));
			  //a_cl    = alpha_1*(beta_2/beta_0);
			  //a_rot = alpha_2*(beta_1/beta_0)*X_ON;

			  a_cl    = (2*Pi*freq/c_snd_z)*sqrt(0.5*(sqrt(1+pow(nu,2))-1)*(1+pow(cchi,2))/((1+pow(nu,2))*(1+pow(sigma*cchi,2))));
			  a_rot   = (2*Pi*freq/c_snd_z)*X_ON*((pow(sigma,2)-1)*chi/(2*sigma))*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/((1+pow(nu,2))*(1+pow(cchi,2))));
			  a_diff  = 0.003*a_cl;

			  //---------Vibrational relaxation-------------------------------------------
			  Tr = pow(T_z/T_o,-1./3.)-1;
			  A1 = (X[0]+X[1])*24*exp(-9.16*Tr);
			  A2 = (X[4]+X[5])*2400;
			  B  = 40400*exp(10*Tr);
			  C  = 0.02*exp(-11.2*Tr);
			  D  = 0.391*exp(8.41*Tr);
			  E  = 9*exp(-19.9*Tr);
			  F  = 60000;
			  G  = 28000*exp(-4.17*Tr);
			  H  = 22000*exp(-7.68*Tr);
			  I  = 15100*exp(-10.4*Tr);
			  J  = 11500*exp(-9.17*Tr);
			  K  = (8.48E08)*exp(9.17*Tr);
			  L  = exp(-7.72*Tr);
			  ZZ = H*X[2]+I*(X[0]+0.5*X[4])+J*(X[1]+0.5*X[5])+K*(X[6]+X[3]);
			  hu = 100*(X[3]+X[6]);
			  f_vib[0] = (P_z/P_o)*(mu_o/mu)*(A1+A2+B*hu*(C+hu)*(D+hu));
			  f_vib[1] = (P_z/P_o)*(mu_o/mu)*(E+F*X[3]+G*X[6]);
			  f_vib[2] = (P_z/P_o)*(mu_o/mu)*ZZ;
			  f_vib[3] = (P_z/P_o)*(mu_o/mu)*(1.2E5)*L;

			  a_vib = 0.;
			  for (m=0; m<4; m++)
			  {
				  C_R        = ((pow(theta[m]/T_z,2))*exp(-theta[m]/T_z))/(pow(1-exp(-theta[m]/T_z),2));
				  A_max      = (X[m]*(Pi/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
				  a_vib_c[m] = (A_max/c_snd_z)*((2*(pow(freq,2))/f_vib[m])/(1+pow(freq/f_vib[m],2)));
				  a_vib      = a_vib + a_vib_c[m];
			  }

			  alpha[ii] = (a_cl + a_rot + a_diff + a_vib)*tweak_abs;
    }
	}
	else { 
	    // load atten. coeff from a text file with columns: z (km AGL) | attn
	    cout << "Loading attenuation coefficients from file " << usrattfile << endl;
	
	    double d1, d2;
	    vector<double> zatt, att;
	    ifstream indata;

	    indata.open(usrattfile.c_str());	
	    if (!indata) {
	        cerr << "file " << usrattfile << " could not be opened.\n";
	        exit(1);
	    }
	
	    indata >> d1 >> d2; // read first line in file
	    while (!indata.eof()) {
	        zatt.push_back(d1*1000.0);
	        att.push_back(d2);
	        indata >> d1 >> d2;
	    }
	    indata.close();
	    
	    //// print data
	    //for (unsigned i=0; i<zatt.size(); i++) {
	    //    cout << i << "  " << zatt[i] << endl;
	    //}

	    // do spline interpolation of attenuation coefficient
	    gsl_interp_accel *acc_;
	    acc_ = gsl_interp_accel_alloc();
	    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, zatt.size());
	    gsl_spline_init(spline, zatt.data(), att.data(), zatt.size());
	    
	    // do spline up to the max height hatt found in the provided atten. file
	    int it = 0;
	    while ((it*dz< zatt[zatt.size()-1]) && (it<n)) {
	        alpha[it] = gsl_spline_eval(spline, it*dz, acc_ );
	        it++;
	    }
	    // if the grid extends above hatt => extend alpha with the last value up to grid height 
      for (int i=it; i< n; i++) {
          alpha[i] = alpha[it-1];
      }
            
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc_);
  }
   
  /*
  // save alpha?
  if (1) {
      FILE *fp = fopen("attn.dat", "w");
      for (int ii=0; ii<n; ii++) {
          fprintf(fp, "%9.4f  %14.6e\n", ii*dz/1000.0, alpha[ii]);
      }
      printf("Attenuation coeff. saved in 'attn.dat'\n");
      fclose(fp);
  }
  */
  
  
  return 0;
}











// this version has   double sourceheight, double receiverheight, as arguments              	
int NCPA::SolveModBB::getModalTraceModESS(\
            int nz, double z_min, double sourceheight, double receiverheight, \
            double dz, SampledProfile *p, \
            double admittance, double freq, double azi, double *diag, \
            double *k_min, double *k_max, bool turnoff_WKB) 
{
  // DV Note: Claus's profile->ceff() computes ceff = sqrt(gamma*R*T) + wind; 
  // this version of getModalTraceModESS computes ceff = sqrt(gamma*P/rho) + wind.
  
  // use atmospherics input for the trace of the matrix.
  // the vector diag can be used to solve the modal problem
  // also returns the bounds on the wavenumber spectrum, [k_min,k_max]
  int    i, top;
  double azi_rad, z_km, z_min_km, dz_km, omega, gamma, bnd_cnd; 
  double cz, windz, ceffmin, ceffmax, ceff_grnd, cefftop; 
  double kk, dkk, k_eff, k_gnd, k_max_full, wkbIntegral, wkbTerm;
  //double rho_factor;
  z_min_km = z_min/1000.0;
  dz_km    = dz/1000.0;
  omega    = 2*Pi*freq;
  
  azi_rad  = p->getPropagationAzimuth()*Pi/180.0;
  if (fabs(azi_rad-azi*Pi/180.0)>1e-12) {
      std::ostringstream es;
      es << "Azimuth conflict; requested was " << azi << " but the actual is "
         << p->getPropagationAzimuth() << " degrees" << endl; 
      throw invalid_argument(es.str());   
  }
  
  double *ceffz;
  ceffz = new double [nz];
  gamma = 1.4;  
  // gamma = 1.371 + 2.46E-04*T - 6.436E-07*pow(T,2) + 5.2E-10*pow(T,3) - 1.796E-13*pow(T,4) + 2.182E-17*pow(T,5);

  z_km     = z_min_km;
  cz       = sqrt(gamma*Pr[0]/rho[0]); // in m/s
  windz    = zw[0]*sin(azi_rad) + mw[0]*cos(azi_rad);
  ceff_grnd = cz + windz;
  ceffmin  = ceff_grnd;  // in m/s; initialize ceffmin
  ceffmax  = ceffmin;    // initialize ceffmax   
  for (i=0; i<nz; i++) {	     
      cz       = sqrt(gamma*Pr[i]/rho[i]); // in m/s
      windz    = zw[i]*sin(azi_rad) + mw[i]*cos(azi_rad);      
      ceffz[i] = cz + windz;
      
      // we neglect rho_factor - it does not make sense to have in this approximation
      //rho_factor is 1/2*rho_0"/rho_0 - 3/4*(rho_0')^2/rho_0^2
      //rho_factor = (-0.75*pow(p->drhodz(z_km)/p->rho(z_km),2)+ 0.5*p->ddrhodzdz(z_km)/p->rho(z_km))*1e-6; // in meters^(-2)!
      //diag[i]    = pow(omega/ceffz[i],2) + rho_factor;
      diag[i] = pow(omega/ceffz[i],2);

      if (ceffz[i] < ceffmin)
       		ceffmin = ceffz[i];  // approximation to find minimum sound speed in problem		   		
      if (ceffz[i] > ceffmax)
		      ceffmax = ceffz[i];

      z_km += dz_km;		
  }
  
  bnd_cnd = (1./(dz*admittance+1))/(pow(dz,2)); // bnd cnd assuming centered fd
  diag[0] = bnd_cnd + diag[0];
  
  
  // use WKB trick for ground to ground propagation.
  if ((fabs(sourceheight)<1.0e-3) && (fabs(receiverheight)<1.0e-3) && (!turnoff_WKB)) {
    //
    // WKB trick for ground to ground propagation. 
    // Cut off lower phasespeed (highest wavenumber) where tunneling 
    // is insignificant (freq. dependent)
    //
    k_max_full = omega/ceffmin; 
    k_gnd      = omega/ceff_grnd;
    //k_gnd      = omega/(p->ceff(z_min_km, azi)*1000.0);
    dkk        = (pow(k_max_full,2)-pow(k_gnd,2))/100.0;
    
    //
    // dv: when ceffmin is at the ground dkk can be very small but non-zero (rounding errors?)
    // in that case we need to skip the next {for loop}; otherwise it will be executed
    // with the small increment dkk
    //
    kk = pow(k_gnd,2); //initialization of kk in case the loop is skipped
    if (dkk >1.0e-10) {  // do this loop only if dkk is nonzero (significantly)
        for (kk = pow(k_gnd,2); kk < pow(k_max_full,2); kk=kk+dkk) {
            wkbIntegral = 0.0;
            wkbTerm     = 1.0;  
            i           = 0;
            z_km        = z_min_km;
            while (wkbTerm > dkk) {
                k_eff       = omega/ceffz[i];
                //wkbTerm     = kk - pow(k_eff,2);
                wkbTerm     = abs(kk - pow(k_eff,2));
                wkbIntegral = wkbIntegral + dz*sqrt(wkbTerm); // dz should be in meters					
                i++;
                z_km += dz_km;
            } 

            if (wkbIntegral >= 10.0) {
                printf("\nWKB fix: new phasevelocity minimum= %6.2f m/s (was %6.2f m/s)\n", \
                				omega/sqrt(kk), omega/k_max_full);
                // printf("WKBIntegral= %12.7f at z = %6.2f km\n", wkbIntegral, z_km);          				
                break;
            }
        }
    }  
    
    *k_max = sqrt(kk);  // use this for ground level 1D Tloss 
						//(uses WKB trick to include only non-vanishing modes at the ground)
    // *k_max = omega/ceffmin; 
  }
  else { // not ground-to-ground propagation
      *k_max = omega/ceffmin; // same as k_max_full
  }  

  top     = nz - ((int) nz/10);
  z_km    = z_min_km + (top+1)*dz_km;
  cz      = sqrt(gamma*Pr[top+1]/rho[top+1]); // in m/s
  windz   = zw[top+1]*sin(azi_rad) + mw[top+1]*cos(azi_rad);
  cefftop = cz + windz;
  *k_min  = omega/cefftop;
   
  // optional save ceff
  if (0) {
      double *target;
      target = new double [p->nz()];
      p->get_ceff(target, p->nz());
      
      FILE *fp = fopen("ceff.nm", "w");
      printf("fp=%p\n", fp);
      
      for (i=0; i<p->nz(); i++) {     
          z_km = p->z(i);
          //printf("%g target[%d] = %g; %g\n", z_km, i, target[i]*1000.0, ceffz[i]);
          fprintf(fp, "%g %15.6e  %15.6e\n", z_km, target[i], ceffz[i]);
      }
      fclose(fp);
      printf("ceff saved in ceff.nm\n");
      delete [] target;
  }
  
  delete [] ceffz;
  return 0;
} // end of getModalTraceModESS



// 20130114 - added global arrays for each atm profiles such as Hgt, zw, mw, T, rho, Pr
// 20121113 - DV modified to recompute ceff=sqrt(gamma*P/rho); 
// this could also be ensured to come out of Claus's p->ceff() - but currently p->ceff gives sqrt(gamma*R*T)
int NCPA::SolveModBB::getModalTraceWMod(int nz, double z_min, double sourceheight, double receiverheight, double dz, SampledProfile *p, double admittance, double freq, double *diag, double *kd, double *md, double *cd, double *k_min, double *k_max, bool turnoff_WKB)
{
  // DV Note: Claus's profile->ceff() computes ceff = sqrt(gamma*R*T) + wind; 
  // this version of getModalTraceWMod computes ceff = sqrt(gamma*P/rho) + wind.  
    
  //     use atmospherics input for the trace of the matrix.
  //     the vector diag can be used to solve the modal problem
  //     also returns the bounds on the wavenumber spectrum, [k_min,k_max]

  int i, top;
  double azi_rad, gamma, omega, bnd_cnd, windz, ceffmin, ceffmax, ceff_grnd, cefftop, cz;
  double z_min_km, z_km, dz_km;
  double kk, dkk, k_eff, k_gnd, k_max_full, wkbIntegral, wkbTerm;
  //double rho_factor; 
  double *ceffz;
  ceffz = new double [nz];

  //FILE *profile= fopen("profile.int", "w");

  gamma    = 1.4;
  z_min_km = z_min/1000.0;
  dz_km    = dz/1000.0;
  omega    = 2*Pi*freq;

  azi_rad = p->getPropagationAzimuth()*Pi/180.0;
  if (fabs(azi_rad-azi*Pi/180.0)>1e-12) {
      std::ostringstream es;
      es << "Azimuth conflict; requested was " << azi << " but the actual is "
         << p->getPropagationAzimuth() << " degrees" << endl; 
      throw invalid_argument(es.str());   
  }

  z_km      = z_min_km; 
  cz        = sqrt(gamma*Pr[0]/rho[0]); // in m/s
  windz     = zw[0]*sin(azi_rad) + mw[0]*cos(azi_rad); 
  
  ceff_grnd = cz + windz;
  ceffmin   = ceff_grnd;  // in m/s; initialize ceffmin
  ceffmax   = ceffmin;    // initialize ceffmax 
  
  for (i=0; i<nz; i++) {
      cz    = sqrt(gamma*Pr[i]/rho[i]); // in m/s
      windz = zw[i]*sin(azi*Pi/180.0) + mw[i]*cos(azi*Pi/180.0);
      
      ceffz[i] = cz + windz;	    
	    
	    // we neglect rho_factor - it does not make sense to have in this approximation
	    //rho_factor is 1/2*rho_0"/rho_0 - 3/4*(rho_0')^2/rho_0^2
	    //rho_factor = (-3/4*pow(p->drhodz(z_km)/p->rho(z_km),2) + 1/2*p->ddrhodzdz(z_km)/p->rho(z_km))*1e-6; // in meters^(-2)!
	    kd[i]      = pow(omega/cz,2); // + rho_factor;
	    md[i]      = pow(windz/cz,2) - 1;
	    cd[i]      = -2*omega*(windz/pow(cz,2));
	    diag[i]    = pow(omega/ceffz[i],2); // + rho_factor;
	    if (ceffz[i] < ceffmin)
	     		ceffmin = ceffz[i];     // approximation to find minimum sound speed in problem
	    if (ceffz[i] > ceffmax)
	     		ceffmax = ceffz[i];
	     
	    z_km += dz_km; 
  }
  //fclose(profile);

  bnd_cnd = (1./(dz*admittance+1))/(pow(dz,2)); // bnd cnd assuming centered fd
  diag[0] = bnd_cnd + diag[0];
  kd[0]   = bnd_cnd + kd[0];


  // use WKB trick for ground to ground propagation.
  if ((fabs(sourceheight)<1.0e-3) && (fabs(receiverheight)<1.0e-3) && (!turnoff_WKB)) {
      //
      // WKB trick for ground to ground propagation. 
      // Cut off lower phasespeed (highest wavenumber) where tunneling 
      // is insignificant (freq. dependent)
      //
      k_max_full = omega/ceffmin;
      k_gnd      = omega/ceff_grnd;
      dkk        = (pow(k_max_full,2)-pow(k_gnd,2))/100.0;

      //
      // dv: when ceffmin is at the ground dkk can be very small but non-zero (rounding errors?)
      // in that case we need to skip the next for loop; otherwise it will be executed
      // with the small increment dkk
      //
      kk = pow(k_gnd,2); //initialization of kk in case the loop is skipped

      if (dkk >1e-10) { // do this loop only if dkk is nonzero (significantly)
		      for (kk = pow(k_gnd,2); kk < pow(k_max_full,2); kk=kk+dkk) {
              i           = 0;
              wkbIntegral = 0.0;
              wkbTerm     = 1.0; 
              z_km = z_min_km; 
              while (wkbTerm > dkk) {
                  k_eff = omega/ceffz[i];
                  wkbTerm = abs(kk - pow(k_eff,2));
                  wkbIntegral = wkbIntegral + dz*sqrt(wkbTerm); // dz should be in meters
                  i++;
                  z_km += dz_km;
              } 
              if (wkbIntegral >= 10.0) {
                  printf("\n WKB fix: new phasevelocity minimum: %6.2f m/s (was %6.2f m/s); \n WKBIntegral= %12.7f at z = %6.2f km\n", omega/sqrt(kk), omega/k_max_full, wkbIntegral, z_km);
                  break;
              }
          }
      }
      
      *k_max = sqrt(kk);   // use this for ground-to-ground 1D Tloss 
											     //(uses WKB trick to include only non-vanishing modes at the ground)
	    // *k_max = k_max_full;
  }
  else { // not ground-to-ground propagation
      *k_max = omega/ceffmin; // same as k_max_full
  } 

  top     = nz - ((int) nz/10);
  z_km    = z_min_km + (top+1)*dz_km;  
  cz      = sqrt(gamma*Pr[top+1]/rho[top+1]); // in m/s
  windz   = zw[top+1]*sin(azi_rad) + mw[top+1]*cos(azi_rad); 
  cefftop = cz + windz;
  *k_min  = omega/cefftop;
  
  if (0) { // disabled DV 20170805
  // check if duct is not formed and modes exist
  if (cefftop<ceff_grnd) {
      printf(" --------------------------------------------------------------\n");
      printf(" ERROR!\n");
      printf(" It appears that ducting conditions are not formed and \n");
      printf(" modes will not exist because the effective sound speed \n");
      printf(" at the top is smaller than the sound speed at the ground. \n");
      printf("  cefftop   = %g m/s\n  ceff_grnd = %g m/s\n  ceffmin   = %g m/s\n", cefftop, ceff_grnd, ceffmin);
      printf(" Suggest increasing maxheight and checking your profile file.\n");
      printf(" --------------------------------------------------------------\n");
      exit(1);
  }
  }
  
  // optional save ceff
  if (0) {
      double *target;
      target = new double [p->nz()];
      p->get_ceff(target, p->nz());
      
      FILE *fp = fopen("ceff.wnm", "w");
      for (i=0; i<p->nz(); i++) {     
          z_km = p->z(i);
          //printf("%g target[%d] = %g; %g\n", z_km, i, target[i]*1000.0, ceffz[i]);
          fprintf(fp, "%8.3f %15.6e  %15.6e\n", z_km, target[i], ceffz[i]);
      }
      fclose(fp);
      printf("ceff saved in ceff.wnm\n");
      delete [] target;
  }
  delete [] ceffz;

  return 0;
} // end getModalTraceWMod





int NCPA::SolveModBB::getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev)
{
  int nev_max, nev_min;
  sturmCount(n,dz,diag,k_max,&nev_max);
  sturmCount(n,dz,diag,k_min,&nev_min);
  *nev = nev_max - nev_min;
  return 0;
}

int NCPA::SolveModBB::sturmCount(int n, double dz, double *diag, double k, int *cnt)
{
  double kk,pot,cup0,cup1,cup2;
  double fd_d_val, fd_o_val;
  int i,pm;

  fd_d_val = -2./pow(dz,2);             // Finite Difference Coefficient on  diagonal
  fd_o_val =  1./pow(dz,2);             // Finite Difference Coefficient off diagonal

  pm   = 0;
  kk   = k*k;
  // Rogers way to start off:
  //cup0 = 1.5*fd_d_val + dd_val(n) - kk  + k/dz     // I don't understand the last division
                                                     // but it doesn't seem to affect the # modes
  cup0 = fd_d_val + diag[n-1] - kk;
  pot  = fd_d_val + diag[n-2] - kk;
  cup1 = cup0*pot;
  if (cup0*cup1 < 0.0) { pm++; } 

  cup0=cup0/fabs(cup1);
  cup1=cup1/fabs(cup1);

  for (i=n-3; i>=0; i--) {
      pot  = fd_d_val + diag[i] - kk;
      cup2 = pot*cup1 - (pow(fd_o_val,2))*cup0;
      if (cup1*cup2 < 0.0) { pm++; }
      cup0=cup1/fabs(cup2);
      cup1=cup2/fabs(cup2);
  }
  *cnt = pm;
  return 0;
}
 
 
int NCPA::SolveModBB::doPerturb2(int nz, double z_min, double dz, int n_modes, double freq, SampledProfile *p, double *k, double **v, double *alpha, complex<double> *k_pert, double *kr, double *ki)
{
  int i, j;
  double absorption, gamma, c_T; //, z;
  double z_km, dz_km;
  double omega = 2*Pi*freq;
  complex<double> I (0.0, 1.0);
  gamma = 1.4;
  // gamma = 1.371 + 2.46E-04*T - 6.436E-07*pow(T,2) + 5.2E-10*pow(T,3) - 1.796E-13*pow(T,4) + 2.182E-17*pow(T,5);
   
  dz_km = dz/1000.0;
  for (j=0; j<n_modes; j++) {
      absorption = 0.0;
      z_km=z_min/1000.0;
      for (i=0; i<nz; i++) {
        c_T = sqrt(gamma*Pr[i]/rho[i]); // in m/s
        absorption = absorption + dz*v[i][j]*v[i][j]*(omega/c_T)*alpha[i]*2;
        z_km = z_km + dz_km;
      }			
      k_pert[j] = sqrt(k[j]*k[j] + I*absorption);
      kr[j] = real(k_pert[j]);
      ki[j] = imag(k_pert[j]);
  }
  return 0;
}


/*
// 20121114 - new doPerturb version implementing a lower absorption above 90km
// via Jelle Assink - and Joel Lonzaga
int NCPA::SolveModNB::doPerturb2(int nz, double z_min, double dz, int n_modes, double freq, \
              SampledProfile *p, double *k, double **v, \
              double *alpha, complex<double> *k_pert, double *kr, double *ki)
{
  int i, j;
  double absorption, gamma, c_T;
  double z_km, dz_km;
  double taper;
  
  double omega      = 2*Pi*freq;
  double rdx_factor = 0.3;
  double rdx_slope  = 0.25E-03;
  double rdx_height = 90.0E03;
  complex<double> I (0.0, 1.0);

  gamma = 1.4;
  dz_km = dz/1000.0;
  for (j=0; j<n_modes; j++) {
      absorption = 0.0;
      z_km=z_min/1000.0;
      for (i=0; i<nz; i++) {
          //gamma = 1.371 + 2.46E-04*T - 6.436E-07*pow(T,2) + 5.2E-10*pow(T,3) - 1.796E-13*pow(T,4) + 2.182E-17*pow(T,5);
          c_T = sqrt(gamma*(p->p(z_km)*100.0)/(p->rho(z_km)*1000.0)); // in m/s
          taper = (1-rdx_factor)/(1.0+exp(rdx_slope*((i+1)*dz-rdx_height)))+rdx_factor;
          absorption = absorption + dz*v[i][j]*v[i][j]*(omega/c_T)*taper*alpha[i]*2;
          z_km = z_km + dz_km;
      }
      k_pert[j] = sqrt(k[j]*k[j] + I*absorption);
      kr[j] = real(k_pert[j]);
      ki[j] = imag(k_pert[j]);
  }
  cout << "!!! Sutherland-Bass absorption reduced by a factor of " << rdx_factor << endl;
  return 0;
}
*/

int NCPA::SolveModBB::doSelectModESS(int nz, int n_modes, double k_min, double k_max, double *k2, double **v, double *k_s, double **v_s, int *select_modes)
{
  int cnt = -1;
  int i, j;
  double k;
  for (j=0; j<n_modes; j++) {
      k = sqrt(k2[j]);
      if ((k >= k_min) && (k <= k_max)) {
          cnt++;
          for (i=0; i<nz; i++) {
              v_s[i][cnt] = v[i][j];
          }
          k_s[cnt] = k;
      }
  }
  *select_modes = cnt + 1; //+1 because cnt started from zero
  return 0;
}


int NCPA::SolveModBB::doSelectWMod(int nz, int n_modes, double k_min, double k_max, double *kH, double **v, double *k_s, double **v_s, int *select_modes)
{
  int cnt = -1;
  int i, j;

  for (j=0; j<n_modes; j++) {
      if ((kH[j] >= k_min) && (kH[j] <= k_max)) {
          cnt++;
          for (i=0; i<nz; i++) {
              v_s[i][cnt] = v[i][j];
          }
          k_s[cnt] = kH[j];
      }
  }
  *select_modes = cnt + 1; //+1 because cnt started from zero
  return 0;
}


// new function - passing a file pointer - 20150720 DV
int NCPA::SolveModBB::writeDispersion_bb_ascii(FILE *fp, int select_modes, double dz, double z_src, double z_rcv, double freq, double *rho, complex<double> *k_pert, double **v_s)
{
  // saves freq #modes rho(@z=src), rho(@z=rcv), k_horiz, modes(@z=src), modes(@z=receiverheight)
  int i;
  int n_zsrc = (int) ceil(z_src/dz);
  int n_zrcv = (int) ceil(z_rcv/dz);

  //fprintf(fp, "%e  %d",freq,select_modes);
  fprintf(fp, "%e  %d %e %e",freq,select_modes, rho[n_zsrc], rho[n_zrcv]);
  for(i=0;i<select_modes;i++) {
      fprintf(fp,"   %.12e   %.12e",real(k_pert[i]),imag(k_pert[i]));
      fprintf(fp,"   %.12e   %.12e",(v_s[n_zsrc][i]),(v_s[n_zrcv][i]));
  }
  fprintf(fp,"\n");
  //fclose(fp); // the file is closed outside of this function after the big frequency loop
  return 0;
}


// Old prototype for writeDispersion_bb_ascii: decommissioned 20150720 and 
// replaced with a new one to allow overwriting the file not appending- DV

int NCPA::SolveModBB::writeDispersion_bb_ascii(string filen, int select_modes, double dz, double z_src, double z_rcv, double freq, double *rho, complex<double> *k_pert, double **v_s)
{
  // saves freq #modes rho(@z=src), rho(@z=rcv), k_horiz, modes(@z=src), modes(@z=receiverheight)
  int i;
  int n_zsrc = (int) ceil(z_src/dz);
  int n_zrcv = (int) ceil(z_rcv/dz);
  char dispf_name[256];

  sprintf(dispf_name,"%s", filen.c_str());
  FILE *fp = fopen(dispf_name,"a");

  //fprintf(fp, "%e  %d",freq,select_modes);
  fprintf(fp, "%e  %d %e %e",freq,select_modes, rho[n_zsrc], rho[n_zrcv]);
  for(i=0;i<select_modes;i++) {
      fprintf(fp,"   %.12e   %.12e",real(k_pert[i]),imag(k_pert[i]));
      fprintf(fp,"   %.12e   %.12e",(v_s[n_zsrc][i]),(v_s[n_zrcv][i]));
  }
  fprintf(fp,"\n");
  fclose(fp);
  return 0;
}


int NCPA::SolveModBB::writeDispersion_bb_bin3(\
					string filen, double freq, int Nfreq, double df, int Nz_grid,  \
					double z_min, int Nz_subgrid, double delZ, int NN, int Nmodes, \
					double dz, double z_src, \
					double *rho, double *kreal, double *kim, double **v_s)
{
  // saves freq #modes, (subgrid details), rho(@z=src), rho(@all z), 
  // k_horiz, modes(@z=src), modes(@ all z)
  int iz;
  int szd = sizeof(double);
  int szi = sizeof(int);
  int n_zsrc = (int) ceil(z_src/dz);
  char dispf_name[256];

  sprintf(dispf_name,"%s_%e_nm.bin", filen.c_str(),freq );
  //FILE *fp = fopen(dispf_name,"ab");
  FILE *fp = fopen(dispf_name,"wb"); // DV 20150720 - changed to 'wb'

  fwrite(&freq,        szd, 1, fp);	
  fwrite(&Nmodes,      szi, 1, fp);	
  fwrite(&Nfreq,       szi, 1, fp);
  fwrite(&df,          szd, 1, fp);
  fwrite(&Nz_grid,     szi, 1, fp);
  fwrite(&dz,          szd, 1, fp);
  fwrite(&Nz_subgrid,  szi, 1, fp);
  fwrite(&delZ,        szd, 1, fp);
  fwrite(&z_min,       szd, 1, fp);
  fwrite(&rho[n_zsrc], szd, 1, fp); 

  fwrite(kreal      , szd, Nmodes, fp);
  fwrite(kim        , szd, Nmodes, fp); 
  fwrite(v_s[n_zsrc], szd, Nmodes, fp); // write all modes at height z_src

  for (iz=0; iz<Nz_grid; iz+=NN) {      // write in total (Nmodes+1) doubles at each z
      fwrite(&rho[iz], szd, 1, fp);     // write rho at all subgrid heights
      fwrite(v_s[iz], szd, Nmodes, fp); // write all modes at all grid heights
  }
  fclose(fp);
  printf("           file %s created\n", dispf_name);
  return 0;
}
