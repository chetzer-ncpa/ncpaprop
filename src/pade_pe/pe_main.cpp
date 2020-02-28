//  Pade PE (PaPE) using the effective sound speed approximation
//  to be used with a complex ground impedance model
// 
//  Author : Jelle Assink
//  Date   : 01 October 2012
//  Doru Velea adapted the code to fit into the NCPA infrasound software architecture
//  Jan. 2013
//
// * Changelog:
// * 201306  : DV modified the get() functions that return strings
// * 20131111: DV modified function load_NthProfile() to abort run if profile files are not of the same format/size;
//
//  
// Compile with:
// g++    -c -o binaryreader.o binaryreader.cpp
// g++ -c -Wall -I../common -I../atmosphere atmlib.cpp
// g++ -c -Wall -I../common -I../atmosphere ProcessOptionsPE.cpp
// g++ -c -Wall -I../common -I../atmosphere pape_main.cpp	
// g++ binaryreader.o atmlib.o ProcessOptionsPE.o pe_main.o -lgsl  -lgslcblas  ../lib/libatmosphere.a  ../lib/libcommon.a -o pape
// cp pape ../bin

#include <iostream>
#include <fstream>
#include <ostream>
#include <complex>
#include <stdlib.h>
#include <stdexcept>
#include <math.h>
#include <time.h>
#include "atmlib.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <string>
#include <vector>
#include <dirent.h>
#include <list>
 
#include "anyoption.h"
#include "ProcessOptionsPE.h"
#include "Atmosphere.h"

using namespace NCPA;
using namespace std;

typedef unsigned int uint;

// PE functions
void buildBCVectors(complex<double>*, complex<double>*, complex<double>*, complex<double>*, complex<double>**, complex<double>**, complex<double>**, complex<double>**, complex<double>**);
void buildQOperatorVectors(complex<double>,complex<double>*,complex<double>*);
void buildAbsorptiveLayer();
void getSqrtPadeCoefficients(complex<double>*,complex<double>*);
//void getStarterField(string,complex<double>*, double, double);
void getStarterField(string, string, complex<double>*, double, double);
void rhsMultiplication(int,complex<double>**,complex<double>**,complex<double>*);
void tridagSolver(int,complex<double>**,complex<double>**,complex<double>**,complex<double>*,complex<double>*);
void marchField(complex<double>*,complex<double>*,double*);
void writeField(complex<double>*,int,int,int,int,int);

// atmosphere, interpolation and allocation functions
void loadG2S_1DAtmosphere(std::ifstream*);
void loadG2S_2DAtmosphere(std::ifstream*);
void copy1DAtmosphere(int);
void getG2SRangeIndex(int, int*);
void loadToyAtmosphere();
void initAtmosphere(ProcessOptionsPE *oPE);
void getImpedance(double,complex<double>*);

int  checkSampling();
void initInterpolation();
void reinitInterpolation();
void doInterpolation();
//void freeInterpolation();
void freeGlobals(int filetype);
//void allocateVectors(double*,double*,complex<double>*,complex<double>*,complex<double>*,complex<double>*,complex<double>*,complex<double>*,complex<double>*);
void deleteVectors(complex<double>*,complex<double>*,complex<double>*,complex<double>*,complex<double>**,complex<double>**,complex<double>**,complex<double>**,complex<double>**);

// common functions
void getCMDParameters(int,char**);
void printUsage();

// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

int getRegionBoundaries(bool flg, double maxrange, double req_profile_step, string prf_ranges_km, int *Nprofiles, vector<double> *R);
void parseReqRanges(std::string str, std::vector<double>& retVal);
//int plotwGNUplot(double freq, bool write_2D_TLoss);
int get_RngDepndProfiles_ascii(int N, string atmosfileorder, uint skiplines, \
    string dirname, string pattern, string windspeed_units, string &atmosfile);
int getFileList(std::string dir, std::list<string> &files, std::string pattern);
void load_1DAtm_ascii(string atmosfile, string atmosfileorder, int skiplines, string wind_units);
int load_2DAtm_ascii(vector<double> Rv, string dirname, string pattern, string atmosfileorder, int skiplines, string wind_units);
void load_NthProfile(int J, string atmosfile, string atmosfileorder, int skiplines, string wind_units, int Nz0);
void get2DAtmRangeIndex(int rr, int *atm_r_index);
//double getMaxheightAtmProfile(ProcessOptionsPE *oPE);


////////////////////////////////////////////////////////////////////
// Global parameters

double  zsrc, zrcv, azi, freq, c0, zmin, zmax, rmax, dz, dr, rng_step;
int     atm_nz, atm_nr, nz, nr, nzrcv, n_pade;

double  *alt, *T, *rho, *pr, *zw, *mw, *atm_rng;
double  **T_2D, **rho_2D, **pr_2D, **zw_2D, **mw_2D;

double  *alt_int, *T_int, *rho_int, *pr_int, *zw_int, *mw_int;
double  *c_int, *abs_sb, *abs_layer;

string wind_units;

complex<double> I (0.0,1.0);
vector<double> Rv(20,0.0);

FILE * fid_tloss1d;
FILE * fid_tloss2d;

gsl_interp_accel *acc_T, *acc_mw, *acc_zw, *acc_rho, *acc_pr;
gsl_spline       *T_fit, *mw_fit, *zw_fit, *rho_fit, *pr_fit;

////////////////////////////////////////////////////////////////////

int main( int nargin, char **argv )
{
  time_t clock_1 = time(NULL); // start timer
   
  // parse options from the command line as well as an options file
  AnyOption *opt = parseInputOptions(nargin, argv); 

  // object to process the options
  ProcessOptionsPE *oPE;
  oPE = new ProcessOptionsPE(opt);
  
  // get parameters; defaults are specified in ProcessOptionsPE
  azi      = oPE->getAzimuth();
  freq     = oPE->getFreq();  
  //zmax     = oPE->getMaxheight();
  zmin     = oPE->getZ_min();  // ground level above MSL
  nz       = oPE->getNz_grid();
  zsrc     = oPE->getSourceheight();
  zrcv     = oPE->getReceiverheight();
  rmax     = oPE->getMaxrange();
  rng_step = oPE->getRngStep();
  n_pade   = oPE->getNpade();
  wind_units = oPE->getWindUnits();
  int plot2d = oPE->getWrite_2D_TLoss();
  
  string starter_type, modstartfile;
  starter_type = oPE->getStarterType(); // gaussian/greene/modal
  modstartfile = oPE->getModalStarterFile(); // filename of pre-computed modal starter
  int filetype = oPE->getFiletype();
  
  
  if (filetype==3) {          // if ascii profiles available in a directory
      string prf_ranges_km;   // string specifying profile ranges
      int Nprofiles = 0;
      prf_ranges_km = oPE->getProfileRanges();
      getRegionBoundaries(oPE->getProfile_ranges_given_flag(), \
                          oPE->getMaxrange(), oPE->getReq_profile_step(), \
                          prf_ranges_km, &Nprofiles, &Rv);
      //cout << "Nprofiles = " << Nprofiles << endl;
      //cout << "Rv size   = " << Rv.size() << endl;
      //for (int i=0; i< (int) Rv.size(); i++) {
      //    printf("Rv[%i] = %g\n", i, Rv[i]);
      //}
  }

  // initialize atmospheric profile; 
  // Note! the global zmax is also set in this call; zmax may be readjusted - see logic below
  initAtmosphere(oPE);
  
  // adjust maxheight (zmax) to the minimum between the --maxheight option and the 
  // max height from the provided atmospheric profile
  zmax = min(zmax,oPE->getMaxheight());
  oPE->setMaxheight(zmax);
  
  dz     = (zmax - zmin)/nz; // note ref. to ground level (zmin)
  nzrcv  = (zrcv/dz);
  c0     = 340.0;
  dr     = (c0/freq)*rng_step;
  if (dr>1000.0) {
      dr = 1000.0;
      printf("Note!! range step reduced to dr = %g from %g m\n", dr, (c0/freq)*rng_step);
  }
  nr     = rmax/dr;

  printf("\n");
  printf("High Angle PE\n");
  printf(" -> Azimuth : %.2f degrees\n", azi);
  printf(" -> Frequency: %.2f Hz\n", freq);
  printf(" -> Source height (AGL): %.2f km\n", zsrc/1000);
  printf(" -> Receiver height(AGL): %.2f km\n", zrcv/1000);
  printf(" -> Ground level (MSL): %.2f km\n", zmin/1000);
  printf(" -> Range step: dr = %.2f m\n", dr);

  int showR = 50.0E3/dr; 
  int cS    = checkSampling();
  if (cS == 1) { return 0; }

  int plotr  = 1;   // save the data in range  at every (plotr*dr)
  if (plotr*dr>1000.0) {
      plotr = (int) floor(1000.0/dr);
  }
  int plotz  = 20;  // save the data in height at every (plotz*dz)
  
  /*

  if (filetype==3) {          // if ascii profiles available in a directory
      string prf_ranges_km;   // string specifying profile ranges
      int Nprofiles = 0;
      prf_ranges_km = oPE->getProfileRanges();
      getRegionBoundaries(oPE->getProfile_ranges_given_flag(), \
                          oPE->getMaxrange(), oPE->getReq_profile_step(), \
                          prf_ranges_km, &Nprofiles, &Rv);
      //cout << "Nprofiles = " << Nprofiles << endl;
      cout << "Rv size   = " << Rv.size() << endl;
      for (int i=0; i< (int) Rv.size(); i++) {
          printf("Rv[%i] = %g\n", i, Rv[i]);
      }
  }
  */

  // interpolate the atm profile(s) to the z-grid defined by zmin, zmax, nz, dz
  initInterpolation();
  doInterpolation();
  
  // select the impedance model
  string grnd_imp_model;
  grnd_imp_model = oPE->getGrnd_imp_model();

  complex<double> impedance, admittance;
  if (!strcmp(grnd_imp_model.c_str(), "rigid")) {
      // 1. option for rigid ground boundary condition
      admittance = 0.0;
  }
  else if (!strcmp(grnd_imp_model.c_str(), "soft")) {
      // 2. option for complex impedance ground boundary
      getImpedance(freq,&impedance);                                               
      //admittance = -1.0*I*(2*PI*freq)*rho_int[0]/(impedance*rho_int[0]*c_int[0])     
  }
  else {
      std::ostringstream es;
      es << "This ground impedance model is not implemented: " << grnd_imp_model
         << endl << "Available models: rigid and soft" << endl;
      throw invalid_argument(es.str());
  }

  // get attenuation
  string usrattfile;
  usrattfile = oPE->getUsrAttFile();
  AtmLibrary *atm_  = new AtmLibrary();
  atm_->getAbsorptionCoefficients(nz,freq,alt_int,T_int,pr_int,c_int,usrattfile, abs_sb);
  delete atm_;
  
  if (oPE->getNoabsorption()) {
      // zero the absorption
      cout << " -> Atmospheric absorption is zero (lossless case)." << endl;
      for (int i=0; i<nz; i++) {
          abs_sb[i] = 0.0;
      }
  }

  buildAbsorptiveLayer();

  complex<double> *Qd, *Qo, *cm_coeff, *cp_coeff, *psi_o, *psi_dr;
  complex<double> **Bo, **Bd, **Cl, **Cd, **Cu;
  Qd       = new complex<double> [ nz ];
  Qo       = new complex<double> [ nz ];
  psi_o    = new complex<double> [ nz ];
  psi_dr   = new complex<double> [ nz ];
  cm_coeff = new complex<double> [ n_pade ];
  cp_coeff = new complex<double> [ n_pade ];

  Bd = new complex<double>* [ nz ];
  Bo = new complex<double>* [ nz ];
  Cl = new complex<double>* [ nz ];
  Cd = new complex<double>* [ nz ];
  Cu = new complex<double>* [ nz ];
  for (int i = 0; i < nz; i++) {
      Bd[ i ]  = new complex<double>[ n_pade ];
      Bo[ i ]  = new complex<double>[ n_pade ];
      Cl[ i ]  = new complex<double>[ n_pade ];
      Cd[ i ]  = new complex<double>[ n_pade ];
      Cu[ i ]  = new complex<double>[ n_pade ];
  }

  //getStarterField(starter_type.c_str(),psi_o, zsrc, zrcv);
  getStarterField(starter_type, modstartfile, psi_o, zsrc, zrcv);

  // save starter field?
  if (0) {
      FILE *f;
      f = fopen("starter_field.dat", "w");
      for (int i=0; i<nz; i++) {
          fprintf(f, "%f %g %g\n", i*dz/1000.0, real(psi_o[i]), imag(psi_o[i]));
      }
      fclose(f);
      printf(" !!!  starter field saved in starter_field.dat\n");
  }
   
  if ( abs(admittance) == 0.0) { printf(" -> Setting up finite-differences (rigid ground) ...\n"); }
  else                         { printf(" -> Setting up finite-differences (complex impedance (%.2f,%.2f)) ...\n", real(impedance), imag(impedance)); }

  cout << " -> Determining Pade coefficients (" << n_pade << ")" << endl;
  getSqrtPadeCoefficients(cm_coeff,cp_coeff);
  cout << " -> Setting up operators ..." << endl;
  buildQOperatorVectors(admittance,Qd,Qo);
  buildBCVectors(cm_coeff,cp_coeff,Qd,Qo,Bd,Bo,Cl,Cd,Cu);

  cout << " -> Marching out in range ..." << endl;
                     fid_tloss1d = fopen("tloss_1d.pe","w");
  if (plot2d == 1) { fid_tloss2d = fopen("tloss_2d.pe","w"); }

  int atm_r0_index = 0;
  int atm_dr_index;
  int rr;
  //
  // big loop: marching out
  //
  for (rr=1; rr<nr; rr++) {
      if (rr % showR == 0) { printf("    -> Range %.f km\n", (rr*dr)/1000); }
      writeField(psi_o,rr,nzrcv,plot2d,plotr,plotz);

      if (filetype!=0 && filetype !=2) {       // if not range-independent 1D atmosphere
          if (filetype==1) {   //g2s env file 
              getG2SRangeIndex(rr,&atm_dr_index); // loads 1D profiles at approx midway between ranges stored in atm_rng[]
          }
          else if (filetype==3) {
              //getG2SRangeIndex(rr,&atm_dr_index);  // loads 1D profiles at approx midway between ranges stored in atm_rng[]
              get2DAtmRangeIndex(rr, &atm_dr_index); // ensures 1D profiles loaded at requested ranges

          }

          if (atm_dr_index != atm_r0_index) {
              printf(" -> using atm. profile #%d from range %g km\n", atm_dr_index, rr*dr/1000.0);
              copy1DAtmosphere(atm_dr_index);
              reinitInterpolation();
              doInterpolation();
              buildQOperatorVectors(admittance,Qd,Qo);
              buildBCVectors(cm_coeff,cp_coeff,Qd,Qo,Bd,Bo,Cl,Cd,Cu);
              atm_r0_index = atm_dr_index;
          }
      }

      for (int i_pade=0; i_pade < n_pade; i_pade++) {
          rhsMultiplication(i_pade,Bd,Bo,psi_o);
          tridagSolver(i_pade,Cl,Cd,Cu,psi_o,psi_dr);
          marchField(psi_o,psi_dr,abs_layer);
      }      
  } // end of big loop

  fclose(fid_tloss1d);
  if (plot2d == 1) { fclose(fid_tloss2d); }

  //// plot?
  //if (oPE->getPlot_flg()) {
  //   plotwGNUplot(freq, plot2d);
  //}
  
  // print run info
  oPE->printParams();
  cout << "Results saved in tloss_1d.pe" << endl;
  if (plot2d == 1) { cout << "Results saved in tloss_2d.pe" << endl; }

  //freeInterpolation(); 
  freeGlobals(filetype);
  deleteVectors(Qd,Qo,psi_o,psi_dr,Bo,Bd,Cl,Cd,Cu);
  delete opt;
  delete oPE;

  time_t clock_2 = time(NULL);
  printf("\nRun-time PE code: %.2f s\n\n", difftime(clock_2,clock_1));

  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Parabolic Equation (PE) functions

void getSqrtPadeCoefficients( complex<double> *cm, complex<double> *cp) {
  // routine to compute Pade square root approximant coefficient for sqrt(1+q) 
  double a_pade, b_pade;
  double k0 = 2*PI*freq/c0;
  for (int i_pade=0; i_pade < n_pade; i_pade++) {
      a_pade     = (2.0/(2.0*n_pade+1.0))*pow(sin((i_pade+1)*PI/(2.0*n_pade+1.0)),2);
      b_pade     = pow(cos((i_pade+1)*PI/(2.0*n_pade+1.0)),2);
      cp[i_pade] = b_pade + I*0.5*k0*dr*a_pade;
      cm[i_pade] = b_pade - I*0.5*k0*dr*a_pade;
  }
}


void getStarterField( string type, string modstartfile, complex<double> *psi_o, double src_height, double rcv_height) {
  double fct = 1.0;
  if (src_height > 1) { // source off the ground
      fct = 2.0;
  }
  //printf("fct = %g\n", fct);
  if (! strcmp(type.c_str(),"gaussian")) {
      cout << " -> Gaussian starter" << endl;

      double k0   = 2*PI*freq/c0;
      double z, A0, A2, A4, B;
      double L, lambda;
      complex<double> I (0.0,1.0);
      A0 = 1.9705;
      A2 = -1.1685;
      A4 = 0.0887;
      B  = 3.0;
      lambda = 2*PI/k0;
      L      = lambda/10;

      for (int i=0; i<nz; i++) {
          z = alt_int[i]; // this 'z' is above ground level AGL
          // if source is on the ground: fct = 1 otherwise fct = 2;
          
          // try factor sqrt(freq)/(2*PI) times eq. 6.100 pg. 366 in Oc Acoust. 
          //psi_o[i] = sqrt(freq)/(fct*PI)*sqrt(k0) * exp(-k0*k0/2*pow(z-zsrc,2)); 
          
          // 6.100 pg. 366 in Oc Acoust. divided by fct
          psi_o[i] = sqrt(k0/fct) * exp(-pow(k0/fct,2)*pow(alt_int[i]-zsrc,2)) / (2 * PI);
          
          //psi_o[i] = sqrt(k0/fac) * exp(-pow(k0/fac,2)*pow(alt_int[i]-zsrc,2)) / (2 * PI); // Jelle's original
          
          //// simply eq. 6.100 pg. 366 in Oc Acoust.  / (4 * PI);
          //psi_o[i] = sqrt(k0) * exp(-k0*k0/2*pow(z-zsrc,2)) / (4 * PI); // eq. 6.100 pg. 366 in Oc Acoust.
          
          ////  eq. 6.100 pg. 366 in Oc Acoust.  / (4 * PI) * fudge factor (freq/0.1)
          //psi_o[i] = sqrt(k0) * exp(-k0*k0/2*pow(z-zsrc,2)) / (4 * PI) * (freq/0.1);          

          // from Salomon's "Computational Atm. Acoustics"
          //psi_o[i] = sqrt(I*k0)*(A0 + A2*(k0*k0*(z-zsrc)*(z-zsrc)) + A4*pow(k0*(z-zsrc),4))*exp(-pow(k0*(z-zsrc),2)/B)/(2*PI*freq);
          
          // standard gaussian starting field  -see e.g. Gilbert, Di; Jasa Express Letter JASA 121(5), 2007 pg EL207
          //psi_o[i] = sqrt(2*PI*I/k0)/2.0*exp(-pow(((z-zsrc)/L),2))/(L*sqrt(PI));
          
          //// integral from 0 to inf of psi = 1/(2*pi)
          //psi_o[i] = (2*sqrt(k0/(2*pi)))*1/(fct*PI)*sqrt(k0) * exp(-k0*k0/2*pow(z-zsrc,2));
      }
  }
  else if (! strcmp(type.c_str(),"greene")) { // ideal for wide-angle PE (away from boundaries) 
                                              // see pages 367-369 in Ocean Acoustics 1994 ed.
      cout << " -> Greene starter" << endl;
      double k0   = 2*PI*freq/c0;
      for (int i=0; i<nz; i++) {
      
         // eq 6.101 pg. 367 in Oc Acoust. divided by fct
         // psi_o[i] = sqrt(freq)/(fct*PI)*sqrt(k0)*(1.4467-0.4201*pow(k0,2)*pow(alt_int[i]-zsrc,2))*exp(-pow(k0*(alt_int[i]-zsrc),2)/3.0512);
          
          // eq 6.101 pg. 367 in Oc Acoust
          psi_o[i] = sqrt(k0)*(1.4467-0.4201*pow(k0,2)*pow(alt_int[i]-zsrc,2))*exp(-pow(k0*(alt_int[i]-zsrc),2)/3.0512);
      }
  }
  else if (! strcmp(type.c_str(),"modal")) {
      //ifstream *starter = new ifstream( "modalstarter.nm", ios_base::in );
      ifstream *starter = new ifstream( modstartfile.c_str(), ios_base::in );
      if (!starter->good()) { cerr << " ERROR: Modal starter file does not exist. Exiting..." << endl; exit(0); } 
      cout << " -> Modal starter loaded from file '" << modstartfile <<  "'" << endl;
      float dummy;
      int n_starter = -1;
      while(!starter->eof() ) {
          *starter >> dummy >> dummy >> dummy; 
          n_starter++; 
      }
      double *x      = new double[ n_starter ];
      double *re_psi = new double [ n_starter ];
      double *im_psi = new double [ n_starter ];
      starter->clear();
      starter->seekg(0, ios::beg);
      for (int i = 0; i < n_starter; i++ ) {
          *starter >> x[i] >> re_psi[i] >> im_psi[i];
      }
      delete starter;

      // respline starter
      gsl_interp_accel* acc_msr = gsl_interp_accel_alloc();
      gsl_interp_accel* acc_msi = gsl_interp_accel_alloc();
      gsl_spline* msr_fit       = gsl_spline_alloc(gsl_interp_cspline, n_starter);
      gsl_spline* msi_fit       = gsl_spline_alloc(gsl_interp_cspline, n_starter);
      gsl_spline_init(msr_fit, x, re_psi, n_starter);
      gsl_spline_init(msi_fit, x, im_psi, n_starter);
      double msr_int, msi_int;
      for (int i=0; i< nz; i++) {
          msr_int  = gsl_spline_eval(msr_fit, alt_int[i]/1000, acc_msr  );
          msi_int  = gsl_spline_eval(msi_fit, alt_int[i]/1000, acc_msi  );
          psi_o[i] = msr_int + I*msi_int;
      }
      gsl_spline_free (msr_fit);
      gsl_interp_accel_free(acc_msr);
      gsl_spline_free (msi_fit);
      gsl_interp_accel_free(acc_msi);

      delete [] x;
      delete [] re_psi;
      delete [] im_psi;
  }
  else {
      std::ostringstream es;
      es << "This starter type is not implemented: " << type << endl
         << "Use 'gaussian', 'greene' or 'modal'. 'gaussian' is the default." << endl;
      throw invalid_argument(es.str());
  }
}

void buildBCVectors(complex<double> *cm_coeff, complex<double> *cp_coeff, complex<double> *Qd, complex<double> *Qo, complex<double> **Bd, complex<double> **Bo, complex<double> **Cl, complex<double> **Cd, complex<double> **Cu) {
  for (int i_pade = 0; i_pade < n_pade; i_pade++) { 
      for (int i=0; i<nz; i++) {
          Bd[i][i_pade]   = 1.0 + cp_coeff[i_pade]*Qd[i];
          Bo[i][i_pade]   =       cp_coeff[i_pade]*Qo[i];
          Cd[i][i_pade]   = 1.0 + cm_coeff[i_pade]*Qd[i];
          Cu[i][i_pade]   =       cm_coeff[i_pade]*Qo[i];
          if (i < (nz-1) ) { 
              Cl[i+1][i_pade] =      cm_coeff[i_pade]*Qo[i];
          }
      }
  }
}

void buildQOperatorVectors(complex<double> alpha,complex<double> *Qd, complex<double> *Qo) {
  // builds the Q operator that has the vertical operator plus omega/c squared
  // in order to expand the square root operator later, the operator is scaled :: q = (Q-k0^2) / k0^2
  // so that the vertical operator becomes sqrt(1+q)
  double omega     = 2*PI*freq;
  double k0        = omega/c0;
  int    i         = 0;
  double wind      = cos((PI/180.)*azi)*mw_int[i] + sin((PI/180.)*azi)*zw_int[i];
  double kk        = pow(omega/(c_int[i]+wind),2) - pow(k0,2);

  complex<double> bndcnd    = (1.0 / ( dz * alpha+ 1.0 ) - 2.0) / pow(dz,2);    // impedance boundary condition
  complex<double> fd_on__dg = bndcnd;
         double  fd_off_dg = 1.0/pow(dz,2);
  Qd[i]      = ( fd_on__dg + kk ) / pow(k0,2);
  Qo[i]      =   fd_off_dg / pow(k0,2);
  fd_on__dg  = -2.0/pow(dz,2);
  for (i=1; i<nz; i++) {
      wind  = cos((PI/180.)*azi)*mw_int[i] + sin((PI/180.)*azi)*zw_int[i];
      kk    = pow(omega/(c_int[i]+wind),2) - pow(k0,2);
      Qd[i] = ( fd_on__dg + kk ) / pow(k0,2);
      if (i < (nz - 1)) { Qo[i] = fd_off_dg / pow(k0,2); }
  }
}

void buildAbsorptiveLayer() {
  abs_layer  = new double [ nz ];

  double mu  = 0.5E-1;
  double z_t = alt_int[nz-1]-1000;
  for (int i=0; i< nz; i++) {
      //abs_layer[i] = mu*exp((alt_int[i]-z_t)/2500);
      //abs_layer[i] = mu*exp((alt_int[i]-z_t)/5000);
      abs_layer[i] = mu*exp((alt_int[i]-z_t)/1000); //original code
  }

  // save absorption in a file
  if (0) {
      FILE *f;
      f = fopen("abs_layer.dat", "w");
      for (int i=0; i< nz; i++) {
          fprintf(f, "%g\n", abs_layer[i]);
      }
      fclose(f);
      printf("abs_layer.dat written\n");
  }
}

void rhsMultiplication(int i_pade, complex<double> **Bd, complex<double> **Bo, complex<double> *psi_o) {
  complex<double> *work = new complex<double>[ nz ];
  work[0] = Bd[0][i_pade]*psi_o[0] + Bo[0][i_pade]*psi_o[1];
  for(int i=1; i<nz-1; i++) {
      work[i] = Bo[i-1][i_pade]*psi_o[i-1] + Bd[i][i_pade]*psi_o[i] + Bo[i][i_pade]*psi_o[i+1];
  }
  work[nz-1] = Bo[nz-2][i_pade]*psi_o[nz-2] + Bd[nz-1][i_pade]*psi_o[nz-1];
  for (int i=0; i<nz-1; i++) {
      psi_o[i] = work[i];
  }
  delete [] work;
}


void tridagSolver(int i_pade, complex<double> **Cl, complex<double> **Cd, complex<double> **Cu, complex<double> *psi_o, complex<double> *psi_dr) {
  complex<double> bet;
  complex<double> *gam = new complex<double>[ nz ];
  if (Cd[0][i_pade] == 0.0) cerr << "Error 1 in tridag" << endl;

  psi_dr[0]=psi_o[0]/(bet=Cd[0][i_pade]);
  for (int j=1;j<nz;j++) {
      gam[j]=Cu[j-1][i_pade]/bet;
      bet=Cd[j][i_pade]-Cl[j][i_pade]*gam[j];
      if (bet == 0.0) cerr << "Error 2 in tridag" << endl;
      psi_dr[j]=(psi_o[j]-Cl[j][i_pade]*psi_dr[j-1])/bet;
  }
  for (int j=(nz-2);j>=0;j--) {
      psi_dr[j] -= gam[j+1]*psi_dr[j+1];
  }

  delete [] gam;
}

void marchField(complex<double> *psi_o,complex<double> *psi_dr,double *abs_layer) {
  double damping;
  double rdx_factor = 0.3;
  double rdx_slope  = 0.25E-03;
  double rdx_height = 90.0E03;
  double taper;
  for (int i=0; i< nz; i++) {
      taper   = (1.0-rdx_factor)/(1.0+exp(rdx_slope*(alt_int[i]-rdx_height)))+rdx_factor;
      //damping = exp(-abs_layer[i]*dr);
      damping = exp(-abs_layer[i]*dr) * exp(-taper*abs_sb[i]*1.0*dr);
      
      // CHH 191029: Rewrote to comply with c++11 syntax:
      //real(psi_o[i]) = real(psi_dr[i])*damping;
      //imag(psi_o[i]) = imag(psi_dr[i])*damping;
      psi_o[ i ].real( psi_dr[ i ].real() * damping );
      psi_o[ i ].imag( psi_dr[ i ].imag() * damping );
  }
}

// Jelle's original
/*
void writeField(complex<double> *psi_o,int rr,int nzrcv,int plot2d,int plotr,int plotz) {
  double  k0   = 2*PI*freq/c0;
  double  R    = rr*dr;
  complex<double> hank = sqrt(2.0/(PI*k0*R))*exp(I*(k0*R - PI/4.0)); // eq 6.4 page 345 in Oc. Acoust.

  if (rr % plotr == 0) {
      fprintf(fid_tloss1d,"%.3f %15.8e %15.8e\n", R/1000, real(psi_o[nzrcv]*hank), imag(psi_o[nzrcv]*hank));
      if (plot2d == 1) {
          for (int i=0; i<nz; i=i+plotz) {
              fprintf(fid_tloss2d,"%.3f %.3f %15.8e %15.8e\n", R/1000, alt_int[i]/1000, real(psi_o[i]*hank), imag(psi_o[i]*hank));
          }
          fprintf(fid_tloss2d,"\n");
      }
  }
}
*/


// DV 20150929 - adjusting to get the modal starter to work and not give answers dependent on frequency

void writeField(complex<double> *psi_o,int rr,int nzrcv,int plot2d,int plotr,int plotz) {
  double  k0   = 2*PI*freq/c0;
  double  R    = rr*dr;
  complex<double> hank = sqrt(2.0/(PI*k0*R))*exp(I*(k0*R - PI/4.0)); // eq 6.4 page 345 in Oc. Acoust.
   //complex<double> hank = sqrt(1.0/(k0*R))*exp(I*(k0*R - PI/4.0)); // without factor in front? : eq 6.4 page 345 in Oc. Acoust.

  if (rr % plotr == 0) {
      fprintf(fid_tloss1d,"%.3f %15.8e %15.8e\n", R/1000, real(psi_o[nzrcv]*hank), imag(psi_o[nzrcv]*hank));
      if (plot2d == 1) {
          for (int i=0; i<nz; i=i+plotz) {
              fprintf(fid_tloss2d,"%.3f %.3f %15.8e %15.8e\n", R/1000, alt_int[i]/1000, real(psi_o[i]*hank), imag(psi_o[i]*hank));
          }
          fprintf(fid_tloss2d,"\n");
      }
  }
}




/*
// DV 20150603 - dividing the previous PE output by fct below for testing purposes
void writeField(complex<double> *psi_o,int rr,int nzrcv,int plot2d,int plotr,int plotz) {
  double  k0   = 2*PI*freq/c0;
  double  R    = rr*dr;
  double fct   = 1.0; // 4*PI;
  complex<double> hank = sqrt(2.0/(PI*k0*R))*exp(I*(k0*R - PI/4.0));
  
  //if (0) {
  //  printf("The pressure field is divided by factor = %g to agree with Modess output\n", fct);
  //}

  if (rr % plotr == 0) {
      fprintf(fid_tloss1d,"%.3f %15.8e %15.8e\n", R/1000, real(psi_o[nzrcv]*hank/fct), imag(psi_o[nzrcv]*hank/fct));
      if (plot2d == 1) {
          for (int i=0; i<nz; i=i+plotz) {
              fprintf(fid_tloss2d,"%.3f %.3f %15.8e %15.8e\n", R/1000, alt_int[i]/1000, real(psi_o[i]*hank/fct), imag(psi_o[i]*hank/fct));
          }
          fprintf(fid_tloss2d,"\n");
      }
  }
}
*/


/////////////////////////////////////////////////////////////////////////////////////////////
// Interpolation functions

void freeGlobals(int filetype) {

  // free global 1D arrays
  gsl_spline_free (pr_fit);
  gsl_interp_accel_free(acc_pr);
  gsl_spline_free (rho_fit);
  gsl_interp_accel_free(acc_rho);
  gsl_spline_free (zw_fit);
  gsl_interp_accel_free(acc_zw);
  gsl_spline_free (mw_fit);
  gsl_interp_accel_free(acc_mw);
  gsl_spline_free (T_fit);
  gsl_interp_accel_free(acc_T);

  delete [] atm_rng;
  delete [] c_int ;
  delete [] abs_layer;
  delete [] abs_sb;
  
  delete [] alt ; delete [] alt_int;
  delete [] zw  ; delete [] zw_int;
  delete [] mw  ; delete [] mw_int;
  delete [] T   ; delete [] T_int;
  delete [] rho ; delete [] rho_int;
  delete [] pr  ; delete [] pr_int;

  if (filetype!=0 && filetype !=2)  {// if 2D arrays exist (this is the range-dep atm.)
      // free the 2D global arrays
      for (int i = 0; i < atm_nz; i++) {
          delete [] T_2D[ i ];
          delete [] rho_2D [ i ];
          delete [] pr_2D  [ i ];
          delete [] zw_2D  [ i ];
          delete [] mw_2D  [ i ];
      }
      
      delete zw_2D;
      delete mw_2D;
      delete T_2D;
      delete rho_2D;
      delete pr_2D;
  }
    
}


void initInterpolation() { 
  acc_T   = gsl_interp_accel_alloc();
  acc_zw  = gsl_interp_accel_alloc();
  acc_mw  = gsl_interp_accel_alloc();
  acc_rho = gsl_interp_accel_alloc();
  acc_pr  = gsl_interp_accel_alloc();
  T_fit   = gsl_spline_alloc(gsl_interp_cspline, atm_nz);
  zw_fit  = gsl_spline_alloc(gsl_interp_cspline, atm_nz);
  mw_fit  = gsl_spline_alloc(gsl_interp_cspline, atm_nz);
  rho_fit = gsl_spline_alloc(gsl_interp_cspline, atm_nz);
  pr_fit  = gsl_spline_alloc(gsl_interp_cspline, atm_nz);
  gsl_spline_init(T_fit  , alt, T  , atm_nz);
  gsl_spline_init(zw_fit , alt, zw , atm_nz);
  gsl_spline_init(mw_fit , alt, mw , atm_nz);
  gsl_spline_init(rho_fit, alt, rho, atm_nz);
  gsl_spline_init(pr_fit , alt, pr , atm_nz);
  alt_int = new double [ nz ];
  T_int   = new double [ nz ];
  rho_int = new double [ nz ];
  pr_int  = new double [ nz ];
  zw_int  = new double [ nz ];
  mw_int  = new double [ nz ];
  c_int   = new double [ nz ];
  abs_sb  = new double [ nz ];
}

void reinitInterpolation() {
  gsl_spline_init(T_fit  , alt, T  , atm_nz);
  gsl_spline_init(zw_fit , alt, zw , atm_nz);
  gsl_spline_init(mw_fit , alt, mw , atm_nz);
  gsl_spline_init(rho_fit, alt, rho, atm_nz);
  gsl_spline_init(pr_fit , alt, pr , atm_nz);
}

void doInterpolation() {
  AtmLibrary *atm_  = new AtmLibrary();
  
  // Note that T_int, rho_int, etc below are computed values Above Ground Level (not MSL)
  for (int i=0; i< nz; i++) {
      alt_int[i] = dz*(i+1); // alt_int is altitude Above Ground Level (AGL)
      T_int[i]   = gsl_spline_eval(T_fit  , alt_int[i] + zmin, acc_T  ); // AGL
      rho_int[i] = gsl_spline_eval(rho_fit, alt_int[i] + zmin, acc_rho);
      pr_int[i]  = gsl_spline_eval(pr_fit , alt_int[i] + zmin, acc_pr );
      zw_int[i]  = gsl_spline_eval(zw_fit , alt_int[i] + zmin, acc_zw );
      mw_int[i]  = gsl_spline_eval(mw_fit , alt_int[i] + zmin, acc_mw );
      c_int[i]   = sqrt(GAMMA*pr_int[i]/rho_int[i]);
  }
  
  // write out interpolated values for check
  char profile_file[40] = "profile_int.dat";
  atm_->writeProfile(profile_file,nz,zmin,alt_int,zw_int,mw_int,T_int,rho_int,pr_int);
  delete atm_;
}


int checkSampling() {
  double z_cnd = (c0/freq)/10;
  double r_cnd = z_cnd;
  if ((dz-z_cnd)>-1.0e-10) {
      printf("WARNING: Altitude sampling is too low! (is %5.2f, should be <= %5.2f)\n", dz, z_cnd);
      return 1;
  } 
  if ((dr - r_cnd)>1.0e-10) {
      printf("WARNING: Range sampling is too low! (is %5.2f, should be <= %5.2f; diff=%g)\n", dr, r_cnd, dr-r_cnd);
      return 1;
  }
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////
// Atmospheric functions
 
void initAtmosphere(ProcessOptionsPE *oPE) {
  int filetype = oPE->getFiletype();
  string atmosfile;
  atmosfile = oPE->getAtmosfile();

  if (filetype==0) { 
      string atmfileord;
      atmfileord = oPE->getAtmosfileorder();
      cout << " -> Loading range-independent (1-D) ASCII profile: " << atmosfile << endl;
      load_1DAtm_ascii(atmosfile, atmfileord, oPE->getSkiplines(), oPE->getWindUnits());      
  }
  else if (filetype==1) {
      ifstream *profile = new ifstream( atmosfile.c_str(), ios_base::in );
      if (!profile->good()) { cerr << " ERROR: G2S file does not exist. Exiting..." << endl; exit(0); }
      cout << " -> Loading range-dependent (2-D) G2S profile ..." << endl;
      loadG2S_2DAtmosphere(profile);                   
      profile->close();
      delete profile;
  }
  else if (filetype==2) {
      cout << " -> Setting up built-in toy atmosphere ..." << endl;
      zmax = oPE->getMaxheight(); // set the global maxheight
      zmin = oPE->getZ_min();
      dz   = (zmax - zmin)/nz;    // set the global dz
      loadToyAtmosphere();
  }
  else if (filetype==3) { // ascii files are given in directory "profiles"
      string dir, atmfileord;
      dir = oPE->getAtm_profile_dir();
      atmfileord = oPE->getAtmosfileorder();
      load_2DAtm_ascii(Rv, dir, "profile", atmfileord, oPE->getSkiplines(), oPE->getWindUnits());
  }
} 
 
void loadG2S_1DAtmosphere(ifstream* profile) {
  AtmLibrary *atm_  = new AtmLibrary();
  atm_->getNumberOfLines( profile, &atm_nz );

  alt = new double [ atm_nz ];
  T   = new double [ atm_nz ];
  rho = new double [ atm_nz ];
  pr  = new double [ atm_nz ];
  zw  = new double [ atm_nz ];
  mw  = new double [ atm_nz ];

  atm_->readG2SAscii(profile, alt, zw, mw, T, rho, pr);
  delete atm_;
}
 

void load_1DAtm_ascii(string atmosfile, string atmosfileorder, int skiplines, string wind_units) {   
  // a little convoluted way to read an ascii profile file
  // using SampledProfile class to allow for any column order
  NCPA::SampledProfile *p;
  bool inMPS = 0;
  if ( strcmp( wind_units.c_str(), "mpersec" ) == 0) {
    inMPS = 1;
  }
  p = new SampledProfile(atmosfile, atmosfileorder.c_str(), skiplines, inMPS);

  atm_nz = p->nz();

  // read/store the 1D profile
  alt = new double [ atm_nz ];
  T   = new double [ atm_nz ];
  rho = new double [ atm_nz ];
  pr  = new double [ atm_nz ];
  zw  = new double [ atm_nz ];
  mw  = new double [ atm_nz ];
  
  p->get_z(alt, atm_nz);
  p->get_u(zw, atm_nz);
  p->get_v(mw, atm_nz);
  p->get_t(T, atm_nz);
  p->get_rho(rho, atm_nz);
  p->get_p(pr, atm_nz);
  delete p;

  // convert to SI units
  double kmps2mps = 1.0;
  //if (!wind_units.compare("kmpersec")) {
      kmps2mps = 1000.0;
  //}
  //cout << "wind units; kmps2mps=" << kmps2mps << endl;
  
  for (int i = 0; i < atm_nz; i++) {
      alt[i] = alt[i]*1000;
      rho[i] = rho[i]*1000;
      pr[i]  = pr[i]*100;
      zw[i]  = zw[i]*kmps2mps;
      mw[i]  = mw[i]*kmps2mps;
  }
  zmax = alt[atm_nz-1]; // set the global zmax
  //printf("zmax= %g\n", zmax);
}

/*
// function to return the maximum height in the provided atmospheric profile text file
// maximum height can then be compared with the value on the option --maxheight
double getMaxheightAtmProfile(ProcessOptionsPE *oPE) {   
  // a little convoluted way to read an ascii profile file
  // using SampledProfile class to allow for any column order
  NCPA::SampledProfile *p;
  string atmosfile;
  oPE->getAtmosfile(atmosfile);
  string atmfileord;
  oPE->getAtmosfileorder(atmfileord);
  
  p = new SampledProfile(atmosfile, atmfileord.c_str(), oPE->getSkiplines());

  double mx = p->z(p->nz()-1)*1000.0; // in meters MSL
  delete p;
  return mx;
}
*/


void loadG2S_2DAtmosphere(ifstream* profile) {
  AtmLibrary *atm_  = new AtmLibrary();
  atm_->getBinaryG2SDimensions( profile, &atm_nz, &atm_nr );

  alt     = new double  [ atm_nz ];
  atm_rng = new double  [ atm_nr ];
  T_2D    = new double* [ atm_nz ];
  rho_2D  = new double* [ atm_nz ];
  pr_2D   = new double* [ atm_nz ];
  zw_2D   = new double* [ atm_nz ];
  mw_2D   = new double* [ atm_nz ];
  for (int i = 0; i < atm_nz; i++) {
      T_2D   [ i ]  = new double[ atm_nr ];
      rho_2D [ i ]  = new double[ atm_nr ];
      pr_2D  [ i ]  = new double[ atm_nr ];
      zw_2D  [ i ]  = new double[ atm_nr ];
      mw_2D  [ i ]  = new double[ atm_nr ];
  }

  atm_->readG2SBinary(profile,alt,atm_rng,zw_2D,mw_2D,T_2D,rho_2D,pr_2D);

  T       = new double [ atm_nz ];
  rho     = new double [ atm_nz ];
  pr      = new double [ atm_nz ];
  zw      = new double [ atm_nz ];
  mw      = new double [ atm_nz ];
  copy1DAtmosphere(0);
  
  zmax = alt[atm_nz-1]; // set the global zmax
  //printf("readG2SBinary: zmax= %g\n", zmax);
}


// DV: 11/11/2013 added the last argument Nz0 i.e. the number of altitudes 
// in the first atm profile loaded.
// The program expects all subsequent atm profiles to have the same format
void load_NthProfile(int J, string atmosfile, string atmosfileorder, int skiplines, string wind_units, int Nz0) {
  // loads the Jth profile into the 2D global arrays 
  // use the SampledProfile object for convenience - to allow for any column order
  int i;
  NCPA::SampledProfile *p;
  p = new SampledProfile(atmosfile, atmosfileorder.c_str(), skiplines);

  atm_nz = p->nz();
  
  if (atm_nz!=Nz0) {
      delete p;
      printf("Atm profiles are expected to have the same format and size. However the file %s has %d z-grid points while the previous profile had %d z-grid points!\n", atmosfile.c_str(), atm_nz, Nz0);
      throw invalid_argument( " ");
  }

  p->get_z(  alt, atm_nz);
  p->get_u(   zw, atm_nz);
  p->get_v(   mw, atm_nz);
  p->get_t(    T, atm_nz);
  p->get_rho(rho, atm_nz);
  p->get_p(   pr, atm_nz); 
  
  // convert to SI units
  double kmps2mps = 1.0;
  //if (!wind_units.compare("kmpersec")) {
      kmps2mps = 1000.0;
  //}
  for (int i = 0; i < atm_nz; i++) {
      alt[i] = alt[i]*1000;
      rho[i] = rho[i]*1000;
      pr[i]  = pr[i]*100;
      zw[i]  = zw[i]*kmps2mps;
      mw[i]  = mw[i]*kmps2mps;
  }
  
  // populate the Nth profile in the 2D atmosphere
  for (i=0; i<atm_nz; i++) {
      T_2D  [i][J] = T  [i];
      rho_2D[i][J] = rho[i];
      pr_2D [i][J] = pr [i];
      zw_2D [i][J] = zw [i];
      mw_2D [i][J] = mw [i];
  }

  delete p;
}
 
/*
void load_NthProfile(int J, string atmosfile, string atmosfileorder, int skiplines) {
  // loads the Jth profile into the 2D global arrays 
  // use the SampledProfile object for convenience - to allow for any column order
  int i;
  NCPA::SampledProfile *p;
  p = new SampledProfile(atmosfile, atmosfileorder.c_str(), skiplines);

  atm_nz = p->nz();

  p->get_z(  alt, atm_nz);
  p->get_u(   zw, atm_nz);
  p->get_v(   mw, atm_nz);
  p->get_t(    T, atm_nz);
  p->get_rho(rho, atm_nz);
  p->get_p(   pr, atm_nz);
  
  // convert to SI units
  for (i = 0; i < atm_nz; i++) {
      alt[i] = alt[i]*1000;
      rho[i] = rho[i]*1000;
      pr[i]  = pr[i]*100;  
  } 
  
  // populate the Nth profile in the 2D atmosphere
  for (i=0; i<atm_nz; i++) {
      T_2D  [i][J] = T  [i];
      rho_2D[i][J] = rho[i];
      pr_2D [i][J] = pr [i];
      zw_2D [i][J] = zw [i];
      mw_2D [i][J] = mw [i];
  }

  delete p;
}
*/


int load_2DAtm_ascii(vector<double> Rv, string dirname, string pattern, string atmosfileorder, int skiplines, string wind_units) {
  //
  // populate 2D atmosphere from the available ascii profiles in directory dirname
  //
  int i;
  string s, atmosfile;
  list<string> files;
  list<string>::iterator it;
  
  // get and sort the files (they should have names such as profile0001.dat, etc)
  getFileList(dirname, files, pattern);
  if (files.size()==0) {
    std::ostringstream es("");
    es << "No files with pattern '" << pattern << "' were found in directory " 
       << dirname << ". Please make sure filenames with that pattern exist in the directory provided." << endl;
    throw std::invalid_argument(es.str());
  }
  
  files.sort();
  if (1) { //print the sorted file list
      cout << "Sorted file list from directory:" << dirname << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }

  it=files.begin();
  atmosfile = dirname + "/" + (*it);
  
  // if the number of files is less than Rv.size-1 then truncate Rv;
  if (files.size()<Rv.size()-1) {
      Rv.erase(Rv.begin()+files.size(), Rv.end()-1);
  }
  
  // get the number of lines in the ASCII profile; needed to properly allocate some globals
  // use the SampledProfile object for convenience - to allow for any column order
  NCPA::SampledProfile *p;
  bool inMPS = 0;
  if ( strcmp( wind_units.c_str(), "mpersec" ) == 0) {
    inMPS = 1;
  }
  p = new SampledProfile(atmosfile, atmosfileorder.c_str(), skiplines, inMPS);

  atm_nz = p->nz();
  atm_nr = (int) Rv.size()-1; // number of ranges at which to ingest a new ASCII profile

  atm_rng = new double  [ atm_nr ];
  alt     = new double  [ atm_nz ];
  T_2D    = new double* [ atm_nz ];
  rho_2D  = new double* [ atm_nz ];
  pr_2D   = new double* [ atm_nz ];
  zw_2D   = new double* [ atm_nz ];
  mw_2D   = new double* [ atm_nz ];
  
  for (i = 0; i < atm_nz; i++) {
      T_2D   [ i ]  = new double[ atm_nr ];
      rho_2D [ i ]  = new double[ atm_nr ];
      pr_2D  [ i ]  = new double[ atm_nr ];
      zw_2D  [ i ]  = new double[ atm_nr ];
      mw_2D  [ i ]  = new double[ atm_nr ];
  }
  
  //printf("atm_nz = %d\n", p->nz());
  for (i=0; i<atm_nr; i++) {
      atm_rng[i] = Rv[i];
      //printf("atm_rng[%d] = %g\n", i, atm_rng[i]);
      //printf("Rv[%d] = %g\n", i, Rv[i]);
  }
  
  T   = new double [ atm_nz ];
  rho = new double [ atm_nz ];
  pr  = new double [ atm_nz ];
  zw  = new double [ atm_nz ];
  mw  = new double [ atm_nz ];

  // populate T_2D,rho_2D, etc. by calling load_NthProfile()
  i = 0;
  for (it=files.begin(); it!=files.end(); ++it) {
      s = (*it);
      atmosfile = dirname + "/" + s;
      if (i<atm_nr) {
          load_NthProfile(i, atmosfile, atmosfileorder, skiplines, wind_units, atm_nz);
          cout << "Atm. file #" << i << ": " << atmosfile << " to be used from "
               << atm_rng[i]/1000.0 << " km" << endl;
          i++;             
      } 
      else {
          break;
      }
  }

  printf(" -> Initializing with atm. profile #0 from range 0 km\n");
  copy1DAtmosphere(0);
  
  zmax = alt[atm_nz-1]; // set the global zmax
  //printf("load_2DAtm_ascii: zmax= %g\n", zmax);
  
  return 0;
}


void copy1DAtmosphere(int index) {
  for (int i=0; i< atm_nz; i++) {
      T[i]   = T_2D[i][index];
      rho[i] = rho_2D[i][index];
      pr[i]  = pr_2D[i][index];
      zw[i]  = zw_2D[i][index];
      mw[i]  = mw_2D[i][index];
  }
}

void getG2SRangeIndex(int rr, int *atm_r_index) {
 //cout << "in getG2SRangeIndex\n";
 // this logic ensures to reload the next 1D profile as soon as we step 
 // into more than half way between the succesive profile distance
 
  double R     = rr*dr;
  double min_r = fabs(atm_rng[0] - R);

  int ii       = 0;
  for (int i=1; i<atm_nr; i++) {
      if (fabs(atm_rng[i] - R) < min_r) { 
          min_r = fabs(atm_rng[i] - R );
          ii    = i;
      }
  }
  *atm_r_index = ii;
}


void get2DAtmRangeIndex(int rr, int *atm_r_index) {
  // this logic ensures to reload the next 1D profile as soon as we marched
  // more than the profile given in atm_rng[i] 
  double R = rr*dr; 
  int ii   = 0;
  for (int i=1; i<atm_nr; i++) {
     if (R>= atm_rng[i]) {
        ii = i;
     }
  }
  *atm_r_index = ii;
}

///// END OF G2S routines /////

void loadToyAtmosphere() {
  atm_nz    = 901;
  alt = new double [ atm_nz ];
  zw  = new double [ atm_nz ];
  mw  = new double [ atm_nz ];
  T   = new double [ atm_nz ];
  rho = new double [ atm_nz ];
  pr  = new double [ atm_nz ];

  // Make temperature, density profiles, based on form from Lingevitch et al.,1999
  double T_o    = 288.2;
  double rho_o  = 1.225;
  double A[8] = { -3.9082017E-02, -1.1526465E-03,  3.2891937E-05, -2.0494958E-07, 
           //  -4.7087295E-02,  1.2506387E-03, -1.5194498E-05,  6.5818877E-08 };
               -4.7087295E-02,  1.2506387E-03, -1.5194498E-05,  6.518877E-08 };
  double B[8] = { -4.9244637E-03, -1.2984142E-06, -1.5701595E-06,  1.5535974E-08,
            //  -2.7221769E-02,  4.2474733E-04, -3.9583181E-06,  1.7295795E-08 }; 
               -2.7221769E-02,  4.247473E-04, -3.958318E-06,  1.7295795E-08 }; 
  double T_nm    = 1.0;
  double T_dnm   = 1.0;
  double rho_nm  = 0.0;
  double rho_dnm = 1.0;
  double dz      = zmax/(atm_nz-1);
  for (int i=0; i<atm_nz; i++) {
  alt[i] = i*dz;
      for (int j=0; j<8; j++) {
          if (j<4) {
              rho_nm  = rho_nm  + A[j]*pow((alt[i]/1000),j+1);
              rho_dnm = rho_dnm + B[j]*pow((alt[i]/1000),j+1);
          }
          else {
              T_nm  = T_nm  + A[j]*pow((alt[i]/1000),j-3);
              T_dnm = T_dnm + B[j]*pow((alt[i]/1000),j-3);
          }
      }
      T[i]   = T_o*(T_nm/T_dnm);
      rho[i] = rho_o*pow(10,rho_nm/rho_dnm);
      pr[i]  = rho[i]*GASCONSTANT*T[i];
      zw[i]  = 0.0;
      mw[i]  = 0.0;
      T_nm    = 1.0;
      T_dnm   = 1.0;
      rho_nm  = 0.0;
      rho_dnm = 1.0;
  }

  // Make wind profiles
  AtmLibrary *atm_  = new AtmLibrary();
  double ampJet_o       = 50.0;                  // create wind fields
  double heightJet_o    = 60.0E03;
  double widthJet_o     = 12.5E03;
  atm_->makeGaussianProfile ( atm_nz, ampJet_o,      heightJet_o, widthJet_o, alt, zw );
  //atm_->makeGaussianProfile ( atm_nz, ampJet_o,      5.0E03 , 1.25E03, alt, mw );
  //char profile_file[40] = "profile_int.dat";
  //atm_->writeProfile(profile_file,901,0,alt,zw,mw,T,rho,pr);
  char profile_file[20]= "toyatm.dat";
  atm_->writeProfile(profile_file, atm_nz,0,alt,zw,mw,T,rho,pr);  // DV 20150401
  delete atm_;
}


void getImpedance(double freq, complex<double>* Z) {
  // GroundImpedance
  // Zg = -P / ( Vz rho c )
  // Note that this is the 'nomalized' (rho c)^-1 impedance. 

  complex<double> I = complex<double>(0.,1.);

  //if ( freq >= 100. ) {
  // Z from fit to Table I Gilbert & Di from Attenborough model
  //  double Z_real = 40.0 * exp( -0.017 * freq ) + 5.0;
  //  double Z_imag = 40.0 * exp( -0.014 * freq ) + 3.0;
  //  *Z = Z_real + I*Z_imag;
  // }
  //else {
  // Waxler Model: (Z/(rho*c)) = |Z0| exp(i phi) (f0 / f)^(1/2)
  // where phi ~ pi/2 and Z0 is a known impedance at f0 > f.
  // JASA 124(5), p 2742-2754, 2008

  // Values of phi from model fits to propagation data:
  // phi = 77.8 deg = 1.3578 rad  Talmadge JASA 124(4), p1956-1962, 2008
  // phi = 75.6 deg = 1.3195 rad  Waxler   JASA 124(5), p2742-2754, 2008

  double Z0  = 26.0;    // 100 Hz Talmadge JASA 124(4), p1956-1962, 2008
  double f0  = 100.0;
  double phi = (PI/180.0)*77.8;
  *Z = Z0 * sqrt(f0/freq) * exp(I*phi);
  //}
}

////////////////////////////////////////////////////////////////////
// Program startup routines

void getCMDParameters(int nargin, char **argv) {
  if ( nargin == 6  && (! strcmp(argv[1],"toy"))) {
   azi    = atof(argv[2]);
   freq   = atof(argv[3]);
   zsrc   = atof(argv[4])*1000;
   zrcv   = atof(argv[5])*1000;
  }
  else if  ( nargin == 7 && ( (! strcmp(argv[1],"g2s_1d")) || (! strcmp(argv[1],"g2s_2d")) )) {
   azi    = atof(argv[3]);
   freq   = atof(argv[4]);
   zsrc   = atof(argv[5])*1000;
   zrcv   = atof(argv[6])*1000;
  }
  else { printUsage(); exit(0); }

  printf("\n");
  printf("High-order PE\n");
  printf(" -> Azimuth : %.2f degrees\n", azi);
  printf(" -> Frequency: %.2f Hz\n", freq);
  printf(" -> Source height: %.2f km\n", zsrc/1000);
  printf(" -> Receiver height: %.2f km\n\n", zrcv/1000);
}

/*void allocateVectors(double *Qd, double *Qo,complex<double> *psi_o,complex<double> *psi_dr,complex<double> *Bo,complex<double> *Bd,complex<double> *Cl,complex<double> *Cd, complex<double> *Cu) 
 {
    Qd     = new double [ nz ];
    Qo     = new double [ nz ];
    psi_o  = new complex<double> [ nz ];
    psi_dr = new complex<double> [ nz ];
    Bo     = new complex<double> [ nz ];
    Bd     = new complex<double> [ nz ];
    Cl     = new complex<double> [ nz ];
    Cd     = new complex<double> [ nz ];
    Cu     = new complex<double> [ nz ];
 }*/

/*
void deleteVectors(complex<double> *Qd,complex<double> *Qo,complex<double> *psi_o,complex<double> *psi_dr,complex<double> **Bo,complex<double> **Bd,complex<double> **Cl,complex<double> **Cd, complex<double> **Cu) {
  delete psi_o ; delete psi_dr;
  delete Qd; delete Qo;

  for (int i = 0; i < n_pade; i++) { 
      delete [] Bd[i]; delete [] Bo[i];
      delete [] Cl[i]; delete [] Cd[i]; delete [] Cu[i];
  }
  delete [] Bd; delete [] Bo;
  delete [] Cl; delete [] Cd; delete [] Cu;
}
*/

void deleteVectors(complex<double> *Qd,complex<double> *Qo,complex<double> *psi_o,complex<double> *psi_dr,complex<double> **Bo,complex<double> **Bd,complex<double> **Cl,complex<double> **Cd, complex<double> **Cu) {
  //delete psi_o ; delete psi_dr;
  //delete Qd; delete Qo;
  
  delete [] psi_o ; delete [] psi_dr;
  delete [] Qd; delete [] Qo;

  for (int i = 0; i < n_pade; i++) { 
      delete [] Bd[i]; delete [] Bo[i];
      delete [] Cl[i]; delete [] Cd[i]; delete [] Cu[i];
  }
  delete [] Bd; delete [] Bo;
  delete [] Cl; delete [] Cd; delete [] Cu;
}

void printUsage() {
  cout << "Usage: pape <toy|g2s_1d|g2s_2d> |<file>| <azimuth> <frequency> <zsrc> <zrcv>" << endl;
}
 
 
// a version of this function is also used in the range-dependent normal mode code
int getRegionBoundaries(bool flg, double maxrange, double req_profile_step, string prf_ranges_km, int *Nprofiles, vector<double> *R) {

// this function populates vector R with the Region boundaries.
// if in the main() the following option is used e.g.:
//  --use_profile_ranges_km 30_100_500_700.5 
// then Rv will be [0, 30, 100, 500, 700.5].
// Rv is constructed from req_profile_step: [step:step:maxrange]

  int i;
  vector<double> Rv;
  
  if (flg) { // if option --use_profile_ranges_km such as 30_100_500_700.5 is provided
      parseReqRanges(prf_ranges_km, Rv);
      
      if (Rv[0]!=0.0) { // if the first element is not zero insert 0 at the beginning
          vector<double>::iterator it;
          it = Rv.begin();
          Rv.insert(it, 0);
      }
      
      if (Rv[Rv.size()-1]>=maxrange) {//if Rv(end) is >= maxrange force the last elem. =maxrange
          Rv[Rv.size()-1] = maxrange;
      }
      else { //add maxrange at end
          Rv.push_back(maxrange);
      }
      (*Nprofiles) = Rv.size()-1;
  }
  else { // if profiles are used at equal intervals
      (*Nprofiles) = (int) ceil(maxrange/req_profile_step); // at least one profile is ensured
      vector<double> Rv1 ((*Nprofiles)+1, 0.0);
      for (i=1; i<(int) Rv1.size(); i++) {
          Rv1[i] = i*req_profile_step;  
      } 
      Rv1[Rv1.size()-1] = maxrange; // last element adjusted to maxrange
      Rv.swap(Rv1);
  }
  //cout << "Nprofiles = " << Nprofiles << endl;
  //cout << "Rv size   = " << Rv.size() << endl;
  //for (i=0; i<Rv.size(); i++) {
  //  printf("Rv[%i] = %g\n", i, Rv[i]);
  //}
  
  (*R).swap(Rv);
  return 0;     
}


void parseReqRanges(std::string str, std::vector<double>& retVal) {
  // parse a string of the form "2_20_34_566_34.35" to extract the numbers into an array

  size_t  ix, it, i, p1; 

  // first pass to find out the number of underscores
  ix = -1;
  it =  0;
  ix = str.find("_");
  //cout << "ix= " << ix << endl;
  if (ix==string::npos) { // no underscore found
      vector<double> x (1,-1);
      retVal.swap(x);
      return; //NULL;
  }
  else { //  if (ix!=string::npos) i.e. underscore found
  it = 1;
      while (ix!=string::npos) {
          ix = str.find("_", ix+1);
          it++;
          //cout << "ix= " << ix << endl;
          //cout << "it= " << it << endl;
      }
  }
  //cout << "it= " << it << endl;

  // second pass to populate x
  vector<double> x (it,100);
  ix = str.find("_");
  //cout << "ix= " << ix << endl;
  if (ix==0) {
      throw invalid_argument("profile_ranges should not start with an underscore.");
      exit(1);
  }
      
  if (ix!=string::npos) {     
    //printf("%s\n", str.substr(0,ix).c_str());
     
      x[0] = atof(str.substr(0,ix).c_str())*1000.0;
      //printf("x[0]=%g\n", x[0]); 

      for (i=1; i<it; i++) {
          p1 = ix+1;
          ix = str.find("_", p1);
          x[i] = atof(str.substr(p1,ix-p1).c_str())*1000.0; // in meters
          //printf("x[%d]=%g\n", i, x[i]);  
      }
  }
  retVal.swap(x);
  return;
}


/*
int plotwGNUplot(double freq, bool write_2D_TLoss) {
    if (1) {   
        // open a pipe to gnuplot and plot the results
        FILE *pipe = popen("gnuplot -persist","w");
        printf("Plotting with gnuplot...\n");             
        fprintf(pipe, "set title 'High-angle PE - Transmission loss; Frequency %g Hz'\n", freq);
        fprintf(pipe, "set xlabel 'Range [km]'\n");
        fprintf(pipe, "set ylabel 'TL [dB]'\n");
        fprintf(pipe, "set grid; show grid;\n");
        fprintf(pipe, "plot './tloss_1d.pe' using 1:(10*log10($2**2 + $3**2)) with lines lt 3 title ''\n"); 
        pclose(pipe);
        
        if (write_2D_TLoss) {
            //To plot a surface plot:
            pipe = popen("gnuplot -persist","w");
            //fprintf(pipe,"set term wxt 1\n");
            fprintf(pipe,"set pm3d map\n");
            fprintf(pipe,"set cbrange [-140:-100]\n");
            fprintf(pipe,"set xlabel 'Range [km]'\n");
            fprintf(pipe,"set ylabel 'Height [km]'\n");
            fprintf(pipe, "set title 'High-angle PE - Transmission loss; Frequency %g Hz'\n", freq);
            fprintf(pipe,"splot './tloss_2d.pe' using 1:2:(20*log10(sqrt($3**2 + $4**2)))\n");
            pclose(pipe);
        }
        //printf("Plotting done.\n");
    }
  return 0;
}
*/



int get_RngDepndProfiles_ascii( int N, string atmosfileorder, uint skiplines, \
   string dirname, string pattern, string windspeed_units, string &atmosfile) {
              
  // get the Nth profile stored in a file named "profile000N" residing in a specified
  // directory. If N is greater than the number of files found then the last file
  // is used repeatedly as necessary             

  list<string> files;
  list<string>::iterator it;
  
  // get and sort the files (they should have names such as profile0001.dat, etc)
  getFileList(dirname, files, pattern);
  files.sort();
  if (1) { //print the sorted file list
      cout << "sorted list contains:" << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }
  
  int j = 1;
  string s;
  if (N<(int)(files.size()+1)) {
      for (it=files.begin(); it!=files.end(); ++it) {
          if (j>=N) {
              cout << "j= " << j << ";  *it= " << *it << endl;
              s = (*it);
              break;
          }
          j++;
      }
  }
  else {
      s = files.back();
  } 

  atmosfile = dirname + "/" + s;
  cout << "Atmosfile name " << atmosfile << endl;
  return 0;
}


int getFileList(std::string dir, std::list<string> &files, std::string pattern) {
  int pos = -1;
  string a;
  DIR *dp;
  struct dirent *dirp;
  if((dp = opendir(dir.c_str())) == NULL) {
      std::ostringstream es;
      es << "Error opening directory:" << dir;
      throw invalid_argument(es.str());
  }

  while ((dirp = readdir(dp)) != NULL) {
      a   = string(dirp->d_name);
      pos = a.find(pattern);
      //cout << "pos = " << pos <<  endl;
      if (pos>=0) {
          //cout << a << " pos=" << pos << endl;
          files.push_back(string(dirp->d_name));
      }
      else {
          //cout << "this file does not fit pattern: " << a << endl;
      }
  }
  closedir(dp);
  return 0;
}


//
// Function to parse the input options (both command lines and in the options file ModessRD.options)
//
AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();

  opt->addUsage( "----------------------------------------------------------------------------" );	
  opt->addUsage( "|                                                                           |" );
  opt->addUsage( "|           High-Angle High-Mach Parabolic Equation  (PaPE)                 |" );
  opt->addUsage( "|                                                                           |" );  
  opt->addUsage( "----------------------------------------------------------------------------" );		
  opt->addUsage( "Usage: " );
  opt->addUsage( "By default the program computes the 1D transmission loss (TL)" );
  opt->addUsage( "at the ground or the specified receiver height and saves the data to:" );
  opt->addUsage( "   file tloss_1d.pe - considering attenuation in the atmosphere" );
	opt->addUsage( "Additionally, if the flag --write_2D_TLoss is given the 2D TL is saved to file tloss_2d.pe" );
	opt->addUsage( "" );
  opt->addUsage( "The options below can be specified in a colon-separated file \"PaPE.options\" or at the command line.  Command-line options override file options." );
  opt->addUsage( "Be sure to precede the options with two minuses: --" );
  opt->addUsage( "" );
  opt->addUsage( " --help -h                Print this message and exit" );
  opt->addUsage( "" );
  opt->addUsage(  " The atmosphere can be specified from one of 4 different sources:");
  opt->addUsage( "    1. An .env file containing the atmospheric specifications at certain ranges:" );
  opt->addUsage( "       use option --g2senvfile <filename>" );
  opt->addUsage( "    2. Several ASCII files stored in a given directory:" );
  opt->addUsage( "       use option --use_1D_profiles_from_dir <mydirname>" );
  opt->addUsage( "    3. A single ASCII file. This will just force a range-independent run." );
  opt->addUsage( "       use option --atmosfile1d <filename>" );
  opt->addUsage( "    4. A built-in NCPA canonical profile." );
  opt->addUsage( "       use option (flag) --ncpatoy" );    
  opt->addUsage( "" );  
  opt->addUsage( "The available options are:" );	
  opt->addUsage( "" );	
  opt->addUsage( "REQUIRED (no default values):" );
  opt->addUsage( " --atmosfileorder         The order of the (z,u,v,w,t,d,p) fields in the file");
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
  opt->addUsage( "    Example: >> ../bin/pape --atmosfileorder zuvwtdp --skiplines 1" );
  opt->addUsage( "                --azimuth 90 --freq 0.1 --use_1D_profiles_from_dir myprofiles_dir" );
  opt->addUsage( "                --use_profile_ranges_km 100_300_500_600_700" );
  opt->addUsage( "" );
  opt->addUsage( " --atmosfile1d  <filename>  Uses an ASCII 1D atmosphere file." ); 
  opt->addUsage( "                          In this case the run will just be range-independent." );   	
  opt->addUsage( "" );	
  opt->addUsage( "OPTIONAL [defaults]:" );
  opt->addUsage( " --maxheight_km           Calculation grid height in km above MSL [150 km]" );
  opt->addUsage( " --zground_km             Height of the ground level above MSL [0 km]" );  
  opt->addUsage( " --Nz_grid                Number of points on the z-grid from ground to maxheight [20000]" );  
  opt->addUsage( " --sourceheight_km        Source height in km Above Ground Level (AGL) [0]" );
  opt->addUsage( " --receiverheight_km      Receiver height in km AGL [0]" );
  opt->addUsage( " --maxrange_km            Maximum horiz. propagation distance from origin [1000 km]" );
  opt->addUsage( " --rng_step               A usually fractional number specifying the range step" );
  opt->addUsage( "                          in wavelengths: e.g. --rng_step 0.1 means a step" );
  opt->addUsage( "                          of 0.1*wavelength [default is 0.1]." );
  opt->addUsage( " --ground_impedance_model Name of the ground impedance models to be employed:" );
  opt->addUsage( "                          [rigid], others TBD" );
	opt->addUsage( " --wind_units             Specify 'kmpersec' if the winds are given" );
	opt->addUsage( "                          in km/s [ mpersec ]" );  
  opt->addUsage( " --n_pade                 Number of Pade coefficients [4]" );  
  opt->addUsage( "" );
  opt->addUsage( " --starter_type           Specifies one of 3 available PE starter" );
  opt->addUsage( "                          fields: gaussian, greene, modal." );
  opt->addUsage( "                          The default is 'gaussian'." );
  opt->addUsage( "                          'modal' requires a precomputed starter field" );
  opt->addUsage( "                          obtained by running Modess with option" ); 
  opt->addUsage( "                          --modal_starter_file." );
  opt->addUsage( "" );	
  opt->addUsage( "" );	
  opt->addUsage( "" );	
  opt->addUsage( "" );	  
  opt->addUsage( " --use_profile_ranges_km" );
  opt->addUsage( "                          e.g. --use_profile_ranges_km  20_50_80.5_300     " );   
  opt->addUsage( "                          The profiles at certain ranges specified by numbers" );
  opt->addUsage( "                          (in km) in a string such as 20_50_80.5_300 are");
  opt->addUsage( "                          requested in the propagation. Note that underscores" );
  opt->addUsage( "                          are necessary to separate the numbers." );
  opt->addUsage( "                          In this example the ranges at which the profiles" );  
  opt->addUsage( "                          are considered are: 0, 20, 50, 80.5, 300 km i.e." );
  opt->addUsage( "                          0 is always the first distance even if not specified." );  
  opt->addUsage( "                          Note also that these are requested ranges;" );
  opt->addUsage( "                          however the left-closest profile available" );
  opt->addUsage( "                          in the .env file will actually be used; " );
  opt->addUsage( "                          for instance we request the profile at 300 km " );
  opt->addUsage( "                          but in the .env file the left-closest profile" );
  opt->addUsage( "                          may be available at 290 km and it is the one used." ); 
  opt->addUsage( "                          This convention may change in the future." );
  opt->addUsage( "                          This option is used in conjunction with" );
  opt->addUsage( "                              --use_1D_profiles_from_dir" );      
  opt->addUsage( "" ); 
  opt->addUsage( " --use_profiles_at_steps_km" );  
  opt->addUsage( "                          e.g. --use_profiles_at_steps_km 100" );
  opt->addUsage( "                          The profiles are requested at equidistant intervals " );
  opt->addUsage( "                          specified by this option [1000]" );
  opt->addUsage( "                          This option is used in conjunction with" );
  opt->addUsage( "                              --use_1D_profiles_from_dir" );
  opt->addUsage( "" );	  
	opt->addUsage( " --use_attn_file          Use it to specify a file name containing user-provided" );
	opt->addUsage( "                          attenuation coefficients to be loaded instead of " );
	opt->addUsage( "                          the default Sutherland-Bass attenuation. " ); 
	opt->addUsage( "                          The text file should contain two columns: " );
	opt->addUsage( "                              height (km AGL) and " );
	opt->addUsage( "                              attenuation coefficients in np/m." );    
  opt->addUsage( "" );                        
  opt->addUsage( "" );  
  opt->addUsage( "FLAGS (no value required):" );
  opt->addUsage( " --ncpatoy                Use built-in NCPA canonical profile" );
  opt->addUsage( " --write_2D_TLoss         Outputs the 2D transmission loss to" );
  opt->addUsage( "                          default file: tloss_2d.pe" );	
  opt->addUsage( " --do_lossless            Computation is done with no atm. absorption" );  
  opt->addUsage( "" );
  opt->addUsage( "" );
  opt->addUsage( " The column order of the output files is as follows (P is complex pressure):" );
  opt->addUsage( "  tloss_1d.pe:           r, 4*PI*Re(P), 4*PI*Im(P)" );
  opt->addUsage( "  tloss_2d.pe:        r, z, 4*PI*Re(P), 4*PI*Im(P)" );   
  opt->addUsage( "" );
  opt->addUsage( "" );
  opt->addUsage( "--------------------------------------------------------------------" );  
  opt->addUsage( "Examples to run with various options (from the 'samples' directory):" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/pape  --ncpatoy --azimuth 90 --freq 0.1 --write_2D_TLoss" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/pape --g2senvfile g2sgcp2011012606L.jordan.env --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.3 --sourceheight_km 0 --receiverheight_km 0 --maxheight_km 180 --starter_type gaussian --n_pade 6 --maxrange_km 500" ); 
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/pape --atmosfile1d NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1 --sourceheight_km 0 --receiverheight_km 0 --maxheight_km 180 --starter_type gaussian --n_pade 4 --maxrange_km 500 --write_2D_TLoss" );  
  opt->addUsage( "" );  
  opt->addUsage( "    ../bin/pape --use_1D_profiles_from_dir ../samples/profiles --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --sourceheight_km 0 --receiverheight_km 0 --maxheight_km 180 --starter_type gaussian --n_pade 6 --maxrange_km 1000  --use_profiles_at_steps_km 20" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/pape --use_1D_profiles_from_dir ../samples/profiles --atmosfileorder zuvwtdp --skiplines 1 --azimuth 90 --freq 0.1 --sourceheight_km 0 --receiverheight_km 0 --maxheight_km 180 --starter_type gaussian --n_pade 6 --maxrange_km 1000  --use_profile_ranges_km 0_20_60_400" );  
  opt->addUsage( "" );
  opt->addUsage( "--------------------------------------------------------------------" );  
  opt->addUsage( "Example to plot 1D TL with gnuplot:" );
  opt->addUsage( " plot './tloss_1d.pe' using 1:(10*log10($2**2 + $3**2))" );  
  opt->addUsage( "" );
  opt->addUsage( "Example to plot 2D TL with gnuplot:" );
  opt->addUsage( " set cbrange [-200:-100]" );
  opt->addUsage( " splot './tloss_2d.pe' using 1:2:(20*log10(sqrt($3**2 + $4**2)))" );  
  opt->addUsage( "" );   
  opt->addUsage( "--------------------------------------------------------------------" );  
  

  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  opt->setFlag( "write_2D_TLoss" );
  opt->setFlag( "ncpatoy" );
  opt->setFlag( "do_lossless" );
  opt->setFlag( "plot" );

  opt->setOption( "atmosfile1d" );
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
  opt->setOption( "rng_step" );
  opt->setOption( "maxheight_km" );
  opt->setOption( "zground_km" );  
  opt->setOption( "stepsize" );
  opt->setOption( "Nz_grid" );
  opt->setOption( "n_pade" );
  opt->setOption( "starter_type" );
  opt->setOption( "modal_starter_file" );
  opt->setOption( "ground_impedance_model" );
  opt->setOption( "Lamb_wave_BC" );
  opt->setOption( "use_profile_ranges_km" ); 
  opt->setOption( "use_profiles_at_steps_km" );
  opt->setOption( "use_attn_file" );

  // Process the command-line arguments
  opt->processFile( "./PaPE.options" );
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
