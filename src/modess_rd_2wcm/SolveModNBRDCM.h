#ifndef _SOLVEMODNBRDCM_H_
#define _SOLVEMODNBRDCM_H_

#include "ProcessOptionsNBRDCM.h"

namespace NCPA {
  class SolveModNBRDCM {
    public:
    
      SolveModNBRDCM(double freq, double azi, string atmosfile, string wind_units, \
            NCPA::SampledProfile *atm_profile, \
            int Nz_grid, double z_min, double maxheight, \
            int Nrng_steps, double maxrange, \
            double sourceheight, double receiverheight, \
            std::string gnd_imp_model, \
            int Lamb_wave_BC, double tol,	int write_2D_TLoss, \
            bool write_phase_speeds, bool write_dispersion, bool write_modes);
      
      SolveModNBRDCM(ProcessOptionsNB *oNB, NCPA::SampledProfile *atm_profile); // constructor 2 
            
      ~SolveModNBRDCM(); //destructor            

      void setParams(double freq, double azi, string atmosfile, string wind_units, \
            NCPA::SampledProfile *atm_profile, \
            int Nz_grid, double z_min, double maxheight, \
            int Nrng_steps, double maxrange, \
            double sourceheight, double receiverheight, \
            std::string gnd_imp_model, \
            int Lamb_wave_BC, double tol, int write_2D_TLoss, \
            bool write_phase_speeds, bool write_dispersion, bool write_modes);		
		
		  void setParams(ProcessOptionsNB *oNB, NCPA::SampledProfile *atm_prof);
		  
      void printParams();

      int computeModes(double kMin, double kMax);	

      //int getAbsorption(int n, double dz, NCPA::SampledProfile *p, double freq, double *alpha);
      int getAbsorption(int n, double dz, NCPA::SampledProfile *p, double freq, string usrattfile, double *alpha);

      int getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, double dz, NCPA::SampledProfile *p, double admittance, double freq, double azi, double *diag, double *k_min, double *k_max, bool turnoff_WKB);		

      int getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev);

      int sturmCount(int n, double dz, double *diag, double k, int *cnt);	

      int doPerturb(int nz, double z_min, double dz, int n_modes, double freq, NCPA::SampledProfile *p, double *k, double **v, double *alpha, std::complex<double> *k_pert);

      int doSelect(int nz, int n_modes, double k_min, double k_max, double *k2, double **v, double *k_s, double **v_s, int *select_modes);

      int getTLoss1D(int select_modes, double dz, int n_r, double dr, double z_src, double z_rcv, complex<double> *k_pert, double **v_s);

      int getTLoss2D(int nz, int select_modes, double dz, int n_r, double dr, double z_src, complex<double> *k_pert, double **v_s);

      int writeDispersion(int select_modes, double dz, double z_src, double z_rcv, double freq, complex<double> *k_pert, double **v_s);

      int writePhaseSpeeds(int select_modes, double freq, std::complex<double> *k_pert);

      int writeEigenFunctions(int nz, int select_modes, double dz, double **v_s);
      
      int getNumberOfModes();
      
      int getMaxNumberOfModes();
      
      int getOptimalNumberOfModes();
      
      complex<double> *getWavenumbers();
      
      double **getWavevectors();
      
      int writeEigenValVecs(string fn, int n_modes);

    private:
      bool   write_2D_TLoss;
      bool   write_phase_speeds;
      bool   write_dispersion;
      bool   write_modes;
      bool   turnoff_WKB;

      int    Nz_grid;
      int    Nrng_steps;
      int    Lamb_wave_BC;
      int    select_modes;
      int    Nmod_max;   // max number of modes for a given k_min, k_max interval
      int    Nmod_optim; // optimal number of modes; automatically determined for an atm profile
      int    skiplines;
           
      double freq;
      double azi;           
      double z_min;
      double maxheight;     
      double maxrange;
      double sourceheight;
      double receiverheight;
      double tol;
      double **v_s;
      double *Hgt, *zw, *mw, *T, *rho, *Pr;

      complex<double> *k_pert;
      
      NCPA::SampledProfile *atm_profile;
      std::string gnd_imp_model;
      std::string atmosfile;
      std::string wind_units;
      std::string usrattfile;
   
  }; // mandatory semicolon here
}

#endif
