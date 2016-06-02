#ifndef _SOLVEMODNB_H_
#define _SOLVEMODNB_H_
#include "ProcessOptionsNBRD.h"

namespace NCPA {
  class SolveModNB {
    public:
    
      SolveModNB(double freq, double azi, string atmosfile, string wind_units, 
            NCPA::SampledProfile *atm_profile, \
            int Nz_grid, double z_min, double maxheight, \
            int Nrng_steps, double maxrange, \
            double sourceheight, double receiverheight, \
            std::string gnd_imp_model, \
            int Lamb_wave_BC, double tol,	int write_2D_TLoss, \
            bool write_phase_speeds, bool write_dispersion, bool write_modes, \
            bool turnoff_WKB);
            
      SolveModNB(ProcessOptionsNB *oNB, NCPA::SampledProfile *atm_profile); // constructor 2
            
      ~SolveModNB(); //destructor            

      void setParams(double freq, double azi, string atmosfile, string wind_units, \
            NCPA::SampledProfile *atm_profile, \
            int Nz_grid, double z_min, double maxheight, \
            int Nrng_steps, double maxrange, \
            double sourceheight, double receiverheight, \
            std::string gnd_imp_model, \
            int Lamb_wave_BC, double tol, int write_2D_TLoss, \
            bool write_phase_speeds, bool write_dispersion, bool write_modes, \
            bool turnoff_WKB);     
            
      void setParams(ProcessOptionsNB *oNB, NCPA::SampledProfile *atm_prof);	
		
      void printParams();

      int computeModes();	

      //int getAbsorption(int n, double dz, NCPA::SampledProfile *p, double freq, double *alpha);
      
      int getAbsorption(int n, double dz, NCPA::SampledProfile *p, double freq, string usrattfile, double *alpha);

      int getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, double dz, NCPA::SampledProfile *p, double admittance, double freq, double azi, double *diag, double *k_min, double *k_max, bool turnoff_WKB);		

      int getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev);

      int sturmCount(int n, double dz, double *diag, double k, int *cnt);	

      int doPerturb(int nz, double z_min, double dz, int n_modes, double freq, NCPA::SampledProfile *p, double *k, double **v, double *alpha, std::complex<double> *k_pert);

      int doSelect(int nz, int n_modes, double k_min, double k_max, double *k2, double **v, double *k_s, double **v_s, int *select_modes);
      
      int getNumberOfModes();
      
      complex<double> *getWavenumbers();
      
      double **getWavevectors();

    private:
      bool   write_2D_TLoss;
      bool   write_phase_speeds;
      bool   write_dispersion;
      bool   write_modes;
      bool   turnoff_WKB;
          
      int    select_modes; 
      int    Nz_grid;
      int    Nrng_steps;
      int    Lamb_wave_BC;
      int    skiplines;
        
      double freq;
      double azi;           
      double z_min;
      double maxheight;     
      double maxrange;
      double sourceheight;
      double receiverheight;
      double tol;	
  
      double *Hgt, *zw, *mw, *T, *rho, *Pr;
      double **v_s;
      complex<double> *k_pert;
      
      NCPA::SampledProfile *atm_profile;
      std::string gnd_imp_model;
      std::string atmosfile;
      std::string wind_units;
      std::string usrattfile;
    
  }; // mandatory semicolon here
}

#endif
