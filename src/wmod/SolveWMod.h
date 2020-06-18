#ifndef _SolveWMod_H_
#define _SolveWMod_H_

#include <complex>
#include "Atmosphere1D.h"
#include "ProcessOptionsNB.h"
#include "parameterset.h"


namespace NCPA {
  class SolveWMod {
    public:
	    
      SolveWMod(ParameterSet *param, Atmosphere1D *atm_profile);
      ~SolveWMod();
      
      void setParams(ParameterSet *param, Atmosphere1D *atm_profile);	
			
      void printParams();
	
      int computeModes();	

      //int getAbsorption(int n, double dz, NCPA::SampledProfile *p, double freq, double *alpha);	
      
      //int getAbsorption(int n, double dz, NCPA::SampledProfile *p, double freq, string usrattfile, double *alpha);      

      int getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, double dz, Atmosphere1D *p, 
            double admittance, double freq, double *diag, double *kd, double *md, double *cd, double *k_min, 
            double *k_max, bool turnoff_WKB);

      int getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev);

      int sturmCount(int n, double dz, double *diag, double k, int *cnt);	

      int doPerturb(int nz, double z_min, double dz, int n_modes, double freq, double *k, double **v, 
            double *alpha, std::complex<double> *k_pert);

      int doSelect(int nz, int n_modes, double k_min, double k_max, double *k2, double **v, double *k_s, double **v_s, int *select_modes);

      int getTLoss1D(int select_modes, double dz, int n_r, double dr, double z_src, double z_rcv, double *rho, complex<double> *k_pert, double **v_s);
      
      int getTLoss1DNx2(double azimuth, int select_modes, double dz, int n_r, double dr, double z_src, double z_rcv, double *rho, complex<double> *k_pert, double **v_s, bool Nx2, int iter);

      int getTLoss2D(int nz, int select_modes, double dz, int n_r, double dr, double z_src, double *rho,  complex<double> *k_pert, double **v_s);

      //int writeDispersion(int select_modes, double dz, double z_src, double z_rcv, double freq, complex<double> *k_pert, double **v_s);
      
      int writeDispersion(int select_modes, double dz, double z_src, double z_rcv, double freq, complex<double> *k_pert, double **v_s, double *rho);

      int writePhaseSpeeds(int select_modes, double freq, std::complex<double> *k_pert);

      int writeEigenFunctions(int nz, int select_modes, double dz, double **v_s);

    private:
            bool   write_2D_TLoss;
            bool   write_phase_speeds;
            bool   write_dispersion;
            bool   write_modes;        
            bool   Nby2Dprop;
            bool   turnoff_WKB;
            bool   wvnum_filter_flg;
            bool   z_min_specified;
          
            int    Nz_grid;
            int    Nrng_steps;
            int    Lamb_wave_BC; 
            int    Naz;
            int    skiplines;
      
            double freq;
            double azi;
            double azi_min;
            double azi_max;
            double azi_step;           
            double z_min;
            double maxheight;     
            double maxrange;
            double sourceheight;
            double receiverheight;
            double tol;
            double *Hgt, *zw, *mw, *T, *rho, *Pr, *c_eff, *alpha;
            double c_min; // for wavenumber filtering option
            double c_max; // for wavenumber filtering option
      
            NCPA::Atmosphere1D *atm_profile;
            std::string gnd_imp_model;
            std::string atmosfile;
            std::string wind_units;
            std::string usrattfile;
      
	}; 
}

#endif
