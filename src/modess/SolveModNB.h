#ifndef _SOLVEMODNB_H_
#define _SOLVEMODNB_H_
#include "ProcessOptionsNB.h"

namespace NCPA {
	class SolveModNB {
	public:           

        
		SolveModNB(ProcessOptionsNB *oNB, NCPA::SampledProfile *atm_profile); // constructor 2
      
            
		void setParams(ProcessOptionsNB *oNB, NCPA::SampledProfile *atm_prof);                      	
		
		void printParams();

		int computeModes();	

		int getAbsorption(int n, double dz, NCPA::SampledProfile *p, double freq, 
			string usrattfile, double *alpha);

		int getModalTrace(int nz, double z_min, double sourceheight, double receiverheight, 
			double dz, NCPA::SampledProfile *p, double admittance, double freq, double azi, 
			double *diag, double *k_min, double *k_max, bool turnoff_WKB, double *c_eff);      

		int getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev);

		int sturmCount(int n, double dz, double *diag, double k, int *cnt);	

		int doPerturb(int nz, double z_min, double dz, int n_modes, double freq, 
			NCPA::SampledProfile *p, double *k, double **v, double *alpha, 
			std::complex<double> *k_pert);

		int doSelect(int nz, int n_modes, double k_min, double k_max, double *k2, double **v, 
			double *k_s, double **v_s, int *select_modes);

		int getTLoss1D(int select_modes, double dz, int n_r, double dr, double z_src, double z_rcv, 
			double *rho, complex<double> *k_pert, double **v_s);     
      

      
		int getTLoss1DNx2(double azimuth, int select_modes, double dz, int n_r, double dr, double z_src, 
			double z_rcv,  double *rho, complex<double> *k_pert, double **v_s, bool Nx2, int iter);      

		int getTLoss2D(int nz, int select_modes, double dz, int n_r, double dr, double z_src, 
			double *rho, complex<double> *k_pert, double **v_s); 
      
		// Modal starter - DV 20151014
		// Modification to apply the sqrt(k0) factor to agree with Jelle's getModalStarter; 
		// this in turn will make 'pape' agree with modess output
		void getModalStarter(int nz, int select_modes, double dz, double freq,  double z_src, 
			double z_rcv, double *rho, complex<double> *k_pert, double **v_s, string modstartfile);

		int writeDispersion(int select_modes, double dz, double z_src, double z_rcv, double freq, 
			complex<double> *k_pert, double **v_s, double *rho);

		int writePhaseSpeeds(int select_modes, double freq, std::complex<double> *k_pert);

		int writeEigenFunctions(int nz, int select_modes, double dz, double **v_s);
      
		int writePhaseAndGroupSpeeds(int nz, double dz, int select_modes, double freq, 
			complex<double> *k_pert, double **v_s, double *c_eff);   

	private:
		bool   write_2D_TLoss;
		bool   write_phase_speeds;
		bool   write_speeds;
		bool   write_dispersion;
		bool   write_modes;        
		bool   Nby2Dprop;
		bool   turnoff_WKB;
		bool   wvnum_filter_flg;
          
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
		double *Hgt, *zw, *mw, *T, *rho, *Pr, *c_eff;
		double c_min; // for wavenumber filtering option
		double c_max; // for wavenumber filtering option
      
		NCPA::SampledProfile *atm_profile;
		std::string gnd_imp_model;
		std::string atmosfile;
		std::string wind_units;
		std::string usrattfile;
		std::string modstartfile; // store the modal starter in this file
		//std::string atmosfileorder;
      
	}; 
}

#endif
