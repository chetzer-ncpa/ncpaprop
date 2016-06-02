#ifndef _SOLVEMODBB_H_
#define _SOLVEMODBB_H_

#include <stdexcept>
#include <ctime>
#include "Atmosphere.h"
#include "anyoption.h"
#include "ModBB_lib.h"


namespace NCPA {
  class SolveModBB {
    public:
      SolveModBB(\
	        string filename, string atmosfile, string wind_units, \
	        NCPA::SampledProfile *atm_profile, \
	        int Nfreq, double f_min, double f_step, double f_max, int Nz_grid, double azi, \
	        double z_min, double maxheight, double sourceheight, double receiverheight,\
	        string gnd_imp_model, int Lamb_wave_BC, \
	        bool out_dispersion, bool out_disp_src2rcv, bool usemodess_flg);

      ~SolveModBB(); //destructor


      void setParams( \
	        string filename1, string atmosfile1, string wind_units1, \
	        NCPA::SampledProfile *atm_profile1, \
	        int Nfreq1, double f_min1, double f_step1, double f_max1, int Nz_grid1, double azi1, \
	        double z_min1, double maxheight1, double sourceheight1, double receiverheight1, \
	        string gnd_imp_model1, int Lamb_wave_BC1, \
	        bool out_dispersion1, bool out_disp_src2rcv1, bool usemodess_flg1);						
		
      void printParams();

      int computeModes(bool usemodess_flg);
      
      int computeModESS();
      int computeWmodes();

      int getAbsorption(int n, double dz, NCPA::SampledProfile *atm_profile, double freq, double *alpha);
								        
      int getModalTraceModESS(int nz, double z_min, \
                    double sourceheight, double receiverheight, \
                    double dz, NCPA::SampledProfile *atm_profile, \
								    double admittance, double freq, double azi, double *diag, \
								    double *k_min, double *k_max, bool turnoff_WKB);
								    
								    
      int getModalTraceWMod(int nz, double z_min, \
                    double sourceheight, double receiverheight, \
                    double dz, NCPA::SampledProfile *atm_profile, \
                    double admittance, double freq, double *diag, \
                    double *kd, double *md, double *cd, \
                    double *k_min, double *k_max, bool turnoff_WKB);								    								        			

      int getNumberOfModes(int n, double dz, double *diag, double k_min, double k_max, int *nev);

      int sturmCount(int n, double dz, double *diag, double k, int *cnt);	

      int doPerturb2( int nz, double z_min, double dz, int n_modes, double freq, \
							        NCPA::SampledProfile *atm_profile, double *k, double **v, \
							        double *alpha, complex<double> *k_pert, double *kr, double *ki);

      int doSelectModESS( int nz, int n_modes, double k_min, double k_max, \
                          double *k2, double **v, double *k_s, double **v_s, \
                          int *select_modes);
                    
      int doSelectWMod( int nz, int n_modes, double k_min, double k_max, \
                        double *kH, double **v, double *k_s, double **v_s, \
                        int *select_modes);                    
			
			// old prototype        
      int writeDispersion_bb_ascii(\
				        std::string filen, int select_modes, double dz, double z_src, double z_rcv, \
				        double freq, double *rho, std::complex<double> *k_pert, double **v_s);
			
			// newer prototype to which we pass a file pointer	(DV 20150720)        
      int writeDispersion_bb_ascii(\
				        FILE *fp, int select_modes, double dz, double z_src, double z_rcv, \
				        double freq, double *rho, std::complex<double> *k_pert, double **v_s);				        
		        
				        
      int writeDispersion_bb_bin3(\
					      string filen, double freq, int Nfreq, double df, int Nz_grid,  \
					      double z_min, int Nz_subgrid, double delZ, int NN, int Nmodes, \
					      double dz, double z_src, \
					      double *rho, double *kreal, double *kim, double **v_s);				        			        	

    private:
      bool   out_dispersion;
      bool   out_disp_src2rcv;
      bool   usemodess_flg;
      
      int    Nfreq;
      int    Nz_grid;  
      int    Nrng_steps;
      int    Lamb_wave_BC;
      			
      double f_min;
      double f_step;
      double f_max;
      double azi;
      double z_min;
      double dz;
      double sourceheight;
      double receiverheight;
      double maxheight;
      double maxrange;
      double tol;
      double *Hgt, *zw, *mw, *T, *rho, *Pr;     

      NCPA::SampledProfile *atm_profile;
      std::string atmosfile; 
      std::string filename;
      std::string gnd_imp_model;
      std::string wind_units;
      std::string subdir;   // directory name using date_time format YYYY_MM_DD_HH_MM_SS
      std::string disp_fn;  // a buffer to store path/filename
		
  }; // mandatory semicolon here
}

#endif
