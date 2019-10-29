#include "Atmosphere.h"
#include <complex>
#include <list>

/*
double **dmatrix(long nr, long nc);
int free_dmatrix(double**v, long nr, long nc);
complex<double> **cmatrix(long nr, long nc);
int free_cmatrix(complex<double>**v, long nr, long nc);
complex<double> ***c3Darray(size_t xlen, size_t ylen, size_t zlen);
void free_c3Darray(complex<double> ***data, size_t xlen, size_t ylen);
*/

int writeDispersion(int select_modes, double dz, double z_src, double freq, std::complex<double> *k_pert, double **v_s);

int writeDispersion_bb(string filen, int select_modes, double dz, double z_src, double freq, complex<double> *k_pert, double **v_s);

int writeDispersion_bb_ascii(string filen, int Nz, double delZ, int NN, int select_modes, double dz, double z_src, double freq, complex<double> *k_pert, double **v_s);

int readDispersion_cbb_ascii(string filen, double *f_vec, int *mode_count, \
								double *prho_zsrc, double *prho_zrcv, complex<double> **re_k, \
								complex<double> **mode_S, complex<double> **mode_R);
					

int count_rows_arbcol(const char *filename);								

int writeDispersion_bb_bin(string filen, int Nz, double delZ, int NN, int select_modes, double dz, double z_src, double freq, complex<double> *k_pert, double **v_s);

int writeDispersion_bb_bin2(string filen, int Nfreq, double df, int Nz_grid, double z_min, \
 														int Nz_subgrid, double delZ, int NN, int select_modes, double dz, \
 														double z_src, double freq, complex<double> *k_pert, double **v_s);
 														
int writeDispersion_bb_bin3(string filen, double freq, int Nfreq, double df, int Nz_grid, double z_min, \
 														int Nz_subgrid, double delZ, int NN, int Nmodes, double dz, \
 														double z_src, double *kreal, double *kim, double **v_s); 														

int SolveModesBB(string filename, int Nfreq, double f_min, double fmax, int Nz_grid, \
		double z_min, double dz, double azi, \
		double sourceheight, NCPA::SampledProfile *prof2, string gnd_imp_model, \
		int Lamb_wave_BC, int flag_2D, \
		bool out_dispersion, bool out_disp_src2grnd);
		
double half_hann(int begin,int end,int i);


void fft_pulse_prop(double t0, int n_freqs, double df, double *f_vec, \
                    double range, complex<double> *dft_vec, \
                    complex<double> *pulse_vec, \
                    int *mode_count, double rho_zsrc, double rho_zrcv, \
                    complex<double> **kc, \
                    complex<double> **mode_S, complex<double> **mode_R);                    								
                    
								
int pulse_prop(	const char *filename,double t0,double RR,  \
								int n_freqs, double f_step, double *f_vec, double scale, \
								int *mode_count, double rho_zsrc, double rho_zrcv, \
								complex<double> **kc, \
								complex<double> **mode_S, complex<double> **mode_R);						
					
int pulse_prop_src2rcv_grid2(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, double f_step, double *f_vec, \
					double f_center, int *mode_count, double rho_zsrc, double rho_zrcv, \
					complex<double> **kc, complex<double> **mode_S, complex<double> **mode_R, \
					int src_flg, string srcfile, int pprop_src2rcv_flg);

int get_source_spectrum( \
								int n_freqs, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile);									 																
									
complex<double> pulse_spec_fit(double scale, double x);		

//int plotXY(double* xData, double* yData, int dataSize, double RR);

//int plotInitialPulse(int src_flg, double freq, double fmax);

int getFile_list(string dir, list<string> &files, string pattern);	

complex<double> ***getPz1z2 (int I1, int I2, double r1, double r2, double dr, string dirn, list<string> files);

// comparison, freq in filename.
bool compare_freq (string first, string second);

int process2DPressure(double R_start_km, double width_km, double height_km, \
											double c_ref, double tmstep, int ntsteps, double f_center, \
											string framefn, string dir_name);

int saveMatrix2bin(const char *filename, double **M, int nr, int nc);

int saveAtm_profile(NCPA::SampledProfile *p, std::string wind_units);

int load_source_pulse_td(std::string srcpulsetdfn, vector<double> &t, vector<double> &tdp );

int load_source_spectrum(string srcspfn, double *freqv, complex<double> *dft_vec, int FFTN );

