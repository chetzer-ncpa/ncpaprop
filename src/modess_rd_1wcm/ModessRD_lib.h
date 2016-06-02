//#include <sstream>
#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <list>

//utility functions
double **dmatrix(long nr, long nc);
int      free_dmatrix(double**v, long nr, long nc);
//int      plotwGNUplot(double freq, bool write_2D_TLoss);

NCPA::SampledProfile * get_RngDepnd_profile(std::string env_file, double R_meters);

NCPA::SampledProfile * get_RngDepndProfiles_ascii(int N, std::string atmosfileorder, int skiplines, std::string dirname, std::string pattern);

int saveSampledProfile(std::string filename, NCPA::SampledProfile *p);

int computeTLoss1D(int Nmodes, double rng, double RR, int n_zsrc, std::complex<double> *k_pert, double **v_s, int *signv2, std::complex<double> *Kintg_atR);

int writeDispersion(int select_modes, int n_zsrc, double freq, std::complex<double> *k_pert, double **v_s, std::string file_stub);

int writeEigenVec(int nz, int select_modes, double dz, double **v_s, std::string file_stub);

int getAcurr(int Nz_grid, int Nm_prev, int Nm, double dz, std::complex<double> *A_prev, double **v_prev, std::complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, std::complex<double> *k_curr, NCPA::SampledProfile *atm_profile, std::complex<double> *A_curr);

int getAcurr_ll(int Nz_grid, int Nm_prev, int Nm, double dz, std::complex<double> *A_prev_ll, double **v_prev, std::complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, std::complex<double> *k_curr, NCPA::SampledProfile *atm_profile, std::complex<double> *A_curr_ll);

void parseReqRanges(std::string str, std::vector<double>& retVal);

int saveTLoss1D (int Nm, int n_zsrc, double sqrtrho_s, int n_zrcv, double sqrtrho_r, int it, std::vector<double> Rv, double rng_step, double **v_prev, std::complex<double> *k_prev, std::complex<double> *A_prev, std::complex<double> *A_prev_ll, const char *wa, double *rng);

int saveTLoss2D (int Nm, int Nz_grid, int n_zsrc, double dz_km, int stepj, \
                  double sqrtrho_s, int it, std::vector<double> Rv, double rng_step, \
                  double *rho_prev, double **v_prev, std::complex<double> *k_prev, \
                  std::complex<double> *A_prev, std::complex<double> *A_prev_ll, \
                  const char *wa, double *prng);

int getRegionBoundaries(bool flg, double maxrange, double req_profile_step, std::string prf_ranges_km, int & Nprofiles, std::vector<double> *R);                 
