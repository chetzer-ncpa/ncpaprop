//#include <sstream>
#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <list>

#include <petscksp.h>

//utility functions
double **dmatrix(long nr, long nc);
int      free_dmatrix(double**v, long nr, long nc);

std::complex<double> **cmatrix(long nr, long nc);
int free_cmatrix(std::complex<double> **v, long nr, long nc);

//int      plotwGNUplot(double freq, bool write_2D_TLoss);

NCPA::SampledProfile * get_RngDepnd_profile(std::string env_file, double R_meters);

NCPA::SampledProfile * get_RngDepnd_profile_ascii(double R);

NCPA::SampledProfile * get_RngDepndProfiles_ascii(int N, std::string atmosfileorder, uint skiplines, std::string dirname, std::string pattern, bool print_flg);

int saveSampledProfile(std::string filename, NCPA::SampledProfile *p);

int getAcurr(int Nz_grid, int Nm, double dz, std::complex<double> *A_prev, double **v_prev, std::complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, std::complex<double> *k_curr, NCPA::SampledProfile *atm_profile, std::complex<double> *A_curr);

int getAcurr_ll(int Nz_grid, int Nm, double dz, std::complex<double> *A_prev_ll, double **v_prev, std::complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, std::complex<double> *k_curr, NCPA::SampledProfile *atm_profile, std::complex<double> *A_curr_ll);

void parseReqRanges(std::string str, std::vector<double>& retVal);

int saveTLoss1D (int Nm, int n_zsrc, double sqrtrho_s, int it, std::vector<double> Rv, double rng_step, double **v_prev, std::complex<double> *k_prev, std::complex<double> *A_prev, std::complex<double> *A_prev_ll, const char *wa, double *rng);

int saveTLoss2D (int Nm, int Nz_grid, int n_zsrc, double dz_km, int stepj, \
                  double sqrtrho_s, int it, std::vector<double> Rv, double rng_step, \
                  double **v_prev, std::complex<double> *k_prev, \
                  std::complex<double> *A_prev, std::complex<double> *A_prev_ll, \
                  const char *wa, double *prng);
                   
int getNmodesNz(std::string fn, int *n_modes, int *nz, double *dz);
               
int readEigenValVecs(std::string fn, std::complex<double> *k_pert, double *rho, double **v, int how_many); 

int getRRmats(std::string fn, double r1, double r2, int n_modes, int nz, double dz, \
              std::complex<double> *k_curr, double *rho_curr, double **v_curr, \
              std::complex<double> *k_next, double *rho_next, double **v_next, \
              std::complex<double> **RRR);
              
int getRRmats_ll(std::string fn, double r1, double r2, int n_modes, int nz, double dz, \
              std::complex<double> *k_curr, double *rho_curr, double **v_curr, \
              std::complex<double> *k_next, double *rho_next, double **v_next, \
              std::complex<double> **RRR);                          

int readRRmat(std::string fn, std::complex<double> **RRR, int n);

int MatCMultiply(std::complex<double> **A, std::complex<double> **B, std::complex<double> **C, int M1, int N1, int M2, int N2);
int MatVecCMultiply(std::complex<double> **A, std::complex<double> *B, std::complex<double> *C, int M1, int N1, int M2);

//int SolveLinSys( complex<double> *avec );
int SolveLinSys(PetscScalar *avec);
int SolveLinSys2(PetscScalar **AA, PetscScalar *bb, PetscScalar *yy, PetscInt N);

int save_printCMatrix(std::complex<double> **A, int nr, int nc, std::string filename, bool flag);
int save_printCVector(std::complex<double> *A, int n, std::string filename, bool flag);

int computePressure_1D( int it, int Nm, int n_z, std::vector<double> Rv, double rng_step, double sqrtrho_z, std::complex<double> *k_curr,  double **v_curr, std::complex<double> *ab_curr, const char *wa, double *prng_curr);

// include lossless version
int computePressure_1D( int it, int Nm, int n_z, std::vector<double> Rv, double rng_step, double sqrtrho_z, std::complex<double> *k_curr,  double **v_curr, std::complex<double> *ab_curr, std::complex<double> *ab_curr_ll, const char *wa, double *prng_curr);

int computePressure_2D( int it, int Nm, int Nz, int stepn, std::vector<double> Rv, double rng_step, double dz, double *rho_curr, std::complex<double> *k_curr,  double **v_curr, std::complex<double> *ab_curr, std::complex<double> *ab_curr_ll, const char *wa, double *prng_curr);

int makeYYYY_MM_DD_subdir(std::string *subdir); 

int getRegionBoundaries(bool flg, double maxrange, double req_profile_step, std::string prf_ranges_km, int *Nprofiles, std::vector<double> *Rv);

int getCeffMinMax(NCPA::SampledProfile *p, double *ceffmin, double *ceffmax); // obtain ceffmin, ceffmax of a profile

int getGlobalceffMinMax(int Nprofiles, std::vector<double> Rv, double azi, std::string atm_profile_dir, std::string atmosfile, std::string atmosfileorder, int skiplines, double *ceffMin, double *ceffMax);


