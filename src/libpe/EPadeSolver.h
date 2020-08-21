#ifndef NCPAPROP_EPADE_PE_H_INCLUDED
#define NCPAPROP_EPADE_PE_H_INCLUDED

#include <vector>

#include "petscksp.h"

#include "AtmosphericTransferFunctionSolver.h"
#include "Atmosphere1D.h"
#include "Atmosphere2D.h"
#include "parameterset.h"

#define NCPAPROP_EPADE_PE_FILENAME_1D "tloss_1d.pe"
#define NCPAPROP_EPADE_PE_FILENAME_2D "tloss_2d.pe"
#define NCPAPROP_EPADE_PE_FILENAME_MULTIPROP "tloss_multiprop.pe"
#define NCPAPROP_EPADE_PE_FILENAME_BROADBAND "tloss_broadband.bin"

namespace NCPA {

	class EPadeSolver : public AtmosphericTransferFunctionSolver {

	public:
		EPadeSolver( NCPA::ParameterSet *param );
		~EPadeSolver();
		int solve();
		void output1DTL( std::string filename, bool append = false );
		void output2DTL( std::string filename );

	protected:

		int epade( int order, double k0, double dr, std::vector<PetscScalar> *P, std::vector<PetscScalar> *Q,
			bool starter = false );
		int get_starter_gaussian( size_t NZ, double *z, double zs, double k0, int ground_index, Vec *psi );
		int get_starter_self( size_t NZ, double *z, double zs, double k0, Mat *qpowers, size_t npade, 
			Vec *psi );
		void absorption_layer( double lambda, double *z, int NZ, double *layer );
		void fill_atm_vector_relative( NCPA::Atmosphere2D *atm, double range, int NZvec, double *zvec, 
			std::string key, double groundheight, double *vec );
		void fill_atm_vector_absolute( NCPA::Atmosphere2D *atm, double range, int NZvec, double *zvec, 
			std::string key, double fill_value, double *vec );

		int make_q_powers( int NZvec, double *zvec, double k0, double h2, std::complex<double> *n, 
			size_t nqp, int boundary_index, Mat *qpowers );
		int make_B_and_C_matrices( Mat *qpowers, int npade, int NZ, 
			std::vector< std::complex< double > > &P, std::vector< std::complex< double > > &Q,
			Mat *B, Mat *C );
		void calculate_atmosphere_parameters( NCPA::Atmosphere2D *atm, int NZvec, double *z_vec, 
			double r, double z_g, bool use_lossless, bool use_top_layer, double freq, bool use_absolute_z,
			double &k0, double &c0, double *c_vec, double *a_vec, std::complex<double> *k_vec, 
			std::complex<double> *n_vec );

		void set_1d_output( bool tf );
		void write_broadband_header( std::string filename, double *az_vec, size_t n_az, 
			double *f_vec, size_t n_f, unsigned int precision_factor );
		void write_broadband_results( std::string filename, double this_az, double this_f, 
			double *r_vec, size_t n_r, double *z_vec, size_t n_z, std::complex< double > **tloss_mat, 
			unsigned int precision_factor );

		double *z = NULL, *z_abs = NULL, *r = NULL, *f = NULL, calc_az;
		std::complex< double > **tl;
		int *zgi_r = NULL;   // ground height index
		double freq;         // current active frequency
		double *azi;
		int NZ, NR, NR_requested, NAz, NF;
		double dz;
		int npade;
		bool use_atm_1d = false, use_atm_2d = false, use_atm_toy = false, use_topo = false;
		bool z_ground_specified = false, lossless = false, top_layer = true;
		bool multiprop = false, write1d = true, write2d = false, calculate_attn = true;
		bool broadband = false;
		double r_max;    // range limits
		double z_max, z_min, z_ground, z_bottom;  // atmosphere profile limits
		double zs, zr;  // source height, receiver height
		double c_underground;
		//double zrcv;
		std::string gnd_imp_model;
		std::string starter;
		std::string attnfile;

		std::vector< double > zt;
		std::vector< int > zti;
		int nzplot;

		double absorption_layer_mu = 0.01;

		//NCPA::Atmosphere1D *atm_profile;
		NCPA::Atmosphere2D *atm_profile_2d;


	};

}







#endif