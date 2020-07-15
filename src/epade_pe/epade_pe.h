#ifndef NCPAPROP_EPADE_PE_H_INCLUDED
#define NCPAPROP_EPADE_PE_H_INCLUDED

#include <vector>

#include "petscksp.h"

#include "Atmosphere1D.h"
#include "Atmosphere2D.h"
#include "parameterset.h"

namespace NCPA {

	class EPadeSolver {

	public:
		EPadeSolver( NCPA::ParameterSet *param );
		~EPadeSolver();
		int computeTLField();
		void output1DTL( std::string filename );
		void output2DTL( std::string filename );

	protected:

		int epade( int order, double k0, double dr, std::vector<PetscScalar> *P, std::vector<PetscScalar> *Q,
			bool starter = false );
		int get_starter_gaussian( size_t NZ, double *z, double zs, double k0, Vec *psi );
		int get_starter_self( size_t NZ, double *z, double zs, double k0, Mat *qpowers, size_t npade, 
			Vec *psi );

		double *z = NULL, *r = NULL, **tl = NULL;
		double freq;
		double azi;
		int NZ, NR;
		int npade;
		bool use_atm_1d = false, use_atm_2d = false, use_atm_toy = false;
		bool z_min_specified = false;
		double r_max;
		double z_max;
		double z_min;
		double zs;
		double zrcv;
		std::string gnd_imp_model;
		std::string starter;

		std::vector< double > zt;
		std::vector< int > zti;
		int nzplot;

		NCPA::Atmosphere1D *atm_profile;
		NCPA::Atmosphere2D *atm_profile_2d;


	};

}







#endif