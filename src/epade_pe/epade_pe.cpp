#include "epade_pe.h"
#include "epade_pe_parameters.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <string>
#include <vector>
#include <cfloat>
#include <fstream>

#include "petscksp.h"

#include "Atmosphere1D.h"
#include "ToyAtmosphere1D.h"
#include "units.h"
#include "util.h"


#ifndef PI
#define PI 3.14159
#endif


NCPA::EPadeSolver::EPadeSolver( NCPA::ParameterSet *param ) {

	// obtain the parameter values from the user's options
	// @todo add units to input scalar quantities
	gnd_imp_model 		= param->getString( "ground_impedence_model" );
	r_max	 			= param->getFloat( "maxrange_km" ) * 1000.0;
  	z_max	 			= param->getFloat( "maxheight_km" ) * 1000.0;      // @todo fix elsewhere that m is required
  	zs			 		= param->getFloat( "sourceheight_km" ) * 1000.0;
  	zrcv		 		= param->getFloat( "receiverheight_km" ) * 1000.0;
  	NZ 					= param->getInteger( "Nz_grid" );    // @todo calculate dynamically
  	NR 					= param->getInteger( "Nrng_steps" );
  	azi 				= param->getFloat( "azimuth" );
	freq 				= param->getFloat( "freq" );
	npade 				= param->getInteger( "npade" );
	starter 			= param->getString( "starter" );

	double dr;
  	if (NR == 0) {
  		dr = 340.0 / freq;
		NR = (int)ceil( r_max / dr );
  	} else {
  		dr = r_max / NR;
  	}
  	r = new double[ NR ];
  	std::memset( r, 0, NR * sizeof(double) );
  	int i;
  	for (i = 0; i < NR; i++) {
  		r[ i ] = ((double)(i+1)) * dr;
  	}

  	lossless 			= param->wasFound( "lossless" );
	use_atm_1d			= param->wasFound( "atmosfile" );
	use_atm_2d			= param->wasFound( "atmosfile2d" );
	use_atm_toy			= param->wasFound( "toy" );
	top_layer			= !(param->wasFound( "disable_top_layer" ));

	if (use_atm_1d) {
		atm_profile = new NCPA::Atmosphere1D( param->getString( "atmosfile" ) );
	} else if (use_atm_toy) {
		atm_profile = new NCPA::ToyAtmosphere1D();
	} else {
		std::cerr << "2-D atmosphere not yet implemented" << std::endl;
		exit(0);
	}

	atm_profile->convert_altitude_units( Units::fromString( "m" ) );
	atm_profile->convert_property_units( "Z0", Units::fromString( "m" ) );
	atm_profile->convert_property_units( "U", Units::fromString( "m/s" ) );
	atm_profile->convert_property_units( "V", Units::fromString( "m/s" ) );
	atm_profile->convert_property_units( "T", Units::fromString( "K" ) );
	atm_profile->convert_property_units( "P", Units::fromString( "Pa" ) );
	atm_profile->convert_property_units( "RHO", Units::fromString( "kg/m3" ) );


	z_max = NCPA::min( z_max, atm_profile->get_maximum_altitude() );
	z_min = atm_profile->get_minimum_altitude();
	z_ground = z_min;
	if (atm_profile->contains_scalar( "Z0" )) {
		z_ground = atm_profile->get( "Z0" );
		if (z_ground < z_min) {
			std::cerr << "Supplied ground height is outside of atmospheric specification." << std::endl;
			exit(0);
		}
	}
	
  
	// fill and convert to SI units
	double dz       = (z_max - z_ground)/(NZ - 1);	// the z-grid spacing
	z = new double[ NZ ];
	z_abs = new double[ NZ ];
	std::memset( z, 0, NZ * sizeof( double ) );
	for (i = 0; i < NZ; i++) {
		z[ i ]     = ((double)i) * dz;
		z_abs[ i ] = z[ i ] + z_ground;
	}
	/*
	atm_profile->resample( dz );
	NZ = atm_profile->nz();
	z = new double[ NZ ];
	atm_profile->get_altitude_vector( z );
	*/
	tl = NCPA::dmatrix( NZ, NR-1 );

	int plotz = 10;
	for (i = 0; i < NZ; i += plotz) {
		zt.push_back( z_abs[ i ] );
		zti.push_back( i );
	}
	nzplot = zt.size();

	//double z_min_km = z_min / 1000.0;
	//double dz_km    = dz / 1000.0;
  
	// Note: the rho, Pr, T, zw, mw are computed wrt ground level i.e.
	// the first value is at the ground level e.g. rho[0] = rho(z_min)
	atm_profile->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO", Units::fromString( "m/s" ) );
	atm_profile->calculate_wind_speed( "_WS_", "U", "V" );
	atm_profile->calculate_wind_direction( "_WD_", "U", "V" );
	atm_profile->calculate_wind_component( "_WC_", "_WS_", "_WD_", azi );
	atm_profile->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );
	atm_profile->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", freq );
}

NCPA::EPadeSolver::~EPadeSolver() {
	delete [] r;
	delete [] z;
	delete [] z_abs;
	NCPA::free_dmatrix( tl, NZ, NR-1 );
	delete atm_profile;
}

void NCPA::EPadeSolver::fill_atm_vector( NCPA::Atmosphere1D *atm, int NZvec, double *zvec, std::string key,
	double groundheight, double *vec ) {

	for (int i = 0; i < NZvec; i++) {
		vec[i] = atm->get( key, zvec[i] + groundheight );
	}
}

int NCPA::EPadeSolver::computeTLField() {
	int i;
	std::complex<double> I( 0.0, 1.0 );
	PetscErrorCode ierr;
	PetscInt Istart, Iend, col[3], *indices;
	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar    value[3];  // for populating tridiagonal matrices
	PetscScalar hank, *contents;
	Mat q, B, C;
	Mat *qpowers;
	Vec psi_o, Bpsi_o; //, psi_temp;
	KSP ksp;
	// PC pc;

	double omega = 2.0 * PI * freq;
	//double gamma = 1.4;
	//double z_top = z_max;
	double dr = r[1] - r[0];
	double h = z[1] - z[0];
	double h2 = h * h;
	//double zrcv = receiverheight;

	std::cout << std::endl << "Infrasound PE code at f = " << freq << " Hz, azi = " 
		<< azi << " deg" << std::endl;

	double *c = new double[ NZ ];
	memset( c, 0, NZ * sizeof(double) );
	double *a_t = new double[ NZ ];
	memset( a_t, 0, NZ * sizeof(double) );
	//atm_profile->get_property_vector( "_CEFF_", c );
	fill_atm_vector( atm_profile, NZ, z, "_CEFF_", z_ground, c );
	double c0 = NCPA::mean( c, NZ );

	// @todo implement lossless version
	if (!lossless) {
		//atm_profile->get_property_vector( "_ALPHA_", a_t );
		fill_atm_vector( atm_profile, NZ, z, "_ALPHA_", z_ground, a_t );
	}
	double *abslayer = new double[ NZ ];
	memset( abslayer, 0, NZ * sizeof(double) );
	if (top_layer) {
		absorption_layer( c0 / freq, z, NZ, abslayer );
	}
	

	
	// check vertical sampling
	double z_cnd = (c0 / freq) / 10.0;
	if (h > z_cnd) {
		std::cerr << "Altitude sampling is too low!  (is " << h << ", should be <= "<< z_cnd << std::endl;
		exit(0);
	}

	double k0 = omega / c0;
	std::complex<double> *k = new std::complex<double>[ NZ ];
	std::complex<double> *n = new std::complex<double>[ NZ ];
	std::memset( k, 0, NZ * sizeof( std::complex< double > ) );
	std::memset( n, 0, NZ * sizeof( std::complex< double > ) );
	// Set up vectors
	indices = new PetscInt[ NZ ];
	for (i = 0; i < NZ; i++) {
		k[ i ] = omega / c[ i ] + (a_t[ i ] + abslayer[ i ]) * I;
		n[ i ] = k[ i ] / k0;
		indices[ i ] = i;
	}

	// Set up matrices
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 3, NULL, &q );CHKERRQ(ierr);
	ierr = MatSetFromOptions( q );CHKERRQ(ierr);

	// populate
	double bnd_cnd = -1.0 / h2;    // @todo add hook for alternate boundary conditions
	double k02 = k0*k0;
	//gsl_complex tempdiag;
	ierr = MatGetOwnershipRange(q,&Istart,&Iend);CHKERRQ(ierr);
	if (Istart==0) FirstBlock=PETSC_TRUE;
    if (Iend==NZ) LastBlock=PETSC_TRUE;
    value[0]=1.0 / h2 / k02; value[2]=1.0 / h2 / k02;
    for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {
    		value[ 1 ] = -2.0/h2/k02 + (n[i]*n[i] - 1);
		    col[0]=i-1; col[1]=i; col[2]=i+1;
		    ierr = MatSetValues(q,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (LastBlock) {
		    i=NZ-1; col[0]=NZ-2; col[1]=NZ-1;
		    value[ 1 ] = -2.0/h2/k02 + (n[i]*n[i] - 1);
		    ierr = MatSetValues(q,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (FirstBlock) {
		    i=0; col[0]=0; col[1]=1; 
		    value[0]=bnd_cnd/k02 + (n[i]*n[i] - 1); value[1]=1.0 / h2 / k02;
		    ierr = MatSetValues(q,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	std::cout << "Finding ePade coefficients..." << std::endl;
	std::vector< std::complex<double> > P, Q;
	epade( npade, k0, dr, &P, &Q );

	// tested to this point against pe_ri_n2_epade.m
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade-1, NULL, &B );CHKERRQ(ierr);
	ierr = MatSetFromOptions( B );CHKERRQ(ierr);
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade+1, NULL, &C );CHKERRQ(ierr);
	ierr = MatSetFromOptions( C );CHKERRQ(ierr);

	ierr = MatGetOwnershipRange(B,&Istart,&Iend);CHKERRQ(ierr);
	value[ 0 ] = 1.0;
	for (i = Istart; i < Iend; i++) {
		ierr = MatSetValues( B, 1, &i, 1, &i, value, INSERT_VALUES );CHKERRQ(ierr);
		ierr = MatSetValues( C, 1, &i, 1, &i, value, INSERT_VALUES );CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	// calculate powers of q
	qpowers = new Mat[ npade+1 ];
	ierr = MatConvert( q, MATSAME, MAT_INITIAL_MATRIX, qpowers );CHKERRQ(ierr);
	for (i = 1; i < npade+1; i++) {
		ierr = MatMatMult( qpowers[i-1], q, MAT_INITIAL_MATRIX, PETSC_DEFAULT, 
			qpowers+i );CHKERRQ(ierr);
	}

	for (i = 1; i < Q.size(); i++) {
		ierr = MatAXPY( C, Q[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	for (i = 1; i < P.size(); i++) {
		ierr = MatAXPY( B, P[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}

	if (starter == "self") {
		get_starter_self( NZ, z, zs, k0, qpowers, npade, &psi_o );
	} else if (starter == "gaussian") {
		get_starter_gaussian( NZ, z, zs, k0, &psi_o );
	} else {
		std::cerr << "Unrecognized starter type: " << starter << std::endl;
		exit(0);
	}

	std::cout << "Marching out field..." << std::endl;
	ierr = VecDuplicate( psi_o, &Bpsi_o );CHKERRQ(ierr);
	contents = new PetscScalar[ NZ ];

	ierr = KSPCreate( PETSC_COMM_SELF, &ksp );CHKERRQ(ierr);
	ierr = KSPSetOperators( ksp, C, C );CHKERRQ(ierr);
	// ierr = KSPGetPC( ksp, &pc );CHKERRQ(ierr);
	// ierr = PCSetType( pc, PCLU );
	ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
	for (PetscInt ir = 0; ir < (NR-1); ir++) {

		//double rr = ((double)ir+1) * dr;
		//r[ ir ] = rr;
		double rr = r[ ir ];
		hank = sqrt( 2.0 / ( PI * k0 * rr ) ) * exp( I * ( k0 * rr - PI/4.0 ) );
		//ierr = VecCopy( psi_o, psi_temp );CHKERRQ(ierr);
		//ierr = VecScale( psi_temp, hank );CHKERRQ(ierr);
		ierr = VecGetValues( psi_o, NZ, indices, contents );
		for (i = 0; i < nzplot; i++) {
			tl[ i ][ ir ] = 20.0 * log10( abs( contents[ zti[i] ] * hank ) );
		}

		if ( fmod( rr, 1.0e5 ) < dr) {
			std::cout << " -> Range " << rr/1000.0 << " km" << std::endl;
		}

		ierr = MatMult( B, psi_o, Bpsi_o );CHKERRQ(ierr);
		ierr = KSPSetOperators( ksp, C, C );CHKERRQ(ierr);  // may not be necessary
		// ierr = KSPGetPC( ksp, &pc );CHKERRQ(ierr);
		// ierr = PCSetType( pc, PCLU );
		ierr = KSPSolve( ksp, Bpsi_o, psi_o );CHKERRQ(ierr);
	}

	ierr = MatDestroy( &q );       CHKERRQ(ierr);
	ierr = MatDestroy( &B );       CHKERRQ(ierr);
	ierr = MatDestroy( &C );       CHKERRQ(ierr);
	ierr = VecDestroy( &psi_o );   CHKERRQ(ierr);
	ierr = VecDestroy( &Bpsi_o );  CHKERRQ(ierr);
	ierr = KSPDestroy( &ksp );     CHKERRQ(ierr);
	for (i = 0; i < npade+1; i++) {    // for some reason throws an error at i == npade
		ierr = MatDestroy( qpowers + i );CHKERRQ(ierr);
	}
	delete [] qpowers;
	delete [] k;
	delete [] n;
	delete [] c;
	delete [] a_t;
	delete [] contents;
	delete [] indices;
	delete [] abslayer;

	return 1;
}

void NCPA::EPadeSolver::absorption_layer( double lambda, double *z, int NZ, double *layer ) {
	double z_t = z[NZ-1] - lambda;
	for (int i = 0; i < NZ; i++) {
		layer[ i ] = absorption_layer_mu * std::exp( (z[i]-z_t) / lambda );
	}
}


int NCPA::EPadeSolver::get_starter_gaussian( size_t NZ, double *z, double zs, double k0, Vec *psi ) {

	double fac = 2.0;
	//double kfac = k0 / fac;
	PetscScalar tempval;
	PetscErrorCode ierr;

	ierr = VecCreate( PETSC_COMM_SELF, psi );CHKERRQ(ierr);
	ierr = VecSetSizes( *psi, PETSC_DECIDE, NZ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( *psi ); CHKERRQ(ierr);

	for (PetscInt i = 0; i < NZ; i++) {
		tempval = -( k0*k0/fac/fac ) * (z[i] - zs) * (z[i] - zs);
		tempval = sqrt( k0/fac ) * exp( tempval );
		ierr = VecSetValues( *psi, 1, &i, &tempval, INSERT_VALUES );CHKERRQ(ierr);
	}
	ierr = VecAssemblyBegin( *psi );CHKERRQ(ierr);
	ierr = VecAssemblyEnd( *psi );CHKERRQ(ierr);
	return 1;
}

int NCPA::EPadeSolver::get_starter_self( size_t NZ, double *z, double zs, double k0, Mat *qpowers, 
	size_t npade, Vec *psi ) {

	Vec rhs, ksi, Bksi, tempvec;
	Mat A, AA, B, C;
	KSP ksp, ksp2;
	PetscScalar I( 0.0, 1.0 ), tempsc;
	PetscInt ii, Istart, Iend;
	PetscErrorCode ierr;
	
	ierr = VecCreate( PETSC_COMM_SELF, &rhs );CHKERRQ(ierr);
	ierr = VecSetSizes( rhs, PETSC_DECIDE, NZ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( rhs );CHKERRQ(ierr);
	ierr = VecSet( rhs, 0.0 );CHKERRQ(ierr);
	
	// find closest index to zs
	PetscInt nzsrc = (PetscInt)find_closest_index( z, NZ, zs );

	double h = z[1] - z[0];
	PetscScalar hinv = 1.0 / h;
	ierr = VecSetValues( rhs, 1, &nzsrc, &hinv, INSERT_VALUES );CHKERRQ(ierr);
	
	// set up identity matrix.  If this works, use it elsewhere
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 3, NULL, &A );CHKERRQ(ierr);
	ierr = MatSetFromOptions( A );CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
	tempsc = 1.0;
    for (ii = Istart; ii < Iend; ii++) {
    	ierr = MatSetValues( A, 1, &ii, 1, &ii, &tempsc, INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAXPY( A, -I, qpowers[0], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);

	// square
	ierr = MatMatMult( A, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA );
	
	// solve first part (Eq. 26)
	ierr = VecDuplicate( rhs, &ksi );CHKERRQ(ierr);
	ierr = VecSet( ksi, 0.0 );CHKERRQ(ierr);
	ierr = KSPCreate( PETSC_COMM_WORLD, &ksp );CHKERRQ(ierr);
	ierr = KSPSetOperators( ksp, AA, AA );CHKERRQ(ierr);
	ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
	ierr = KSPSolve( ksp, rhs, ksi );CHKERRQ(ierr);
	
	// get starter
	std::cout << "Finding ePade starter coefficients..." << std::endl;
	double r_ref = 2 * PI / k0;
	std::vector<PetscScalar> P, Q;
	epade( npade, k0, r_ref, &P, &Q, true );

	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade-1, NULL, &B );CHKERRQ(ierr);
	ierr = MatSetFromOptions( B );CHKERRQ(ierr);
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade+1, NULL, &C );CHKERRQ(ierr);
	ierr = MatSetFromOptions( C );CHKERRQ(ierr);

	ierr = MatGetOwnershipRange(B,&Istart,&Iend);CHKERRQ(ierr);
	PetscScalar value = 1.0;
	for (ii = Istart; ii < Iend; ii++) {
		ierr = MatSetValues( B, 1, &ii, 1, &ii, &value, INSERT_VALUES );CHKERRQ(ierr);
		ierr = MatSetValues( C, 1, &ii, 1, &ii, &value, INSERT_VALUES );CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    for (ii = 1; ii < Q.size(); ii++) {
		ierr = MatAXPY( C, Q[ ii ], qpowers[ ii-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	for (ii = 1; ii < P.size(); ii++) {
		ierr = MatAXPY( B, P[ ii ], qpowers[ ii-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}

	PetscScalar hank_inv = pow( sqrt( 2.0 / ( PI * k0 * r_ref ) ) * exp( I * (k0 * r_ref - PI/4.0 ) ),
		-1.0 );

	// Original Matlab: psi = AA * ( C \ (B * ksi) ) / hank
	// compute product of B and ksi
	ierr = VecDuplicate( ksi, &Bksi );
	ierr = VecDuplicate( ksi, &tempvec );
	ierr = VecDuplicate( ksi, psi );
	ierr = MatMult( B, ksi, Bksi );

	// solve for tempvec = C \ Bksi
	ierr = KSPCreate( PETSC_COMM_WORLD, &ksp2 );CHKERRQ(ierr);
	ierr = KSPSetOperators( ksp2, C, C );CHKERRQ(ierr);
	ierr = KSPSetFromOptions( ksp2 );CHKERRQ(ierr);
	ierr = KSPSolve( ksp2, Bksi, tempvec );CHKERRQ(ierr);
	
	// multiply and scale
	ierr = MatMult( AA, tempvec, *psi );CHKERRQ(ierr);
	ierr = VecScale( *psi, hank_inv );CHKERRQ(ierr);

	// clean up
	ierr = MatDestroy( &A );CHKERRQ(ierr);
	ierr = MatDestroy( &AA );CHKERRQ(ierr);
	ierr = MatDestroy( &B );CHKERRQ(ierr);
	ierr = MatDestroy( &C );CHKERRQ(ierr);
	ierr = VecDestroy( &rhs );CHKERRQ(ierr);
	ierr = VecDestroy( &ksi );CHKERRQ(ierr);
	ierr = VecDestroy( &Bksi );CHKERRQ(ierr);
	ierr = VecDestroy( &tempvec );CHKERRQ(ierr);
	ierr = KSPDestroy( &ksp );CHKERRQ(ierr);
	ierr = KSPDestroy( &ksp2 );CHKERRQ(ierr);

	return 1;
}

int NCPA::EPadeSolver::epade( int order, double k0, double dr, std::vector<PetscScalar> *P, 
		std::vector<PetscScalar> *Q, bool starter ) {

	int M = order, N = 2*order;
	int L = N - 1 - M;
	std::complex<double> j( 0.0, 1.0 );
	double delta = k0 * dr;
	int m;

	PetscErrorCode ierr;
	PetscInt Istart, Iend, ii, jj, pM = M, pL = L, *indices;
	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar tempsc, *contents;
	Mat A;
	Vec x, y;
	KSP ksp;
	//PC pc;

	// Create the temporary coefficient vector c
	// using complex<double> to make use of the built-in calculations
	std::complex<double> *c = new std::complex<double>[ N ];
	std::memset( c, 0, N*sizeof(std::complex<double>));
	c[ 0 ] = 1.0;
	c[ 1 ] = j * 0.5 * delta;
	for (m = 1; m <= (2*M - 2); m++) {
		int idx = m+1;
		double dm = (double)m;
		c[ idx ] = -((2.0*dm - 1.0) / (2.0*dm + 2.0)) * c[idx-1]
				   - (delta*delta / (4.0*dm*(dm+1.0))) * c[idx-2];
	}

	// apply modification for starter calculation
	if (starter) {
		std::complex<double> *d = new std::complex<double>[ N ];
		std::memset( d, 0, N*sizeof(std::complex<double>) );
		d[0] = 1.0;
		for (m = 1; m <= (2*M - 1); m++) {
			d[ m ] = (j*delta*c[m-1] - (((double)(1 + 2*(m-1))) * d[m-1])) / (2.0*m);
		}
		std::memcpy( c, d, N*sizeof(std::complex<double>) );
		delete [] d;
	}

	// Create and populate matrix system A
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, N-1, N-1, M+1, NULL, &A );CHKERRQ(ierr);
	ierr = MatSetFromOptions( A );CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
	if (Istart==0) FirstBlock=PETSC_TRUE;
    if (Iend==(N-1)) LastBlock=PETSC_TRUE;
    tempsc = 0.0;
    for (ii = Istart; ii < Iend; ii++) {
    	ierr = MatSetValues( A, 1, &ii, 1, &ii, &tempsc, INSERT_VALUES);CHKERRQ(ierr);
    }
    for (ii = Istart; ii < min(pM,Iend-1); ii++) {
    	for (jj = ii; jj < (N-1); jj++) {
    		ierr = MatSetValues(A,1,&jj,1,&ii,c+(jj-ii),INSERT_VALUES);CHKERRQ(ierr);
    	}
    }
    tempsc = -1.0;
    for (ii = Istart; ii < min(pL,Iend); ii++) {
    	jj = ii+pM;
    	ierr = MatSetValues(A,1,&ii,1,&jj,&tempsc,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

/*
    cout << "A Matrix Contents:" << endl;
    ierr = MatView( A, PETSC_VIEWER_STDOUT_SELF );CHKERRQ(ierr);
    cout << endl;
*/

	// setup right-side vector
	ierr = VecCreate( PETSC_COMM_SELF, &x );CHKERRQ(ierr);
	ierr = VecSetSizes( x, PETSC_DECIDE, N-1 );CHKERRQ(ierr);
	ierr = VecSetFromOptions( x ); CHKERRQ(ierr);
	ierr = VecDuplicate( x, &y );CHKERRQ(ierr);

	indices = new PetscInt[ N-1 ];
	for (ii = 0; ii < (N-1); ii++) {
		tempsc = -c[ ii+1 ];
		ierr = VecSetValues( y, 1, &ii, &tempsc, INSERT_VALUES );CHKERRQ(ierr);
		indices[ ii ] = ii;
	}
	ierr = VecAssemblyBegin( y );CHKERRQ(ierr);
	ierr = VecAssemblyEnd( y );CHKERRQ(ierr);

/*
	cout << "y Vector contents:" << endl;
	ierr = VecView( y, PETSC_VIEWER_STDOUT_SELF );CHKERRQ(ierr);
	cout << endl;
*/

	ierr = VecSet( x, 0.0 );CHKERRQ(ierr);

	// solve
	ierr = KSPCreate( PETSC_COMM_WORLD, &ksp );CHKERRQ(ierr);
	ierr = KSPSetOperators( ksp, A, A );CHKERRQ(ierr);
	//ierr = KSPGetPC( ksp, &pc );CHKERRQ(ierr);
	//ierr = PCSetType( pc, PCLU );CHKERRQ(ierr);
	ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
	ierr = KSPSolve( ksp, y, x );CHKERRQ(ierr);
/*
	cout << "x Vector contents:" << endl;
	ierr = VecView( x, PETSC_VIEWER_STDOUT_SELF );CHKERRQ(ierr);
	cout << endl;
*/

	// populate P and Q vectors as complex<double> instead of gsl_complex
	Q->push_back( 1.0 );
	contents = new PetscScalar[ N-1 ];
	ierr = VecGetValues( x, N-1, indices, contents );

	for (ii = 0; ii < M; ii++) {
		Q->push_back( contents[ ii ] );
	}
	P->push_back( c[ 0 ] );
	for (ii = M; ii < (M+L); ii++) {
		P->push_back( contents[ ii ] );
	}
	delete [] contents;
	delete [] indices;
	delete [] c;
	
	// clean up memory
	ierr = KSPDestroy( &ksp );CHKERRQ(ierr);
	ierr = VecDestroy( &x );CHKERRQ(ierr);
	ierr = VecDestroy( &y );CHKERRQ(ierr);
	ierr = MatDestroy( &A );CHKERRQ(ierr);
	return 0;
}

void NCPA::EPadeSolver::output1DTL( std::string filename ) {
	std::ofstream out_1d( "tloss_1d.pe", std::ofstream::out | std::ofstream::trunc );
	for (int i = 0; i < (NR-1); i++) {
		out_1d << r[ i ]/1000.0 << " " << tl[ 0 ][ i ] << std::endl;
	}
	out_1d.close();
}

void NCPA::EPadeSolver::output2DTL( std::string filename ) {
	std::ofstream out_2d( "tloss_2d.pe", std::ofstream::out | std::ofstream::trunc );
	for (int i = 0; i < (NR-1); i++) {
		for (int j = 0; j < nzplot; j++) {
			out_2d << r[ i ]/1000.0 << " " << zt[ j ]/1000.0 << " " << tl[ j ][ i ] << " 0.0" << std::endl;
		}
		out_2d << std::endl;
	}
	out_2d.close();
}