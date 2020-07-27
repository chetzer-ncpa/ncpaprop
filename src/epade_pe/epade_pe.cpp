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
#include "Atmosphere2D.h"
#include "ToyAtmosphere1D.h"
#include "StratifiedAtmosphere2D.h"
#include "ProfileSeriesAtmosphere2D.h"
#include "units.h"
#include "util.h"


#ifndef PI
#define PI 3.14159
#endif


NCPA::EPadeSolver::EPadeSolver( NCPA::ParameterSet *param ) {

	c_underground		= 5000.0;

	// obtain the parameter values from the user's options
	// @todo add units to input scalar quantities
	gnd_imp_model 		= param->getString( "ground_impedence_model" );
	r_max	 			= param->getFloat( "maxrange_km" ) * 1000.0;
  	z_max	 			= param->getFloat( "maxheight_km" ) * 1000.0;      // @todo fix elsewhere that m is required
  	zs			 		= param->getFloat( "sourceheight_km" ) * 1000.0;
  	NR 					= param->getInteger( "Nrng_steps" );
  	freq 				= param->getFloat( "freq" );
	npade 				= param->getInteger( "npade" );
	starter 			= param->getString( "starter" );

	// flags
	lossless 			= param->wasFound( "lossless" );
	size_t n_atm_selected = 0;
	use_atm_1d			= param->wasFound( "atmosfile" );
	use_atm_2d			= param->wasFound( "atmosfile2d" );
	use_atm_toy			= param->wasFound( "toy" );
	top_layer			= !(param->wasFound( "disable_top_layer" ));
	use_topo			= param->wasFound( "topo" );
	write2d 			= param->wasFound( "write_2d_tloss" );
	multiprop 			= param->wasFound( "multiprop" );

	// Handle differences based on single vs multiprop
	double min_az, max_az, step_az;
	if (multiprop) {
		if (use_atm_2d) {
			std::cerr << "Range-dependent 2-D atmosphere incompatible with multiple azimuth propagation"
					  << std::endl;
			exit(0);
		}
		if (write2d) {
			std::cout << "Multi-azimuth propagation requested, disabling 2-D output" << std::endl;
			write2d = false;
		}
		if (use_topo) {
			std::cout << "Multi-azimuth propagation requested, disabling topography flag" << std::endl;
		}
		min_az 			= param->getFloat( "azimuth_start" );
		max_az 			= param->getFloat( "azimuth_end" );
		step_az 		= param->getFloat( "azimuth_step" );
		NAz 			= (int) ((max_az - min_az)/step_az) + 1;
	} else {
		NAz 			= 1;
		min_az 			= param->getFloat( "azimuth" );
		max_az 			= min_az;
		step_az 		= 0;
	}
	azi = new double[ NAz ];
	memset( azi, 0, NAz * sizeof( double ) );
	for (int i = 0; i < NAz; i++) {
		azi[ i ] = min_az + i * step_az;
	}

	double dr;
  	if (NR == 0) {
  		dr = 340.0 / freq;
		NR = (int)ceil( r_max / dr );
  	} else {
  		dr = r_max / NR;
  	}
  	r = new double[ NR ];
  	std::memset( r, 0, NR * sizeof(double) );
  	zgi_r = new int[ NR ];
  	std::memset( zgi_r, 0, NR * sizeof( int ) );
  	int i;
  	for (i = 0; i < NR; i++) {
  		r[ i ] = ((double)(i+1)) * dr;
  	}

  	

	//NCPA::Atmosphere1D *atm_profile_1d;
	if (use_atm_1d) {
		atm_profile_2d = new NCPA::StratifiedAtmosphere2D( param->getString( "atmosfile" ) );
	} else if (use_atm_toy) {
		NCPA::Atmosphere1D *tempatm = new NCPA::ToyAtmosphere1D();
		atm_profile_2d = new NCPA::StratifiedAtmosphere2D( tempatm );
		delete tempatm;
	} else if (use_atm_2d) {
		atm_profile_2d = new NCPA::ProfileSeriesAtmosphere2D( param->getString( "atmosfile2d" ) );
		atm_profile_2d->convert_range_units( NCPA::Units::fromString( "m" ) );
	} else {
		std::cerr << "Unknown atmosphere option selected" << std::endl;
		exit(0);
	}

	z_min = atm_profile_2d->get_minimum_altitude( 0.0 );
	z_ground = z_min;
	if (param->wasFound("groundheight_km")) {
		z_ground = param->getFloat( "groundheight_km" ) * 1000.0;
		z_ground_specified = true;
	}

	// set units
	atm_profile_2d->convert_altitude_units( Units::fromString( "m" ) );
	atm_profile_2d->convert_property_units( "Z0", Units::fromString( "m" ) );
	atm_profile_2d->convert_property_units( "U", Units::fromString( "m/s" ) );
	atm_profile_2d->convert_property_units( "V", Units::fromString( "m/s" ) );
	atm_profile_2d->convert_property_units( "T", Units::fromString( "K" ) );
	atm_profile_2d->convert_property_units( "P", Units::fromString( "Pa" ) );
	atm_profile_2d->convert_property_units( "RHO", Units::fromString( "kg/m3" ) );

	// calculate derived quantities
	atm_profile_2d->calculate_sound_speed_from_pressure_and_density( "_C0_", "P", "RHO", Units::fromString( "m/s" ) );
	atm_profile_2d->calculate_wind_speed( "_WS_", "U", "V" );
	atm_profile_2d->calculate_wind_direction( "_WD_", "U", "V" );
	// atm_profile_2d->calculate_wind_component( "_WC_", "_WS_", "_WD_", azi[0] );
	// atm_profile_2d->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );
	if (param->wasFound( "attnfile" ) ) {
		atm_profile_2d->read_attenuation_from_file( "_ALPHA_", param->getString( "attnfile" ) );
	} else {
		atm_profile_2d->calculate_attenuation( "_ALPHA_", "T", "P", "RHO", freq );
	}

	// calculate/check z resolution
	dz = 				param->getFloat( "dz_m" );
	double c0 = atm_profile_2d->get( 0.0, "_C0_", z_ground );
	double lambda0 = c0 / freq;
  	if (dz <= 0.0) {
  		dz = lambda0 / 20.0;
  		double nearestpow10 = std::pow( 10.0, std::floor( std::log10( dz ) ) );
  		double factor = std::floor( dz / nearestpow10 );
  		dz = nearestpow10 * factor;
  		std::cout << "Setting dz to " << dz << " m" << std::endl;
  	}
  	if (dz > (c0 / freq / 10.0) ) {
  		std::cerr << "Altitude resolution is too coarse.  Must be <= " << lambda0 / 10.0 << " meters." 
  			<< std::endl;
  		exit(0);
  	}
}

NCPA::EPadeSolver::~EPadeSolver() {
	delete [] r;
	delete [] z;
	delete [] z_abs;
	delete [] zgi_r;
	delete [] azi;
	NCPA::free_cmatrix( tl, NZ, NR-1 );
	delete atm_profile_2d;
}

int NCPA::EPadeSolver::computeTLField() {
	int i;
	std::complex<double> I( 0.0, 1.0 );
	PetscErrorCode ierr;
	PetscInt Istart, Iend, col[3], *indices;
	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar    value[3];  // for populating tridiagonal matrices
	PetscScalar hank, *contents;
	Mat B, C;   // , q;
	Mat *qpowers, *qpowers_starter;
	Vec psi_o, Bpsi_o; //, psi_temp;
	KSP ksp;
	// PC pc;

	// set up z grid for flat ground.  When we add terrain we will need to move this inside
	// the range loop
	//int profile_index = atm_profile_2d->get_profile_index( 0.0 );
	int profile_index;
	double minlimit, maxlimit;
	atm_profile_2d->get_maximum_altitude_limits( minlimit, maxlimit );
	z_max = NCPA::min( z_max, minlimit );    // lowest valid top value
	int ground_index = 0;

	// truncate multiprop file if needed
	if (multiprop) {
		std::ofstream ofs( "tloss_multiprop.pe", std::ofstream::out | std::ofstream::trunc );
		ofs.close();
	}

	if (use_topo) {
		z_bottom = -5000.0;    // make this eventually depend on frequency
		z_bottom -= fmod( z_bottom, dz );
		z_ground = atm_profile_2d->get( 0.0, "Z0" );
		NZ = ((int)std::floor((z_max - z_bottom) / dz)) + 1;
		z = new double[ NZ ];
		z_abs = new double[ NZ ];
		std::memset( z, 0, NZ * sizeof( double ) );
		std::memset( z_abs, 0, NZ * sizeof( double ) );
		indices = new PetscInt[ NZ ];
		std::memset( indices, 0, NZ*sizeof(PetscInt) );
		for (i = 0; i < NZ; i++) {
			z[ i ] = ((double)i) * dz + z_bottom;
			z_abs[ i ] = z[ i ];
			indices[ i ] = i;
		}
		zs = NCPA::max( zs, z_ground );
		ground_index = NCPA::find_closest_index( z, NZ, z_ground );
		if ( z[ ground_index ] < z_ground ) {
			ground_index++;
		}
		// if (z[ ground_index ] > z_ground && ground_index > 0) {
		// 	ground_index--;
		// }

	} else {
		atm_profile_2d->get_minimum_altitude_limits( minlimit, z_min );
		//z_min = atm_profile_2d->get_highest_minimum_altitude();
		if ( (!z_ground_specified) && atm_profile_2d->contains_scalar( 0.0, "Z0" )) {
			z_ground = atm_profile_2d->get( 0.0, "Z0" );
		}
		if (z_ground < z_min) {
			std::cerr << "Supplied ground height is outside of atmospheric specification." << std::endl;
			exit(0);
		}
	  	z_bottom = z_min;
		// fill and convert to SI units
		//double dz       = (z_max - z_ground)/(NZ - 1);	// the z-grid spacing
		NZ = ((int)std::floor((z_max - z_ground) / dz)) + 1;
		z = new double[ NZ ];
		z_abs = new double[ NZ ];
		std::memset( z, 0, NZ * sizeof( double ) );
		std::memset( z_abs, 0, NZ * sizeof( double ) );
		indices = new PetscInt[ NZ ];
		std::memset( indices, 0, NZ*sizeof(PetscInt) );
		for (i = 0; i < NZ; i++) {
			z[ i ]     = ((double)i) * dz;
			z_abs[ i ] = z[ i ] + z_ground;
			indices[ i ] = i;
		}
		zs = NCPA::max( zs-z_ground+dz, dz );
	}
	tl = NCPA::cmatrix( NZ, NR-1 );
	
	/*
	int plotz = 10;
	for (i = ((int)fmod((double)ground_index,(double)plotz)); i < NZ; i += plotz) {
		zt.push_back( z_abs[ i ] );
		zti.push_back( i );
	}
	nzplot = zt.size();
	*/

	// constants for now
	double omega = 2.0 * PI * freq;
	double dr = r[1] - r[0];
	//double h = z[1] - z[0];
	double h = dz;
	double h2 = h * h;

	// set up for source atmosphere
	double k0 = 0.0, c0 = 0.0;
	double *c = new double[ NZ ];
	double *a_t = new double[ NZ ];
	std::complex<double> *k = new std::complex<double>[ NZ ];
	std::complex<double> *n = new std::complex<double>[ NZ ];
	

	for (int azind = 0; azind < NAz; azind++) {
		std::cout << "Infrasound PE code at f = " << freq << " Hz, azi = " 
			<< azi[ azind ] << " deg" << std::endl;

		profile_index = -1;
		atm_profile_2d->calculate_wind_component( "_WC_", "_WS_", "_WD_", azi[ azind ] );
		atm_profile_2d->calculate_effective_sound_speed( "_CEFF_", "_C0_", "_WC_" );
		calc_az = azi[ azind ];

		//std::cout << "Using atmosphere index " << profile_index << std::endl;
		calculate_atmosphere_parameters( atm_profile_2d, NZ, z, 0.0, z_ground, lossless, top_layer, freq, use_topo,
			k0, c0, c, a_t, k, n );


		// calculate q matrices
		qpowers = new Mat[ npade+1 ];
		//qpowers_starter = new Mat[ npade+1 ];
		make_q_powers( NZ, z, k0, h2, n, npade+1, 0, qpowers );


		if (starter == "self") {
			qpowers_starter = new Mat[ npade+1 ];
			make_q_powers( NZ, z, k0, h2, n, npade+1, ground_index, qpowers_starter );
			//get_starter_self( NZ, z, zs, k0, qpowers, npade, &psi_o );
			get_starter_self( NZ, z, zs, k0, qpowers_starter, npade, &psi_o );
		} else if (starter == "gaussian") {
			qpowers_starter = qpowers;
			get_starter_gaussian( NZ, z, zs, k0, ground_index, &psi_o );
		} else {
			std::cerr << "Unrecognized starter type: " << starter << std::endl;
			exit(0);
		}

		std::cout << "Finding ePade coefficients..." << std::endl;
		std::vector< std::complex<double> > P, Q;
		epade( npade, k0, dr, &P, &Q );
		make_B_and_C_matrices( qpowers_starter, npade, NZ, P, Q, &B, &C );

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
			// check for atmosphere change
			if (atm_profile_2d->get_profile_index( rr ) != profile_index) {
			//if (rr > 20000.0) {

				profile_index = atm_profile_2d->get_profile_index( rr );
				//z_ground = atm_profile_2d->get( rr, "Z0" );
				calculate_atmosphere_parameters( atm_profile_2d, NZ, z, rr, z_ground, lossless, top_layer, freq, 
					use_topo, k0, c0, c, a_t, k, n );
				for (i = 0; i < npade+1; i++) {
					ierr = MatDestroy( qpowers + i );CHKERRQ(ierr);
				}
				make_q_powers( NZ, z, k0, h2, n, npade+1, 0, qpowers );
				epade( npade, k0, dr, &P, &Q );
				ierr = MatZeroEntries( B );CHKERRQ(ierr);
				ierr = MatZeroEntries( C );CHKERRQ(ierr);
				make_B_and_C_matrices( qpowers, npade, NZ, P, Q, &B, &C );
				std::cout << "Switching to atmosphere index " << profile_index 
					<< " at range = " << rr/1000.0 << " km" << std::endl;
			}



			hank = sqrt( 2.0 / ( PI * k0 * rr ) ) * exp( I * ( k0 * rr - PI/4.0 ) );
			//ierr = VecCopy( psi_o, psi_temp );CHKERRQ(ierr);
			//ierr = VecScale( psi_temp, hank );CHKERRQ(ierr);
			ierr = VecGetValues( psi_o, NZ, indices, contents );
			// for (i = 0; i < nzplot; i++) {
			// 	tl[ i ][ ir ] = 20.0 * log10( abs( contents[ zti[i] ] * hank ) );
			// }
			for (i = 0; i < NZ; i++) {
				//tl[ i ][ ir ] = 20.0 * log10( abs( contents[ i ] * hank ) );
				tl[ i ][ ir ] = contents[ i ] * hank;
			}

			if (use_topo) {
				double z0g = atm_profile_2d->get( rr, "Z0" );
				zgi_r[ ir ] = NCPA::find_closest_index( z, NZ, z0g );
				if ( z[ zgi_r[ ir ] ] < z0g ) {
					zgi_r[ ir ]++;
				}
			} else {
				zgi_r[ ir ] = 0.0;
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
		std::cout << "Stopped at range " << r[ NR-1 ]/1000.0 << " km" << std::endl;

		atm_profile_2d->remove_property( "_CEFF_" );
		atm_profile_2d->remove_property( "_WC_" );

		if (multiprop) {
			std::cout << "Writing 1-D output to tloss_multiprop.pe" << std::endl;
			output1DTL( "tloss_multiprop.pe", true );
		} else { 
			std::cout << "Writing 1-D output to tloss_1d.pe" << std::endl;
			output1DTL( "tloss_1d.pe" );
			if (write2d) {
				std::cout << "Writing 2-D output to tloss_2d.pe" << std::endl;
				output2DTL( "tloss_2d.pe" );
			}
		}
		std::cout << std::endl;

		for (i = 0; i < npade+1; i++) {
			ierr = MatDestroy( qpowers + i );CHKERRQ(ierr);
			if (starter == "self") {
				ierr = MatDestroy( qpowers_starter + i );CHKERRQ(ierr);
			}
		}
		delete [] qpowers;
		if (starter == "self") {
			delete [] qpowers_starter;
		}
	}

	//ierr = MatDestroy( &q );       CHKERRQ(ierr);
	ierr = MatDestroy( &B );       CHKERRQ(ierr);
	ierr = MatDestroy( &C );       CHKERRQ(ierr);
	ierr = VecDestroy( &psi_o );   CHKERRQ(ierr);
	ierr = VecDestroy( &Bpsi_o );  CHKERRQ(ierr);
	ierr = KSPDestroy( &ksp );     CHKERRQ(ierr);
	// for (i = 0; i < npade+1; i++) {
	// 	ierr = MatDestroy( qpowers + i );CHKERRQ(ierr);
	// 	ierr = MatDestroy( qpowers_starter + i );CHKERRQ(ierr);
	// }
	
	delete [] k;
	delete [] n;
	delete [] c;
	delete [] a_t;
	delete [] contents;
	delete [] indices;

	return 1;
}

void NCPA::EPadeSolver::calculate_atmosphere_parameters( NCPA::Atmosphere2D *atm, int NZvec, double *z_vec, 
	double r, double z_g, bool use_lossless, bool use_top_layer, double freq, bool absolute, 
	double &k0, double &c0, double *c_vec, double *a_vec, std::complex<double> *k_vec, 
	std::complex<double> *n_vec ) {

	std::complex<double> I( 0.0, 1.0 );

	std::memset( c_vec, 0, NZvec * sizeof(double) );
	std::memset( a_vec, 0, NZvec * sizeof(double) );
	std::memset( k_vec, 0, NZvec * sizeof( std::complex< double > ) );
	std::memset( n_vec, 0, NZvec * sizeof( std::complex< double > ) );

	// z_vec is relative to ground
	if (absolute) {
		fill_atm_vector_absolute( atm, r, NZvec, z_vec, "_CEFF_", c_underground, c_vec );
		//int zeroind = find_closest_index( z_vec, NZvec, 0.0 );
		//c0 = NCPA::mean( c_vec+zeroind, NZvec-zeroind );
	} else {
		fill_atm_vector_relative( atm, r, NZvec, z_vec, "_CEFF_", z_g, c_vec );
		//c0 = NCPA::mean( c_vec, NZvec );
		//c0 = c_vec[ 0 ];
	}
	c0 = atm->get( r, "_CEFF_", z_g );

	if (!use_lossless) {
		if (absolute) {
			fill_atm_vector_absolute( atm, r, NZvec, z_vec, "_ALPHA_", 0.0, a_vec );
		} else {
			fill_atm_vector_relative( atm, r, NZvec, z_vec, "_ALPHA_", z_g, a_vec );
		}
	}
	double *abslayer = new double[ NZvec ];
	memset( abslayer, 0, NZvec * sizeof(double) );
	if (use_top_layer) {
		absorption_layer( c0 / freq, z_vec, NZvec, abslayer );
	}
	
	k0 = 2.0 * PI * freq / c0;
	
	// Set up vectors
	//indices = new PetscInt[ NZ ];
	for (int i = 0; i < NZ; i++) {
		// if (absolute) {
		// 	if (z_vec[i] < z_g) {
		// 		k_vec[ i ] = 0.0;
		// 	} else {
		// 		k_vec[ i ] = 2.0 * PI * freq / c_vec[ i ] + (a_vec[ i ] + abslayer[ i ]) * I;
		// 	}
		// } else {
			k_vec[ i ] = 2.0 * PI * freq / c_vec[ i ] + (a_vec[ i ] + abslayer[ i ]) * I;
		// }
		n_vec[ i ] = k_vec[ i ] / k0;
	}
}

void NCPA::EPadeSolver::fill_atm_vector_relative( NCPA::Atmosphere2D *atm, double range, int NZvec, double *zvec, 
	std::string key, double groundheight, double *vec ) {

	for (int i = 0; i < NZvec; i++) {
		vec[i] = atm->get( range, key, zvec[i] + groundheight );
	}
}

void NCPA::EPadeSolver::fill_atm_vector_absolute( NCPA::Atmosphere2D *atm, double range, int NZvec, double *zvec, 
	std::string key, double fill_value, double *vec ) {

	double zmin = atm->get( range, "Z0" );
	// double bound_val = atm->get( range, key, zmin );
	// double tempval;
	for (int i = 0; i < NZvec; i++) {
		// if (zvec[i] < (zmin - 500.0)) {
		if (zvec[i] < zmin) {
			vec[i] = fill_value;
		// } else if (zvec[i] < zmin) {
		// 	double factor = 0.5 - 0.5 * std::cos( PI * (zmin - zvec[i]) / 500.0  );
		// 	tempval = (fill_value - bound_val) * factor + bound_val;
		// 	vec[i] = (fill_value - bound_val) * factor + bound_val;
		} else {
			vec[i] = atm->get( range, key, zvec[i] );
		}
	}
}



int NCPA::EPadeSolver::make_B_and_C_matrices( Mat *qpowers, int npade, int NZ, 
	std::vector< std::complex< double > > &P, std::vector< std::complex< double > > &Q,
	Mat *B, Mat *C ) {

	PetscErrorCode ierr;
	PetscInt Istart, Iend, i;
	PetscScalar value;

	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade-1, NULL, B );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *B );CHKERRQ(ierr);
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZ, NZ, 2*npade+1, NULL, C );CHKERRQ(ierr);
	ierr = MatSetFromOptions( *C );CHKERRQ(ierr);

	ierr = MatGetOwnershipRange(*B,&Istart,&Iend);CHKERRQ(ierr);
	value = 1.0;
	for (i = Istart; i < Iend; i++) {
		ierr = MatSetValues( *B, 1, &i, 1, &i, &value, INSERT_VALUES );CHKERRQ(ierr);
	}
	ierr = MatGetOwnershipRange( *C, &Istart, &Iend );CHKERRQ(ierr);
	for (i = Istart; i < Iend; i++) {
		ierr = MatSetValues( *C, 1, &i, 1, &i, &value, INSERT_VALUES );CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(*C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*C,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	for (i = 1; i < Q.size(); i++) {
		ierr = MatAXPY( *C, Q[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	for (i = 1; i < P.size(); i++) {
		ierr = MatAXPY( *B, P[ i ], qpowers[ i-1 ], DIFFERENT_NONZERO_PATTERN );CHKERRQ(ierr);
	}
	return 1;
}

int NCPA::EPadeSolver::make_q_powers( int NZvec, double *zvec, double k0, double h2, 
	std::complex<double> *n, size_t nqp, int boundary_index, Mat *qpowers ) {

	Mat q;
	PetscInt Istart, Iend, col[3];
	PetscBool FirstBlock = PETSC_FALSE, LastBlock = PETSC_FALSE;
	PetscErrorCode ierr;
	PetscScalar value[3];
	PetscInt i;

	// Set up matrices
	ierr = MatCreateSeqAIJ( PETSC_COMM_SELF, NZvec, NZvec, 3, NULL, &q );CHKERRQ(ierr);
	ierr = MatSetFromOptions( q );CHKERRQ(ierr);

	// populate
	double bnd_cnd = -1.0 / h2;    // @todo add hook for alternate boundary conditions
	//double bnd_cnd = -2.0 / h2;      // pressure release condition
	double k02 = k0*k0;
	
	ierr = MatGetOwnershipRange(q,&Istart,&Iend);CHKERRQ(ierr);
	if (Istart==0) FirstBlock=PETSC_TRUE;
    if (Iend==NZ) LastBlock=PETSC_TRUE;
    value[0]=1.0 / h2 / k02; value[2]=1.0 / h2 / k02;
    for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {
    		if (i < boundary_index)  {
    			value[ 0 ] = 0.0;
    			value[ 1 ] = 0.0;  
    			value[ 2 ] = 0.0;
    		} else if (i == boundary_index) {
    			value[ 0 ] = 0.0;
    			value[ 1 ] = bnd_cnd/k02 + (n[i]*n[i] - 1);
    			value[ 2 ] = 1.0 / h2 / k02;
    		} else {
    			value[ 0 ] = 1.0 / h2 / k02;
    			value[ 1 ] = -2.0/h2/k02 + (n[i]*n[i] - 1);
    			value[ 2 ] = 1.0 / h2 / k02;
    		}
		    col[0]=i-1; col[1]=i; col[2]=i+1;
		    ierr = MatSetValues(q,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (LastBlock) {
		    i=NZ-1; col[0]=NZ-2; col[1]=NZ-1;
		    value[ 0 ] = 1.0 / h2 / k02;
		    value[ 1 ] = -2.0/h2/k02 + (n[i]*n[i] - 1);
		    ierr = MatSetValues(q,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (FirstBlock) {
		    i=0; col[0]=0; col[1]=1; 
		    if (i < boundary_index)  {
    			value[ 0 ] = 0.0;
    			value[ 1 ] = 0.0;
    		} else {
    			value[ 0 ] = bnd_cnd/k02 + (n[i]*n[i] - 1);
    			value[ 1 ] = 1.0 / h2 / k02;
    		}
		    //value[0]=bnd_cnd/k02 + (n[i]*n[i] - 1); 
		    //value[1]=1.0 / h2 / k02;
		    ierr = MatSetValues(q,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(q,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    // calculate powers of q
	//qpowers = new Mat[ nqp ];
	ierr = MatConvert( q, MATSAME, MAT_INITIAL_MATRIX, qpowers );CHKERRQ(ierr);
	for (i = 1; i < nqp; i++) {
		ierr = MatMatMult( qpowers[i-1], qpowers[0], MAT_INITIAL_MATRIX, PETSC_DEFAULT, 
			qpowers+i );CHKERRQ(ierr);
	}
	ierr = MatDestroy( &q );CHKERRQ(ierr);
	return 1;
}

void NCPA::EPadeSolver::absorption_layer( double lambda, double *z, int NZ, double *layer ) {
	double z_t = z[NZ-1] - lambda;
	for (int i = 0; i < NZ; i++) {
		layer[ i ] = absorption_layer_mu * std::exp( (z[i]-z_t) / lambda );
	}
}


int NCPA::EPadeSolver::get_starter_gaussian( size_t NZ, double *z, double zs, double k0, int ground_index,
	Vec *psi ) {

	double fac = 2.0;
	//double kfac = k0 / fac;
	PetscScalar tempval;
	PetscErrorCode ierr;

	ierr = VecCreate( PETSC_COMM_SELF, psi );CHKERRQ(ierr);
	ierr = VecSetSizes( *psi, PETSC_DECIDE, NZ );CHKERRQ(ierr);
	ierr = VecSetFromOptions( *psi ); CHKERRQ(ierr);
	ierr = VecSet( *psi, 0.0 );

	for (PetscInt i = 0; i < NZ; i++) {
		//if (z[i] >= zg) {
			tempval = -( k0*k0/fac/fac ) * (z[i] - zs) * (z[i] - zs);
			tempval = sqrt( k0/fac ) * exp( tempval );
			ierr = VecSetValues( *psi, 1, &i, &tempval, INSERT_VALUES );CHKERRQ(ierr);
		//}
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
	PetscScalar I( 0.0, 1.0 ), tempsc, zeroval = 0.0;
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
	P->clear();
	Q->clear();

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

void NCPA::EPadeSolver::output1DTL( std::string filename, bool append ) {
	std::ofstream out_1d;
	if (append) {
		out_1d.open( filename, std::ofstream::out | std::ofstream::app );
		out_1d << std::endl;
	} else {
		out_1d.open( filename, std::ofstream::out | std::ofstream::trunc );
	}
	for (int i = 0; i < (NR-1); i++) {
		//out_1d << calc_az << " " << r[ i ]/1000.0 << " " << tl[ zgi_r[ i ] ][ i ] << std::endl;
		out_1d << calc_az << " " << r[ i ]/1000.0 << " " << tl[ zgi_r[ i ] ][ i ].real()
		       << " " << tl[ zgi_r[ i ] ][ i ].imag() << std::endl;
	}
	out_1d.close();
}

void NCPA::EPadeSolver::output2DTL( std::string filename ) {
	std::ofstream out_2d( filename, std::ofstream::out | std::ofstream::trunc );
	int zplot_int = 10;
	for (int i = 0; i < (NR-1); i++) {
		for (int j = 0; j < NZ; j += zplot_int) {
			//out_2d << r[ i ]/1000.0 << " " << z[ j ]/1000.0 << " " << tl[ j ][ i ] << " 0.0" << std::endl;
			out_2d << r[ i ]/1000.0 << " " << z[ j ]/1000.0 << " " << tl[ j ][ i ].real() 
			<< " " << tl[ j ][ i ].imag() << std::endl;
		}
		out_2d << std::endl;
	}
	out_2d.close();
}