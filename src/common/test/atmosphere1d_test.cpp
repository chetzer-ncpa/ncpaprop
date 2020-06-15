/* Compile with

g++ -o atmosphere1d_test -ggdb -I.. -I../../atmosphere ../../atmosphere/Atmosphere1D.cpp ../../atmosphere/AtmosphericProperty1D.cpp ../util.cpp ../units.cpp atmosphere1d_test.cpp -lgsl -lgslcblas

*/
#include "Atmosphere1D.h"
#include "AtmosphericProperty1D.h"
#include "util.h"
#include "units.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>


using namespace std;
using namespace NCPA;


int main( int argc, char **argv ) {

	string filename;

	if (argc > 1) {
		filename = argv[ 1 ];
	} else {
		cerr << "No input file specified!" << endl;
		return 0;
	}

	double dummy_z[5] = { 0, 10, 20, 30, 40 };
	double dummy_t[5] = { 5, 5, 3, 5, 3 };
	AtmosphericProperty1D *T = new AtmosphericProperty1D( 5, dummy_z, UNITS_DISTANCE_METERS, dummy_t, UNITS_TEMPERATURE_CELSIUS );
	AtmosphericProperty1D Tcopy( *T );
	exit( 0 );
	
	//ifstream ifs( filename );
	Atmosphere1D *atm = new Atmosphere1D( filename );

	// start playing with atmosphere
	atm->calculate_sound_speed_from_temperature( "C_T", "T", Units::fromString( "m/s" ) );
	atm->calculate_sound_speed_from_pressure_and_density( "C_P", "P", "RHO", Units::fromString( "m/s" ) );

	// get wind vector quantities from components
	atm->calculate_wind_speed( "WSPD", "U", "V" );
	atm->calculate_wind_direction( "WDIR", "U", "V" );

	// calculate wind component and effective sound speed in four directions
	atm->calculate_wind_component( "WC_90", "WSPD", "WDIR", 90.0 );
	atm->calculate_wind_component( "WC_45", "WSPD", "WDIR", 45.0 );
	atm->calculate_wind_component( "WC_00", "WSPD", "WDIR", 0.0 );
	atm->calculate_wind_component( "WC_270", "WSPD", "WDIR", -90.0 );

	// calculate effective sound speeds in four directions
	atm->calculate_effective_sound_speed( "C_90", "C_P", "WC_90" );
	atm->calculate_effective_sound_speed( "C_45", "C_P", "WC_45" );
	atm->calculate_effective_sound_speed( "C_00", "C_P", "WC_00" );
	atm->calculate_effective_sound_speed( "C_270", "C_P", "WC_270" );


	size_t nz = atm->get_basis_length();
	double *z = new double[ nz ], *u = new double[ nz ], *t = new double[ nz ], *c_t = new double[ nz ], *c_p = new double[ nz ];
	units_t z_units, u_units, t_units, c_t_units, c_p_units;
	atm->convert_altitude_units( Units::fromString( "m" ) );
	atm->get_altitude_vector( z, &z_units );
	atm->get_property_vector( "U", u, &u_units );
	atm->get_property_vector( "T", t, &t_units );
	atm->get_property_vector( "C_T", c_t, &c_t_units );
	atm->get_property_vector( "C_P", c_p, &c_p_units );

	cout << "C_P has units " << Units::toString( c_p_units ) << endl;

	// check keys
	if (atm->contains_vector( "U" )) {
		cout << "U found (this is good)" << endl;
	} else {
		cout << "U not found (this is bad)" << endl;
	}
	if (atm->contains_vector( "X" )) {
		cout << "X found (this is bad)" << endl;
	} else {
		cout << "X not found (this is good)" << endl;
	}
	if (atm->contains_scalar( "Z0" )) {
		cout << "Z0 found (this is good)" << endl;
	} else {
		cout << "Z0 not found (this is bad)" << endl;
	}
	if (atm->contains_scalar( "Z1" )) {
		cout << "Z1 found (this is bad)" << endl;
	} else {
		cout << "Z1 not found (this is good)" << endl;
	}


	// convert z to meters
	//Units::convert( z, nz, z_units, Units::fromString("m"), z );
	ofstream os( "atm_out_m.dat" );
	atm->print_atmosphere( "Z", os );
	//unsigned int i;
	//for (i = 0; i < nz; i += 5) {
	//	os << z[ i ] << "  " << t[ i ] << "  " << u[ i ] << "  " << c_t[ i ] << "  " << c_p[ i ] << endl;
	//}
	os.close();

	atm->convert_altitude_units( Units::fromString( "km" ) );
	//atm->get_altitude_vector( z, &z_units );
	os = ofstream( "atm_out_km.dat" );
	atm->print_atmosphere( "Z", os );
	//for (i = 0; i < nz; i += 5) {
	//	os << z[ i ] << "  " << t[ i ] << "  " << u[ i ] << "  " << c_t[ i ] << "  " << c_p[ i ] << endl;
	//}
	os.close();

	atm->resample( 0.5 );
	os = ofstream( "atm_out_resampled.dat" );
	atm->print_atmosphere( "Z", os );
	os.close();

	delete [] z, u, t;
	delete atm;

	return 1;
}