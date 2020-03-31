/* Compile with

g++ -o atmosphere1d_test -ggdb -I.. -I../../atmosphere ../../atmosphere/Atmosphere1D.cpp ../../atmosphere/Atmosphere.cpp ../util.cpp ../units.cpp atmosphere1d_test.cpp

*/
#include "Atmosphere1D.h"
#include "util.h"
#include "units.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>


using namespace std;
using namespace NCPA;


int main( int argc, char **argv ) {

	string filename, line;
	ostringstream oss;  // for exceptions
	unsigned int i;

	if (argc > 1) {
		filename = argv[ 1 ];
	} else {
		cerr << "No input file specified!" << endl;
		return 0;
	}

	vector< string > atmlines, headerlines;

	// read in the header.  Descriptive lines start with "#%"
	// descriptive line format is column number,key,units string
	// column 0 is assumed to be altitude, so the key string is ignored, but
	// the other fields are still required
	ifstream infile( filename.c_str() );
	getline( infile, line );

	while ( infile.good() ) {
		// lines will either be comments (# ), field descriptions (#% ), or
		// field contents
		line = deblank( line );
		if (line[ 0 ] == '#') {
			// check second character
			if (line.size() > 1 && line[ 1 ] == '%') {
				headerlines.push_back( line.substr( 2 ) );
			} // otherwise it's a regular comment and can be ignored

		} else if (line.size() == 0) {
			// skip empty lines
		} else {
			atmlines.push_back( line );
		}

		getline( infile, line );
	}
	infile.close();
	cout << "Found " << headerlines.size() << " header lines" << endl;
	cout << "Found " << atmlines.size() << " data lines" << endl;

	// parse them out
	size_t nfields = headerlines.size();
	if (nfields == 0) {
		cerr << "No descriptive fields found" << endl;
		return 0;
	}
	vector< string > keys( nfields ), units( nfields ), fields;
	for (vector< string >::const_iterator it = headerlines.begin(); it != headerlines.end(); ++it) {
		fields = split( *it, "," );
		if ( fields.size() != 3 ) {
			oss << "Error parsing descriptive line:" << endl << line << endl 
				<< "Can't split into three comma-separated fields" << endl;
			throw invalid_argument( oss.str() );
		}

		// process fields
		// field line needs to be parseable as an integer
		unsigned int col;
		try {
			col = (unsigned int)stoi( fields[ 0 ] );
		} catch ( invalid_argument &e ) {
			oss << "Error parsing descriptive line:" << endl << line << endl 
				<< "First field not parseable as a positive integer" << endl;
			throw invalid_argument( oss.str() );
		}

		if (keys.size() < col) {
			keys.reserve( 2 * col );
			units.reserve( 2 * col );
		}
		
		cout << "keys[ " << col-1 << " ] = " << fields[ 1 ] << endl;
		keys[ col-1 ] = deblank(fields[ 1 ]);
		cout << "units[ " << col-1 << " ] = " << fields[ 2 ] << endl;
		units[ col-1 ] = deblank(fields[ 2 ]);
	}

	int nlines = atmlines.size();
	int ncols = keys.size();

	vector< double * > columns( ncols );
	for ( i = 0; i < ncols; i++ ) {
		double *col = new double[ nlines ];
		columns[ i ] = col;
	}


	// step through the descriptive lines
	size_t row = 0;
	for (vector< string >::const_iterator it = atmlines.cbegin(); it != atmlines.cend(); ++it ) {
		fields.clear();
		fields = split( *it, " \t," );
		if (fields.size() != ncols ) {
			oss << "Error parsing data line:" << endl << *it << endl
				<< ncols << " columns expected, " << fields.size() << " columns found.";
			for ( i = 0; i < ncols; i++ ) {
				delete [] columns[ i ];
			}
			throw invalid_argument( oss.str() );
		}
		for ( i = 0; i < ncols; i++ ) {
			try {
				double *thiscol = columns[ i ];
				thiscol[ row ] = stof( fields[ i ] );
			} catch (invalid_argument &e) {
				oss << "Error parsing data line:" << endl << *it << endl
					<< "Can't parse field " << fields[ i ] << " as a double";
				for ( i = 0; i < ncols; i++ ) {
					delete [] columns[ i ];
				}
				throw invalid_argument( oss.str() );
			}
		}
		row++;
	}

	Atmosphere1D *atm = new Atmosphere1D( nlines, columns[ 0 ], Units::fromString( units[ 0 ] ) );
	for (i = 1; i < ncols; i++) {
		atm->add_property( keys[ i ], nlines, columns[ i ], Units::fromString( units[ i ] ) );
	}
	for ( i = 0; i < ncols; i++ ) {
		delete [] columns[ i ];
	}

	// start playing with atmosphere
	atm->calculate_sound_speed_from_temperature( "C_T", "T" );
	atm->calculate_sound_speed_from_pressure_and_density( "C_P", "P", "RHO" );
	size_t nz = atm->get_basis_length();
	double *z = new double[ nz ], *u = new double[ nz ], *t = new double[ nz ], *c_t = new double[ nz ], *c_p = new double[ nz ];
	units_t z_units, u_units, t_units, c_t_units, c_p_units;
	atm->convert_altitude_units( Units::fromString( "m" ) );
	atm->get_altitude_basis( z, &z_units );
	atm->get_property_basis( "U", u, &u_units );
	atm->get_property_basis( "T", t, &t_units );
	atm->get_property_basis( "C_T", c_t, &c_t_units );
	atm->get_property_basis( "C_P", c_p, &c_p_units );

	// convert z to meters
	//Units::convert( z, nz, z_units, Units::fromString("m"), z );

	ofstream os( "atm_out_m.dat" );
	for (i = 0; i < nz; i += 5) {
		os << z[ i ] << "  " << t[ i ] << "  " << u[ i ] << "  " << c_t[ i ] << "  " << c_p[ i ] << endl;
	}
	os.close();

	atm->revert_altitude_units();
	atm->get_altitude_basis( z, &z_units );
	os = ofstream( "atm_out_km.dat" );
	for (i = 0; i < nz; i += 5) {
		os << z[ i ] << "  " << t[ i ] << "  " << u[ i ] << "  " << c_t[ i ] << "  " << c_p[ i ] << endl;
	}
	os.close();


	delete [] z, u, t;
	delete atm;

	return 1;
}