/* Compile with

g++ -o atmosphere1d_test -ggdb -I.. -I../../atmosphere ../../atmosphere/Atmosphere1D.cpp ../../atmosphere/AtmosphericProperty1D.cpp ../util.cpp ../units.cpp atmosphere1d_test.cpp -lgsl -lgslcblas

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
#include <algorithm>


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

	vector< string > atmlines, headerlines, scalarlines;

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

	// hold contents
	vector< string > keys, fields;
	vector< unsigned int > column_numbers;
	vector< double > values;
	vector< units_t > units;

	for (vector< string >::const_iterator it = headerlines.begin(); it != headerlines.end(); ++it) {
		fields = split( *it, "," );
		if ( fields.size() != 3 && fields.size() != 4 ) {
			oss << "Error parsing descriptive line:" << endl << line << endl 
				<< "Must be formatted as:" << endl
				<< "column,key,units[,value]" << endl
				<< "Use column=0 and specify value for scalar quantities." << endl;
			throw invalid_argument( oss.str() );
		}

		// process fields
		// field line needs to be parseable as an integer
		unsigned int col;
		try {
			col = (unsigned int)stoi( fields[ 0 ] );
		} catch ( invalid_argument &e ) {
			oss << "Error parsing descriptive line:" << endl << line << endl 
				<< "First field not parseable as an integer" << endl;
			throw invalid_argument( oss.str() );
		}

/*
		if (keys.size() < col) {
			keys.reserve( 2 * col );
			units.reserve( 2 * col );
			column_numbers.reserve( 2 * col );
			values.reserve( 2 * col );
		}
*/

		double tempval = 0.0;
		units_t tempunits = UNITS_NONE;
		try {
			tempunits = Units::fromString( deblank(fields[ 2 ]) );
		} catch (out_of_range& oor) {
			throw invalid_argument( oor.what() );
		}
		if (fields.size() == 4) {
			tempval = stof( deblank(fields[ 3 ]) );
		}

		// add to header vectors
		column_numbers.push_back( col );
		keys.push_back( deblank(fields[ 1 ]) );
		units.push_back( tempunits );
		values.push_back( tempval );  // this will be ignored for vector quantities

		cout << "Column " << col << ": " << deblank(fields[1]) << " [ " << Units::toStr( tempunits ) << " ] = " << tempval << endl;
	}

	// check for uniqueness in keys
	vector< string > tempkeys( keys );
	sort( tempkeys.begin(), tempkeys.end() );
	vector< string >::iterator uit = unique( tempkeys.begin(), tempkeys.end() );
	if (uit != tempkeys.end()) {
		cout << "Keys are not unique" << endl;
		exit( 1 );
	}


	int nlines = atmlines.size();
	int ncols = 0;
	units_t depunits = UNITS_NONE;
	for (i = 0; i < column_numbers.size(); i++) {
		ncols = column_numbers[ i ] > ncols ? column_numbers[ i ] : ncols;
		if (column_numbers[ i ] == 1) {
			depunits = units[ i ];
		}
	}

	vector< double * > columns( ncols );
	for ( i = 0; i < ncols; i++ ) {
		double *col = new double[ nlines ];
		columns[ i ] = col;
	}


	// step through the data lines
	size_t row = 0;
	vector< string >::const_iterator it;
	for (it = atmlines.cbegin(); it != atmlines.cend(); ++it ) {
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

	Atmosphere1D *atm = new Atmosphere1D( nlines, columns[ 0 ], depunits );
	for (i = 0; i < keys.size(); i++) {
		if (column_numbers[ i ] == 0) {
			atm->add_property( keys[ i ], values[ i ], units[ i ] );
		} else {
			atm->add_property( keys[ i ], nlines, columns[ column_numbers[ i ] - 1 ], units[ i ] );
		}
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
	atm->get_altitude_vector( z, &z_units );
	atm->get_property_vector( "U", u, &u_units );
	atm->get_property_vector( "T", t, &t_units );
	atm->get_property_vector( "C_T", c_t, &c_t_units );
	atm->get_property_vector( "C_P", c_p, &c_p_units );

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
	for (i = 0; i < nz; i += 5) {
		os << z[ i ] << "  " << t[ i ] << "  " << u[ i ] << "  " << c_t[ i ] << "  " << c_p[ i ] << endl;
	}
	os.close();

	atm->revert_altitude_units();
	atm->get_altitude_vector( z, &z_units );
	os = ofstream( "atm_out_km.dat" );
	for (i = 0; i < nz; i += 5) {
		os << z[ i ] << "  " << t[ i ] << "  " << u[ i ] << "  " << c_t[ i ] << "  " << c_p[ i ] << endl;
	}
	os.close();


	delete [] z, u, t;
	delete atm;

	return 1;
}