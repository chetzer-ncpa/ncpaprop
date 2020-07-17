#include "ProfileSeriesAtmosphere2D.h"
#include "Atmosphere1D.h"
#include "units.h"
#include <fstream>


NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D() : Atmosphere2D() { }

NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D( const std::string &filename ) : Atmosphere2D() {

	std::ifstream infile( filename );
	double range;
	std::string atmfile;
	infile >> range >> atmfile;
	Atmosphere1D *tempatm;
	set_insert_range_units( NCPA::Units::fromString( "km" ) );
	while (infile.good()) {
		tempatm = new Atmosphere1D( atmfile );
		insert_profile( tempatm, range );
		infile >> range >> atmfile;
	}
	infile.close();
	sort_profiles();
}

NCPA::ProfileSeriesAtmosphere2D::~ProfileSeriesAtmosphere2D() { }

//NCPA::ProfileSeriesAtmosphere2D::ProfileSeriesAtmosphere2D( const ProfileSeriesAtmosphere2D &atm ) : Atmosphere2D( atm ) { }