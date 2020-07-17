#include "StratifiedAtmosphere2D.h"
#include "Atmosphere1D.h"
#include "units.h"


NCPA::StratifiedAtmosphere2D::StratifiedAtmosphere2D( const Atmosphere1D *atm ) : Atmosphere2D() {
	Atmosphere1D *tempatm = new Atmosphere1D( *atm );
	set_insert_range_units( NCPA::Units::fromString( "km" ) );
	insert_profile( tempatm, 0.0 );
	sorted_ = true;
}

NCPA::StratifiedAtmosphere2D::StratifiedAtmosphere2D( const std::string &filename ) : Atmosphere2D() {
	Atmosphere1D *tempatm = new Atmosphere1D( filename );
	set_insert_range_units( NCPA::Units::fromString( "km" ) );
	insert_profile( tempatm, 0.0 );
	sorted_ = true;
}

NCPA::StratifiedAtmosphere2D::~StratifiedAtmosphere2D() { }

//NCPA::StratifiedAtmosphere2D::StratifiedAtmosphere2D( const StratifiedAtmosphere2D &atm ) : Atmosphere2D( atm ) { }