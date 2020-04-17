#include "Atmosphere.h"
#include "Atmosphere1D.h"
#include "AtmosphericProperty1D.h"
#include "units.h"
#include "util.h"
#include <map>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#define GAMMA_FOR_C 1.4
#define R_FOR_C 287.0

#ifndef PI
#define PI 3.14159
#endif


NCPA::Atmosphere1D::Atmosphere1D( size_t n_altitude_points, double *altitude_points, units_t altitude_units ) {
	contents_.clear();
	nz_ = n_altitude_points;
	z_ = new double[ nz_ ];
	std::memcpy( z_, altitude_points, nz_ * sizeof(double) );
	z_units_.push( altitude_units );
}

NCPA::Atmosphere1D::Atmosphere1D( std::istream& in ) {

	if ( ! in.good() ) {
		throw std::runtime_error( "Atmosphere1D - Input stream not in good state" );
	}

	std::string line;
	std::vector< std::string > atmlines, headerlines, scalarlines;
	std::ostringstream oss;     // for exceptions
	size_t i;                   // repeated index variable

	std::getline( in, line );
	while ( in.good() ) {
		// lines will either be comments (# ), field descriptions (#% ), or
		// field contents
		line = NCPA::deblank( line );
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

		getline( in, line );
	}
	//in.close();
	//cout << "Found " << headerlines.size() << " header lines" << endl;
	//cout << "Found " << atmlines.size() << " data lines" << endl;

	// parse them out
	size_t nfields = headerlines.size();
	if (nfields == 0) {
		throw std::runtime_error( "Atmosphere1D - No descriptive fields found." );
	}

	// hold contents
	std::vector< std::string > keys, fields;
	std::vector< unsigned int > column_numbers;
	std::vector< double > values;
	std::vector< NCPA::units_t > units;

	for (std::vector< std::string >::const_iterator it = headerlines.begin(); it != headerlines.end(); ++it) {
		fields = NCPA::split( NCPA::deblank( *it ), " ," );
		if ( fields.size() != 3 && fields.size() != 4 ) {
			oss << "Atmosphere1D - Error parsing descriptive line:" << std::endl << line << std::endl 
				<< "Must be formatted as:" << std::endl
				<< "column,key,units[,value]" << std::endl
				<< "Use column=0 and specify value for scalar quantities." << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		// process fields
		// field line needs to be parseable as an integer
		unsigned int col;
		try {
			col = (unsigned int)std::stoi( fields[ 0 ] );
		} catch ( std::invalid_argument &e ) {
			oss << "Atmosphere1D - Error parsing descriptive line:" << std::endl << line << std::endl 
				<< "First field not parseable as an integer" << std::endl;
			throw std::invalid_argument( oss.str() );
		}

		double tempval = 0.0;
		NCPA::units_t tempunits = NCPA::UNITS_NONE;
		try {
			tempunits = NCPA::Units::fromString( NCPA::deblank( fields[ 2 ] ) );
		} catch (std::out_of_range& oor) {
			throw std::invalid_argument( oor.what() );
		}
		if (fields.size() == 4) {
			tempval = std::stof( deblank(fields[ 3 ]) );
		}

		// add to header vectors
		column_numbers.push_back( col );
		keys.push_back( deblank( fields[ 1 ] ) );
		units.push_back( tempunits );
		values.push_back( tempval );  // this will be ignored for vector quantities

		//cout << "Column " << col << ": " << deblank(fields[1]) << " [ " << Units::toStr( tempunits ) << " ] = " << tempval << endl;
	}

	// check for uniqueness in keys
	std::vector< std::string > tempkeys( keys );
	std::sort( tempkeys.begin(), tempkeys.end() );
	std::vector< std::string >::iterator uit = std::unique( tempkeys.begin(), tempkeys.end() );
	if (uit != tempkeys.end()) {
		throw std::invalid_argument( "Atmosphere1D: key values are not unique" );
	}


	int nlines = atmlines.size();
	int ncols = 0;
	NCPA::units_t depunits = NCPA::UNITS_NONE;
	for (i = 0; i < column_numbers.size(); i++) {
		ncols = column_numbers[ i ] > ncols ? column_numbers[ i ] : ncols;
		if (column_numbers[ i ] == 1) {
			depunits = units[ i ];
		}
	}

	std::vector< double * > columns( ncols );
	for ( i = 0; i < ncols; i++ ) {
		double *col = new double[ nlines ];
		columns[ i ] = col;
	}


	// step through the data lines
	size_t row = 0;
	std::vector< std::string >::const_iterator it;
	for (it = atmlines.cbegin(); it != atmlines.cend(); ++it ) {
		fields.clear();
		fields = NCPA::split( NCPA::deblank( *it ), " \t," );
		if (fields.size() != ncols ) {
			oss << "Error parsing data line:" << std::endl << *it << std::endl
				<< ncols << " columns expected, " << fields.size() << " columns found.";
			for ( i = 0; i < ncols; i++ ) {
				delete [] columns[ i ];
			}
			throw std::invalid_argument( oss.str() );
		}
		for ( i = 0; i < ncols; i++ ) {
			try {
				double *thiscol = columns[ i ];
				thiscol[ row ] = std::stof( fields[ i ] );
			} catch (std::invalid_argument &e) {
				oss << "Error parsing data line:" << std::endl << *it << std::endl
					<< "Can't parse field " << fields[ i ] << " as a double";
				for ( i = 0; i < ncols; i++ ) {
					delete [] columns[ i ];
				}
				throw std::invalid_argument( oss.str() );
			}
		}
		row++;
	}

	contents_.clear();
	nz_ = nlines;
	z_ = new double[ nz_ ];
	std::memcpy( z_, columns[ 0 ], nz_ * sizeof(double) );
	z_units_.push( depunits );

	for (i = 0; i < keys.size(); i++) {
		if (column_numbers[ i ] == 0) {
			this->add_property( keys[ i ], values[ i ], units[ i ] );
		} else if (column_numbers[ i ] > 1) {
			this->add_property( keys[ i ], nlines, columns[ column_numbers[ i ] - 1 ], units[ i ] );
		}
	}

	for ( i = 0; i < ncols; i++ ) {
		delete [] columns[ i ];
	}
}


NCPA::Atmosphere1D::~Atmosphere1D() {
	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it=contents_.begin(); it != contents_.end(); ++it) {
		delete it->second;
	}
	contents_.clear();
	while (! z_units_.empty()) {
		z_units_.pop();
	}
}

size_t NCPA::Atmosphere1D::get_basis_length() const {
	return nz_;
}

void NCPA::Atmosphere1D::calculate_wind_speed( std::string new_key, std::string we_wind_speed_key, 
	std::string sn_wind_speed_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *u_prop = contents_.at( we_wind_speed_key );
	NCPA::AtmosphericProperty1D *v_prop = contents_.at( sn_wind_speed_key );
	units_t u_units, v_units;
	double *u = new double[ nz_ ];
	double *v = new double[ nz_ ];
	u_prop->get_vector( u, &u_units );
	v_prop->get_vector( v, &v_units );
	if (u_units != v_units) {
		throw std::runtime_error( "Units mismatch between wind speed properties" );
	}

	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( u[ i ] * u[ i ] + v[ i ] * v[ i ] );
	}
	delete [] u;
	delete [] v;
	add_property( new_key, nz_, c, u_units );
	delete [] c;
}

void NCPA::Atmosphere1D::calculate_wind_direction( std::string new_key, std::string we_wind_speed_key, 
	std::string sn_wind_speed_key, units_t direction_units ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *u_prop = contents_.at( we_wind_speed_key );
	NCPA::AtmosphericProperty1D *v_prop = contents_.at( sn_wind_speed_key );
	units_t u_units, v_units;
	double *u = new double[ nz_ ];
	double *v = new double[ nz_ ];
	u_prop->get_vector( u, &u_units );
	v_prop->get_vector( v, &v_units );
	if (u_units != v_units) {
		throw std::runtime_error( "Units mismatch between wind speed properties" );
	}

	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = PI/2 - std::atan2( v[ i ], u[ i ] );
	}
	delete [] u;
	delete [] v;
	NCPA::Units::convert( c, nz_, NCPA::UNITS_ANGLE_RADIANS, NCPA::UNITS_ANGLE_DEGREES, c );
	NCPA::Units::convert( c, nz_, NCPA::UNITS_DIRECTION_DEGREES_CLOCKWISE_FROM_NORTH, direction_units, c );
	add_property( new_key, nz_, c, direction_units );
	delete [] c;
}

void NCPA::Atmosphere1D::calculate_wind_component( std::string new_key, std::string wind_speed_key, 
	std::string wind_direction_key, double azimuth ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If it's there we
		// erase it so we can create a new one
		c_prop = contents_.at( new_key );
		contents_.erase( new_key );
		delete c_prop;
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *sp_prop = contents_.at( wind_speed_key );
	NCPA::AtmosphericProperty1D *dir_prop = contents_.at( wind_direction_key );
	//NCPA::ScalarWithUnits *prop_dir = scalar_contents_.at( prop_direction_key );
	units_t sp_units, dir_units;
	double *sp = new double[ nz_ ];
	double *dir = new double[ nz_ ];
	sp_prop->get_vector( sp, &sp_units );
	dir_prop->get_vector( dir, &dir_units );
	//double propaz_rad = NCPA::Units::convert( prop_dir->get(), prop_dir->get_units(), NCPA::UNITS_ANGLE_RADIANS );

	// convert wind direction to radians for trig
	NCPA::Units::convert( dir, nz_, NCPA::UNITS_ANGLE_DEGREES, NCPA::UNITS_ANGLE_RADIANS, dir );
	double az_rad = NCPA::Units::convert( azimuth, NCPA::UNITS_ANGLE_DEGREES, NCPA::UNITS_ANGLE_RADIANS );

	double *wc = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		wc[ i ] = sp[ i ] * std::cos( dir[ i ] - az_rad );
	}
	delete [] sp;
	delete [] dir;
	add_property( new_key, nz_, wc, sp_units );
	delete [] wc;
}

void NCPA::Atmosphere1D::calculate_effective_sound_speed( std::string new_key, std::string sound_speed_key, 
	std::string wind_component_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If it's there we
		// erase it so we can create a new one
		c_prop = contents_.at( new_key );
		contents_.erase( new_key );
		delete c_prop;
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *c0_prop = contents_.at( sound_speed_key );
	NCPA::AtmosphericProperty1D *ws_prop = contents_.at( wind_component_key );
	units_t c0_units, ws_units;
	double *c0 = new double[ nz_ ];
	double *ws = new double[ nz_ ];
	c0_prop->get_vector( c0, &c0_units );
	ws_prop->get_vector( ws, &ws_units );

	// make sure units are consistent
	if (c0_units != ws_units) {
		NCPA::Units::convert( ws, nz_, ws_units, c0_units, ws );
		ws_units = c0_units;
	}

	double *ceff = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		ceff[ i ] = ws[ i ] + c0[ i ];
	}
	delete [] ws;
	delete [] c0;
	add_property( new_key, nz_, ceff, c0_units );
	delete [] ceff;
}

void NCPA::Atmosphere1D::get_altitude_vector( double *buffer, units_t *buffer_units ) const {
	*buffer_units = z_units_.top();
	std::memcpy( buffer, z_, nz_* sizeof( double ) );
}

void NCPA::Atmosphere1D::get_property_vector( std::string key, double *buffer, units_t *buffer_units ) const {
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );
	prop->get_vector( buffer, buffer_units );
}

NCPA::units_t NCPA::Atmosphere1D::get_property_units( std::string key ) {
	try {
		return contents_.at( key )->get_units();
	} catch (std::out_of_range& oor) {}
	try {
		return scalar_contents_.at( key )->get_units();
	} catch (std::out_of_range& oor) {
		throw std::out_of_range( "No vector or scalar quantity found with key " + key );
	}
}


double NCPA::Atmosphere1D::get_minimum_altitude() const {
	return z_[0];
	//return NCPA::Units::convert( z_[0], z_units_, altitude_units );
}

double NCPA::Atmosphere1D::get_maximum_altitude() const {
	return z_[ nz_ - 1 ];
	//return NCPA::Units::convert( z_[ nz_-1 ], z_units_, altitude_units );
}

void NCPA::Atmosphere1D::add_property( std::string key, size_t n_points, double *quantity_points, units_t quantity_units ) {

	// see if we already have one with that key
	NCPA::AtmosphericProperty1D *prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		prop = contents_.at( key );
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	if (n_points != nz_) {
		std::ostringstream oss;
		oss << "Array length " << n_points << " for key " << key << " does not match number of altitude points " << nz_;
		throw std::runtime_error( oss.str() );
	}

	prop = new AtmosphericProperty1D( n_points, z_, z_units_.top(), quantity_points, quantity_units );
	contents_[ key ] = prop;
}

void NCPA::Atmosphere1D::add_property( std::string key, double value, NCPA::units_t units ) {
	NCPA::ScalarWithUnits *scalar;
	try {
		scalar = scalar_contents_.at( key );
		throw std::runtime_error( "Requested key " + key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	scalar = new ScalarWithUnits( value, units );
	scalar_contents_[ key ] = scalar;
}

double NCPA::Atmosphere1D::get( std::string key, double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop;
	try {
		prop = contents_.at( key );
	} catch (std::out_of_range& oor) {
		throw std::out_of_range( "No vector quantity \"" + key + "\" found" );
	}
	return prop->get( altitude );
}

double NCPA::Atmosphere1D::get( std::string key ) const {
	NCPA::ScalarWithUnits *prop;
	try {
		prop = scalar_contents_.at( key );
	} catch (std::out_of_range& oor) {
		throw std::out_of_range( "No scalar quantity \"" + key + "\" found" );
	}
	return prop->get();
}

bool NCPA::Atmosphere1D::contains_vector( std::string key ) const {
	NCPA::AtmosphericProperty1D *prop;
	try {
		prop = contents_.at( key );
	} catch (std::out_of_range& oor) {
		return false;
	}
	return true;
}

bool NCPA::Atmosphere1D::contains_scalar( std::string key ) const {
	NCPA::ScalarWithUnits *prop;
	try {
		prop = scalar_contents_.at( key );
	} catch (std::out_of_range& oor) {
		return false;
	}
	return true;
}

double NCPA::Atmosphere1D::get_first_derivative( std::string key, double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_first_derivative( altitude );
}


double NCPA::Atmosphere1D::get_second_derivative( std::string key, double altitude ) const {
	
	// if it's not there it'll throw an out_of_range exception
	NCPA::AtmosphericProperty1D *prop = contents_.at( key );  
	return prop->get_second_derivative( altitude );
}


void NCPA::Atmosphere1D::calculate_sound_speed_from_temperature( std::string new_key, std::string temperature_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *t_prop = contents_.at( temperature_key );
	t_prop->convert_units( NCPA::UNITS_TEMPERATURE_KELVIN );
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( GAMMA_FOR_C * R_FOR_C * t_prop->get( z_[i] ) );
	}
	t_prop->revert_units();

	add_property( new_key, nz_, c, NCPA::UNITS_SPEED_METERS_PER_SECOND );
}

void NCPA::Atmosphere1D::calculate_sound_speed_from_pressure_and_density( std::string new_key, 
			std::string pressure_key, std::string density_key ) {

	NCPA::AtmosphericProperty1D *c_prop;
	try {
		// see if the requested key exists in the map.  If not, it throws
		// an out_of_range exception that we catch and ignore; if it is then
		// we throw an exception
		c_prop = contents_.at( new_key );
		throw std::runtime_error( "Requested key " + new_key + " already exists in atmosphere" );
	} catch (const std::out_of_range& oor) { }

	NCPA::AtmosphericProperty1D *p_prop = contents_.at( pressure_key );
	p_prop->convert_units( NCPA::UNITS_PRESSURE_PASCALS );
	NCPA::AtmosphericProperty1D *r_prop = contents_.at( density_key );
	r_prop->convert_units( NCPA::UNITS_DENSITY_KILOGRAMS_PER_CUBIC_METER );
	double *c = new double[ nz_ ];
	for (size_t i = 0; i < nz_; i++) {
		c[ i ] = std::sqrt( GAMMA_FOR_C * p_prop->get( z_[i] ) / r_prop->get( z_[i] ) );
	}

	p_prop->revert_units();
	r_prop->revert_units();
	add_property( new_key, nz_, c, NCPA::UNITS_SPEED_METERS_PER_SECOND );
}

void NCPA::Atmosphere1D::convert_altitude_units( units_t new_units ) {
	// will throw out_of_range and leave original units unchanged if there's an error
	// if there's no change in units, don't bother with the calculation, just push another
	// one onto the stack so reversion can happen properly
	if (new_units != z_units_.top()) {
		do_units_conversion_( nz_, z_, z_units_.top(), new_units );
	}
	z_units_.push( new_units );

	// update all contents
	for ( std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it = contents_.begin();
						it != contents_.end(); ++it ) {
		(*it).second->convert_altitude_units( new_units );
	}
}

void NCPA::Atmosphere1D::revert_altitude_units() {
	if (z_units_.size() < 2) {
		return;
	}

	NCPA::units_t current_units, last_units;
	current_units = z_units_.top();
	z_units_.pop();
	last_units = z_units_.top();
	if (current_units != last_units) {
		try {
			do_units_conversion_( nz_, z_, current_units, last_units );
		} catch (std::out_of_range &oor) {
			z_units_.push( current_units );
			throw;
		}
	}

	// update all contents
	for ( std::map< std::string, NCPA::AtmosphericProperty1D * >::iterator it = contents_.begin();
						it != contents_.end(); ++it ) {
		(*it).second->revert_altitude_units();
	}
}

void NCPA::Atmosphere1D::convert_property_units( std::string key, units_t new_units ) {
	// if it's not there it'll throw an out_of_range exception
	try {
		contents_.at( key )->convert_units( new_units );
		return;
	} catch (std::out_of_range& oor) {}
	try {
		scalar_contents_.at( key )->convert_units( new_units );
	} catch (std::out_of_range& oor) {
		throw std::out_of_range( "No vector or scalar quantity found with key " + key );
	}
}

void NCPA::Atmosphere1D::revert_property_units( std::string key ) {
	try {
		contents_.at( key )->revert_units();
		return;
	} catch (std::out_of_range& oor) {}
	try {
		scalar_contents_.at( key )->revert_units();
	} catch (std::out_of_range& oor) {
		throw std::out_of_range( "No vector or scalar quantity found with key " + key );
	}
}



void NCPA::Atmosphere1D::do_units_conversion_( size_t n_points, double *inplace, 
			NCPA::units_t fromUnits, NCPA::units_t toUnits ) {

	// try to convert
	double *units_buffer = new double[ n_points ];
	std::memset( units_buffer, 0, n_points * sizeof( double ) );
	
	// throws out_of_range if conversion is undefined
	NCPA::Units::convert( inplace, n_points, fromUnits, toUnits, units_buffer );

	// successful, so record the units change
	std::memcpy( inplace, units_buffer, n_points * sizeof( double ) );
	delete [] units_buffer;
}

std::vector< std::string > NCPA::Atmosphere1D::get_keys() const {
	std::vector< std::string > keys;
	for (std::map< std::string, NCPA::AtmosphericProperty1D * >::const_iterator it = contents_.begin();
			it != contents_.end(); ++it) {
		keys.push_back( (*it).first );
	}
	return keys;
}

void NCPA::Atmosphere1D::print_atmosphere( const std::vector< std::string >& columnorder, std::string altitude_key, 
	std::ostream& os ) {

	// check columnorder variable for key validity
	std::vector< std::string >::const_iterator vit;
	for (vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit) {
		if (! contains_vector( *vit ) ) {
			throw std::invalid_argument( "No vector quantity exists with key " + *vit );
		}
	}

	// first we do the header.  That contains column descriptions as well as scalar values
	// scalars first
	for (auto mit = scalar_contents_.cbegin(); mit != scalar_contents_.cend(); ++mit ) {
		os  << "#% 0, " << (*mit).first << ", " 
			<< NCPA::Units::toStr( (*mit).second->get_units() ) << ", "
			<< (*mit).second->get() << std::endl;
	}

	// Now column descriptors.  Altitude first
	os  << "#% 1, " << altitude_key << ", " 
		<< NCPA::Units::toStr( z_units_.top() ) << std::endl;
	unsigned int column = 2;
	for ( vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit ) {
		os  << "#% " << column << ", "
			<< *vit << ", "
			<< NCPA::Units::toStr( get_property_units( *vit ) ) << std::endl;
		column++;
	}

	// Now columns
	os.setf( std::ios::scientific, 	std::ios::floatfield );
	os.setf( std::ios::right, 		std::ios::adjustfield );
	os.precision( 6 );
	os.width( 9 );
	os.fill( ' ' );
	for ( size_t i = 0; i < nz_; i++) {
		os << z_[ i ];
		for (vit = columnorder.cbegin(); vit != columnorder.cend(); ++vit ) {
			os << " " << get( *vit, z_[ i ] );
		}
		os << std::endl;
	}
	os.flush();
}

void NCPA::Atmosphere1D::print_atmosphere( std::string altitude_key, std::ostream& os ) {
	print_atmosphere( get_keys(), altitude_key, os );
}