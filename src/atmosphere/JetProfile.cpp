#include "JetProfile.h"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <cstring>


#define GAM 1.4
#define PI 3.14159
#define R 287.0


NCPA::JetProfile::JetProfile() {
	
	_jets.clear();
	T0 = 288.2;
	rho0 = 1.225e-3;
	eps_z = 0.001;
	
	
    // set up with default values.  This will eventually be deprecated, but is
    // good enough for testing for now.
	/*
    v_strat = 0.06;
    v_jet = 0.0;
    v_noct = 0.0;
    z_strat = 60;
    z_jet = 10;
    z_noct = 0.3;
    strat_jet_width = 17.5;
    jet_width = 2;
    noct_width = 0.6;

    azimuth = 0.0 / 180.0 * PI;
    jet_azimuth = 0.0 / 180.0 * PI;
    noct_azimuth = 0.0 / 180.0 * PI;
*/
    
}

NCPA::JetProfile::JetProfile(std::string filename) {
	// Start parsing file
	std::ifstream paramfile;
	paramfile.exceptions( std::ifstream::failbit | std::ifstream::badbit );
	paramfile.open( filename.c_str(), std::ifstream::in );
	if (!paramfile.good()) {
		std::string error = "Can't open parameter file " + filename + "!";
		throw new std::ios_base::failure( error );
	}
	
	std::string line;
	int njets = 0;
	eps_z = 0.01;
	T0 = 288.2;
	rho0 = 1.225;
	
	// REad in the line and get rid of leading whitespace
	//std::getline( paramfile, line );
	NCPA::safe_getline(paramfile, line);
	line = NCPA::deblank( line );
	while (line.length() == 0 || line[0] == '#')    // ignore comments
		NCPA::safe_getline(paramfile, line);
		//std::getline( paramfile, line );
	
	njets = std::atoi( line.c_str() );    // first non-commented line needs to be the number of jets described
	if (njets == 0) {
		std::string error = "Error reading number of jets from parameter file " + filename + "!";
		throw new std::ios_base::failure( error );
	}
	char cline[ 4096 ];
	for (int i = 0; i < njets; i++) {
		double height = 0.0, speed = 0.0, width = 0.0, azimuth = 0.0;
		for (int k = 0; k < 4; k++) {
			
			// read in the four quantities
			//std::getline( paramfile, line );
			NCPA::safe_getline( paramfile, line );
			line = NCPA::deblank( line );
			while (line.length() == 0 || line[0] == '#')
				NCPA::safe_getline( paramfile, line );
				//std::getline( paramfile, line );
			
			// put it in a char* to use strtok()
			std::strncpy( cline, line.c_str(), line.length() );
			cline[ line.length() ] = '\0';
			
			char *key = std::strtok( cline, " ,:" );
			char *val = std::strtok( NULL, " ,:" );
			if (std::strcmp(key,"height") == 0) {
				height = std::atof( val );
			} else if (std::strcmp(key,"speed") == 0) {
				speed = std::atof( val );
			} else if (std::strcmp(key,"width") == 0) {
				width = std::atof( val );
			} else if (std::strcmp(key,"azimuth") == 0) {
				azimuth = NCPA::deg2rad(std::atof( val ));
			} else {
				std::string error = "Error: Unrecognized quantity '";
				error = error + key;
				throw new std::ios_base::failure( error );
			}
		}
		this->addJet( height, speed, width, azimuth );
	}
}

// put a new jet into the mix
void NCPA::JetProfile::addJet( double height, double speed, double width, double azimuth ) {
	_jetStruct jet;
	jet.height = height;
	jet.speed = speed;
	jet.width = width;
	jet.azimuth = azimuth;
	_jets.push_back( jet );
}

double NCPA::JetProfile::z0() {
    return 0.0;
}

double NCPA::JetProfile::t(double z) {
    double poly_A, poly_B;

    poly_A = 1.0 + A5 * z + A6 * pow(z, 2) + A7 * pow(z, 3) + A8 * pow(z, 4);
    poly_B = 1.0 + B5 * z + B6 * pow(z, 2) + B7 * pow(z, 3) + B8 * pow(z, 4);

    return T0 * poly_A / poly_B; //define temp function and derivative
}

double NCPA::JetProfile::u(double z) {
    double answer = 0.0, cup;
    
    for (unsigned int i = 0; i < _jets.size(); i++) {
	    _jetStruct jet = _jets[ i ];
	    cup = (z - jet.height) / jet.width;
	    cup *= cup;
	    answer += jet.speed * std::sin( jet.azimuth ) * std::exp( -cup );
    }

    return answer;
}

double NCPA::JetProfile::v(double z) {
    double answer = 0.0, cup;
    
    for (unsigned int i = 0; i < _jets.size(); i++) {
	    _jetStruct jet = _jets[ i ];
	    cup = (z - jet.height) / jet.width;
	    cup *= cup;
	    answer += jet.speed * std::cos( jet.azimuth ) * std::exp( -cup );
    }

    return answer;

    /*
    cup = (z - z_noct) / noct_width;
    cup = cup*cup;
    answer += v_noct * cos(noct_azimuth - azimuth) * exp(-cup);
    */

}

double NCPA::JetProfile::w(double z) {
    return 0.0;
}

double NCPA::JetProfile::c0(double z) {
    return 1.0e-3 * sqrt(GAM * R * this->t(z));
}

double NCPA::JetProfile::ceff(double z, double phi) {
    return this->c0(z) + this->wcomponent(z, phi);
}

double NCPA::JetProfile::rho(double z) {
    double polyA, polyB;

    polyA = A1 * z + A2 * pow(z, 2) + A3 * pow(z, 3) + A4 * pow(z, 4);
    polyB = 1 + B1 * z + B2 * pow(z, 2) + B3 * pow(z, 3) + B4 * pow(z, 4);

    return rho0 * pow(10, polyA / polyB);
}

double NCPA::JetProfile::p( double z ) {
	
	// compute via polynomial fit per Sutherland and Bass (2004)
	double P0, P1, P2, P3, P4, P5;
	if (z > 90) {
		P0=-1.09914e+2;
		P1=4.7143;
		P2=-8.2123e-2;
		P3=6.622e-4;
		P4=-2.5593e-6;
		P5=3.84825e-9;
	} else if (z > 83) {
		P0=-1.1279e+1;
		P1=3.1982e-1;
		P2=-5.8093e-3;
		P3=2.2331e-5;
		P4=0.;
		P5=0.;
	} else if (z > 50) {
		P0=0.76338;
		P1=-2.58298e-1;
		P2=3.76139e-3;
		P3=-4.20887e-5;
		P4=1.602e-7;
		P5=-1.92509e-10;
	} else if (z > 28) {
		P0=2.18198;
		P1=-4.11497e-1;
		P2=1.33665e-2;
		P3=-3.59519e-4;
		P4=5.10097e-6;
		P5=-2.89056e-8;
	} else if (z > 18) {
		P0=0.984143;
		P1=-0.269769;
		P2=8.52275e-3;
		P3=-3.96203e-4;
		P4=1.01465e-5;
		P5=-1.02643e-7;
	} else if (z > 11) {
		P0=-7.99108e-2;
		P1=-8.10464e-2;
		P2=-5.55224e-3;
		P3=3.1117e-4;
		P4=-1.66878e-5;
		P5=3.832e-7;
	} else {
		P0=1.68716e-2;
		P1=-1.14252e-1;
		P2=-1.36123e-3;
		P3=7.36241e-5;
		P4=-1.08003e-5;
		P5=3.30464e-7;
	}
	double B=(P0 + P1*z + P2*z*z + P3*z*z*z + P4*z*z*z*z + P5*z*z*z*z*z);
	double P=(std::exp(B))*10.e+2;
	return P;

}