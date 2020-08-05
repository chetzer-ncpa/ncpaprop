#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <complex>
#include <cstring>
#include <ctime>

#include "Atmosphere1D.h"
//#include "SolveModNB.h"
#include "modes.h"
#include "util.h"
#include "modess_parameters.h"

//#ifndef Pi
//#define Pi 3.141592653589793
//#endif

#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

/*
* This is Jelle Assink's Normal Modes code effective sound speed- moved to C++
* Contains some of Roger Waxler's original design also.
* @version 2.0
* @date 2020-04-21
* @authors Jelle Assink; Roger Waxler; Claus Hetzer; Doru Velea;
* 
* Changelog:
* 20130326: DV added turnoff_WKB flag
* 201305  : DV added use_attn_file option (to allow atten. coeff loaded from a text file)
* 201306  : DV modified the get() functions that return strings
* 20130827: DV added the computation of modal group velocities
* 20200421: CH Retrofitted to use new ParameterSet and Atmosphere1D libraries
*/

//
// main
//
int main( int argc, char **argv ) {
  
	// object to process the options
	ParameterSet *param = new ParameterSet();
	configure_modess_parameter_set( param );
	param->parseCommandLine( argc, argv );

	// check for help text
	if (param->wasFound( "help" ) || param->wasFound("h") ) {
		param->printUsage( cout );
		return 1;
	}

	// See if an options file was specified
	string paramFile = param->getString( "paramfile" );
	param->parseFile( paramFile );

	// parse command line again, to override file options
	param->parseCommandLine( argc, argv );

	// see if we want a parameter summary
	if (param->wasFound( "printparams" ) ) {
		param->printParameters();
	}

	// run parameter checks
	if (! param->validate() ) {
		cout << "Parameter validation failed:" << endl;
		param->printFailedTests( cout );
		return 0;
	}
	
	// set up to measure the duration of this run
	time_t tm1 = time(NULL);

	// open the file
	string atmosfile = 			param->getString( "atmosfile" );
	Atmosphere1D *atm_profile = new Atmosphere1D( atmosfile );
	 
	// get solver object
	//SolveModNB *a = new SolveModNB( param, atm_profile );
    ESSModeSolver *a = new ESSModeSolver( param, atm_profile );
                                         
	//   					 
	// compute modes - main action happens here
	//				 
	a->solve();
	a->printParams();

	// save atm. profile if requested
	if (param->getBool( "write_atm_profile" ) ) {
		ofstream ofs( "atm_profile.nm" );
		atm_profile->print_atmosphere( "Z", ofs );
		ofs.close();
	}
  
  	// clean up dynamic memory
	delete a;	 
	delete atm_profile;
	delete param;
	
	cout << endl << " ... main() is done." << endl;
	time_t tm2 = time(NULL);
	cout << "Run duration: " << difftime(tm2,tm1) << " seconds." << endl;	

	return (0);
} // end of main();
