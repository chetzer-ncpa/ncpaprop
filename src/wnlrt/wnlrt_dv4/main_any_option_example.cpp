#include <iostream>
//#include <cmath>
//#include <fstream>
//#include <vector>
//#include <stdexcept>
//#include <sstream>
//#include <complex>

//#include "Atmosphere.h"
#include "anyoption.h"
//#include "ProcessOptionsNB.h"
//#include "SolveModNB.h"
//#include "Modess_lib.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

//using namespace NCPA;
using namespace std;


// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

//
// main
//
//int main( int argc, char **argv ) {
  
  // Physical values are usually in SI (System International) units 
  // unless mentioned otherwise.
  
  cout << "In main" << endl;

  // parse options from the command line as well as an options file
  AnyOption *opt = parseInputOptions( argc, argv ); 

  // object to process the options
  //ProcessOptionsNB *oNB = new ProcessOptionsNB(opt);



  // set up to measure the duration of this run
  time_t tm1 = time(NULL);

  // obtain the parameter values from the user's options
  //atmosfile      = oNB->getAtmosfile();


  delete opt;
  //delete oNB;
  cout << "\n ... main() is done." << endl;
  time_t tm2 = time(NULL);
  cout << "Run duration: " << difftime(tm2,tm1) << " seconds." << endl;	

  return (0);
} // end of main();


AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();

  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "|                             Help                                         |" );
  opt->addUsage( "----------------------------------------------------------------------------" );	
  opt->addUsage( "Usage: " );
  opt->addUsage( "Insert Help text here" );
  opt->addUsage( "" );  

  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  opt->setFlag( "write_2D_TLoss" );

  opt->setOption( "atmosfile" );
  opt->setOption( "atmosfileorder" );
  opt->setOption( "wind_units" );


  // Process the command-line arguments
  opt->processFile( "../samples/Modess.options" );
  opt->processCommandArgs( argc, argv );

  if( ! opt->hasOptions()) { // print usage if no options
		  opt->printUsage();
		  delete opt;
		  exit( 1 );
  }

  // Check to see if help text was requested
  if ( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) {
	  opt->printUsage();
	  exit( 1 );
  }

  return opt;
}


 
