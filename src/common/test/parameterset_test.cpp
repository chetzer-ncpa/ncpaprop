#include "parameterset.h"
#include <iostream>

using namespace std;
using namespace NCPA;


int main( int argc, char **argv ) {

	// Create parameter set object
	ParameterSet *ps = new ParameterSet();
	ParameterTest *test = NULL;

	// set up expected commands
	ps->setStrict( false );
	ps->setComments( "#%" );

	ps->addParameter( new FlagParameter( "help" ) );
	ps->addParameter( new FlagParameter( "h" ) );
	ps->addUsageLine( "  --help, -h       Prints help text" );

	ps->addParameter( new IntegerParameter( "testval" ) );
	ps->addUsageLine( "  --testval        Test integer value" );
	//test = ps->addTest( "testval", PARAMETER_TEST_INTEGER_GREATER_THAN );
	//test->addIntegerParameter( 10 );
	ps->addTest( new IntegerGreaterThanTest( "testval", 10 ) );

	ps->addParameter( new FloatParameter( "testfloat" ) );
	ps->addUsageLine( "  --testfloat      Test double value" );

	// full suite of tests. for each known test type
	//ps->addParameter( new FlagParameter( "flagtwo" ) );
	//ps->addParameter( new FlagParameter( "flagthree" ) );
	//test = ps->addTest( "oneflag", PARAMETER_TEST_RADIO_BUTTON );
	//test->addStringParameter( "flagone" );
	//test->addStringParameter( "flagtwo" );
	//test->addStringParameter( "flagthree" );
	string flagList[ 3 ] = { "flagone", "flagtwo", "flagthree" };
	for (unsigned int i = 0; i < 3; i++) {
		ps->addParameter( new FlagParameter( flagList[ i ] ) );
	}
	ps->addTest( new RadioButtonTest( "oneflag", 3, flagList ) );



	// positive integer
	ps->addParameter( new IntegerParameter( "int_positive" ) );
	ps->addTest( new IntegerGreaterThanTest( "int_positive", 0 ) );

	ps->addParameter( new IntegerParameter( "int_not_negative"));
	ps->addTest( new IntegerGreaterThanOrEqualToTest( "int_not_negative", 0 ) );

	
	

	// first, check command line to see if help was requested
	unsigned int ncomm = ps->parseCommandLine( argc, argv );
	if (ps->getParameter("help")->wasFound() || ps->getParameter("h")->wasFound()) {
		ps->printUsage();
		exit(0);
	}

	// read in command line 
	try {
		ncomm = ps->parseFile( "test.options" );
		if (ncomm == 0) {
			cout << "No file options read" << endl;
		}
		ncomm = ps->parseCommandLine( argc, argv );
		if (ncomm == 0) {
			cout << "No command line options read" << endl;
		}
	} catch (invalid_argument &err) {
		cout << "Error parsing arguments: " << err.what() << endl;
		exit( 1 );
	}

	cout << "Parameters:" << endl << endl;
	ps->printParameters( cout );

	cout << endl << "Validation:" << endl;
	if (ps->validate()) {
		cout << "All tests passed!" << endl;
	} else {
		ps->printFailedTests( cout );
	}
}