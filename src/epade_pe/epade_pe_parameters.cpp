#include "parameterset.h"
#include "epade_pe_parameters.h"
#include <string>
#include <iostream>

void NCPA::configure_epade_pe_parameter_set( NCPA::ParameterSet *ps ) {

	NCPA::ParameterTest *test = NULL;

	// general configuration
	// set up expected commands
	ps->setStrict( true );
	ps->setComments( "#%" );
	ps->setDelimiters( "=: " );

	// Add header instructions
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addHeaderTextVerbatim("|                             NCPA Infrasound                              |");
	ps->addHeaderTextVerbatim("|                         ePade Parabolic Equation                         |");
	ps->addHeaderTextVerbatim("|         Single Frequency - Effective Sound Speed Approximation           |");
	ps->addHeaderTextVerbatim("----------------------------------------------------------------------------");
	ps->addBlankHeaderLine();
	ps->addHeaderText("By default the program computes both the 1D (at the ground) and 2D transmission loss (TL) and saves the data to 2 files:" );
	ps->setHeaderIndent( 4 );
	ps->addHeaderText("file tloss_1d.pe - transmission loss at the ground" );
	ps->addHeaderText("file tloss_2d.pe - full 2-D transmission loss field" );
	ps->resetHeaderIndent();
	ps->addBlankHeaderLine();
	ps->addHeaderText("The options below can be specified in a colon-separated file \"epade_pe.param\" or at the command line. Command-line options override file options.");

	// Parameter descriptions
	ps->addParameter( new FlagParameter( "help" ) );
	ps->addParameter( new FlagParameter( "h" ) );
	ps->addParameterDescription( "Options Control", "--help", "Prints help test" );

	ps->addParameter( new StringParameter( "paramfile", "epade_pe.param") );
	ps->addParameterDescription( "Options Control", "--paramfile", "Parameter file name [epade_pe.param]" );

	ps->addParameter( new FlagParameter( "printparams" ) );
	ps->addParameterDescription( "Options Control", "--printparams", "Print parameter summary to screen" );

	// Atmosphere
	std::string atmosphere_types[ 3 ] = { "atmosfile", "atmosfile2d", "toy" };
	ps->addParameter( new NCPA::StringParameter( atmosphere_types[ 0 ] ) );
	ps->addParameter( new NCPA::StringParameter( atmosphere_types[ 1 ] ) );
	ps->addParameter( new NCPA::FlagParameter( atmosphere_types[ 2 ] ) );
	ps->addTest( new NCPA::RadioButtonTest( "atmosphere_type", 3, atmosphere_types ) );
	ps->addParameterDescription( "Atmosphere", "--atmosfile", "1-D atmospheric profile filename" );
	ps->addParameterDescription( "Atmosphere", "--atmosfile2d", "2-D atmospheric profile filename" );
	ps->addParameterDescription( "Atmosphere", "--toy", "Use NCPA toy atmosphere" );


	// Required parameters
	ps->addParameter( new NCPA::FloatParameter( "freq" ) );
	ps->addTest( new NCPA::RequiredTest( "freq" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "freq", 0.0 ) );
	ps->addParameterDescription( "Required Parameters", "--freq", "Frequency of analysis (Hz)" );

	ps->addParameter( new NCPA::StringParameter( "starter", "self" ) );
	// ps->addTest( new NCPA::RequiredTest( "starter" ) );
	// test = ps->addTest( new NCPA::StringSetTest( "starter" ) );
	// test->addStringParameter( "self" );
	// test->addStringParameter( "gaussian" );
	// ps->addParameterDescription( "Required Parameters", "--starter", "Starter type: one of { self, gaussian }" );

	ps->addParameter( new NCPA::FloatParameter( "maxrange_km" ) );
	ps->addTest( new NCPA::RequiredTest( "maxrange_km" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxrange_km", 0.01 ) );
	ps->addParameterDescription( "Required Parameters", "--maxrange_km", "Maximum range in km to use for modeling" );

	// Modes of operation
	std::string modes_of_operation[ 2 ] = { "singleprop", "multiprop" };
	for (unsigned int i = 0; i < 2; i++) {
		std::string tmpStr( modes_of_operation[ i ] );
		ps->addParameter( new NCPA::FlagParameter( tmpStr ) );
	}
	ps->addTest( new NCPA::RadioButtonTest( "operation_mode", 2, modes_of_operation ) );
	ps->addParameterDescription( "Modes of Operation", "--singleprop", "Single azimuth propagation.  Requires --azimuth" );
	
	// for single propagation, must specify azimuth
	ps->addParameter( new NCPA::FloatParameter( "azimuth" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth", "singleprop" ) );
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth", "Propagation azimuth ( degrees clockwise from north, [0,360) )" );
	ps->resetParameterIndent();

	// for multiple propagation, must specify azimuth start, end, and step

	ps->addParameterDescription( "Modes of Operation", "--multiprop", "Multiple azimuth propagation.  Requires --azimuth_start, --azimuth_end, and --azimuth_step, disables --atmosfile2d and --topo" );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_start" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_start", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_start", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_start", "multiprop" ) );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_end" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_end", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_end", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_end", "multiprop" ) );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_step" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_step", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth_step", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_step", "multiprop" ) );
	ps->setParameterIndent( 2 * DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_start", "Starting azimuth, in degrees CW from North [0,360)" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_end", "Ending azimuth, in degrees CW from North [0,360)" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_step", "Azimuth step, in degrees CW from North [0,360)" );
	ps->resetParameterIndent();



	// optional parameters
	ps->addParameter( new NCPA::IntegerParameter( "npade", 4 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "npade", 3 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--npade", "Number of Pade coefficients to use [4]" );

	ps->addParameter( new NCPA::FloatParameter( "dz_m", -1.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--dz_m", "Altitude resolution in meters [automatic]" );

	ps->addParameter( new NCPA::FloatParameter( "maxheight_km", 150.0 ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxheight_km", 0.1 ));
	ps->addParameterDescription( "Optional Parameters [default]", "--maxheight_km", "Maximum height of analysis in km [150.0, or max height of atmosphere]" );

	/*
	ps->addParameter( new NCPA::IntegerParameter( "Nz_grid", 5001 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "Nz_grid", 10 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nz_grid", "Number of vertical grid points to use [5001]" );
	*/


	ps->addParameter( new NCPA::FloatParameter( "sourceheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--sourceheight_km", "Source height in km [ground]" );

	ps->addParameter( new NCPA::FloatParameter( "receiverheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--receiverheight_km", "Receiver height in km [ground]" );

	ps->addParameter( new NCPA::FloatParameter( "groundheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--groundheight_km", "Ground height in km [Z0 parameter in profile, or 0.0]" );

	ps->addParameter( new NCPA::IntegerParameter( "Nrng_steps", 0 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "Nrng_steps", 0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nrng_steps", "Number of range steps to use [decide internally]" );

	ps->addParameter( new NCPA::StringParameter( "ground_impedence_model", "rigid" ) );
	test = ps->addTest( new NCPA::StringSetTest( "ground_impedence_model" ) );
	test->addStringParameter( "rigid" );
	ps->addParameterDescription( "Optional Parameters [default]", "--ground_impedence_model", "Impedence model to use.  Currently only \"rigid\" is supported. [rigid]" );

	ps->addParameter( new NCPA::StringParameter( "attnfile", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--attnfile", "File name containing attenuation, to override default Sutherland/Bass [n/a]. Columns are #n# Height(km) Attenuation(np/m)" );


	// Setup flags
	//ps->addUsageLine( "Flags:" );
	ps->addParameter( new NCPA::FlagParameter( "write_2d_tloss" ) );
	ps->addParameterDescription( "Flags", "--write_2d_tloss", "Output 2-D transmission loss to tloss_2d.pe" );
	ps->addParameter( new NCPA::FlagParameter( "write_atm_profile" ) );
	ps->addParameterDescription( "Flags", "--write_atm_profile", "Output atmospheric profile to atm_profile.pe" );
	ps->addParameter( new NCPA::FlagParameter( "lossless" ) );
	ps->addParameterDescription( "Flags", "--lossless", "Ignore atmospheric attenuation" );
	ps->addParameter( new NCPA::FlagParameter( "topo" ) );
	ps->addParameterDescription( "Flags", "--topo", "Use topography.  Requires presence of 'Z0' parameter in atmospheric files" );
	ps->addParameter( new NCPA::FlagParameter( "disable_top_layer" ) );
	ps->addParameter( new NCPA::FlagParameter( "broadband" ) );
	//ps->addParameterDescription( "Flags", "--broadband", "Calculate at multiple frequencies" );


	// Footer with file formats and sample commands
	ps->addBlankFooterLine();
	ps->addFooterText("OUTPUT Files:  Format description (column order):");
	ps->addFooterTextVerbatim("  tloss_1d.pe:                 az (deg), r (km), TL (real), TL (imag)");
	ps->addFooterTextVerbatim("  tloss_multiplot.pe           az (deg), r (km), TL (real), TL (imag)");
	ps->addFooterTextVerbatim("  tloss_2d.pe:                 r, z, TL (real), TL (imag)");
	ps->addFooterTextVerbatim("  atm_profile.pe:              z,u,v,w,t,d,p,c,c_eff");
	ps->addBlankFooterLine();
	ps->addFooterText("Examples (run from 'samples' directory):");
	ps->setFooterIndent( 4 );
	ps->setFooterHangingIndent( 4 );
	ps->setCommandMode( true );
	// ps->addFooterText("../bin/ePape --singleprop --starter gaussian --toy --freq 0.1 --azimuth 90 --maxrange_km 1000" );
	ps->addFooterText("../bin/ePape --singleprop --toy --freq 0.1 --azimuth 90 --maxrange_km 1000" );
	ps->addBlankFooterLine();
	// ps->addFooterText("../bin/ePape --singleprop --starter self --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.1 --azimuth 90 --maxrange_km 1000" );
	ps->addFooterText("../bin/ePape --singleprop --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.1 --azimuth 90 --maxrange_km 1000" );
	ps->addBlankFooterLine();
	//ps->addFooterText("../bin/ePape --singleprop --starter self --atmosfile2d atmosphere_2d_summary.dat --freq 0.5 --azimuth 90 --maxrange_km 1000 --lossless" );
	ps->addFooterText("../bin/ePape --singleprop --atmosfile2d atmosphere_2d_summary.dat --freq 0.5 --azimuth 90 --maxrange_km 1000 --lossless" );
	ps->addBlankFooterLine();
	// ps->addFooterText("../bin/ePape --singleprop --starter self --atmosfile2d atmosphere_2d_summary.dat --freq 1.0 --azimuth 90 --maxrange_km 500 --topo");
	ps->addFooterText("../bin/ePape --singleprop --atmosfile2d atmosphere_2d_summary.dat --freq 1.0 --azimuth 90 --maxrange_km 500 --topo");
	ps->addBlankFooterLine();
	// ps->addFooterText("../bin/ePape --multiprop --starter self --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 --maxrange_km 1000");
	ps->addFooterText("../bin/ePape --multiprop --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 2 --maxrange_km 1000");
	ps->setFooterHangingIndent( 0 );
	ps->setCommandMode( false );
	ps->resetFooterIndent();
}
