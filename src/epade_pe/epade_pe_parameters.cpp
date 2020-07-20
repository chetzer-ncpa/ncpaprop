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
	std::string modes_of_operation[ 3 ] = { "atmosfile", "atmosfile2d", "toy" };
	ps->addParameter( new NCPA::StringParameter( modes_of_operation[ 0 ] ) );
	ps->addParameter( new NCPA::StringParameter( modes_of_operation[ 1 ] ) );
	ps->addParameter( new NCPA::FlagParameter( modes_of_operation[ 2 ] ) );
	ps->addTest( new NCPA::RadioButtonTest( "atmosphere_type", 3, modes_of_operation ) );
	ps->addParameterDescription( "Atmosphere", "--atmosfile", "1-D atmospheric profile filename" );
	ps->addParameterDescription( "Atmosphere", "--atmosfile2d", "2-D atmospheric profile filename" );
	ps->addParameterDescription( "Atmosphere", "--toy", "Use NCPA toy atmosphere" );


	// Required parameters
	ps->addParameter( new NCPA::FloatParameter( "freq" ) );
	ps->addTest( new NCPA::RequiredTest( "freq" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "freq", 0.0 ) );
	ps->addParameterDescription( "Required Parameters", "--freq", "Frequency of analysis (Hz)" );

	ps->addParameter( new NCPA::StringParameter( "starter" ) );
	ps->addTest( new NCPA::RequiredTest( "starter" ) );
	test = ps->addTest( new NCPA::StringSetTest( "starter" ) );
	test->addStringParameter( "self" );
	test->addStringParameter( "gaussian" );
	ps->addParameterDescription( "Required Parameters", "--starter", "Starter type: one of { self, gaussian }" );

	ps->addParameter( new NCPA::FloatParameter( "azimuth" ) );
	ps->addTest( new NCPA::RequiredTest( "azimuth" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanOrEqualToTest( "azimuth", 360.0 ) );
	ps->addParameterDescription( "Required Parameters", "--azimuth", "Propagation azimuth (degrees clockwise from north)" );

	ps->addParameter( new NCPA::FloatParameter( "maxrange_km" ) );
	ps->addTest( new NCPA::RequiredTest( "maxrange_km" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "maxrange_km", 0.01 ) );
	ps->addParameterDescription( "Required Parameters", "--maxrange_km", "Maximum range in km to use for modeling" );

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
/*
	ps->addParameter( new NCPA::FloatParameter( "receiverheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--receiverheight_km", "Receiver height in km [ground]" );
*/
	ps->addParameter( new NCPA::FloatParameter( "groundheight_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--groundheight_km", "Ground height in km [0.0, or Z0 parameter in profile]" );

	ps->addParameter( new NCPA::IntegerParameter( "Nrng_steps", 0 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "Nrng_steps", 0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nrng_steps", "Number of range steps to use [decide internally]" );

	ps->addParameter( new NCPA::StringParameter( "ground_impedence_model", "rigid" ) );
	test = ps->addTest( new NCPA::StringSetTest( "ground_impedence_model" ) );
	test->addStringParameter( "rigid" );
	ps->addParameterDescription( "Optional Parameters [default]", "--ground_impedence_model", "Impedence model to use.  Currently only \"rigid\" is supported. [rigid]" );


/*
	ps->addParameter( new NCPA::StringParameter( "use_attn_file", "" ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--use_attn_file", "File name containing attenuation, to override default Sutherland/Bass [n/a]. Columns are #n# Height(km) Attenuation(np/m)" );
*/

	// Setup flags
	//ps->addUsageLine( "Flags:" );
	ps->addParameter( new NCPA::FlagParameter( "write_atm_profile" ) );
	ps->addParameterDescription( "Flags", "--write_atm_profile", "Output atmospheric profile to atm_profile.nm" );
	ps->addParameter( new NCPA::FlagParameter( "lossless" ) );
	ps->addParameterDescription( "Flags", "--lossless", "Ignore atmospheric attenuation" );
	ps->addParameter( new NCPA::FlagParameter( "topo" ) );
	ps->addParameterDescription( "Flags", "--topo", "Use topography.  Requires presence of Z0 parameter in atmospheric files" );
	ps->addParameter( new NCPA::FlagParameter( "disable_top_layer" ) );


	// Footer with file formats and sample commands
	ps->addBlankFooterLine();
	ps->addFooterText("OUTPUT Files:  Format description (column order):");
	ps->addFooterTextVerbatim("  tloss_1d.pe:                 r (km), TL (dB)");
	ps->addFooterTextVerbatim("  tloss_2d.nm:                 r, z, TL (dB), 0.0");
	ps->addFooterTextVerbatim("  atm_profile.nm               z,u,v,w,t,d,p,c,c_eff");
	ps->addBlankFooterLine();
	ps->addFooterText("Examples (run from 'samples' directory):");
	ps->setFooterIndent( 4 );
	ps->setFooterHangingIndent( 4 );
	ps->setCommandMode( true );
	//ps->addFooterText("../bin/Modess --singleprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1" );
	ps->addFooterText("../bin/ePape --starter gaussian --toy --freq 0.1 --azimuth 90 --maxrange_km 1000" );
	ps->addBlankFooterLine();
	//ps->addFooterText("../bin/Modess --singleprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --freq 0.1 --write_2D_TLoss" );
	ps->addFooterText("../bin/ePape --starter self --atmosfile NCPA_canonical_profile_trimmed.dat --freq 0.1 --azimuth 90 --maxrange_km 1000" );
	ps->addBlankFooterLine();
	ps->addFooterText("../bin/ePape --starter self --atmosfile2d atmosphere_2d_summary.dat --freq 1.0 --azimuth 90 --maxrange_km 1000 --lossless" );
	//ps->addFooterText("../bin/Modess --Nby2Dprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --freq 0.1 --azimuth_start 0 --azimuth_end 360 --azimuth_step 1" );
	//ps->addFooterText("../bin/Modess --Nby2Dprop --atmosfile NCPA_canonical_profile_zuvwtdp.dat --freq 0.1 --azimuth_start 0 --azimuth_end 360 --azimuth_step 1" );
	ps->setFooterHangingIndent( 0 );
	ps->setCommandMode( false );
	ps->resetFooterIndent();
}
