#include "parameterset.h"
#include "modess_parameters.h"
#include <string>


void NCPA::configure_modess_parameter_set( NCPA::ParameterSet *ps ) {

	NCPA::ParameterTest *test = NULL;

	// general configuration
	// set up expected commands
	ps->setStrict( true );
	ps->setComments( "#%" );
	ps->setDelimiters( "=: " );

	// @todo: Add general descriptive text here
	ps->addParameter( new FlagParameter( "help" ) );
	ps->addParameterDescription( "Options Control", "--help", "Prints help test" );
	//ps->addUsageLine( "  --help, -h              Prints help text" );
	//ps->addUsageLine( "" );
	ps->addParameter( new StringParameter( "paramfile", "modess.param") );
	//ps->addUsageLine( "  --paramfile             Parameter file [modess.param]" );
	ps->addParameterDescription( "Options Control", "--paramfile", "Parameter file name [modess.param]" );

	// Required parameters
	//ps->addUsageLine( "Required Parameters:" );
	ps->addParameter( new NCPA::StringParameter( "atmosfile" ) );
	ps->addTest( new NCPA::RequiredTest( "atmosfile" ) );
	//ps->addUsageLine( "  --atmosfile             Atmospheric profile file" );
	ps->addParameterDescription( "Required Parameters", "--atmosfile", "Atmospheric profile filename" );

	ps->addParameter( new NCPA::StringParameter( "atmosfileorder" ) );
	ps->addTest( new NCPA::RequiredTest( "atmosfileorder" ) );
	/*
	ps->addUsageLine( "  --atmosfileorder        The order of the (z,u,v,w,t,d,p) fields in --atmosfile." );
	ps->addUsageLine( "                          The units assumed in the ASCII file are z[km], t [kelvin]," );
	ps->addUsageLine( "                          d [g/cm^3], p [hectoPa]" );
	ps->addUsageLine( "                          The wind speeds are in m/s by default;" );
	ps->addUsageLine( "                          however if the winds are given in km/s then use" ); 
	ps->addUsageLine( "                          option --wind_units kmpersec" );
	*/
	ps->addParameterDescription( "Required Parameters", "--atmosfileorder", "The order of the (z,u,v,w,t,d,p) fields in --atmosfile. The units assumed in the ASCII file are z[km], t [kelvin], d [g/cm^3], p [hectoPa]. The wind speeds are in m/s by default; however if the winds are given in km/s then use option --wind_units kmpersec");

	ps->addParameter( new NCPA::FloatParameter( "freq" ) );
	ps->addTest( new NCPA::RequiredTest( "freq" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "freq", 0.0 ) );
	//ps->addUsageLine( "  --freq                  Frequency of analysis (Hz)" );
	ps->addParameterDescription( "Required Parameters", "--freq", "Frequency of analysis (Hz)" );

	// optional parameters
	//ps->addUsageLine( "Optional Parameters [default value]:" );
	ps->addParameter( new NCPA::IntegerParameter( "skiplines", 0 ) );
	ps->addTest( new NCPA::IntegerGreaterThanOrEqualToTest( "skiplines", 0 ) );
	//ps->addUsageLine( "  --skiplines             Number of header lines to skip in --atmosfile [0]" );
	ps->addParameterDescription( "Optional Parameters [default]", "--skiplines", "Number of header lines to skip in --atmosfile [0]" );


	ps->addParameter( new NCPA::FloatParameter( "maxheight_km", 150.0 ) );
	//ps->addUsageLine( "  --maxheight_km          Maximum height of analysis in km [150.0]" );
	ps->addParameterDescription( "Optional Parameters [default]", "--maxheight_km", "Maximum height of analysis in km [150.0]" );


	ps->addParameter( new NCPA::FloatParameter( "zground_km", 0.0 ) );
	//ps->addUsageLine( "  --zground_km            Ground height [take from Z0 parameter in atmosfile, or 0.0]" );
	ps->addParameterDescription( "Optional Parameters [default]", "--zground_km", "Ground height [take from Z0 parameter in atmosfile, or 0.0]" );


	ps->addParameter( new NCPA::IntegerParameter( "Nz_grid", 20000 ) );
	ps->addTest( new NCPA::IntegerGreaterThanTest( "Nz_grid", 0 ) );
	//ps->addUsageLine( "  --Nz_grid               Number of vertical grid points to use" );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nz_grid", "Number of vertical grid points to use [2000]" );

	ps->addParameter( new NCPA::FloatParameter( "sourceheight_km", 0.0 ) );
	//ps->addUsageLine( "  --sourceheight_km       Source height in km [0.0]" );
	ps->addParameterDescription( "Optional Parameters [default]", "--sourceheight_km", "Source height in km [0.0]" );


	ps->addParameter( new NCPA::FloatParameter( "receiverheight_km", 0.0 ) );
	//ps->addUsageLine( "  --receiverheight_km     Receiver height in km [0.0]" );
	ps->addParameterDescription( "Optional Parameters [default]", "--receiverheight_km", "Receiver height in km [0.0]" );

	ps->addParameter( new NCPA::FloatParameter( "maxrange_km", 1000.0 ) );
	//ps->addUsageLine( "  --maxrange_km           Maximum range in km to use for modeling [1000.0]" );
	ps->addTest( new NCPA::FloatGreaterThanTest( "maxrange_km", 0.0 ) );
	ps->addParameterDescription( "Optional Parameters [default]", "--maxrange_km", "Maximum range in km to use for modeling [1000.0]" );

	ps->addParameter( new NCPA::IntegerParameter( "Nrng_steps", 1000 ) );
	ps->addTest( new NCPA::IntegerGreaterThanTest( "Nrng_steps", 0 ) );
	//ps->addUsageLine( "  --Nrng_steps            Number of range steps to use [1000]" );
	ps->addParameterDescription( "Optional Parameters [default]", "--Nrng_steps", "Number of range steps to use [1000]" );

	ps->addParameter( new NCPA::StringParameter( "ground_impedence_model", "rigid" ) );
	test = ps->addTest( new NCPA::StringSetTest( "ground_impedence_model" ) );
	test->addStringParameter( "rigid" );
	//ps->addUsageLine( "  --ground_impedence_model" );
	//ps->addUsageLine( "                          Impedence model to use.  Currently only \"rigid\" is supported." );
	ps->addParameterDescription( "Optional Parameters [default]", "--ground_impedence_model", "Impedence model to use.  Currently only \"rigid\" is supported. [rigid]" );

	ps->addParameter( new NCPA::StringParameter( "wind_units", "mpersec" ) );
	test = ps->addTest( new NCPA::StringSetTest( "wind_units" ) );
	test->addStringParameter( "mpersec" );
	test->addStringParameter( "kmpersec" );
	//ps->addUsageLine( "  --wind_units            Units wind speed is presented in, kmpersec or mpersec [mpersec]" );
	ps->addParameterDescription( "Optional Parameters [default]", "--wind_units", "Units wind speed is presented in, kmpersec or mpersec [mpersec]" );


	ps->addParameter( new NCPA::StringParameter( "use_attn_file" ) );
	//ps->addUsageLine( "  --use_attn_file         File name containing attenuation, to override default Sutherland/Bass." );
	//ps->addUsageLine( "                          Columns are Height(km)  Attenuation (np/m)" );
	ps->addParameterDescription( "Optional Parameters [default]", "--use_attn_file", "File name containing attenuation, to override default Sutherland/Bass. Columns are #n# Height(km) Attenuation(np/m)" );

	ps->addParameter( new NCPA::StringParameter( "modal_starter_file" ) );
	ps->addUsageLine( "  --modal_starter_file    Filename to output a starter file for the pape module" );
	ps->addParameterDescription( "Optional Parameters [default]", "--modal_starter_file", "Filename to output a starter file for the pape module" );
	


	// Modes of operation
	std::string modes_of_operation[ 2 ] = { "singleprop", "Nby2Dprop" };
	for (unsigned int i = 0; i < 3; i++) {
		std::string tmpStr( modes_of_operation[ i ] );
		ps->addParameter( new NCPA::FlagParameter( tmpStr ) );
	}
	ps->addTest( new NCPA::RadioButtonTest( "operation_mode", 2, modes_of_operation ) );
	//ps->addUsageLine( "Modes of Operation - Specify One:" );
	//ps->addUsageLine( "  --singleprop            Single azimuth propagation" );
	//ps->addUsageLine( "    REQUIRED:" );
	ps->addParameterDescription( "Modes of Operation", "--singleprop", "Single azimuth propagation.  Requires --azimuth" );
	
	// for single propagation, must specify azimuth
	ps->addParameter( new NCPA::FloatParameter( "azimuth" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth", "singleprop" ) );
	//ps->addUsageLine( "    --azimuth             Azimuth of propagation, in degrees CW from North [0,360)" );
	//ps->addUsageLine( "" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth", "Azimuth of propagation, in degrees CW from North [0,360)", 2*DEFAULT_PARAMETER_INDENT );

	//ps->addUsageLine( "  --Nby2Dprop             Multiple azimuth propagation" );
	ps->addParameterDescription( "Modes of Operation", "--Nby2Dprop", "Multiple azimuth propagation.  Requires --azimuth_start, --azimuth_end, and --azimuth_step" );

	// for multiple propagation, must specify azimuth start, end, and step
	ps->addParameter( new NCPA::FloatParameter( "azimuth_start" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_start", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth_start", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_start", "Nby2Dprop" ) );
	//ps->addUsageLine( "    --azimuth_start       Starting azimuth, in degrees CW from North [0,360)" );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_end" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_end", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth_end", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_end", "Nby2Dprop" ) );
	//ps->addUsageLine( "    --azimuth_end         Ending azimuth, in degrees CW from North [0,360)" );
	ps->addParameter( new NCPA::FloatParameter( "azimuth_step" ) );
	ps->addTest( new NCPA::FloatGreaterThanOrEqualToTest( "azimuth_step", 0.0 ) );
	ps->addTest( new NCPA::FloatLessThanTest( "azimuth_step", 360.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "azimuth_step", "Nby2Dprop" ) );
	//ps->addUsageLine( "    --azimuth_step        Azimuth step, in degrees CW from North [0,360)" );
	//ps->addUsageLine( "" );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_start", "Starting azimuth, in degrees CW from North [0,360)", 2*DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_end", "Ending azimuth, in degrees CW from North [0,360)", 2*DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Modes of Operation", "--azimuth_step", "Azimuth step, in degrees CW from North [0,360)", 2*DEFAULT_PARAMETER_INDENT );



	// Setup flags
	//ps->addUsageLine( "Flags:" );
	ps->addParameter( new NCPA::FlagParameter( "write_2D_TLoss" ) );
	ps->addParameterDescription( "Flags", "--write_2D_TLoss", "Output 2-D transmission loss to tloss2D.nm" );

	ps->addParameter( new NCPA::FlagParameter( "write_phase_speeds" ) );
	ps->addParameterDescription( "Flags", "--write_phase_speeds", "Output phase speeds to phasespeeds.nm" );

	ps->addParameter( new NCPA::FlagParameter( "write_speeds" ) );
	ps->addParameterDescription( "Flags", "--write_speeds", "Output both the phase speeds and the group speeds to speeds.nm" );

	ps->addParameter( new NCPA::FlagParameter( "write_modes" ) );
	ps->addParameterDescription( "Flags", "--write_modes", "Output modes to mode_###.nm.  Also implies --write_speeds" );

	ps->addParameter( new NCPA::FlagParameter( "write_dispersion" ) );
	ps->addParameterDescription( "Flags", "--write_dispersion", "Output dispersion to dispersion_FFF.nm" );

	ps->addParameter( new NCPA::FlagParameter( "write_atm_profile" ) );
	ps->addParameterDescription( "Flags", "--write_atm_profile", "Output atmospheric profile to atm_profile.nm" );

	ps->addParameter( new NCPA::FlagParameter( "Lamb_wave_BC" ) );
	ps->addParameterDescription( "Flags", "--Lamb_wave_BC", "Use admittance = -1/2*dln(rho)/dz" );

	ps->addParameter( new NCPA::FlagParameter( "turnoff_WKB" ) );
	ps->addParameterDescription( "Flags", "--turnoff_WKB", "Turn off the WKB least phase speed estimation" );

	ps->addParameter( new NCPA::FlagParameter( "wvnum_filter" ) );
	//ps->addUsageLine( "  --wvnum_filter          Use wavenumber filter by phase speed.  Requires:" );
	ps->addParameterDescription( "Flags", "--wvnum_filter", "Use wavenumber filter by phase speed.  Requires --c_min and --c_max" );

	ps->addParameter( new NCPA::FloatParameter( "c_min" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "c_min", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "c_min", "wvnum_filter" ) );
	//ps->addUsageLine( "    --c_min               Minimum phase speed to keep" );
	ps->addParameter( new NCPA::FloatParameter( "c_max" ) );
	ps->addTest( new NCPA::FloatGreaterThanTest( "c_max", 0.0 ) );
	ps->addTest( new NCPA::RequiredIfOtherIsPresentTest( "c_max", "wvnum_filter" ) );
	//ps->addUsageLine( "    --c_max               Maximum phase speed to keep" );
	ps->addParameterDescription( "Flags", "--c_min", "Minimum phase speed to keep", 2*DEFAULT_PARAMETER_INDENT );
	ps->addParameterDescription( "Flags", "--c_max", "Maximum phase speed to keep", 2*DEFAULT_PARAMETER_INDENT );

}