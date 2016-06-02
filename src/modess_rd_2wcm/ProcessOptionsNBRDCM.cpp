#include <sstream>
#include <stdexcept>
#include "anyoption.h"
#include "ProcessOptionsNBRDCM.h"

#include "Atmosphere.h"
#include "ModessRDCM_lib.h"

//#include <vector>

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

// constructor
NCPA::ProcessOptionsNB::ProcessOptionsNB(AnyOption *opt) {

  // defaults
  z_min            = 0.0;         // meters 
  maxrange         = 1.0E6;       // meters
  maxheight        = 150000.0;    // meters
  sourceheight     = z_min;       // meters
  receiverheight   = z_min;       // meters
  Nz_grid          = 20000;       // number of points on the z-grid
  Nrng_steps       = 1000;        // number of range steps	
  skiplines        = 0;           // skiplines in "atmosfile"
  Lamb_wave_BC     = 0;           // 1 to enforce the Lamb wave BC
  gnd_imp_model    = "rigid";     // rigid ground
  wind_units       = "mpersec";   // m/s
  usrattfile       = "";             // user-provided attenuation filename
  req_profile_step = maxrange;    // specifies the range step to request a new profile 
  tol              = 1.0E-8;      // tolerance for Slepc calculations
  
  write_2D_TLoss     = opt->getFlag( "write_2D_TLoss");
  write_phase_speeds = opt->getFlag( "write_phase_speeds" );   
  write_modes        = opt->getFlag( "write_modes" );
  write_dispersion   = opt->getFlag( "write_dispersion" );
  turnoff_WKB        = opt->getFlag( "turnoff_WKB" ); // if ==1 turns off the WKB least phase speed approx
  plot_flg           = opt->getFlag( "plot");         // flag to plot results with gnuplot  
  
  
  // Parse arguments based on file type selected and set defaults
  // Declare and populate variables
  // First, file type
  enum AtmosphericFileType { G2SENVFILE, ATMOSFILE, ASCIIDIR, JETFILE, SLICEFILE };   // Expand this enum as we put in more file types
  //AtmosphericFileType filetype;
  int numTypesDeclared = 0;
  
  if ( opt->getValue( "use_1D_profiles_from_dir" ) != NULL ) {
      atm_profile_dir = opt->getValue( "use_1D_profiles_from_dir" );
      filetype = ASCIIDIR;
      numTypesDeclared++;
  }

  if ( opt->getValue( "jetfile" ) != NULL ) {
      atmosfile = opt->getValue( "jetfile" );
      filetype = JETFILE;
      numTypesDeclared++;
  }
  
  //if ( opt->getValue( "atmosfile" ) != NULL ) {
  //    atmosfile = opt->getValue( "atmosfile" );
  //    filetype = ATMOSFILE;
  //    numTypesDeclared++;
  //    //cout << "Atm file: " << atmosfile << endl;
  //}
  
  if (opt->getValue( "g2senvfile" ) != NULL ) {
      atmosfile = opt->getValue( "g2senvfile" );
      filetype = G2SENVFILE;
      numTypesDeclared++;
  } 

  if (opt->getValue( "slicefile" ) != NULL ) {
      atmosfile = opt->getValue( "slicefile" );
      filetype = SLICEFILE;
      numTypesDeclared++;
      cout << "Slice file: " << atmosfile << endl;
  }

  // Make sure only one file type was requested
  if (numTypesDeclared != 1) {
      delete opt;
      throw invalid_argument( "Choose only one source for the atmospheric profiles: either option --use_1D_profiles_from_dir or --g2senvfile.");
  }
  
  if (opt->getValue( "atmosfileorder" ) != NULL) {
      atmosfileorder = opt->getValue( "atmosfileorder" );
  } 
  else {
      delete opt;
      throw invalid_argument( "Option --atmosfileorder is required for ATMOSFILE files!" );
  }  
  
  // set skiplines to skip in the atmosfile
  if (opt->getValue( "skiplines" ) != NULL) {
      skiplines = atoi( opt->getValue( "skiplines" ) );
      //cout << "skiplines = " << skiplines << endl;
  } else {
      delete opt;
      throw invalid_argument( "Option --skiplines is required!" );
  }

  if (opt->getValue( "azimuth" ) != NULL) {
      azi = atof( opt->getValue("azimuth") );
      //cout << "azimuth = " << azi << endl;
  } else {
      delete opt;
      throw invalid_argument( "Option --azimuth is required!" );
  }	
	
  if ( opt->getValue( "freq" ) != NULL ) {
      freq = atof(opt->getValue( "freq" ));
      if (freq < 0) {
          //cout << "freq = " << freq << endl;
          delete opt;
          throw invalid_argument("Frequency must be positive.");
      }
  }
  else {
      delete opt;
      throw invalid_argument( "Option --freq is required!" );
  }	                 

  if (opt->getValue( "maxrange_km" ) != NULL) {
      maxrange = atof( opt->getValue( "maxrange_km" ))*1000.0;
      if (maxrange < 1000) {
      //cout << "maxrange = " << maxrange << " meters" << endl; // in meters
      delete opt;
      throw invalid_argument("maxrange is too short (< 1 km).");
      }
  }

  if (opt->getValue( "maxheight_km" ) != NULL) {
      maxheight = atof( opt->getValue( "maxheight_km" ))*1000.0;
      if (maxheight < 1000) {
          //cout << "maxheight = " << maxheight << " meters" << endl; // in meters
          delete opt;
          throw invalid_argument("maxheight is too short (< 1 km).");
      }  
  }

  if (opt->getValue( "sourceheight_km" ) != NULL) {
      sourceheight = atof( opt->getValue( "sourceheight_km" ))*1000.0;
      if (sourceheight > maxheight) {
          //cout << "sourceheight = " << sourceheight << " meters" << endl; // in meters
          delete opt;
          throw invalid_argument("sourceheight cannot be higher than maxheight.");
      }
  }
  
  if (opt->getValue( "receiverheight_km" ) != NULL) {
      receiverheight = atof( opt->getValue( "receiverheight_km" ))*1000.0;
      if (receiverheight > maxheight || receiverheight < 0) {
          //cout << "receiverheight = " << receiverheight << " meters" << endl; // in meters
          delete opt;
          throw invalid_argument("receiverheight cannot be higher than maxheight or less than 0.");
      }
  }
  
  if (opt->getValue( "zground_km" ) != NULL) {
      z_min = atof( opt->getValue( "zground_km" ))*1000.0;
      if (z_min > maxheight || z_min < 0) {
          //cout << "z_min = " << z_min << " meters" << endl; // in meters
          delete opt;
          throw invalid_argument("z_min cannot be higher than maxheight or less than 0.");
      }
  }        

  if ( opt->getValue( "Nz_grid" ) != NULL ) {
      Nz_grid = atoi(opt->getValue( "Nz_grid" ));
      if (Nz_grid < 10) {
          //cout << "Nz_grid = " << Nz_grid << endl;
          delete opt;
          throw invalid_argument("Nz_grid value is too low (typically 20000)");
      }
  }	
        
  if ( opt->getValue( "Nrng_steps" ) != NULL ) {
      Nrng_steps = atoi(opt->getValue( "Nrng_steps" ));
      if (Nrng_steps < 1) {
          // cout << "Nrng_steps = " << Nrng_steps << endl;
          delete opt;
          throw invalid_argument("Nrng_steps should be at least 1.");
      }
  }	

  //string gnd_imp_model ("rigid");
  if ( opt->getValue( "ground_impedance_model" ) != NULL ) {
      gnd_imp_model.assign(opt->getValue( "ground_impedance_model" ));
      if (gnd_imp_model.compare("rigid")!=0) {
          //cout << "ground_impedance_model = " << gnd_imp_model << endl;
          std::ostringstream es;
          es << "This ground impedance model is not implemented yet: " << gnd_imp_model;
          delete opt;
          throw invalid_argument(es.str());
      }   
  }	
  

  //int Lamb_wave_BC = 0; // for rigid ground: if ==1 then admittance = -1/2*dln(rho)/dz
  if ( opt->getValue( "Lamb_wave_BC" ) != NULL ) {
      Lamb_wave_BC = atoi(opt->getValue( "Lamb_wave_BC" ));
      //cout << "Lamb_wave_BC = " << Lamb_wave_BC << endl;
  }	
  
   
  if ( opt->getValue( "use_profile_ranges_km" ) != NULL ) {
      profile_ranges_given = 1;
      prf_ranges_km = opt->getValue( "use_profile_ranges_km" );
      //cout << "prf_ranges_km= " << prf_ranges_km << endl;
  }
  
  if (opt->getValue( "use_profiles_at_steps_km" ) != NULL) {
      req_profile_step = atof( opt->getValue( "use_profiles_at_steps_km" ))*1000.0;
  }
  
  if ( opt->getValue( "wind_units" ) != NULL ) {
      wind_units = opt->getValue( "wind_units" );
  }
  
  if ( opt->getValue( "use_attn_file" ) != NULL ) {
      usrattfile = opt->getValue( "use_attn_file" );
      cout << "User-provided attenuation file: " << usrattfile << endl;
  }  
  
} // ------- End of processOptions ---------------------------------------------

// utility to print the parameters to the screen
void NCPA::ProcessOptionsNB::printParams() {
  printf("\nRange Dependent Two-Way Coupled Modes run info:\n");
  printf("                   freq : %g\n", freq);
  printf("                azimuth : %g\n", azi);
  printf("                Nz_grid : %d\n", Nz_grid);
  printf("      z_min (meters MSL): %g\n", z_min);
  printf("      maxheight_km (MSL): %g\n", maxheight/1000.0);
  printf("   sourceheight_km (AGL): %g\n", sourceheight/1000.0);
  printf(" receiverheight_km (AGL): %g\n", receiverheight/1000.0);   
  printf("             Nrng_steps : %d\n", Nrng_steps);
  printf("            maxrange_km : %g\n", maxrange/1000.0);
  printf("          gnd_imp_model : %s\n", gnd_imp_model.c_str());
  printf("Lamb wave boundary cond : %d\n", Lamb_wave_BC);
  printf("  SLEPc tolerance param : %g\n", tol);
  printf("    write_2D_TLoss flag : %d\n", write_2D_TLoss);
  printf("             wind_units : %s\n", wind_units.c_str());  
  if (filetype==2) {
  printf("atmospheric profile dir : %s\n", atm_profile_dir.c_str());
  } else {
  printf("    atmospheric profile : %s\n", atmosfile.c_str());
  }    
}


// "get" functions
double NCPA::ProcessOptionsNB::getZ_min() {
  return z_min;
}

int    NCPA::ProcessOptionsNB::getFiletype() {
  return filetype;
}

std::string   NCPA::ProcessOptionsNB::getAtmosfile() {
  return atmosfile;
}

std::string   NCPA::ProcessOptionsNB::getAtmosfileorder() {
  return atmosfileorder;
}

std::string   NCPA::ProcessOptionsNB::getAtm_profile_dir() {
  return atm_profile_dir;
}

std::string   NCPA::ProcessOptionsNB::getUsrAttFile() {
  return usrattfile;
}

std::string   NCPA::ProcessOptionsNB::getWindUnits() {
  return wind_units;
}

int    NCPA::ProcessOptionsNB::getSkiplines() {
  return skiplines;
}

double NCPA::ProcessOptionsNB::getFreq() {
  return freq;
}

double NCPA::ProcessOptionsNB::getAzimuth() {
  return azi;
}

double NCPA::ProcessOptionsNB::getMaxrange() {
  return maxrange;
}

double NCPA::ProcessOptionsNB::getMaxheight() {
  return maxheight;
}

double NCPA::ProcessOptionsNB::getSourceheight() {
  return sourceheight;
}

double NCPA::ProcessOptionsNB::getReceiverheight() {
  return receiverheight;
}

double NCPA::ProcessOptionsNB::getSlepcTolerance() {
  return tol;
}

int    NCPA::ProcessOptionsNB::getNrng_steps() {
  return Nrng_steps;
}

int    NCPA::ProcessOptionsNB::getNz_grid() {
  return Nz_grid;
}

std::string   NCPA::ProcessOptionsNB::getGnd_imp_model() {
  return gnd_imp_model;
}

int   NCPA::ProcessOptionsNB::getLamb_wave_BC() {
  return Lamb_wave_BC;
}

bool   NCPA::ProcessOptionsNB::getWrite_2D_TLoss() {
  return write_2D_TLoss;
}

bool   NCPA::ProcessOptionsNB::getWrite_phase_speeds() {
  return write_phase_speeds;
}

bool   NCPA::ProcessOptionsNB::getWrite_modes() {
  return write_modes;
}

bool   NCPA::ProcessOptionsNB::getWrite_dispersion() {
  return write_dispersion;
}

std::string   NCPA::ProcessOptionsNB::getProfileRanges() {
  return prf_ranges_km;
}

bool   NCPA::ProcessOptionsNB::getProfile_ranges_given_flag() {
  return profile_ranges_given;
}

double NCPA::ProcessOptionsNB::getReq_profile_step() {
  return req_profile_step;
}

bool   NCPA::ProcessOptionsNB::getTurnoff_WKB() {
  return turnoff_WKB;
}

bool   NCPA::ProcessOptionsNB::getPlot_flg() {
  return plot_flg;
}

