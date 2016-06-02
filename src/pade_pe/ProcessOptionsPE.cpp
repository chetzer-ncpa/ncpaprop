#include <sstream>
#include <stdexcept>
#include "anyoption.h"
#include "ProcessOptionsPE.h"
#include "Atmosphere.h"


#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

// constructor
NCPA::ProcessOptionsPE::ProcessOptionsPE(AnyOption *opt) {

  // defaults
  z_min            = 0.0;         // meters 
  maxrange         = 1.0E6;       // meters
  maxheight        = 150000.0;    // meters
  sourceheight     = z_min;       // meters
  receiverheight   = z_min;       // meters
  rng_step         = 1.0/10.0;    // wavelengths
  Nz_grid          = 20000;       // number of points on the z-grid	
  skiplines        = 0;           // skiplines in "atmosfile"
  grnd_imp_model   = "rigid";     // rigid ground
  wind_units       = "mpersec";   // m/s
  req_profile_step = maxrange;    // specifies the range step to request a new profile 
  starter_type     = "gaussian";  // starter field type
  usrattfile       = "";          // user-provided attenuation filename
  n_pade           = 4;
  do_lossless      = 0;           // flag; if=1 => no atmospheric absorption
  ncpatoy          = 0;

  plot_flg         = opt->getFlag( "plot"); // flag to plot results with gnuplot  
  
  // Parse arguments based on file type selected and set defaults
  // Declare and populate variables
  // First, file type
  enum AtmosphericFileType { ATMOSFILE, G2SENVFILE, NCPACANONICAL, ASCIIDIR, SLICEFILE};   // Expand this enum as we put in more file types
  //AtmosphericFileType filetype;
  int numTypesDeclared = 0;
  
  if ( opt->getValue( "use_1D_profiles_from_dir" ) != NULL ) {
      atm_profile_dir = opt->getValue( "use_1D_profiles_from_dir" );
      filetype = ASCIIDIR;
      numTypesDeclared++;
  }

  //if ( opt->getValue( "jetfile" ) != NULL ) {
  //    atmosfile = opt->getValue( "jetfile" );
  //    filetype = JETFILE;
  //    numTypesDeclared++;
  //}
  
  if ( opt->getValue( "atmosfile1d" ) != NULL ) {
      atmosfile = opt->getValue( "atmosfile1d" );
      filetype = ATMOSFILE;
      numTypesDeclared++;
      //cout << "Atm file: " << atmosfile << endl;
  }

  if (opt->getValue( "slicefile" ) != NULL ) {
      atmosfile = opt->getValue( "slicefile" );
      filetype = SLICEFILE;
      numTypesDeclared++;
      cout << "Slice file: " << atmosfile << endl;
  }
  
  if (opt->getValue( "g2senvfile" ) != NULL ) {
      atmosfile = opt->getValue( "g2senvfile" );
      filetype = G2SENVFILE;
      numTypesDeclared++;
  }
  
  if (opt->getValue( "ncpatoy" ) != NULL ) {
      //atmosfile = opt->getValue( "ncpatoy" );
      ncpatoy = 1;
      filetype = NCPACANONICAL;
      numTypesDeclared++;
      //cout << "Toy file: " << atmosfile << endl;
  }  

  // Make sure only one file type was requested
  if (numTypesDeclared != 1) {
      delete opt;
      throw invalid_argument( "Choose only one source for the atmospheric profiles: options --use_1D_profiles_from_dir or --g2senvfile or --atmosfile1d or --ncpatoy.");
  }
  
  if (filetype!=2) {
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
  
  if (opt->getValue( "rng_step" ) != NULL) {
      rng_step = atof( opt->getValue("rng_step") );
      //cout << "rng_step = " << rng_step << endl;
  }     

  if ( opt->getValue( "Nz_grid" ) != NULL ) {
      Nz_grid = atoi(opt->getValue( "Nz_grid" ));
      if (Nz_grid < 10) {
          //cout << "Nz_grid = " << Nz_grid << endl;
          delete opt;
          throw invalid_argument("Nz_grid value is too low (typically 20000)");
      }
  }	    
  
  if (opt->getValue( "starter_type" ) != NULL) {
      starter_type = opt->getValue( "starter_type" );
  }
  
  
  modstartfile = "none";
  if (! strcmp(starter_type.c_str(), "modal")) {
    if (opt->getValue( "modal_starter_file" ) != NULL) {
      modstartfile = opt->getValue( "modal_starter_file" );
    }
    else {
      delete opt;
      throw invalid_argument("Please provide the modal starter file name. Use option --modal_starter_file");
    }
  }
  
  if (opt->getValue( "n_pade" ) != NULL) {
      n_pade = atoi(opt->getValue( "n_pade" ));
  }  

  //string grnd_imp_model ("rigid");
  if ( opt->getValue( "ground_impedance_model" ) != NULL ) {
      //cout << opt->getValue( "ground_impedance_model" ) << endl;
      //grnd_imp_model.assign(opt->getValue( "ground_impedance_model" ));
      grnd_imp_model = opt->getValue( "ground_impedance_model" );
      
      //if (grnd_imp_model.compare("rigid")!=0) {
      //    //cout << "ground_impedance_model = " << grnd_imp_model << endl;
      //    std::ostringstream es;
      //    es << "This ground impedance model is not implemented yet: " << grnd_imp_model;
      //    delete opt;
      //    throw invalid_argument(es.str());
      //}   
  }	
  
  if ( opt->getValue( "use_profile_ranges_km" ) != NULL ) {
      profile_ranges_given = 1;
      prf_ranges_km = opt->getValue( "use_profile_ranges_km" );
      //cout << "prf_ranges_km= " << prf_ranges_km << endl;
  }

  if (opt->getValue( "use_profiles_at_steps_km" ) != NULL) {
      profile_ranges_given = 0;
      req_profile_step = atof( opt->getValue( "use_profiles_at_steps_km" ))*1000.0;
  }
  
  if ( opt->getValue( "use_attn_file" ) != NULL ) {
      usrattfile = opt->getValue( "use_attn_file" );
      cout << "User-provided attenuation file: " << usrattfile << endl;
  }  
  
  if ( opt->getValue( "wind_units" ) != NULL ) {
      wind_units = opt->getValue( "wind_units" );
      if (!( !wind_units.compare("mpersec") ||  !wind_units.compare("kmpersec") ) ) {
        delete opt;
        throw invalid_argument("Bad wind units: they can only be mpersec or kmpersec.");
      }
  }

  
  do_lossless    = opt->getFlag("do_lossless");
  write_2D_TLoss = opt->getFlag("write_2D_TLoss");
  
} // ------- End of processOptions -----------------------------------


// utility to print the parameters to the screen
void NCPA::ProcessOptionsPE::printParams() {
  printf("\n High-Angle PE run info:\n");
  printf("                   freq : %g\n", freq);
  printf("                azimuth : %g\n", azi);
  printf("                Nz_grid : %d\n", Nz_grid);
  printf("      z_min (meters MSL): %g\n", z_min);
  printf("      maxheight_km (MSL): %g\n", maxheight/1000.0);
  printf("   sourceheight_km (AGL): %g\n", sourceheight/1000.0);
  printf(" receiverheight_km (AGL): %g\n", receiverheight/1000.0);    
  printf("            maxrange_km : %g\n", maxrange/1000.0);  
  printf("        PE starter_type : %s\n", starter_type.c_str());
  printf("  N Pade coeffs (n_pade): %d\n", n_pade); 
  printf("         grnd_imp_model : %s\n", grnd_imp_model.c_str());
  if (do_lossless) {
  printf("         atm absorption : %s\n", "not considered");
  } else {
  printf("         atm absorption : %s\n", "yes");
  }
  printf("    write_2D_TLoss flag : %d\n", write_2D_TLoss);
  if (ncpatoy) {
  printf("use NCPA canonical prof : %d\n", ncpatoy);
  } else if (filetype==3) {
  printf("atmospheric profile dir : %s\n", atm_profile_dir.c_str());
  } else {
  printf("    atmospheric profile : %s\n", atmosfile.c_str());
  }
  printf("             wind_units : %s\n", wind_units.c_str());
  if (!usrattfile.empty()) {
  printf("  User attenuation file : %s\n", usrattfile.c_str());
  }
}

// set functions
void  NCPA::ProcessOptionsPE::setMaxheight(double newmax) {
  maxheight = newmax;
}

// "get" functions
int   NCPA::ProcessOptionsPE::getFiletype() {
  return filetype;
}

std::string   NCPA::ProcessOptionsPE::getAtmosfile() {
  return atmosfile;
}

std::string   NCPA::ProcessOptionsPE::getWindUnits() {
  return wind_units;
}

std::string   NCPA::ProcessOptionsPE::getAtmosfileorder() {
  return atmosfileorder;
}

std::string   NCPA::ProcessOptionsPE::getAtm_profile_dir() {
  return atm_profile_dir;
}

std::string   NCPA::ProcessOptionsPE::getUsrAttFile() {
  return usrattfile;
}

std::string   NCPA::ProcessOptionsPE::getStarterType() {
  return starter_type;
}

std::string   NCPA::ProcessOptionsPE::getModalStarterFile() {
  return modstartfile;
}

int   NCPA::ProcessOptionsPE::getSkiplines() {
  return skiplines;
}

double NCPA::ProcessOptionsPE::getFreq() {
  return freq;
}

double NCPA::ProcessOptionsPE::getAzimuth() {
  return azi;
}

double NCPA::ProcessOptionsPE::getMaxrange() {
  return maxrange;
}

double NCPA::ProcessOptionsPE::getMaxheight() {
  return maxheight;
}

double NCPA::ProcessOptionsPE::getZ_min() { // ground level above MSL
  return z_min;
}

double NCPA::ProcessOptionsPE::getSourceheight() {
  return sourceheight;
}

double NCPA::ProcessOptionsPE::getReceiverheight() {
  return receiverheight;
}

double NCPA::ProcessOptionsPE::getRngStep() {
  return rng_step;
}

int   NCPA::ProcessOptionsPE::getNz_grid() {
  return Nz_grid;
}

int   NCPA::ProcessOptionsPE::getNpade() {
  return n_pade;
}

std::string   NCPA::ProcessOptionsPE::getGrnd_imp_model() {
  return grnd_imp_model;
}

bool   NCPA::ProcessOptionsPE::getNoabsorption() {
  return do_lossless;
}

bool   NCPA::ProcessOptionsPE::getWrite_2D_TLoss() {
  return write_2D_TLoss;
}

std::string   NCPA::ProcessOptionsPE::getProfileRanges() {
  return prf_ranges_km;
}

bool   NCPA::ProcessOptionsPE::getProfile_ranges_given_flag() {
  return profile_ranges_given;
}

double NCPA::ProcessOptionsPE::getReq_profile_step() {
  return req_profile_step;
}

bool   NCPA::ProcessOptionsPE::getPlot_flg() {
  return plot_flg;
}

