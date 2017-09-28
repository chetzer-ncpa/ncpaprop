#include <sstream>
#include <stdexcept>
#include "anyoption.h"
#include "ProcessOptionsNB.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

// constructor
NCPA::ProcessOptionsNB::ProcessOptionsNB(AnyOption *opt) {

  // defaults; may be changed by the command line options
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
  usrattfile       = "";          // user-provided attenuation filename
  tol              = 1.0E-08;     // tolerance for Slepc calculations
  c_min            = 0.0;         // minimum sound speed requested by user to do wavenumber filtering
  c_max            = 0.0;         // maximum sound speed requested by user to do wavenumber filtering
  
  write_2D_TLoss     = opt->getFlag( "write_2D_TLoss");
  write_phase_speeds = opt->getFlag( "write_phase_speeds" );
  write_speeds       = opt->getFlag( "write_speeds" ); // save both phase and group speeds 
  write_modes        = opt->getFlag( "write_modes" );
  write_dispersion   = opt->getFlag( "write_dispersion" );
  write_atm_profile  = opt->getFlag( "write_atm_profile" );
  Nby2Dprop          = opt->getFlag( "Nby2Dprop");    //(N by 2D) propagation flag
  turnoff_WKB        = opt->getFlag( "turnoff_WKB" ); // if ==1 turns off the WKB least phase speed approx
  plot_flg           = opt->getFlag( "plot"); // flag to plot results with gnuplot
  
  // Parse arguments based on file type selected and set defaults
  // Declare and populate variables
  // First, file type
  enum AtmosphericFileType { ATMOSFILE, JETFILE, SLICEFILE };   // Expand this enum as we put in more file types
  AtmosphericFileType filetype;
  int numTypesDeclared = 0;

  if ( opt->getValue( "jetfile" ) != NULL ) {
      atmosfile = opt->getValue( "jetfile" );
      filetype = JETFILE;
      numTypesDeclared++;
  }
  
  if ( opt->getValue( "atmosfile" ) != NULL ) {
      atmosfile = opt->getValue( "atmosfile" );
      filetype = ATMOSFILE;
      numTypesDeclared++;
      //cout << "Atm file: " << atmosfile << endl;
  }

  if (opt->getValue( "slicefile" ) != NULL ) {
      atmosfile = opt->getValue( "slicefile" );
      filetype = SLICEFILE;
      numTypesDeclared++;
  }

  // Make sure only one file type was requested
  if (numTypesDeclared != 1) {
      delete opt;
      throw invalid_argument( "One and only one atmospheric file type must be specified!");
  }
  
  if (opt->getValue( "atmosfileorder" ) != NULL) {
      atmosfileorder = opt->getValue( "atmosfileorder" );
  } 
  else {
      delete opt;
      throw invalid_argument( "Option --atmosfileorder is required for ATMOSFILE files!" );
  }  
  
  // set skiplines to skip in the atmosfile
  //cout << "opt->getValue( ""skiplines"" ) = " << opt->getValue( "skiplines" ) << endl;
  if (opt->getValue( "skiplines" ) != NULL) {
      skiplines = atoi( opt->getValue( "skiplines" ) );
      //cout << "skiplines = " << skiplines << endl;
  } else {
      delete opt;
      throw invalid_argument( "Option --skiplines is required!" );
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

  if (!Nby2Dprop) {
      if (opt->getValue( "azimuth" ) != NULL) {
          azi = atof( opt->getValue("azimuth") );
          //cout << "azimuth = " << azi << endl;
      } else {
          delete opt;
          throw invalid_argument( "Option --azimuth is required (unless you intend to specify flag --Nby2Dprop!)" );
      }	
  }
  else { // (N by 2D) requested
      if (opt->getValue( "azimuth_start" ) != NULL) {
        azi_min = atof( opt->getValue("azimuth_start") );
        //cout << "azimuth_start = " << azi << endl;
      } else {
        delete opt;
        throw invalid_argument( "Option --azimuth_start is required!" );
      }	
      
      if (opt->getValue( "azimuth_end" ) != NULL) {
        azi_max = atof( opt->getValue("azimuth_end") );
        //cout << "azimuth_end = " << azi << endl;
      } else {
        delete opt;
        throw invalid_argument( "Option --azimuth_end is required!" );
      }	
      
      if (azi_min > azi_max) {
        throw invalid_argument( "Option --azimuth_start has to be less than azimuth_end!" );
      }
      
      if (opt->getValue( "azimuth_step" ) != NULL) {
        azi_step = atof( opt->getValue("azimuth_step") );
        if (azi_step==0) {
          throw invalid_argument( "azimuth_step must be non-zero numerical value - in degrees" );
        }
        //cout << "azimuth_step = " << azi << endl;
      } else {
        delete opt;
        throw invalid_argument( "Option --azimuth_step is required!" );
      }	           
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
  
  if ( opt->getValue( "wind_units" ) != NULL ) {
      wind_units = opt->getValue( "wind_units" );
  }
  
  if ( opt->getValue( "use_attn_file" ) != NULL ) {
      usrattfile = opt->getValue( "use_attn_file" );
      cout << "User-provided attenuation file: " << usrattfile << endl;
  }
  
  
  modal_starter_file = ""; // default value for modal_starter_file name
  if ( opt->getValue( "modal_starter_file" ) != NULL ) {
      modal_starter_file = opt->getValue( "modal_starter_file" );
  }

  // wavenumber filtering option
  wvnum_filter_flg = 0;

  if ( opt->getValue( "wvnum_filter" ) != NULL ) {
      wvnum_filter_flg = 1;
              //cout << "wvnum_filter_flg = " << wvnum_filter_flg << endl;
      if ( opt->getValue( "c_min" ) != NULL ) {
          c_min = atof(opt->getValue( "c_min" ));  
      }
      else {
          delete opt;
          throw invalid_argument( "For wavenumber filtering provide a minimum sound speed c_min.");
      }

      if ( opt->getValue( "c_max" ) != NULL ) {
          c_max = atof(opt->getValue( "c_max" ));  
      }
      else {
          delete opt;
          throw invalid_argument( "For wavenumber filtering provide a maximum sound speed c_max.");
      }
  }
  

}


// "get" functions
double NCPA::ProcessOptionsNB::getZ_min() {
  return z_min;
}

std::string   NCPA::ProcessOptionsNB::getAtmosfile() {
  return atmosfile;
}

std::string   NCPA::ProcessOptionsNB::getWindUnits() {
  return wind_units;
}

std::string   NCPA::ProcessOptionsNB::getAtmosfileorder() {
  return atmosfileorder;
}

std::string   NCPA::ProcessOptionsNB::getUsrAttFile() {
  return usrattfile;
}

std::string   NCPA::ProcessOptionsNB::getModalStarterFile() {
  return modal_starter_file;
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

double NCPA::ProcessOptionsNB::getAzimuthStart() {
  return azi_min;
}

double NCPA::ProcessOptionsNB::getAzimuthEnd() {
  return azi_max;
}

double NCPA::ProcessOptionsNB::getAzimuthStep() {
  return azi_step;
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

double NCPA::ProcessOptionsNB::getC_min() {
  return c_min;
}

double NCPA::ProcessOptionsNB::getC_max() {
  return c_max;
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

int    NCPA::ProcessOptionsNB::getLamb_wave_BC() {
  return Lamb_wave_BC;
}

bool   NCPA::ProcessOptionsNB::getWrite_2D_TLoss() {
  return write_2D_TLoss;
}

bool   NCPA::ProcessOptionsNB::getWrite_phase_speeds() {
  return write_phase_speeds;
}

bool   NCPA::ProcessOptionsNB::getWrite_speeds() {
  return write_speeds;
}

bool   NCPA::ProcessOptionsNB::getWrite_modes() {
  return write_modes;
}

bool   NCPA::ProcessOptionsNB::getWrite_dispersion() {
  return write_dispersion;
}

bool   NCPA::ProcessOptionsNB::getNby2Dprop() {
  return Nby2Dprop;
}

bool   NCPA::ProcessOptionsNB::getWriteAtmProfile() {
  return write_atm_profile ;
}

bool   NCPA::ProcessOptionsNB::getTurnoff_WKB() {
  return turnoff_WKB;
}

bool   NCPA::ProcessOptionsNB::getPlot_flg() {
  return plot_flg;
}

bool   NCPA::ProcessOptionsNB::getWvnum_filter_flg() {
  return wvnum_filter_flg;
}

