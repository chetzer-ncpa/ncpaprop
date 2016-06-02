#include <sstream>
#include <stdexcept>
#include "anyoption.h"
#include "ProcessOptionsNRT.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;

// constructor
NCPA::ProcessOptionsNRT::ProcessOptionsNRT(AnyOption *opt) {

  // defaults; may be changed by the command line options
  

 //z_min            = 0.0;            // meters 
  maxrange         = 1.0E6;          // meters
 // maxheight        = 150000.0;       // meters
  sourceheight     = z_min;          // meters
  receiverheight   = z_min;          // meters
 // Nz_grid          = 20000;          // number of points on the z-grid
 // Nrng_steps       = 1000;           // number of range steps	
  skiplines        = 0;              // skiplines in "atmosfile"
 // Lamb_wave_BC     = 0;              // 1 to enforce the Lamb wave BC
 // gnd_imp_model    = "rigid";        // rigid ground
//  wind_units       = "mpersec";      // m/s
 // usrattfile       = "";             // user-provided attenuation filename
 // tol              = 1.0E-08;        // tolerance for Slepc calculations
  

  write_atm_profile  = opt->getFlag( "write_atm_profile" );

  //plot_flg           = opt->getFlag( "plot"); // flag to plot results with gnuplot
  
  
  
  // default values
  src_z = 0;
  rcv_x = 0;
  rcv_y = 0;
  rcv_z = 0;
  range = 1.0E6;
  theta = 45;
  azi = 0;
  dth = 1;
  daz = 10;
  tol = 1000;
  
  wftype     = "Nwave";
  wfAmpl     = 1000.0;  // Pa
  wfDuration = 0.5;     // seconds

  findEigenray = opt->getFlag( "findeigenray");
  shootRay     = opt->getFlag( "shootray" );
  
  // Parse command line arguments 

  if  ((opt->getValue( "eigenrayfile" ) == NULL) &&  \
      (opt->getValue( "findeigenray" ) == NULL) &&  \
      (opt->getValue( "shootray" ) == NULL) ) {
      delete opt;
      throw invalid_argument("Please specify one of 3 options: eigenrayfile, findeigenray or shootray");
  }

  if (opt->getValue( "eigenrayfile" ) != NULL) {
    eigenfile = opt->getValue("eigenrayfile");
    findEigenray = false;
  } 
  
  if (opt->getValue( "ampl" ) != NULL) { 
    wfAmpl = atof(opt->getValue( "ampl" ));
  }
  
  if ( (opt->getValue( "duration" ) != NULL) ) { 
    wfDuration = atof(opt->getValue( "duration" ));
  }  

  if (opt->getValue( "waveform" ) != NULL) {
    wftype = opt->getValue( "waveform" );
    //if (wftype.compare("Nwave") || wftype.compare("nwave")) 
    if (opt->getValue( "ampl" ) != NULL) { 
      wfAmpl = atof(opt->getValue( "ampl" ));
    }
    else {
      delete opt;
      throw invalid_argument( "Option --ampl (related to waveform) is required!" );
    }
    
    //cout << "wftype = " << wftype << endl;
    //cout << "wftype.compare(""Nwave"") = " << wftype.compare("Nwave") << endl;
    
    if ((wftype.compare("Nwave")==0)) {
      if ( (opt->getValue( "duration" ) != NULL) ) { 
        wfDuration = atof(opt->getValue( "duration" ));
      }
      else {
        delete opt;
        throw invalid_argument( "Option --duration (related to waveform) is required!" );
      }
    }
     
  } else {
     cout << "--> Assuming default N wave as the source waveform" << endl;
     cout << "--> amplitude = " << wfAmpl << " Pa;  duration = " << wfDuration \
          << " sec" << endl;
  }

  if (findEigenray || shootRay) {

    if ( opt->getValue( "atmosfile" ) != NULL ) {
        atmosfile = opt->getValue( "atmosfile" );
    } else {
        delete opt;
        throw invalid_argument( "Option --atmosfile is required" );
    }	

    //if (opt->getValue( "atmosfileorder" ) != NULL) {
    //   atmosfileorder = opt->getValue( "atmosfileorder" );
    //} 
    //else {
    //    delete opt;
    //    throw invalid_argument( "Option --atmosfileorder is required!" );
    //}  
    
    if (opt->getValue( "src_z" ) != NULL) {
      src_z = atof( opt->getValue("src_z") )*1000;
    } else {
        delete opt;
        throw invalid_argument( "Option --src_z is required" );
    }	
    
        
    if ( opt->getValue( "inclin" ) != NULL ) {
      //cout << opt->getValue( "elevation" )  << endl;
      theta = atof(opt->getValue( "inclin" ));
      if ( (theta < -90) || (theta>90) ) {
        delete opt;
        throw invalid_argument("Inclination angle must be between -90 and 90 degrees");
      }
    }

    if (opt->getValue( "azimuth" ) != NULL) {
        azi = atof( opt->getValue("azimuth") );
        azi = 90 - azi; // back to trig angle
        //cout << "azimuth = " << azi << endl;
    } 
    //else {
    //    delete opt;
    //    throw invalid_argument( "Option --azimuth is required" );
    //}	
  }
   
  if (findEigenray) {
    if (opt->getValue( "rcv_x" ) != NULL) {
        rcv_x = atof( opt->getValue("rcv_x") )*1000;
    } else {
        delete opt;
        throw invalid_argument( "Option --rcv_x is required" );
    }	
    
    if (opt->getValue( "rcv_y" ) != NULL) {
        rcv_y = atof( opt->getValue("rcv_y") )*1000;
    } else {
        delete opt;
        throw invalid_argument( "Option --rcv_y is required" );
    }	   

    if (opt->getValue( "rcv_z" ) != NULL) {
        rcv_z = atof( opt->getValue("rcv_z") )*1000;
    } else {
        delete opt;
        throw invalid_argument( "Option --rcv_z is required" );
    } 
      
    if (opt->getValue( "dth" ) != NULL) {
        dth = atof( opt->getValue("dth") );
    }     
    
    if (opt->getValue( "daz" ) != NULL) {
        daz = atof( opt->getValue("daz") );
    } 
 
    if (opt->getValue( "tol" ) != NULL) {
        tol = atof( opt->getValue("tol") )*1000;
    }  
   
  }
  
   
  
   if (shootRay) {
      if (opt->getValue( "range" ) != NULL) {
          range = atof( opt->getValue("range") );
      } else {
          delete opt;
          throw invalid_argument( "Option --range is required" );
      }	   

      if (opt->getValue( "rcv_z" ) != NULL) {
          rcv_z = atof( opt->getValue("rcv_z") );
      } else {
          delete opt;
          throw invalid_argument( "Option --rcv_z is required" );
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


}


// "get" functions


double NCPA::ProcessOptionsNRT::getInclination() {
  return theta;
}

double NCPA::ProcessOptionsNRT::getAzimuth() {
  return azi;
}

double NCPA::ProcessOptionsNRT::getSrcz() {
  return src_z;
}

double NCPA::ProcessOptionsNRT::getRcvx() {
  return rcv_x;
}

double NCPA::ProcessOptionsNRT::getRcvy() {
  return rcv_y;
}

double NCPA::ProcessOptionsNRT::getRcvz() {
  return rcv_z;
}

double NCPA::ProcessOptionsNRT::getRange() {
  return range;
}

double NCPA::ProcessOptionsNRT::getDth() {
  return dth;
}

double NCPA::ProcessOptionsNRT::getDaz() {
  return daz;
}

double NCPA::ProcessOptionsNRT::getTol() {
  return tol; // tolerance for the eigeray to hit the target [km]
}


bool NCPA::ProcessOptionsNRT::getFindEigenray() {
  return findEigenray;
}

std::string   NCPA::ProcessOptionsNRT::getAtmosfile() {
  return atmosfile;
}

std::string   NCPA::ProcessOptionsNRT::getEigenfile() {
  return eigenfile;
}

std::string   NCPA::ProcessOptionsNRT::getWftype() {
  return wftype;
}

double NCPA::ProcessOptionsNRT::getWfampl() {
  return wfAmpl;
}

double NCPA::ProcessOptionsNRT::getWfduration() {
  return wfDuration;
}


















std::string   NCPA::ProcessOptionsNRT::getAtmosfileorder() {
  return atmosfileorder;
}


int    NCPA::ProcessOptionsNRT::getSkiplines() {
  return skiplines;
}

double NCPA::ProcessOptionsNRT::getFreq() {
  return freq;
}



double NCPA::ProcessOptionsNRT::getAzimuthStart() {
  return azi_min;
}

double NCPA::ProcessOptionsNRT::getAzimuthEnd() {
  return azi_max;
}

double NCPA::ProcessOptionsNRT::getAzimuthStep() {
  return azi_step;
}

double NCPA::ProcessOptionsNRT::getMaxrange() {
  return maxrange;
}

double NCPA::ProcessOptionsNRT::getMaxheight() {
  return maxheight;
}

double NCPA::ProcessOptionsNRT::getSourceheight() {
  return sourceheight;
}

double NCPA::ProcessOptionsNRT::getReceiverheight() {
  return receiverheight;
}



std::string   NCPA::ProcessOptionsNRT::getGnd_imp_model() {
  return gnd_imp_model;
}


