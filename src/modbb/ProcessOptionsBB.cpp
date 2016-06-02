#include "anyoption.h"
#include "SolveModBB.h"
#include "ProcessOptionsBB.h"


#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;


// constructor
NCPA::ProcessOptionsBB::ProcessOptionsBB(AnyOption *opt)
{
  // handle options to write dispersion file(s)
  
  // defaults/initializations for some parameters
  azi            = 0.0;        // degrees
  maxrange       = 1.0E6;      // meters
  maxheight      = 150000;     // meters
  z_min          = 0.0;        // meters
  sourceheight   = z_min;      // meters
  receiverheight = z_min;      // meters
  Nz_grid        = 20000;      // number of points on the z-grid
  Nrng_steps     = 1000;       // number of range steps		
  skiplines      = 0;          // skiplines in "atmosfile"
  gnd_imp_model  = "rigid";    // rigid ground
  wind_units     = "mpersec";  // m/s
  Nfreq          = 1;          // a default number of frequencies
  srcfile        = "";
  f_center       = -1;  // must have negative initialization value; 
		                    // used later to decide if f_center is reset by input 
		                    // option used for built-in pulse
		                    
  NFFT           = -1;  // number of fft points; this initial negative value
                        // is necessary to signal the code later that in fact
                        // NFFT was not requested as an option; 
                        // NFFT will then change to 4*fmax/f_step 
                        // see function 'pulse_prop_src2rcv_grid2()' where this 
                        // change happens         

  plot_flg       = opt->getFlag( "plot"); // flag to plot results with gnuplot
  
  // dispersion: source to receiver (1D)
  w_disp_src2rcv_flg = 0;
  if ( opt->getValue( "out_disp_src2rcv_file" ) != NULL ) {
      w_disp_src2rcv_flg = 1;
      out_disp_file.assign(opt->getValue( "out_disp_src2rcv_file" ));
      cout << "out_disp_src2rcv_file = " << out_disp_file << endl;
  }

  // full 2D dispersion 
  w_disp_flg = 0;
  if ( opt->getValue( "out_dispersion_files" ) != NULL ) {
      w_disp_flg = 1;
      out_disp_file.assign(opt->getValue( "out_dispersion_files" ));
      cout << "out_dispersion_files = " << out_disp_file << endl;
  }
  
  // if computing dispersion ask which mode computation to be used: Modess or WMod
  if (((int)w_disp_src2rcv_flg + (int)w_disp_flg)>0) { 
      usemodess_flg = 0;
      if ( opt->getValue( "use_modess" ) != NULL ) {
              usemodess_flg = 1;
              //cout << "usemodess_flg = " << usemodess_flg << endl;
      }
      
      usewmod_flg = 0;
      if ( opt->getValue( "use_wmod" ) != NULL ) {
                usewmod_flg = 1;
                //cout << "usewmod_flg = " << usewmod_flg << endl;
      }

      if (((int)usemodess_flg+(int)usewmod_flg)!=1) {
          delete opt;
          throw invalid_argument( "Select one and only option: either: --use_modess or --use_wmod!" );
      }
  }
  
  if ( opt->getValue( "wind_units" ) != NULL ) {
      wind_units = opt->getValue( "wind_units" );
  }  

  // handle options for pulse propagation
  pprop_grid_flg = 0;
  if ( opt->getValue( "pulse_prop_grid" ) != NULL ) {
      pprop_grid_flg  = 1;
      pprop_grid_dirname.assign(opt->getValue( "pulse_prop_grid" ));
      //cout << "pulse_prop_grid dispersion files = " << pprop_grid_dirname << endl;
  }	

  pprop_s2r_flg = 0;
  if ( opt->getValue( "pulse_prop_src2rcv" ) != NULL ) {
      pprop_s2r_flg  = 1;
      pprop_s2r_disp_file.assign(opt->getValue( "pulse_prop_src2rcv" ));
      cout << "--> Source-to-ground dispersion file = " << pprop_s2r_disp_file << endl;
  }	

  pprop_s2r_grid_flg = 0;
  if ( opt->getValue( "pulse_prop_src2rcv_grid" ) != NULL ) {
      pprop_s2r_grid_flg  = 1;
      pprop_s2r_disp_file.assign(opt->getValue( "pulse_prop_src2rcv_grid" ));
      cout << "--> Source-to-ground dispersion file = " << pprop_s2r_disp_file << endl;
  }	

  // only one of these options allowed at a time
  int only1 = (int) (w_disp_flg + w_disp_src2rcv_flg + \
              pprop_grid_flg + pprop_s2r_flg + pprop_s2r_grid_flg);	
						
  //cout << "w_disp_flg          = " << w_disp_flg << endl;
  //cout << "w_disp_src2rcv_flg  = " << w_disp_src2rcv_flg << endl;
  //cout << "pprop_grid_flg      = " << pprop_grid_flg << endl;
  //cout << "pprop_s2r_flg       = " << pprop_s2r_flg << endl;
  //cout << "pprop_s2r_grid_flg  = " << pprop_s2r_grid_flg << endl;									

  if (only1==0) {
      cout << "\n\nChoose one and only one option and make sure you provide the file/dir name (if necessary): " << endl;
      cout << "--out_dispersion_files\n--out_disp_src2rcv_file\n--pulse_prop_src2rcv\n--pulse_prop_src2rcv_grid\n--pulse_prop_grid\n ...exiting." << endl;
      exit(1);
  }
  else if (only1>1) {
      cout << "\n\nAmbiguous request. Choose only one option: " << endl;
      cout << "--out_dispersion_files\n--out_disp_src2rcv_file\n--pulse_prop_src2rcv\n--pulse_prop_src2rcv_grid\n--pulse_prop_grid\n ...exiting." << endl;
      exit(1);
  }

  max_cel = 300.0; //default celerity
  if ( opt->getValue( "max_celerity" ) != NULL ) {
      max_cel = atof(opt->getValue( "max_celerity" ));
      cout << "max_celerity = " << max_cel << " m/s" << endl;
  }
  else {
      cout << "Using default max_celerity = " << max_cel << " m/s" << endl;
  }
  
  
  // what kind of source?
  src_flg = 0;
  if (pprop_grid_flg || pprop_s2r_flg || pprop_s2r_grid_flg) {
  
   		if ( opt->getValue( "get_impulse_resp" ) != NULL ) {
                src_flg = 0;
                cout << "impulse response requested ..." << endl;
		  }  

   		if ( opt->getValue( "use_builtin_pulse" ) != NULL ) {
                src_flg = 1;
                cout << "Using builtin_pulse ..." << endl;
		  } 
		  
   		if ( opt->getValue( "src_spectrum_file" ) != NULL ) {
                srcfile = opt->getValue( "src_spectrum_file" );
                src_flg = 2;
                cout << "source spectrum file =  " << srcfile << endl;
		  } 

   		if ( opt->getValue( "src_waveform_file" ) != NULL ) {
                srcfile = opt->getValue( "src_waveform_file" );
                src_flg = 3;
                cout << "source waveform file =  " << srcfile << endl;
		  }  
  }

  RR = 0.0;
  if (pprop_s2r_flg) {
      //make sure range_R_km is provided
      if ( opt->getValue( "range_R_km" ) != NULL ) {
                RR = atof(opt->getValue( "range_R_km" ))*1000.0;
                cout << "RR = " << RR << " meters" << endl;
      }
      else {
          delete opt;
          throw invalid_argument( "Option --range_R_km is required!" );
      }  		
  }

  R_start = RR; // R_start should always have the same value as RR initially
  R_end = DR = 0.0;
  if (pprop_s2r_grid_flg) {
		  //make sure R_start, R_end and DR are provided
   		if ( opt->getValue( "R_start_km" ) != NULL ) {
          R_start = atof(opt->getValue( "R_start_km" ))*1000.0;
          cout << "R_start_km = " << R_start/1000.0 << " km" << endl;
          if (R_start<=0) {
              delete opt;
              throw invalid_argument( "Option --R_start_km should be greater than 0.0!" );
          }                
		  }
		  else {
          delete opt;
          throw invalid_argument( "Option --R_start_km is required!" );
		  }
		
		  if ( opt->getValue( "R_end_km" ) != NULL ) {
          R_end = atof(opt->getValue( "R_end_km" ))*1000.0;
          cout << "R_end_km   = " << R_end/1000.0 << " km" << endl;
          if (R_end<=0) {
              delete opt;
              throw invalid_argument( "Option --R_end_km should be greater than 0.0!" );
          }             
		  }
		  else {
          delete opt;
          throw invalid_argument( "Option --R_end_km is required!" );
		  }
		
   		if ( opt->getValue( "DR_km" ) != NULL ) {
          DR = atof(opt->getValue( "DR_km" ))*1000.0;
          cout << "DR_km      = " << DR/1000.0 << " km" << endl;
          if (DR<=0) {
              delete opt;
              throw invalid_argument( "Option --DR_km should be greater than 0.0!" );
          }                
		  }
		  else {
          delete opt;
          throw invalid_argument( "Option --DR_km is required!" );
		  }  				
  }
  
  // the center frequency for the built-in pulse; handle it here
  if ( opt->getValue( "f_center" ) != NULL ) {
            f_center = atof(opt->getValue( "f_center" ));
            cout << "f_center     = " << f_center << " Hz" << endl;
            if (f_center < 0.0) {
                delete opt;
                throw invalid_argument( "f_center should be strictly positive (and <= f_max/5)!" );
            }
  }
  
  // the number of FFT points
  if ( opt->getValue( "nfft" ) != NULL ) {
            NFFT = (int) round(atof(opt->getValue( "nfft" )));
            cout << "NFFT     = " << NFFT << endl;   
            if (NFFT < 0.0) {
                delete opt;
                throw invalid_argument( "NFFT should be strictly a positive integer" );
            }
  }  

  //
  // logic to handle either computing the dispersion and modal values 
  // or to propagate a pulse
  //
  if (w_disp_src2rcv_flg || w_disp_flg) { //write dispersion	

      z_min = 0.0; //default value is at the ground

      // Parse arguments based on file type selected and set defaults
      // Declare and populate variables
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

      if ( opt->getValue( "slicefile" ) != NULL ) {
                atmosfile = opt->getValue( "slicefile" );
                filetype = SLICEFILE;
                numTypesDeclared++;
      }

      // Make sure only one file type was requested
      if (numTypesDeclared == 0) {
          delete opt;
          throw invalid_argument( "Please specify one atmospheric file type!");
      }
      else if (numTypesDeclared > 1) {
          delete opt;
          throw invalid_argument( "Ambiguous request: please specify only one atmospheric file type!");
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
      }
      else {
          delete opt;
          throw invalid_argument( "Option --skiplines is required for ATMOSFILE files!" );
      }		
		
      if (opt->getValue( "azimuth" ) != NULL) {
          azi = atof( opt->getValue("azimuth") );
          //cout << "azimuth = " << azi << endl;
      } else {
          delete opt;
          throw invalid_argument( "Option --azimuth is required!" );
      }	

      if (opt->getValue( "maxrange_km" ) != NULL) {
          maxrange = atof( opt->getValue( "maxrange_km" ))*1000.0;
          //cout << "maxrange = " << maxrange << " meters" << endl; // in meters
      }

      if (opt->getValue( "maxheight_km" ) != NULL) {
          maxheight = atof( opt->getValue( "maxheight_km" ))*1000.0;
          //cout << "maxheight = " << maxheight << " meters" << endl; // in meters
          if (maxheight > 200000.0) {
              delete opt;
              throw invalid_argument( "maxheight cannot exceed 200 km!" );
          }
      }    

      if (opt->getValue( "sourceheight_km" ) != NULL) {
          sourceheight = atof( opt->getValue( "sourceheight_km" ))*1000.0;
          //cout << "sourceheight = " << sourceheight << " meters" << endl; // in meters
          if (sourceheight > maxheight) {
              delete opt;
              throw invalid_argument( "sourceheight cannot exceed maxheight!" );
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

      // Number of points on the z-grid
      if ( opt->getValue( "Nz_grid" ) != NULL ) {
                Nz_grid = atoi(opt->getValue( "Nz_grid" ));
                //cout << "Nz_grid = " << Nz_grid << endl;
      }	
            
      if ( opt->getValue( "Nrng_steps" ) != NULL ) {
                Nrng_steps = atoi(opt->getValue( "Nrng_steps" ));
                //cout << "Nrng_steps = " << Nrng_steps << endl;
      }	
        

      //string gnd_imp_model ("rigid");
      if ( opt->getValue( "ground_impedance_model" ) != NULL ) {
                gnd_imp_model.assign(opt->getValue( "ground_impedance_model" ));
                //cout << "ground_impedance_model = " << gnd_imp_model << endl;
      }	

      //int Lamb_wave_BC = 0; // for rigid ground: if ==1 then admittance = -1/2*dln(rho)/dz
      if ( opt->getValue( "Lamb_wave_BC" ) != NULL ) {
                Lamb_wave_BC = atoi(opt->getValue( "Lamb_wave_BC" ));
                //cout << "Lamb_wave_BC = " << Lamb_wave_BC << endl;
      }	

      // f_step = frequency step 
      if ( opt->getValue( "f_step" ) != NULL ) {
                f_step = atof(opt->getValue( "f_step" ));
                //cout << "f_step = " << f_step << endl;
                if (f_step < 1.0e-16) {
                    delete opt;
                    throw invalid_argument(" f_step must be strictly positive.");
                }
      } 
      else {
          delete opt;
          throw invalid_argument( "Option --f_step is required!" );
      } 					

      //double f_min; minimum frequency
      if ( opt->getValue( "f_min" ) != NULL ) {
                f_min = atof(opt->getValue( "f_min" ));
                //cout << "f_min = " << f_min << endl;
      }
      else { 
          f_min = f_step;
      }

      //double f_max = 0;	// Nyquist frequency// set it ot f_step
      if ( opt->getValue( "f_max" ) != NULL ) {
          f_max = atof(opt->getValue( "f_max" ));
          //cout << "f_max = " << f_max << endl;
          if (f_max < f_min) {
              cout << "f_min = " << f_min << endl;
              cout << "f_max = " << f_max << endl;
              delete opt;
              throw invalid_argument("f_max must be greater than f_min");
          }
      }
      else {
          delete opt;
          throw invalid_argument( "Option --f_max is required!" );
      } 				
  } // End of dispersion options

	else if (pprop_s2r_flg || pprop_s2r_grid_flg) { //pulse propagation
      //make sure we have the waveform output file name
      if ( opt->getValue( "waveform_out_file" ) != NULL ) {
          waveform_out_file.assign(opt->getValue( "waveform_out_file"));
          // cout << "waveform_out_file = " << waveform_out_file << endl;
      }
      else {
          delete opt;
          throw invalid_argument( "Option --waveform_out_file is required!" );
      }	
  }
	else if (pprop_grid_flg) {	
		
		  // default values; most will be overidden by input options below
		  R_start_km = 1.0;
      width_km   = 50.0;
		  height_km  = maxheight/1000.0;
		  tmstep     = 30.0;
		  ntsteps    = 2;
		  frame_file_stub  = "Pressure2D";	
      //max_cel    = 300.0; // already defined	

		  //make sure R_start, R_end and DR are provided
   		if ( opt->getValue( "R_start_km" ) != NULL ) {
                R_start_km = atof(opt->getValue( "R_start_km" ));
                cout << "R_start_km   = " << R_start_km << " km" << endl;
		  }
		  else {
          delete opt;
          throw invalid_argument( "Option --R_start_km is required!" );
		  }
		
		  if ( opt->getValue( "width_km" ) != NULL ) {
                width_km = atof(opt->getValue( "width_km" ));
                cout << "width_km     = " << width_km << " km" << endl;
		  }
		  else {
          delete opt;
          throw invalid_argument( "Option --width_km (2D window width) is required!" );
		  }
		
		  if ( opt->getValue( "height_km" ) != NULL ) {
                height_km = atof(opt->getValue( "height_km" ));
                cout << "height_km    = " << height_km << " km" << endl;
		  }
		  else {
          cout << "default height_km = " << height_km << " km" << endl;
		  }		
		
   		if ( opt->getValue( "tmstep" ) != NULL ) {
                tmstep = atof(opt->getValue( "tmstep" ));
                cout << "tmstep       = " << tmstep << " seconds" << endl;
		  }
		  else {
          delete opt;
          throw invalid_argument( "Option --tmstep (time step (seconds)) is required!" );
		  }
		
   		if ( opt->getValue( "ntsteps" ) != NULL ) {
                ntsteps = atoi(opt->getValue( "ntsteps" ));
                cout << "ntsteps      = " << ntsteps << endl;
		  }
		  else {
          delete opt;
          throw invalid_argument( "Option --ntsteps (number of time steps) is required!" );
		  }	
		  
		  if ( opt->getValue( "frame_file_stub" ) != NULL ) {
                frame_file_stub.assign(opt->getValue( "frame_file_stub" ));
                cout << "file name stub for the pressure field: frame_file_stub = " \
                     << frame_file_stub << endl;
		  }
		  else {
          cout << "Default file name stub for the pressure field: frame_file_stub = " \
               << frame_file_stub << endl;
      }	
  }
}

//get functions
bool NCPA::ProcessOptionsBB::getW_disp_src2rcv_flg() {
  return w_disp_src2rcv_flg;
}

std::string   NCPA::ProcessOptionsBB::getOut_disp_file() {
  return out_disp_file;
}

bool   NCPA::ProcessOptionsBB::getW_disp_flg() {
  return w_disp_flg;
}

std::string   NCPA::ProcessOptionsBB::getPprop_grid_dirname() {
  return pprop_grid_dirname;
}

bool   NCPA::ProcessOptionsBB::getPprop_grid_flg() {
  return pprop_grid_flg;
}

std::string   NCPA::ProcessOptionsBB::getPprop_s2r_disp_file() {
  return pprop_s2r_disp_file;
}

bool   NCPA::ProcessOptionsBB::getPprop_s2r_flg() {
  return pprop_s2r_flg;
}

bool   NCPA::ProcessOptionsBB::getPprop_s2r_grid_flg() {
  return pprop_s2r_grid_flg;
}

bool   NCPA::ProcessOptionsBB::getUsemodess_flg() {
  return usemodess_flg;
}

double NCPA::ProcessOptionsBB::getZ_min() {
  return z_min;
}

double NCPA::ProcessOptionsBB::getMax_celerity() {
  return max_cel;
}

double NCPA::ProcessOptionsBB::getRR() {
  return RR;
}

double NCPA::ProcessOptionsBB::getR_start() {
  return R_start;
}

double NCPA::ProcessOptionsBB::getR_end() {
  return R_end;
}

double NCPA::ProcessOptionsBB::getDR() {
  return DR;
}

std::string   NCPA::ProcessOptionsBB::getAtmosfile() {
  return atmosfile;
}

std::string   NCPA::ProcessOptionsBB::getWindUnits() {
  return wind_units;
}

std::string   NCPA::ProcessOptionsBB::getAtmosfileorder() {
  return atmosfileorder;
}

int    NCPA::ProcessOptionsBB::getSkiplines() {
  return skiplines;
}

double NCPA::ProcessOptionsBB::getAzimuth() {
  return azi;
}

double NCPA::ProcessOptionsBB::getMaxrange() {
  return maxrange;
}

double NCPA::ProcessOptionsBB::getMaxheight() {
  return maxheight;
}

double NCPA::ProcessOptionsBB::getSourceheight() {
  return sourceheight;
}

double NCPA::ProcessOptionsBB::getReceiverheight() {
  return receiverheight;
}

int    NCPA::ProcessOptionsBB::getNrng_steps() {
  return Nrng_steps;
}

int    NCPA::ProcessOptionsBB::getNz_grid() {
  return Nz_grid;
}

std::string   NCPA::ProcessOptionsBB::getGnd_imp_model() {
  return gnd_imp_model;
}

int    NCPA::ProcessOptionsBB::getLamb_wave_BC() {
  return Lamb_wave_BC;
}

int    NCPA::ProcessOptionsBB::getNfreq() {
  return Nfreq;
}

double NCPA::ProcessOptionsBB::getF_min() {
  return f_min;
}

double NCPA::ProcessOptionsBB::getF_step() {
  return f_step;
}

double NCPA::ProcessOptionsBB::getF_max() {
  return f_max;
}

double NCPA::ProcessOptionsBB::getF_center() {
  return f_center;
}

std::string NCPA::ProcessOptionsBB::getWaveform_out_file() {
  return waveform_out_file;
}

double NCPA::ProcessOptionsBB::getR_start_km() {
  return R_start_km;
}

double NCPA::ProcessOptionsBB::getWidth_km() {
  return width_km;
}

double NCPA::ProcessOptionsBB::getHeight_km() {
  return height_km;
}

double NCPA::ProcessOptionsBB::getTmstep() {
  return tmstep;
}

int    NCPA::ProcessOptionsBB::getNtsteps() {
  return ntsteps;
}

std::string   NCPA::ProcessOptionsBB::getFrame_file_stub() {
  return frame_file_stub;
}

std::string   NCPA::ProcessOptionsBB::getSrcfile() {
  return srcfile;
}

int    NCPA::ProcessOptionsBB::getSrc_flg() {
  return src_flg;
}

bool   NCPA::ProcessOptionsBB::getPlot_flg() {
  return plot_flg;
}

int    NCPA::ProcessOptionsBB::getNFFT() {
  return NFFT;
}

	
			




