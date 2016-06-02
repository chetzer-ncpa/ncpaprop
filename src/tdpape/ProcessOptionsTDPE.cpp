#include <stdexcept>
#include "anyoption.h"
#include "ProcessOptionsTDPE.h"


#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000

using namespace NCPA;
using namespace std;


// constructor
NCPA::ProcessOptionsTDPE::ProcessOptionsTDPE(AnyOption *opt)
{

  // defaults/initializations for some parameters
  Nrng_steps     = 1000;       // number of range steps		
  Nfreq          = 1;          // a default number of frequencies
  srcfile        = "";
  f_center       = -1;  // must have negative initialization value; 
		                    // used later to decide if f_center is reset by input option
		                    // used for built-in pulse
  NFFT           = -1;  // number of fft points; this initial negative value
                        // is necessary to signal the code later that in fact
                        // NFFT was not requested as an option; 
                        // NFFT will then change to 4*fmax/f_step 
                        // see function 'pulse_prop_src2rcv_grid2()' where this 
                        // change happens 		                   

  plot_flg       = opt->getFlag( "plot"); // flag to plot results with gnuplot
  

  // handle options for pulse propagation

  pprop_s2r_flg = 0;
  if ( opt->getValue( "pulse_prop_src2rcv" ) != NULL ) {
      pprop_s2r_flg  = 1;
      pape_out_dir.assign(opt->getValue( "pulse_prop_src2rcv" ));
      cout << "Pre-computed pape output directory = " << pape_out_dir << endl;
  }

  pprop_s2r_grid_flg = 0;
  if ( opt->getValue( "pulse_prop_src2rcv_grid" ) != NULL ) {
      pprop_s2r_grid_flg  = 1;
      pape_out_dir.assign(opt->getValue( "pulse_prop_src2rcv_grid" ));
      cout << "Pre-computed pape output directory = " << pape_out_dir << endl;
  }


  // only one of these options allowed at a time
  int only1 = (int) (pprop_s2r_flg + pprop_s2r_grid_flg);	
  //cout << "pprop_s2r_flg       = " << pprop_s2r_flg << endl;
  //cout << "pprop_s2r_grid_flg  = " << pprop_s2r_grid_flg << endl;									

  if (only1==0) {
      cout << "\n\nPlease choose one (and only one) option: " << endl;
      cout << "--pulse_prop_src2rcv  - for propagation to a single receiver\n--pulse_prop_src2rcv_grid - for propagation to several equispaced receivers\n- ...exiting." << endl;
      delete opt;      
      exit(1);
  }
  else if (only1>1) {
      cout << "\n\nAmbiguous request. Choose only one option: " << endl;
      cout << "--pulse_prop_src2rcv\n--pulse_prop_src2rcv_grid\n ...exiting." << endl;
      delete opt;
      exit(1);
  }

  max_cel = 340.0; //default celerity
  if ( opt->getValue( "max_celerity" ) != NULL ) {
      max_cel = atof(opt->getValue( "max_celerity" ));
      cout << "max_celerity = " << max_cel << " m/s" << endl;
  }
  else {
      cout << "Using default max_celerity = " << max_cel << " m/s" << endl;
  }
  
  
  // what kind of source?
  src_flg = 0;
  if (pprop_s2r_flg || pprop_s2r_grid_flg) {
      int onlyone = 0;
   		if ( opt->getValue( "get_impulse_resp" ) != NULL ) {
                src_flg = 0;
                cout << "impulse response requested ..." << endl;
                onlyone++;
		  }  

   		if ( opt->getValue( "use_builtin_pulse" ) != NULL ) {
                src_flg = 1;
                cout << "Using builtin_pulse ..." << endl;
                onlyone++;
		  } 
		  
   		if ( opt->getValue( "src_spectrum_file" ) != NULL ) {
                srcfile = opt->getValue( "src_spectrum_file" );
                src_flg = 2;
                cout << "source spectrum file =  " << srcfile << endl;
                onlyone++;
		  } 

   		if ( opt->getValue( "src_waveform_file" ) != NULL ) {
                srcfile = opt->getValue( "src_waveform_file" );
                src_flg = 3;
                cout << "source waveform file =  " << srcfile << endl;
                onlyone++;
		  }
		  
		  if (onlyone>1) {
        delete opt;
        throw invalid_argument( "\nYou are requesting more than one source type. Please revise the command." );
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
            NFFT = (int) (atof(opt->getValue( "nfft" )));
            cout << "NFFT     = " << NFFT << endl;   
            if (NFFT < 0.0) {
                delete opt;
                throw invalid_argument( "NFFT should be strictly a positive integer" );
            }
  }  

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

//get functions




std::string   NCPA::ProcessOptionsTDPE::getPape_out_dir() {
  return pape_out_dir;
}

bool   NCPA::ProcessOptionsTDPE::getPprop_s2r_flg() {
  return pprop_s2r_flg;
}

bool   NCPA::ProcessOptionsTDPE::getPprop_s2r_grid_flg() {
  return pprop_s2r_grid_flg;
}

double NCPA::ProcessOptionsTDPE::getMax_celerity() {
  return max_cel;
}

double NCPA::ProcessOptionsTDPE::getRR() {
  return RR;
}

double NCPA::ProcessOptionsTDPE::getR_start() {
  return R_start;
}

double NCPA::ProcessOptionsTDPE::getR_end() {
  return R_end;
}

double NCPA::ProcessOptionsTDPE::getDR() {
  return DR;
}

double NCPA::ProcessOptionsTDPE::getMaxrange() {
  return maxrange;
}

int    NCPA::ProcessOptionsTDPE::getNrng_steps() {
  return Nrng_steps;
}


int    NCPA::ProcessOptionsTDPE::getNfreq() {
  return Nfreq;
}

double NCPA::ProcessOptionsTDPE::getF_min() {
  return f_min;
}

double NCPA::ProcessOptionsTDPE::getF_step() {
  return f_step;
}

double NCPA::ProcessOptionsTDPE::getF_max() {
  return f_max;
}

double NCPA::ProcessOptionsTDPE::getF_center() {
  return f_center;
}

std::string NCPA::ProcessOptionsTDPE::getWaveform_out_file() {
  return waveform_out_file;
}

double NCPA::ProcessOptionsTDPE::getR_start_km() {
  return R_start_km;
}

double NCPA::ProcessOptionsTDPE::getTmstep() {
  return tmstep;
}

int    NCPA::ProcessOptionsTDPE::getNtsteps() {
  return ntsteps;
}

std::string   NCPA::ProcessOptionsTDPE::getSrcfile() {
  return srcfile;
}

int    NCPA::ProcessOptionsTDPE::getSrc_flg() {
  return src_flg;
}

bool   NCPA::ProcessOptionsTDPE::getPlot_flg() {
  return plot_flg;
}

int    NCPA::ProcessOptionsTDPE::getNFFT() {
  return NFFT;
}

	
			




