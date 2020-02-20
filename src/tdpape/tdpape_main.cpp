
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <vector>

#include <dirent.h>
#include <list>
//#include <stdio.h>
//#include <math.h>
#include <complex>
//#include <ctime>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <fftw3.h>

#include "anyoption.h"
#include "ProcessOptionsTDPE.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif

#define FFTN 16*1024
#define MAX_MODES 4000

// to compile
// g++ -o tdpape1 anyoption.cpp ProcessOptionsTDPE.cpp tdpape1.cpp -lgsl -lgslcblas -lfftw3

// to run
//./tdpape1 --pulse_prop_src2rcv /home/doru/infra/Cpp_code/DoruV/NCPA_Infra_20131114_no_gnuplot/samples/PapeBB_examp3_0.5Hz_canonic --range_R_km 280 --waveform_out_file mywf.dat

// or
// ./tdpape1 --pulse_prop_src2rcv_grid /home/doru/infra/Cpp_code/DoruV/NCPA_Infra_20131114_no_gnuplot/samples/PapeBB_examp3_0.5Hz_canonic  --waveform_out_file mywf.dat --R_start_km 200 --R_end_km 230 --DR_km 10



using namespace NCPA;
using namespace std;

// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

int getFile_list(string dir, std::list<string> &files, string pattern);

vector<double> getFreq_vector(string pape_out_dir, string filepattern);

// comparison, freq in filename.
bool compare_freq (string first, string second);

complex<double> getInterp(string filename, double x);

vector < complex<double> > interpDir(string pape_output_dir, string filepattern, double R);

int pulse_prop_src2rcv_grid2(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, double f_step, double *f_vec, \
					double f_center, vector< complex<double> > PP,  \
					int src_flg, string srcfile, int pprop_src2rcv_flg);
					
int pulse_prop_src2rcv_grid3(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, double f_step, double *f_vec, \
					double f_center, \
					string pape_out_dir, string filepattern, \
					int src_flg, string srcfile, int pprop_src2rcv_flg);

//20151020 DV: added NFFT as argument
int pulse_prop_src2rcv_grid4(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int NFFT, int n_freqs, double f_step, double *f_vec, \
					double f_center, \
					string pape_out_dir, string filepattern, \
					int src_flg, string srcfile, int pprop_src2rcv_flg);								
					
void fft_pulse_prop(\
          double t0, int n_freqs, double df, double *f_vec, \
          vector< complex<double> > PP, \
          complex<double> *dft_vec, complex<double> *pulse_vec);
                    					
int get_source_spectrum( \
					int n_freqs, double f_step, double *f_vec, double f_center, \
					complex<double> *dft_vec, complex<double> *pulse_vec, \
					complex<double> *arg_vec, int src_flg, string srcfile);
					
// DV 20151019: Added NFFT as argument; normalizes the builtin pulse
// the source spectrum is loaded/computed in this function
int get_source_spectrum( \
								int n_freqs, int NFFT, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile);					
								
double half_hann(int begin,int end,int i);

complex<double> pulse_spec_fit(double scale, double x);

int load_source_pulse_td(string srcpulsetdfn, vector<double> &t, vector<double> &tdp );

int load_source_spectrum(string srcspfn, double *freqv, complex<double> *dft_vec, int fftn);



								
int main ( int argc, char **argv )  {


  // parse options from the command line as well as an options file
  AnyOption *opt = parseInputOptions( argc, argv ); 

  // object to process the options
  ProcessOptionsTDPE *oTDPE = new ProcessOptionsTDPE(opt);
  

  int src_flg, NFFT;
  double max_cel, R_start, DR, R_end, f_center; 
  string waveform_out_file, pape_output_dir;
  string src_file;

  // some flags defining pulse propagation to single receiver or a series
  // or equispaced receivers
  bool   pprop_s2r_flg;
  //bool   pprop_s2r_grid_flg;    // CHH 191029 Unused
  bool   plt_flg;


  // get propagation flags
  pprop_s2r_flg        = oTDPE->getPprop_s2r_flg();
  //pprop_s2r_grid_flg   = oTDPE->getPprop_s2r_grid_flg();
  
  
 
//  f_center = 0.099;
  NFFT     = oTDPE->getNFFT();
  max_cel  = oTDPE->getMax_celerity();     
  R_start  = oTDPE->getR_start();
  R_end    = oTDPE->getR_end();
  DR       = oTDPE->getDR();
  f_center = oTDPE->getF_center();
  src_flg  = oTDPE->getSrc_flg();
  src_file = oTDPE->getSrcfile();
  plt_flg  = oTDPE->getPlot_flg();    
  pape_output_dir   = oTDPE->getPape_out_dir();
  waveform_out_file = oTDPE->getWaveform_out_file();
  
  cout << "max_cel = " << max_cel << endl;
  cout << "R_start = " << R_start << endl;
  cout << "R_end   = " << R_end << endl;
  cout << "DR      = " << DR    << endl;
  cout << "src_flg = " << src_flg << endl;
  
  // interpolate pape output
  //string filepattern ("tloss");
  string filepattern ("papeTL");
  
  //vector< complex<double> > P;
  //P = interpDir(pape_output_dir, filepattern, R_start);



/*
  list<string> files;
  list<string>::iterator it;
  string filesep = "/";

  //string pape_output_dir = "/home/doru/infra/Cpp_code/DoruV/NCPA_Infra_20131114_no_gnuplot/samples/PapeBB_examp3_0.5Hz_canonic/"; //".";

  // get and sort the files (they have the frequency in the name)
  getFile_list(pape_output_dir, files, filepattern);
  files.sort();
 
  if (0) { //print the sorted file list
      cout << "sorted list off files with pattern " << filepattern << " contains:" << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }
  */


  // read PP
  vector< complex<double> > PP;
  if (0) {
    double dat1, dat2;   
    complex<double> c;    
    complex<double> I (0.0, 1.0);
    ifstream indata;
      
	  indata.open( "Ppe.dat" ); 	// opens the file
	  if(!indata)  				        // file couldn't be opened
	  { 
      std::ostringstream es;
      es << "file << " << "Ppe.dat" << " could not be opened.";
      throw invalid_argument(es.str());
	  }	
	  indata >> dat1 >> dat2;  // read the first line 
	  //cout << dat1 << "  " << dat2 << endl;	

	  int ii = 1;
	  while ( !indata.eof() ) 		// keep reading until end-of-file
	  {
	    //c = dat1 + I*dat2;
	    c = dat2 - I*dat1; // apparently we need -i*imag(PE) to agree with normal modes Pressure
	    PP.push_back(c);
		  indata >> dat1 >> dat2;
		  ii = ii + 1;
		  //printf("it=%d\n", ii);
		  //cout << c << endl;
	  }
	  indata.close();
	  printf("Interpolated values loaded from 'Ppe.dat'\n");

  } 
  
  //cout << PP[0] << endl;
  
  
  vector<double> fv;
  double f_step;    // CHH 191029: fmax unused
  fv = getFreq_vector(pape_output_dir, filepattern);
  int Nfreq = fv.size();
  f_step = fv[1]-fv[0];

  // print frequency vector
  if (0) {
  for (int i=0; i<Nfreq; i++) {
    cout << fv[i] << endl;
  }
  }

	// all set to propagate the pulse				
	
  // This call uses FFTN from #define FFTN
  // pulse_prop_src2rcv_grid3( waveform_out_file.c_str(), max_cel, 
  //                           R_start, DR, R_end, Nfreq, f_step, fv.data(), 
  //  		               f_center, 
  //			       pape_output_dir, filepattern, 
  //			       src_flg, src_file, pprop_s2r_flg);

	// This call uses NFFT as an argument			                    
  pulse_prop_src2rcv_grid4( waveform_out_file.c_str(), max_cel, \
								            R_start, DR, R_end, NFFT, Nfreq, \
								            f_step, fv.data(), f_center, \
				                    pape_output_dir, filepattern, \
				                    src_flg, src_file, pprop_s2r_flg);				                    										          
  
  // (gnu)plot results if requested; calls a bash script
  // @todo get rid of this and its associated options
  if (plt_flg) {
    char buffer [256];
    sprintf(buffer, "%s %s", "./xPlotMods.script plotpulse ", waveform_out_file.c_str());
    printf ("Executing command %s\n", buffer);
    system (buffer);
    //printf ("The value returned was: %d.\n",i);
  }									          						                    
                    
  //delete[] f_vec;
  delete oTDPE;
  delete opt;                 
   
  return 0;
} // end of "main"


AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();

  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "|                             NCPA Infrasound                              |" );  	
  opt->addUsage( "|                 tdpape: Time Domain Parabolic Equation                   |" );
  opt->addUsage( "|                   Synthesizes waveforms based on pape                    |" );
  opt->addUsage( "|                                   2015                                   |" );  
  opt->addUsage( "----------------------------------------------------------------------------" );	
  opt->addUsage( "Usage: " );
  opt->addUsage( "" );
  opt->addUsage( "The options below can be specified at the command line or in a colon-separated" );
  opt->addUsage( "file \"tdpape.options\". Command-line options override file options." ); 
  opt->addUsage( "Be sure to precede all options with two minuses (--)." );
  opt->addUsage( "" );  
  opt->addUsage( " --help -h                Print this message and exit" );
  opt->addUsage( "" );
  opt->addUsage( "" );
  opt->addUsage( "To propagate a pulse (waveform), 2 steps must be completed:");
  opt->addUsage( " 1. Pre-computed single-frequency output files from pape must be available" );
  opt->addUsage( "    and stored in the same directory. If N frequencies were computed" );
  opt->addUsage( "    the template for the file names is" );
  opt->addUsage( "    <#n>papeTL<freq>; e.g. 004_papeTL_0.87776" );
  opt->addUsage( "    Please refer to script xrun_papeBB.sh for an example on how to" );
  opt->addUsage( "    compute single-frequency files in batch mode. " );  
  opt->addUsage( "" );  
  opt->addUsage( " 2. Perform pulse propagation for 2 scenarios:");
  opt->addUsage( "    a. source-to-one-receiver at one range (see option --pulse_prop_src2rcv)");	
  opt->addUsage( "    b. source-to-several-receivers at equally spaced ranges (i.e. on a grid)" );
  opt->addUsage( "       (see option --pulse_prop_src2rcv_grid)");
  opt->addUsage( "" );
  opt->addUsage( " Output:  " );
  opt->addUsage( "    The text output file has the column order:" ); 
  opt->addUsage( "    for option --pulse_prop_src2rcv: | Time (seconds) | Waveform | " );
  opt->addUsage( "    for option --pulse_prop_src2rcv_grid:" );
  opt->addUsage( "    | Range [km] | Time [seconds] | Waveform |" );
  opt->addUsage( "    Here 'Range' changes to reflect all receiver ranges on the defined grid." );
  opt->addUsage( "" );
  opt->addUsage( "" );
  opt->addUsage( " Four types of sources are available:" );
  opt->addUsage( "       delta function              -> see option --get_impulse_resp" );  
  opt->addUsage( "       built-in pulse              -> see option --use_builtin_pulse" );
  opt->addUsage( "       user-provided spectrum file -> see option --src_spectrum_file" );
  opt->addUsage( "       user-provided waveform file -> see option --src_waveform_file" );
  opt->addUsage( "" );
  opt->addUsage( "" ); 
  opt->addUsage( "Options:" );
  opt->addUsage( " --pulse_prop_src2rcv <directory name of pre-computed pape files> ");
  opt->addUsage( "                    Propagate pulse from source to 1 receiver");
  opt->addUsage( "                    at a distance specified by option --range_R_km; " );
  opt->addUsage( " --range_R_km       Propagate pulse to this range [km]" );
  opt->addUsage( " --waveform_out_file <waveform filename>   Name of the waveform output file." );
  opt->addUsage( "" );
  opt->addUsage( " --pulse_prop_src2rcv_grid <directory name of pre-computed pape files>");
  opt->addUsage( "                    Propagate pulse from source to array of ");
  opt->addUsage( "                    horizontally equally-spaced receivers" );  
  opt->addUsage( "" );
  opt->addUsage(" REQUIRED additional options:" );
  opt->addUsage( " --R_start_km       Propagation from this range to R_end_km in DR_km steps." );
  opt->addUsage( " --R_end_km         Pulse is propagated from R_start_km to this range." );
  opt->addUsage( " --DR_km            Range step to propagate from R_start_km to R_end_km." );
  //opt->addUsage( " --waveform_out_file <waveform filename> ");
  //opt->addUsage( "                    Name of the waveform output file." );
  opt->addUsage( "" );
  opt->addUsage( " OPTIONAL [defaults]:" );
  opt->addUsage( " --max_celerity     Maximum celerity [300 m/s]." );
  opt->addUsage( " --f_center         The center frequency of the pulse; must be <= [f_max/5]." );
  opt->addUsage( " --nfft             Number of points used in the FFT computation. ");
  opt->addUsage( "                    Defaults to [4*f_max/f_step]." );	
  opt->addUsage( "" );
  opt->addUsage( "SOURCE TYPE options: Use one of the following 4 options to specify the source:" );
  opt->addUsage( " --get_impulse_resp       Flag to use a (band-limited) delta function as source" );
  opt->addUsage( "                          and to output the impulse response." );
  opt->addUsage( "                          (this is the default)." );  
  
  opt->addUsage( " --use_builtin_pulse      Flag to request the use of the built-in source pulse." );
  opt->addUsage( "                          Note: Use --f_center to request the central frequency" ); 
  opt->addUsage( "                          of the pulse. f_center is restricted to a maximum" );
  opt->addUsage( "                          value of fmax/5 where fmax is the maximum frequency" ); 
  opt->addUsage( "                          defined by the dispersion file." );
  opt->addUsage( "                          The input waveform and spectrum are also saved for the" );
  opt->addUsage( "                          user's reference such that:" ); 
  opt->addUsage( "                          The built-in source spectrum is outputted in the file." );
  opt->addUsage( "                           'source_spectrum_input.dat' with the format." );
  opt->addUsage( "                             | Freq (Hz) | Re(S) | Imag(S) |." );
  opt->addUsage( "                          The input source waveform is outputted in the file." );
  opt->addUsage( "                           'source_waveform_input.dat' with the format." );
  opt->addUsage( "                             | Time [s] | Amplitude |." ); 
  
    
  opt->addUsage( " --src_spectrum_file      Specify the file name of the source spectrum");
  opt->addUsage( "                          at positive frequencies. The file must have 3 columns" );
  opt->addUsage( "                             | Freq | Real(Spectrum) | Imag(Spectrum) |" );
  opt->addUsage( " --src_waveform_file      Specify the file name of the user-provided " );
  opt->addUsage( "                          source waveform. The file must have 2 columns" );
  opt->addUsage( "                             |Time | Amplitude |" ); 
  opt->addUsage( "   If none of the source type options are specified the delta function source");
  opt->addUsage( "   is the default i.e. the output is the impulse response." );  
  opt->addUsage( "" );
  opt->addUsage( " QUICK-START EXAMPLES (run from the 'samples' directory):" );
  opt->addUsage( " (Assume that the pre-computed single frequency files reside in myTDPape_dir.) " );
  opt->addUsage( "" );  
  opt->addUsage( " Example 1: Pulse propagation to a point on the ground at range_R_km" ); 
  opt->addUsage( "            and output the impulse response:" );
  opt->addUsage( "" );
  //opt->addUsage( "   a. Compute dispersion file that will be used to compute the pressure pulse at 1 receiver. Assume that we want to end up with a pulse having a spectrum with a maximum frequency of f_max=0.5 Hz. Also assume that we want the pulse represented on a time record of T=512 seconds. The number of positive frequencies necessary for the calculation is T*f_max = 256 i.e.256 frequencies between 0 and 0.5 Hz. Thus we know f_max=0.5 Hz and f_step=f_max/256=0.001953125 Hz. The corresponding run command is:" );

  opt->addUsage( "    ../bin/tdpape --pulse_prop_src2rcv myTDPape_dir --range_R_km 240 --waveform_out_file mywavf.dat --max_celerity 320 --get_impulse_resp" );
  opt->addUsage( "" );
  opt->addUsage( " Example 2: Pulse propagation to a point on the ground at range_R_km" );
  opt->addUsage( "            and employ the built-in source pulse:" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/tdpape --pulse_prop_src2rcv myTDPape_dir --range_R_km 240 --waveform_out_file mywavef.dat --max_celerity 320 --use_builtin_pulse" );
  opt->addUsage( "" );
  opt->addUsage( " Example 3: Pulse propagation to several points on the ground 20 km apart" );
  opt->addUsage( "            and employ the user-provided source waveform:" );
  opt->addUsage( "" );  
  opt->addUsage( "    ../bin/tdpape --pulse_prop_src2rcv_grid  myTDPape_dir  --R_start_km 200 --R_end_km 240 --DR_km 20 --waveform_out_file mywavef.dat --max_celerity 320 --src_waveform_file source_waveform_input_example.dat" );  
  opt->addUsage( "" );   
  opt->addUsage( " Example 4: Pulse propagation to a point on the ground at range_R_km" );
  opt->addUsage( "            and employ the user-provided source spectrum:" );
  opt->addUsage( "" );
  opt->addUsage( "   ../bin/tdpape --pulse_prop_src2rcv myTDPape_dir --range_R_km 240 --waveform_out_file mywavf.dat --max_celerity 320 --src_spectrum_file source_spectrum_example.dat" );    
  opt->addUsage( "" );
  opt->addUsage( "" ); 
  opt->addUsage( "" );
  


  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  opt->setFlag( "get_impulse_resp" );
  opt->setFlag( "use_builtin_pulse" );
  opt->setFlag( "plot");
  
  
  opt->setOption( "maxrange_km" );
  opt->setOption( "stepsize" );
  opt->setOption( "Nz_grid" );
  opt->setOption( "Nrng_steps" );
  opt->setOption( "out_TL_2D" );
  opt->setOption( "ground_impedance_model" );
  opt->setOption( "Lamb_wave_BC" );
  opt->setOption( "f_min" );
  opt->setOption( "f_step" );	
  opt->setOption( "f_max" );
  opt->setOption( "f_center" );
  opt->setOption( "pulse_prop_grid" );
  opt->setOption( "pulse_prop_src2rcv" );
  opt->setOption( "pulse_prop_src2rcv_grid" );
  opt->setOption( "R_start_km" );
  opt->setOption( "R_end_km" );
  opt->setOption( "DR_km" );
  opt->setOption( "max_celerity" );
  opt->setOption( "range_R_km" );
  //opt->setOption( "out_dispersion_files" );
  opt->setOption( "out_disp_src2rcv_file" );
  opt->setOption( "waveform_out_file" );
  opt->setOption( "width_km" );
  opt->setOption( "height_km" );
  opt->setOption( "tmstep" );
  opt->setOption( "ntsteps" );
  opt->setOption( "frame_file_stub" );
  opt->setOption( "disp_pape_output_dirame" );
  opt->setOption( "src_spectrum_file" );
  opt->setOption( "src_waveform_file" );
  opt->setOption( "nfft" );

  // Process the command-line arguments
  opt->processFile( "./tdpape.options" );
  opt->processCommandArgs( argc, argv );

  if( ! opt->hasOptions()) { // print usage if no options
      opt->printUsage();
      delete opt;
      exit( 1 );
  }

  // Check to see if help text was requested
  if ( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) {
      opt->printUsage();
      exit( 1 );
  }

  return opt;
}
















// -----------------------------------------------------------------


complex<double> getInterp(string filename, double x) {
  double dat1, dat2, dat3, xintrp, xintrp2;
  vector <double> xx, p1, p2;
  gsl_interp_accel *acc_;
  gsl_spline *splin; 
	//char line[512]; // to store the headerline text
	
	ifstream indata;
    
	indata.open( filename.c_str() ); 	// opens the file
	if(!indata)  				// file couldn't be opened
	{ 
		//cerr << "Error: File " << filename << " could not be opened!" << endl;
    //    exit(1);
        
	  std::ostringstream es;
    es << "Error: File " << filename << " could not be opened!";
    throw invalid_argument(es.str());        
	}	
	// indata.getline(line, 512);

	indata >> dat1 >> dat2 >> dat3;  // read the first line 
	//cout << dat1 << "  " << dat2 << endl;	

	int ii = 1;
	while ( !indata.eof() ) 		// keep reading until end-of-file
	{
	  // in SI units
	  xx.push_back(dat1);
	  p1.push_back(dat2);
	  p2.push_back(dat3);

		indata >> dat1 >> dat2 >> dat3;
		ii = ii + 1;
		//printf("it=%d\n", ii);
	}
	indata.close();

	// abort if the the value to interpolate at is outside the available range
	if (dat1<=x) {
    std::ostringstream es;
    es << "Error: cannot interpolate at requested range R = " << x 
       << " km. Maximum range is Rmax = " << dat1 << " km.";
    throw invalid_argument(es.str());
	}

  // do interpolation
  acc_  = gsl_interp_accel_alloc();
  splin = gsl_spline_alloc(gsl_interp_cspline, xx.size());
  //printf("splin[%d]=%p\n", j, splin[j]);
  
  // interpolate the real part
  gsl_spline_init(splin, xx.data(), p1.data(), xx.size());
  xintrp  = gsl_spline_eval(splin, x, acc_ );
  
  // interpolate the imaginary part
  gsl_spline_init(splin, xx.data(), p2.data(), xx.size());
  xintrp2 = gsl_spline_eval(splin, x, acc_ );
  //cout << "x interp = " << xintrp << endl;
  //cout << "x interp2 = " << xintrp2 << endl;
       
  complex<double> pp (xintrp, xintrp2);
  //cout << " " << pp << endl; 
  
  delete acc_;
  delete splin;

  return pp;
}





vector < complex<double> > interpDir(string pape_output_dir, string filepattern, double R) {

// interpolates data at range R from single-frequency files residing in 
// directory pape_output_dir. The files are found based on a pattern, 
// say, "*papeTL*"

  vector< complex<double> > P;
  list<string> files;
  list<string>::iterator it;
  string filesep = "/";
  //string pattern ("tloss");
  //string pape_output_dir = "/home/doru/infra/Cpp_code/DoruV/NCPA_Infra_20131114_no_gnuplot/samples/PapeBB_examp3_0.5Hz_canonic/"; //".";

  // get and sort the files (they have the frequency in the name)
  getFile_list(pape_output_dir, files, filepattern);
  files.sort();

  if (1) { //print the sorted file list
      cout << "Sorted list of files with pattern " << filepattern << " contains:" << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }
  
  
  if (1)  {
    complex<double> pp;
    string fullname;

    if (1) { //assemble interpolated values from the sorted file list
      //cout << "sorted list of files with pattern " << filepattern << " contains:" << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
        fullname = pape_output_dir + filesep + (*it);
        pp = getInterp(fullname, R); // interpolate here
        P.push_back(pp);
        cout << *it << "  " << pp << endl;
      }  
    }
    cout << endl;
    
    
    // show interpolated values
    if (0) {
      for (unsigned int j=0; j<P.size(); j++) {
        printf("%3d  P = %g  %g \n", j, real(P[j]), imag(P[j]));
      } 
    }

    // save P
    if (0) {
        FILE *fp = fopen("Ppe.dat", "w");
        for (unsigned int ii=0; ii<P.size(); ii++) {
            fprintf(fp, "%18.12e  %18.12e\n", real(P[ii]), imag(P[ii]));
        }
        printf("Interpolated values saved in 'Ppe.dat'\n");
        fclose(fp);
    }
    
  }
  
  return P;

}











int getFile_list(string dir, list<string> &files, string pattern)
{
  int pos = -1;
  string a;
  DIR *dp;
  struct dirent *dirp;
  if((dp = opendir(dir.c_str())) == NULL) {
      std::ostringstream es;
      es << "Error opening directory:" << dir;
      throw invalid_argument(es.str());
  }

  while ((dirp = readdir(dp)) != NULL) {
      a = string(dirp->d_name);
      pos = a.find(pattern);
      if (pos>=0) {
          //cout << a << " pos=" << pos << endl;
          files.push_back(string(dirp->d_name));
      }
      else {
          //cout << "this file does not fit pattern: " << a << endl;
      }
  }
  closedir(dp);
  return 0;
}



vector<double> getFreq_vector(string pape_out_dir, string filepattern) {

    vector<double> fv;
    int Nfreq;
    double fmax, f_step;
    size_t pos1;
    string a,s;
    
    list<string> files;
    list<string>::iterator it;
    string filesep = "/";

    //string pape_output_dir = "/home/doru/infra/Cpp_code/DoruV/NCPA_Infra_20131114_no_gnuplot/samples/PapeBB_examp3_0.5Hz_canonic/"; //".";

    // get and sort the files (they have the frequency in the name)
    getFile_list(pape_out_dir, files, filepattern);
    
    if (files.size()==0) {
      std::ostringstream es("");
      es << "No files with pattern '" << filepattern << "' were found in directory '" 
         << pape_out_dir << "'. Please make sure filenames with that pattern exist in the directory provided." << endl;
      throw std::invalid_argument(es.str());
    }

    files.sort();
    Nfreq = files.size();

    s      = files.back();
    pos1   = s.rfind("_");
    a      = s.substr(pos1+1, s.length()-pos1-1);
    fmax   = strtod(a.c_str(), NULL);
    f_step = fmax/Nfreq;
    //cout << "pos1= " << pos1 << " len s= " << s.length() << "  a = " << a << endl;
    //cout << "fmax= " << fmax  << "  N = " << N << endl;

    // the frequency vector
    for (int i=0; i<Nfreq; i++) {
      fv.push_back((i)*f_step); // should it be (i+1) here?
    }
    
    return fv;
}
  



// comparison, freq in filename.
bool compare_freq (string first, string second)
{
  string a, b;
  string patt1 ("_");
  string patt2 ("_nm.bin");
  size_t pos1, pos2;
  double x, y;
    
  pos2 = first.find(patt2);
  pos1 = first.rfind(patt1, pos2-1) + patt1.length();
  
  a = first.substr(pos1, pos2-pos1);
  b = second.substr(pos1, pos2-pos1);
  
  x = strtod(a.c_str(), NULL);
  y = strtod(b.c_str(), NULL);

  if (x<y) {return true; }
  else     {return false;}
}


double half_hann(int begin,int end,int i) {
  double answer;
  if(i<begin) answer=1.0;
  else if((begin<=i) && (i<=end)){
      answer=0.5*(cos(Pi*((i-begin))/((int)(end-begin)))+1.0);
  }
  else answer=0.0;
  return answer;
}

complex<double> pulse_spec_fit(double scale, double x) { 
  double fnorm = ((sqrt(4.0+0.25)-0.5)*0.5);
  //double fnorm = 1.0;
  complex<double> answer, I;
  I = complex<double> (0.0, 1.0);
  // fit to effective source spectrum for the propane cannon
  // scale down the frequency to get infrasound

  x=scale*fnorm*x;

  // fit to fourier transform of signal at 10 m
  answer=0.2421405*x*exp((-x-x*x)/2)*exp(I*(-0.5*Pi));
  /* calibrate (times 10 divided by 2 for pressure doubling) and time shift by 0.01 sec */
  /* multiply by 4 pi to get a delta func source strength */
  answer=4.0*Pi*5.0*answer*exp(I*2.0*Pi*x*(0.01));

  return answer;
}


// the source spectrum is loaded/computed in this function
// this version uses FFTN from #define FFTN; 
// the newer version further below has NFFT as argument
int get_source_spectrum( \
								int n_freqs, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile) 
{
  int i;
  double fmx, scale = 1.0;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  fftw_plan p;

  fmx = ((double)FFTN)*f_step; // max frequency
  
  if (src_flg==0) { // use spectrum of a delta function => to get impulse response;
      for(i=0; i<n_freqs; i++) {
          dft_vec[i] = 1.0 + I*0.0; // for all positive freqencies
      }
  }
  else if (src_flg==1) { //use the built-in Roger pulse synthesised from a given spectrum          
      // pulse center frequency f_center should be set <= f_max/5
      
      // f_center is initialized with negative value if the option --f_center was not used
      // if that's the case redefine f_center to be =f_max/5;
      if (f_center < 0) {
          f_center = f_vec[n_freqs-1]/5.0;
      }
      
      // get the scale for the built-in pulse
      if ((f_center > 0) & (f_center <= f_vec[n_freqs-1]/5.0+1.0E-16)) {	
          scale = 1/f_center;
      }
      else if (f_center==0) {
          std::ostringstream es;
          es << endl << "Found an invalid frequency f_center = " << f_center << " Hz." << endl
             << "It appears that the files with the pattern *papeTL* are not of assumed format. They should be 3-column text files similar to 'tloss_1d.pe' output from pape.";
          throw invalid_argument(es.str());
      }
      else {
          std::ostringstream es;
          es << endl << "f_center = " << f_center << " Hz is too large." << endl
             << "For the built-in pulse f_center should be set smaller than f_max/5 = " 
             << f_vec[n_freqs-1]/5.0;
          throw invalid_argument(es.str());
      }

      f=fopen("source_spectrum.dat","w"); // save the source spectrum
      for(i=0; i<n_freqs; i++) {
          //dft_vec[i]=f_step*pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale);
          dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale); // source spectrum
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The source spectrum is saved in file 'source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;
  }
  else if (src_flg==2) { //use custom source spectrum from a file
      cout << "Loading source spectrum from file " << srcfile << endl;
      double *freqv;
      freqv = new double [n_freqs];
      load_source_spectrum(srcfile, freqv, dft_vec, n_freqs);
      
      if (fabs(freqv[0]-f_vec[0])>1.0e-10) { // compare frequencies
                std::ostringstream es;
          es << "The frequencies from the source spectrum file " << srcfile
             << " do not appear to match the ones from the single-frequency pape output files in the provided directory." << endl;
          throw invalid_argument(es.str());
      }
      
      
      
      delete [] freqv;
  }

  if ((src_flg<3)) { // common block only if the source spectrum (at positive frequencies) is available
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; the left zero pad above was unnecessary
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec has src spectrum in it
      for( ; i<FFTN; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of FFTN points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( FFTN, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
  }
  else if (src_flg==3) { // use custom source pulse (time domain) from a file
      double dt;
      vector<double> t, tdp;
      load_source_pulse_td(srcfile, t, tdp);
      dt = t[1]-t[0];
      for (i=0; i<(int)t.size(); i++) {
          pulse_vec[i] = tdp[i]; // automatic conversion to complex<double>
          //cout << i << "  " << t[i] << " " << tdp[i] << endl;
      }
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with FFTN points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( FFTN, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt;
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
  }
  else {
      cerr << "Unknown option " << src_flg << endl;
      exit(1);
  }
  
  // save initial pulse: pulse_vec
  if (1) {
      f = fopen("source_waveform.dat","w");
      if (src_flg==1) {
        for(i=0;i<(FFTN*f_step*5*scale);i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 because we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, 2*real(pulse_vec[i])); 
        }
      fclose(f);
      }
      else { // if not the builtin pulse
        double fct = 2;
        if (src_flg==3) fct = 1;
        
        for(i=0;i<FFTN/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 if we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, fct*real(pulse_vec[i]));
        }    
      }
      cout << "Initial source waveform saved in file 'source_waveform.dat' " << endl
           << " with format: | Time (s) | Amplitude |" << endl;
           //<< " with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
  }

  //// save arg_vec to file
  //printf("in pulse_prop: saving arg_vec to file arg_vec.dat; n_freqs=%d\n", n_freqs);
  //fmx = ((double)FFTN)*f_step;
  //f   = fopen("arg_vec.dat","w");
  //for(i=0;i<FFTN; i++) {
  //    fprintf(f,"%e   %e   %e\n", 1.0*i/fmx, real(arg_vec[i]), imag(arg_vec[i]));
  //}
  //fclose(f);
  //

  return 0;
}


// DV 20151021: Added some code to save fft and ifftw results to compare with Matlab
// we have established the relationship between FFTW_BACKWARD and Matlab ifft 
//      FFTW_BACKWARD(x,NFFT) = NFFT*ifft_matlab(x, NFFT)
// and  FFTW_FORWARD(x,NFFT)  = fft_matlab(x, NFFT)
// Thus spectrum = FFTW_BACKWARD(pulse,NFFT)*dt
//      pulse    = FFTW_FORWARD(spectrum,NFFT)*df
// where df = f_max/f_step and dt = 1/(NFFT*df)
//
// Added NFFT as argument; normalizes the builtin pulse
// the source spectrum is loaded/computed in this function
int get_source_spectrum( \
								int n_freqs, int NFFT, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile) 
{
  int i;
  double dt, fmx, scale;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  fftw_plan p;

  fmx = ((double)NFFT)*f_step; // max frequency; the resulting time interval is dt = 1/fmx
  dt  = 1.0/fmx;
  
  if (src_flg==0) { // use spectrum of a delta function => to get the impulse response;
      for(i=0; i<n_freqs; i++) {
          dft_vec[i] = 1.0 + I*0.0; // 'ones' for all positive frequencies
          dft_vec[i] = dft_vec[i]*exp(I*2.0*Pi*f_vec[i]*(double)n_freqs/2.0); // add delay - DV 20150720
      }
  }
  else if (src_flg==1) { //use the built-in Roger pulse synthesised from a given spectrum          
          
      // make sure the dispersion file contains what we need to generate a good pulse
      // pulse center frequency f_center should be set <= f_max/5
      
      // f_center is initialized with negative value if the option --f_center was not used
      // if that's the case redefine f_center to be = f_max/5;
      if (f_center < 0) {
          f_center = f_vec[n_freqs-1]/5.0;
      }
      
      // get the scale for the built-in pulse
      if ((f_center > 0) & (f_center <= f_vec[n_freqs-1]/5.0+1.0E-16)) {	
          scale = 1/f_center;
      }
      else {
          std::ostringstream es;
          es << endl << "f_center = " << f_center << " Hz is too large." << endl
             << "For the built-in pulse f_center should be set smaller than f_max/5 = " 
             << f_vec[n_freqs-1]/5.0;
          throw invalid_argument(es.str());
      }
      
      f=fopen("builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          //dft_vec[i]=f_step*pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale);
          //dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale); // source spectrum
          dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale)*exp(I*2.0*Pi*f_vec[i]*5.0*scale/4.0); // DV 20150720 - time delay of '5.0*scale/4.0' added; source spectrum
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The built-in source spectrum is saved in file 'builtin_source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;


      // The built-in source spectrum 'dft_vec' gives a time domain source pulse that 
      // is not normalized. Next we will normalize it to have a max positive
      // amplitude of one. So 1. we go into time domain; 2. normalize the pulse
      // then 3. come back to frequency domain with a 'normalized' dft_vec[].

      // Step 1 (of 3): from dft_vec to time domain pulse_vec
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; no need for the left zero pad  if fmin=f_step
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec has src spectrum in it
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // temporarily save the non-normalized pulse
      if (1) {
        f = fopen("unnormalized_pulse_input.dat","w");
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "Unnormalized Initial source waveform saved in file " 
             << "'unnormalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }

      // Step 2 (of 3) - normalize the builtin pulse to an amplitude of 1
      double pmax = 0.0; // will store the max(abs) value
      int    imax = 0;   // index of max
      for (i=0; i<NFFT; i++) {
        if ( pmax<abs(real(pulse_vec[i])) ) {
          pmax = abs(real(pulse_vec[i]));
          imax = i;
        };
      }
      
      // the pulse maximum should not be zero at this point; if it is, something went wrong
      if (pmax==0.0) {
        std::ostringstream es;
        es << endl << "The built-in pulse appears to be all zeros." << endl
           << "Check the parameters defining the built-in source spectrum "
           << " such as f_max, f_center." << endl;
        throw std::invalid_argument(es.str());
      }
      
      // normalize the time domain pulse; the positive max amplitude will be =1.
      // first redefine pmax to have the sign of real(pulse_vec[i])
      pmax = real(pulse_vec[imax]);
      for (i=0; i<NFFT; i++) {
        pulse_vec[i] = real(pulse_vec[i])/pmax; 
      }

      // temporarily save the normalized pulse
      if (0) {
        f = fopen("normalized_pulse_input.dat","w");
          //for(i=0;i<(NFFT*f_step*5*scale);i++) {
          for(i=0;i<NFFT;i++) {
            	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "Normalized Initial source waveform saved in file "
             << "'normalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }
          
      // Step 3 (of 3) - go from the normalized pulse back to the frequency domain
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      //
      // we have established the relationship between FFTW_BACKWARD and Matlab ifft 
      // FFTW_BACKWARD(x,NFFT) = NFFT*ifft_matlab(x, NFFT)
      //

      // temporary save of pure FFTW_BACKWARD(normalized pulse)
      if (0) {
              f = fopen("FFTW_BACKWARD_of_norm_pulse.dat","w");
          //for(i=0;i<(NFFT*f_step*5*scale);i++) {
          for(i=0;i<NFFT;i++) {
            	fprintf(f,"%14.9f %15.6e %15.6e\n", i*f_step, real(arg_vec[i]), imag(arg_vec[i]));
        }      
        fclose(f);
        cout << "FFTW_BACKWARD (normalized pulse) time domain waveform saved in file "
             << "'FFTW_BACKWARD_of_norm_pulse.dat' " << endl
             << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;  
             
      }
      
      
      // temporary test showing FFTW_FORWARD( FFTW_BACKWARD(normalized_pulse_vec)).
      // Our choice here:
      // FFTW_BACKWARD(pulse) takes us to frequency domain
      // FFTW_FORWARD() take us to time domain
      if (0) {
      // here arg_vec = FFTW_BACKWARD(normalized_pulse_vec)
      //
      // perform FFTW_FORWARD to obtain the pulse at the source (time domain) 
      // of NFFT points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p); 
      
      // now pulse_vec stores the result of 
      // FFTW_FORWARD( FFTW_BACKWARD(normalized_pulse_vec)). Call this P2.
      // P2 is in fact equal to P2 = NFFT*normalized_pulse_vec i.e.
      // FFTW_FORWARD( FFTW_BACKWARD(normalized_pulse_vec))= NFFT*normalized_pulse_vec 

        f = fopen("FFTW_FORWARD_arg_vec.dat","w");
          //for(i=0;i<(NFFT*f_step*5*scale);i++) {
          for(i=0;i<NFFT;i++) {
            	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "FFTW_FORWARD_arg_vec.dat time domain waveform saved in file "
             << "'FFTW_FORWARD_arg_vec.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;      
      
      }    

      
      // We just computed arg_vec = FFTW_BACKWARD(normalized pulse_vec))
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt; // multiply by the dt factor
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
                 
      //f=fopen("norm_builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      f=fopen("source_spectrum_input.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The built-in source spectrum of the built-in normalized source pulse"
           << "is saved in file 'source_spectrum_input.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;           
       
  }
  else if (src_flg==2) { //use custom source spectrum from a file
      cout << "Loading source spectrum from file " << srcfile << endl;
      double *freqv;
      freqv = new double [n_freqs];
      load_source_spectrum(srcfile, freqv, dft_vec, n_freqs);
      if (fabs(freqv[0]-f_vec[0])>1.0e-10) { // compare frequencies
          cerr << "The frequencies from the source spectrum file " << srcfile
               << " do not appear to match the ones from the provided dispersion file" 
               << endl << " ... aborting." << endl;
               exit(1);
      }
      delete [] freqv;
  }

  // obtain the pulse at the source: either from builtin or user provided
  if ((src_flg<3)) { // common block only if the source spectrum (at positive frequencies) 
                     // is available: either the builtin spectrum or loaded from a file
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; the left zero pad
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec contains src spectrum times df
    
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      // The multiplication my df to complete the Fouruier integral should have
      // been done above.
      // Note though that since arg_vec contains only the positive freq spectrum
      // we'll have to double the result when we convert to the time domain
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      
  }
  else if (src_flg==3) { // use custom source pulse (time domain) from a file
      double dt;
      vector<double> t, tdp;
      load_source_pulse_td(srcfile, t, tdp);
      dt = t[1]-t[0];
      for (i=0; i<(int)t.size(); i++) {
          pulse_vec[i] = tdp[i]; // automatic conversion to complex<double>
          //cout << i << "  " << t[i] << " " << tdp[i] << endl;
      }
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt;
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
  }
  else { // the code should never get here
      cerr << "Unknown option " << src_flg << endl;
      exit(1);
  }
  
  // save initial pulse: pulse_vec
  if (1) {
      f = fopen("source_waveform_input.dat","w");
      if (src_flg==1) { // the builtin pulse
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 because we only had the spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, 2*real(pulse_vec[i])); 
        }
      fclose(f);
      }
      else { // if not the builtin pulse
        double fct = 2;
        if (src_flg==3) fct = 1;  // use custom source pulse (time domain) from a file
        
        for(i=0;i<NFFT/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 if we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, fct*real(pulse_vec[i]));
        }    
      }
      cout << "Initial source waveform saved in file 'source_waveform_input.dat' " << endl
           << " with format: | Time (s) | Amplitude |" << endl;
           //<< " with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
  }

  //// save arg_vec to file
  //printf("in pulse_prop: saving arg_vec to file arg_vec.dat; n_freqs=%d\n", n_freqs);
  //fmx = ((double)FFTN)*f_step;
  //f   = fopen("arg_vec.dat","w");
  //for(i=0;i<FFTN; i++) {
  //    fprintf(f,"%e   %e   %e\n", 1.0*i/fmx, real(arg_vec[i]), imag(arg_vec[i]));
  //}
  //fclose(f);
  //

  return 0;
}



/*
// DV 20151020: Added NFFT as argument; normalizes the builtin pulse
// the source spectrum is loaded/computed in this function
int get_source_spectrum( \
								int n_freqs, int NFFT, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile) 
{
  int i;
  double fmx, scale;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  fftw_plan p;

  fmx = ((double)NFFT)*f_step; // max frequency
  
  if (src_flg==0) { // use spectrum of a delta function => to get the impulse response;
      for(i=0; i<n_freqs; i++) {
          dft_vec[i] = 1.0 + I*0.0; // 'ones' for all positive frequencies
          dft_vec[i] = dft_vec[i]*exp(I*2.0*Pi*f_vec[i]*(double)n_freqs/2.0); // add delay - DV 20150720
      }
  }
  else if (src_flg==1) { //use the built-in Roger pulse synthesised from a given spectrum          
          
      // make sure the dispersion file contains what we need to generate a good pulse
      // pulse center frequency f_center should be set <= f_max/5
      
      // f_center is initialized with negative value if the option --f_center was not used
      // if that's the case redefine f_center to be = f_max/5;
      if (f_center < 0) {
          f_center = f_vec[n_freqs-1]/5.0;
      }
      
      // get the scale for the built-in pulse
      if ((f_center > 0) & (f_center <= f_vec[n_freqs-1]/5.0+1.0E-16)) {	
          scale = 1/f_center;
      }
      else {
          std::ostringstream es;
          es << endl << "f_center = " << f_center << " Hz is too large." << endl
             << "For the built-in pulse f_center should be set smaller than f_max/5 = " 
             << f_vec[n_freqs-1]/5.0;
          throw invalid_argument(es.str());
      }
      
      f=fopen("builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          //dft_vec[i]=f_step*pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale);
          //dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale); // source spectrum
          dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale)*exp(I*2.0*Pi*f_vec[i]*5.0*scale/4.0); // DV 20150720 - time delay of '5.0*scale/4.0' added; source spectrum
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The built-in source spectrum is saved in file 'builtin_source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;


      // The built-in source spectrum 'dft_vec' gives a time domain source pulse that 
      // is not normalized. Next we will normalize it to have a max positive
      // amplitude of one. So 1. we go into time domain; 2. normalize the pulse
      // then 3. come back to frequency domain with a 'normalized' dft_vec[].

      // Step 1 (of 3): from dft_vec to time domain pulse_vec
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; no need for the left zero pad  if fmin=f_step
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec has src spectrum in it
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // temporarily save the non-normalized pulse
      if (1) {
        f = fopen("unnormalized_pulse_input.dat","w");
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "Unnormalized Initial source waveform saved in file " 
             << "'unnormalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }

      // Step 2 (of 3) - normalize the builtin pulse to an amplitude of 1
      double pmax = 0.0; // will store the max(abs) value
      int    imax = 0;   // index of max
      for (i=0; i<NFFT; i++) {
        if ( pmax<abs(real(pulse_vec[i])) ) {
          pmax = abs(real(pulse_vec[i]));
          imax = i;
        };
      }
      
      // the pulse maximum should not be zero at this point; if it is, something went wrong
      if (pmax==0.0) {
        std::ostringstream es;
        es << endl << "The built-in pulse appears to be all zeros." << endl
           << "Check the parameters defining the built-in source spectrum "
           << " such as f_max, f_center." << endl;
        throw std::invalid_argument(es.str());
      }
      
      // normalize the time domain pulse; the positive max amplitude will be =1.
      // first redefine pmax to have the sign of real(pulse_vec[i])
      pmax = real(pulse_vec[imax]);
      for (i=0; i<NFFT; i++) {
        pulse_vec[i] = real(pulse_vec[i])/pmax; 
      }

      // temporarily save the normalized pulse
      if (1) {
        f = fopen("normalized_pulse_input.dat","w");
          for(i=0;i<(NFFT*f_step*5*scale);i++) {
            	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "Normalized Initial source waveform saved in file "
             << "'normalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }
          
      // Step 3 (of 3) - go from the normalized pulse back to the frequency domain
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      //double dt = f_vec[n_freqs-1]/n_freqs;
      //double dt = double (n_freqs)/((double)NFFT);
      
      // We just computed arg_vec = FFTW_BACKWARD(normalized pulse_vec))
      // Since FFTW has the property: FFTW_FORWARD(FFTW_BACKWARD(x)) = NFFT*x then
      // FFTW_FORWARD(arg_vec) = NFFT*normalized_pulse. 
      // Thus if we divide by arg_vec by NFFT below we'll reobtain the
      // normalized_pulse when we compute FFTW_FORWARD(arg_vec)
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]/((double)NFFT); // divide by NFFT in FFTW
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
                 
      f=fopen("norm_builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The built-in source spectrum is saved in file 'norm_builtin_source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;           
               
  }
  else if (src_flg==2) { //use custom source spectrum from a file
      cout << "Loading source spectrum from file " << srcfile << endl;
      double *freqv;
      freqv = new double [n_freqs];
      load_source_spectrum(srcfile, freqv, dft_vec, n_freqs);
      if (fabs(freqv[0]-f_vec[0])>1.0e-10) { // compare frequencies
          cerr << "The frequencies from the source spectrum file " << srcfile
               << " do not appear to match the ones from the provided dispersion file" 
               << endl << " ... aborting." << endl;
               exit(1);
      }
      delete [] freqv;
  }

  if ((src_flg<3)) { // common block only if the source spectrum (at positive frequencies) 
                     //is available: either the builtin spectrum or loaded from a file
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; the left zero pad
      int i0 = i;
      //for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec has src spectrum in it
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]; // arg_vec has src spectrum in it      
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      // Note though that since arg_vec contains only the positive freq spectrum
      // we'll have to double the result when we convert to the time domain
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      
  }
  else if (src_flg==3) { // use custom source pulse (time domain) from a file
      double dt;
      vector<double> t, tdp;
      load_source_pulse_td(srcfile, t, tdp);
      dt = t[1]-t[0];
      for (i=0; i<(int)t.size(); i++) {
          pulse_vec[i] = tdp[i]; // automatic conversion to complex<double>
          //cout << i << "  " << t[i] << " " << tdp[i] << endl;
      }
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt;
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
  }
  else { // the code should never get here
      cerr << "Unknown option " << src_flg << endl;
      exit(1);
  }
  
  // save initial pulse: pulse_vec
  if (1) {
      f = fopen("source_waveform_input.dat","w");
      if (src_flg==1) {
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 because we only had the spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, 2*real(pulse_vec[i])); 
        }
      fclose(f);
      }
      else { // if not the builtin pulse
        double fct = 2;
        if (src_flg==3) fct = 1;  // use custom source pulse (time domain) from a file
        
        for(i=0;i<NFFT/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 if we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, fct*real(pulse_vec[i]));
        }    
      }
      cout << "Initial source waveform saved in file 'source_waveform_input.dat' " << endl
           << " with format: | Time (s) | Amplitude |" << endl;
           //<< " with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
  }

  //// save arg_vec to file
  //printf("in pulse_prop: saving arg_vec to file arg_vec.dat; n_freqs=%d\n", n_freqs);
  //fmx = ((double)FFTN)*f_step;
  //f   = fopen("arg_vec.dat","w");
  //for(i=0;i<FFTN; i++) {
  //    fprintf(f,"%e   %e   %e\n", 1.0*i/fmx, real(arg_vec[i]), imag(arg_vec[i]));
  //}
  //fclose(f);
  //

  return 0;
}
*/


// 'fft_pulse_prop'
void fft_pulse_prop(double t0, int n_freqs, double df, double *f_vec, \
                    vector< complex<double> > PP, complex<double> *dft_vec, \
                    complex<double> *pulse_vec)
{
  int i,i0,smooth_space;    // CHH 191029: j unused
  complex<double> cup,t_phase,*arg_vec;
  complex<double> I (0.0, 1.0);  
  fftw_plan p;

  if(FFTN < n_freqs){
      throw invalid_argument("fft too short (i.e. FFTN < n_freqs), exiting.");
  }
  
  arg_vec = new complex<double> [FFTN];
  
  
  // old left zero-pad; sometimes with a C compiler it gives i=2 even if fvec[0]=df
  // so I rewrote it below; DV
  //for(i=0;i*df<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum;
  
  // new left zero pad up to f_min present in the spectrum;
  for(i=0;(fabs(i*df-f_vec[0])>1.0e-10) & (i*df<f_vec[0]); i++) arg_vec[i]=0.0;
  //printf("in fft_pulse_prop:f_vec[0] = %f; df = %f; i=%d\n", f_vec[0], df, i);

  i0 = i;
  for(i=0;i<n_freqs;i++){
      t_phase=exp(-I*2.0*Pi*f_vec[i]*t0); // corresponding to reduced time t0

      cup=PP[i]*t_phase;
      
      // up to a factor ((delayed Fourier pressure component)*source_spectrum*df) 
      // see eq. 8.9 pg 480 in Comp. Oc. Acoust. 1994 ed.)
      cup=cup*dft_vec[i]*df;
      arg_vec[i0+i] = cup; // note sqrt_rho_ratio
  }

  smooth_space=(int)floor(0.1*n_freqs); // smoothly zero out on right; as in RW (July 2012)

  for(i=n_freqs-smooth_space;i<n_freqs;i++){
      // arg_vec[i]=arg_vec[i]*half_hann(n_freqs-smooth_space,n_freqs-df,i);
      arg_vec[i]=arg_vec[i]*half_hann(n_freqs-smooth_space,n_freqs-1,i); // changed df to 1 as df doesn't make sense here
  }
  for(;i<FFTN;i++) arg_vec[i]=0.0; // right zero pad
  
  // save arg_vec for debugging purposes
  if (1) {
      FILE *fp = fopen("arg_vec.dat", "w");
      for (i=0; i<FFTN; i++) {
          fprintf(fp, "%12.6f %15.6e %15.6e\n", i*df, real(arg_vec[i]), imag(arg_vec[i]));
      }
      fclose(fp);
      printf("arg_vec saved in file 'arg_vec.dat' with format: i*df | Re(arg_vec) | Im(arg_vec)\n");
  }
   
  // 
  //perform fft of arg_vec to obtain the propagated time domain pulse (synthesis step)
  //
  p=fftw_plan_dft_1d(FFTN, reinterpret_cast<fftw_complex*>(arg_vec), \
												   reinterpret_cast<fftw_complex*>(pulse_vec), \
												   FFTW_FORWARD,FFTW_ESTIMATE); 												 
  fftw_execute(p); 
  fftw_destroy_plan(p);
  delete [] arg_vec;
}  // end of debug version of 'fft_pulse_prop'




int pulse_prop_src2rcv_grid2(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, double f_step, double *f_vec, \
					double f_center, vector< complex<double> > PP,  \
					int src_flg, string srcfile, int pprop_src2rcv_flg) 
{
  int i,n;
  double rr, tskip, fmx, t0;	
  complex<double> cup,*dft_vec,*pulse_vec,*arg_vec;
  //complex<double> I = complex<double> (0.0, 1.0);  // CHH 191029: Unused
  FILE *f;

  dft_vec   = new complex<double> [n_freqs];
  pulse_vec = new complex<double> [FFTN];
  arg_vec   = new complex<double> [FFTN];

  fmx = ((double)FFTN)*f_step; // max frequency
  
  
  get_source_spectrum(n_freqs,  f_step, f_vec,  f_center, \
								      dft_vec, pulse_vec, arg_vec, src_flg,  srcfile);
								      
	
	if (pprop_src2rcv_flg) { // propagation to one receiver at distance RR from source
	    cout << "Doing pulse propagation source-to-receiver at one range: " 
	         << R_start/1000.0 << " km" << endl;
	         
	    t0 = R_start/max_cel;
	    //
	    // fft propagation: note that here 'pulse_vec' is rewritten with the propagated pulse
	    //							

	    fft_pulse_prop(t0, n_freqs, f_step, f_vec, \
                    PP, dft_vec, pulse_vec);
                    
                    
                    
	    // save propagated pulse to file						   															
      f = fopen(filename,"w");
      for(i=0; i<FFTN; i++) {
          //fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, real(pulse_vec[i]), imag(pulse_vec[i]));
          fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx+t0, real(pulse_vec[i]));
      }
      fclose(f);
      //cout << "Propagated pulse saved in file: '" << filename << "'" << endl
      //     << "' with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;      
      cout << "Propagated pulse saved in file: '" << filename << "'" << endl
           << "with format: | Time [s] | Waveform |" << endl;      


	}
	else { // propagation to several receivers
	    cout << "Doing pulse propag source-to-receiver on grid ..." << endl;					      
      printf("----------------------------------------------\n");
      printf("max_celerity     t0             R\n");
      printf("    m/s          sec            km\n");
      printf("----------------------------------------------\n");

      f=fopen(filename,"w");
      tskip = 0.0;
      for(n=0; n<=(int)(floor((R_end-R_start)/DR)); n++) {
          rr= R_start + DR*n;
          t0=tskip+rr/max_cel;
          printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

	        fft_pulse_prop(t0, n_freqs, f_step, f_vec, \
                        PP, dft_vec, pulse_vec);
							              
          for(i=0;i<FFTN;i++){
              fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]));
          }
          fprintf(f,"\n");
      }
      fclose(f);
      printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
      printf("Time array length = %d; delta_t = %g s\n", FFTN, 1.0/fmx);
      printf("Propagation results saved in file: %s\n", filename);
      printf("with columns: | R [km] | Time [s] | Waveform(R,t) |\n");
  }

  delete [] dft_vec;
  delete [] arg_vec;
  delete [] pulse_vec;

  return 0;
}


int pulse_prop_src2rcv_grid3(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, double f_step, double *f_vec, \
					double f_center, \
					string pape_out_dir, string filepattern, \
					int src_flg, string srcfile, int pprop_src2rcv_flg) 
{
  int i,n;
  double rr, tskip, fmx, t0;	
  complex<double> cup,*dft_vec,*pulse_vec,*arg_vec;
  //complex<double> I = complex<double> (0.0, 1.0);  // CHH 191029: Unused
  vector< complex<double> > PP;
  FILE *f;

  dft_vec   = new complex<double> [n_freqs];
  pulse_vec = new complex<double> [FFTN];
  arg_vec   = new complex<double> [FFTN];

  fmx = ((double)FFTN)*f_step; // max frequency
  
  
  get_source_spectrum(n_freqs,  f_step, f_vec,  f_center, \
								      dft_vec, pulse_vec, arg_vec, src_flg,  srcfile);
								      
	
	if (pprop_src2rcv_flg) { // propagation to one receiver at distance RR from source
	    cout << "Doing pulse propagation source-to-receiver at one range: " 
	         << R_start/1000.0 << " km" << endl;
	         
	    t0 = R_start/max_cel;
	    
	    
	    // interpolate pape output
      PP = interpDir(pape_out_dir, filepattern, R_start/1000.0);
      
      // show PP values
      //for (int j=0; j<PP.size(); j++) {
      //  printf("%3d  PP = %g  %g \n", j, real(PP[j]), imag(PP[j]));
      //}

	    
	    //
	    // fft propagation: note that here 'pulse_vec' is rewritten with the propagated pulse
	    //							
	    fft_pulse_prop(t0, n_freqs, f_step, f_vec, \
                    PP, dft_vec, pulse_vec);         
             
	    // save propagated pulse to file
	    // DV 20150609 - factor of two necessary because we only used 
	    // the positive freq. spectrum in the fft							   															
      f = fopen(filename,"w");
      for(i=0; i<FFTN; i++) {
          //fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, real(pulse_vec[i]), imag(pulse_vec[i]));
          fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, 2.0*real(pulse_vec[i]), 2.0*imag(pulse_vec[i]));
      }
      fclose(f);
      cout << "Propagated pulse saved in file: '" << filename << "'" << endl
           << "' with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
      //cout << "Propagated pulse saved in file: '" << filename << "'" << endl
      //     << "with format: | Time (s) | Waveform |" << endl << endl;            
	}
	else { // propagation to several receivers
	    cout << "Doing pulse propag source-to-receiver on grid ..." << endl;					      
      printf("----------------------------------------------\n");
      printf("max_celerity     t0             R\n");
      printf("    m/s          sec            km\n");
      printf("----------------------------------------------\n");

      f=fopen(filename,"w");
      tskip = 0.0;
      for(n=0; n<=(int)(floor((R_end-R_start)/DR)); n++) {
          rr= R_start + DR*n;
          t0=tskip+rr/max_cel;
          printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

          // interpolate pape output
          PP = interpDir(pape_out_dir, filepattern, rr/1000.0);
          
          // fft
	        fft_pulse_prop(t0, n_freqs, f_step, f_vec, \
                        PP, dft_vec, pulse_vec);
							              
          for(i=0;i<FFTN;i++){
              //fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]));
              fprintf(f,"%10.3f %12.6f %15.6e %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i])); // output the imaginary part which is in fact what modess calls the real part !
          }
          fprintf(f,"\n");
      }
      fclose(f);
      printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
      printf("Time array length = %d; delta_t = %g s\n", FFTN, 1.0/fmx);
      printf("Propagation results saved in file: %s\n", filename);
      printf("with columns: R (km) | time (s) | pulse(R,t) |\n");
  }

  delete [] dft_vec;
  delete [] arg_vec;
  delete [] pulse_vec;

  return 0;
}


//20151020 DV: added NFFT as argument
int pulse_prop_src2rcv_grid4(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int NFFT, int n_freqs, double f_step, double *f_vec, \
					double f_center, \
					string pape_out_dir, string filepattern, \
					int src_flg, string srcfile, int pprop_src2rcv_flg) 
{
  int i,n;
  double rr, tskip, fmx, t0;	
  complex<double> cup,*dft_vec,*pulse_vec,*arg_vec;
  //complex<double> I = complex<double> (0.0, 1.0);  // CHH 191029: Unused
  vector< complex<double> > PP;
  FILE *f;
  
  // if the option --nfft was not passed to the main program then the 
  // NFFT default was -1. Otherwise it has the requested positive value. 
  // Here we make sure that whatever the current value of NFFT is 
  // it's not less than 4*n_freqs
  if (NFFT<n_freqs) {
    NFFT = 4*n_freqs;
  }

  dft_vec   = new complex<double> [n_freqs];
  pulse_vec = new complex<double> [FFTN];
  arg_vec   = new complex<double> [FFTN];

  fmx = ((double)FFTN)*f_step; // max frequency
  
  

								      
  get_source_spectrum(n_freqs,  NFFT, f_step, f_vec,  f_center,
		      dft_vec, pulse_vec, arg_vec, src_flg,  srcfile);								      
								      
	
	if (pprop_src2rcv_flg) { // propagation to one receiver at distance RR from source
	    cout << "Doing pulse propagation source-to-receiver at one range: " 
	         << R_start/1000.0 << " km" << endl;
	         
	    t0 = R_start/max_cel;
	    
	    
	    // interpolate pape output
      PP = interpDir(pape_out_dir, filepattern, R_start/1000.0);
      
      // show PP values
      //for (int j=0; j<PP.size(); j++) {
      //  printf("%3d  PP = %g  %g \n", j, real(PP[j]), imag(PP[j]));
      //}

	    
	    //
	    // fft propagation: note that here 'pulse_vec' is rewritten with the propagated pulse
	    //							
	    fft_pulse_prop(t0, n_freqs, f_step, f_vec, \
                    PP, dft_vec, pulse_vec);         
             
	    // save propagated pulse to file
	    // DV 20150609 - factor of two necessary because we only used 
	    // the positive freq. spectrum in the fft							   															
      f = fopen(filename,"w");
      for(i=0; i<FFTN; i++) {
          //fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, real(pulse_vec[i]), imag(pulse_vec[i]));
          fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, 2.0*real(pulse_vec[i]), 2.0*imag(pulse_vec[i]));
      }
      fclose(f);
      cout << "Propagated pulse saved in file: '" << filename << "'" << endl
           << "' with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
      //cout << "Propagated pulse saved in file: '" << filename << "'" << endl
      //     << "with format: | Time (s) | Waveform |" << endl << endl;            
	}
	else { // propagation to several receivers
	    cout << "Doing pulse propag source-to-receiver on grid ..." << endl;					      
      printf("----------------------------------------------\n");
      printf("max_celerity     t0             R\n");
      printf("    m/s          sec            km\n");
      printf("----------------------------------------------\n");

      f=fopen(filename,"w");
      tskip = 0.0;
      for(n=0; n<=(int)(floor((R_end-R_start)/DR)); n++) {
          rr= R_start + DR*n;
          t0=tskip+rr/max_cel;
          printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

          // interpolate pape output
          PP = interpDir(pape_out_dir, filepattern, rr/1000.0);
          
          // fft
	        fft_pulse_prop(t0, n_freqs, f_step, f_vec, \
                        PP, dft_vec, pulse_vec);
							              
          for(i=0;i<FFTN;i++){
              //fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]));
              fprintf(f,"%10.3f %12.6f %15.6e %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i])); // output the imaginary part which is in fact what modess calls the real part !
          }
          fprintf(f,"\n");
      }
      fclose(f);
      printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
      printf("Time array length = %d; delta_t = %g s\n", FFTN, 1.0/fmx);
      printf("Propagation results saved in file: %s\n", filename);
      printf("with columns: R (km) | time (s) | pulse(R,t) |\n");
  }

  delete [] dft_vec;
  delete [] arg_vec;
  delete [] pulse_vec;

  return 0;
}




int load_source_pulse_td(string srcpulsetdfn, vector<double> &t, vector<double> &tdp ) {
  // load time domain source pulse: time (s) | amplitude
  cout << "Loading time domain source pulse from file " << srcpulsetdfn << endl;

  double d1, d2; //, d3;
  ifstream indata;

  indata.open(srcpulsetdfn.c_str());	
  if (!indata) {
      cerr << "file " << srcpulsetdfn << " could not be opened.\n";
      exit(1);
  }

  //indata >> d1 >> d2 >> d3; // read first line in file
  indata >> d1 >> d2; // read first line in file
  while (!indata.eof()) {
      t.push_back(d1);
      tdp.push_back(d2);
      //indata >> d1 >> d2 >> d3;
      indata >> d1 >> d2;
  }
  indata.close();

  //// print data
  //for (unsigned i=0; i<t.size(); i++) {
  //    cout << i << "  " << t[i] << endl;
  //}

  return 0;
}


int load_source_spectrum(string srcspfn, double *freqv, complex<double> *dft_vec, int fftn) {
  // load source spectrum
  //cout << "Loading source spectrum from file " << srcspfn << endl;
  
  int i;
  double d1, d2, d3;
  complex<double> I (0.0, 1.0);
  ifstream indata;

  indata.open(srcspfn.c_str());	
  
  if (!indata) {
      std::ostringstream es("");
      es << "File " << srcspfn << " could not be opened." << endl;
      throw std::invalid_argument(es.str());
  } 

  i = 0;
  indata >> d1 >> d2 >> d3; // read first line in file
  while (!indata.eof()) {
      freqv[i] = d1;
      dft_vec[i] = d2 + I*d3;
      indata >> d1 >> d2 >> d3;

      if (i>(fftn-1)) {
        std::ostringstream es("");
        es << "The source spectrum file " << srcspfn << " appears to have more frequencies"
           << " than the number of frequencies derived from single-frequency pape files in"
           << " the directory provided. The two numbers need to match. Found " 
           << fftn << " frequencies in the pre-computed pape output directory." << endl;
        throw std::invalid_argument(es.str());
      }      

      i++;
  }
  indata.close();

  return 0;
}	  
