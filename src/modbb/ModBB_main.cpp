#include <complex>
#include "Atmosphere.h"
#include "anyoption.h"
#include "ModBB_lib.h"
#include "SolveModBB.h"
#include "ProcessOptionsBB.h"


/*
 * This is the Normal Modes Broadband code 
 * based on either the Effective Sound Speed Approximation (Modess) or
 * on the Wide-Angle High-Mach modal code (WMod)
 * @authors Doru Velea; Jelle Assink; Roger Waxler; Claus Hetzer; 
 */

using namespace NCPA;
using namespace std;

#ifndef Pi
#define Pi 3.141592653589793
#endif
#define MAX_MODES 4000


// Function to parse the options from the command line/config file
AnyOption *parseInputOptions( int argc, char **argv );

//
// main
//
int main( int argc, char **argv ) {

  // Read options from the command line as well as an options file
  AnyOption *opt = parseInputOptions( argc, argv ); 

  //
  // Physical values are usually in SI (Systeme International) units 
  // unless mentioned otherwise.
  //

  // some flags defining the kind of output:
  // flags related to writing dispersion file(s)
  bool   w_disp_src2rcv_flg;
  bool   w_disp_flg;
  //bool   usemodess_flg;

  // flags related to pulse propagation
  bool   pprop_grid_flg;
  bool   pprop_s2r_flg;
  bool   pprop_s2r_grid_flg;

  // object to process the options
  ProcessOptionsBB *oBB;
  oBB = new ProcessOptionsBB(opt);

  // get dispersion / propagation flags
  w_disp_src2rcv_flg   = oBB->getW_disp_src2rcv_flg();
  w_disp_flg           = oBB->getW_disp_flg();
  pprop_grid_flg       = oBB->getPprop_grid_flg();
  pprop_s2r_flg        = oBB->getPprop_s2r_flg();
  pprop_s2r_grid_flg   = oBB->getPprop_s2r_grid_flg();

  // set the timer to measure running time
  std::time_t tm1 = std::time(NULL);

  //
  // logic to handle either computing the dispersion and modal values 
  // or to propagate a pulse
  //
  if (w_disp_src2rcv_flg || w_disp_flg) { //compute and write dispersion
      //
      // declare some variables to use in this bloc; 
      //		
      string   atmosfile     ;	// stores the atmospheric profile name 
      string   atmosfileorder;  // column order in atmosfile
      string   wind_units    ;  // default 'mpersec'
      string   out_disp_file ;  // dispersion filename
      string   gnd_imp_model ;  // default is "rigid"; corresponding to admittance = 0;

      SampledProfile  *atm_profile;  // the atmospheric profile object
  
      // create the atmopsheric profile object and set the azimuth
      atmosfileorder = oBB->getAtmosfileorder();
      //atm_profile = new SampledProfile( atmosfile, atmosfileorder.c_str(), oBB->getSkiplines() );
      atm_profile = new SampledProfile( oBB->getAtmosfile(), \
                        oBB->getAtmosfileorder().c_str(), oBB->getSkiplines() );
      atm_profile->setPropagationAzimuth(oBB->getAzimuth());
      //cout << "oBB->getAzimuth() = " << oBB->getAzimuth() << endl;
      // cout << "1  atm_profile->getPropagationAzimuth() = " << atm_profile->getPropagationAzimuth() << endl;
									        
      SolveModBB *a = new SolveModBB(oBB, atm_profile); 									        
				         	 					 
      // solve for modes/ loop over frequencies (broadband) - main action happens here					 
      //a->computeModes(usemodess_flg);
      a->computeModes(oBB->getUsemodess_flg());
      a->printParams();
      
      // save atm. profile if requested
      if (1) {
          saveAtm_profile(atm_profile, wind_units);
      }	

      delete a;
      cout << "Finished computing and saving dispersion data" << endl;	
								
  } // end computing dispersion data
  else if (pprop_s2r_flg || pprop_s2r_grid_flg) { //pulse propagation source-to-receiver
      bool zero_attn_flg;
      int    Nfreq, NFFT, src_flg, *mode_count;
      double f_center, f_step, max_cel, R_start, R_end, DR;
      double rho_zsrc, rho_zrcv;
      //double f_center, scale;   
      double *f_vec;
      double **re_k, **im_k, **mode_S, **mode_R;

      string pprop_s2r_disp_file; 
      string waveform_out_file;
      string src_file;
      
      NFFT                = oBB->getNFFT();
      waveform_out_file   = oBB->getWaveform_out_file();
      pprop_s2r_disp_file = oBB->getPprop_s2r_disp_file();	
      //cout << "pprop_s2r_disp_file = " << pprop_s2r_disp_file << endl;	

      // get number of (positive) frequencies from the dispersion file
      Nfreq = count_rows_arbcol(pprop_s2r_disp_file.c_str());
      Nfreq = Nfreq + 1; // count the DC frequency

      f_vec      = new double [Nfreq];
      mode_count = new int [Nfreq];
      re_k       = dmatrix(Nfreq, MAX_MODES);
      im_k       = dmatrix(Nfreq, MAX_MODES);
      mode_S     = dmatrix(Nfreq, MAX_MODES);
      mode_R     = dmatrix(Nfreq, MAX_MODES);
      	
      // read the dispersion file
      readDispersion_bb_ascii(pprop_s2r_disp_file, f_vec, mode_count, \
                              &rho_zsrc, &rho_zrcv, re_k, im_k, mode_S, mode_R);     
											        															
      f_step   = f_vec[1] - f_vec[0];
      max_cel  = oBB->getMax_celerity();     
  	  R_start  = oBB->getR_start();
      R_end    = oBB->getR_end();
      DR       = oBB->getDR();
      f_center = oBB->getF_center();
      src_flg  = oBB->getSrc_flg();
      src_file = oBB->getSrcfile();
      zero_attn_flg = oBB->getZeroAttn_flg();
            	
      pulse_prop_src2rcv_grid2( waveform_out_file.c_str(), max_cel, \
								R_start, DR, R_end, Nfreq, NFFT, f_step, f_vec, \
								f_center, mode_count, rho_zsrc, rho_zrcv, \
								re_k, im_k, mode_S, mode_R, \
								src_flg, src_file, pprop_s2r_flg, zero_attn_flg);
      
      //// temporary plotInitialPulse
      //if (oBB->getPlot_flg())  {
          //printf("plotting initial_pulse.dat\n");
          //plotInitialPulse(oBB->getSrc_flg(), oBB->getF_center(), f_vec[Nfreq-1]);
      //}
      if (zero_attn_flg) {
        printf("The attenuation was set to sero with the option use_zero_attn flag = %d.\n", zero_attn_flg);
      }				

      // clean up			
      delete[] f_vec;
      delete[] mode_count;
      free_dmatrix(re_k,   Nfreq, MAX_MODES);
      free_dmatrix(im_k,   Nfreq, MAX_MODES);								
      free_dmatrix(mode_S, Nfreq, MAX_MODES);
      free_dmatrix(mode_R, Nfreq, MAX_MODES);
  }
  else if (pprop_grid_flg) {
      // compute the 2D pressure in the frame of width_km by height_km, 
      // starting at R_start_km, for ntsteps every tmstep seconds.
      int ntsteps;
      double R_start_km, width_km, height_km, max_cel, tmstep, f_center;
      string frame_file_stub;
      string pprop_grid_dirname;
      pprop_grid_dirname = oBB->getPprop_grid_dirname();

      cout << "Computing 2D broadband field from directory "
           << pprop_grid_dirname << endl;
		
      R_start_km = oBB->getR_start_km();
      width_km   = oBB->getWidth_km();
      height_km  = oBB->getHeight_km();
      max_cel    = oBB->getMax_celerity();	
      tmstep     = oBB->getTmstep();
      ntsteps    = oBB->getNtsteps();
      f_center   = oBB->getF_center();

      frame_file_stub = oBB->getFrame_file_stub();
			
      process2DPressure(R_start_km, width_km, height_km, \
                        max_cel, tmstep, ntsteps, f_center, \
                        frame_file_stub, pprop_grid_dirname, \
                        oBB->getSrc_flg(), oBB->getSrcfile(), \
                        oBB->getZeroAttn_flg());
  }
  else {
      cout << "Nothing done! Reserved for future option here" << endl;
  }

  delete(oBB); // delete object processing options

  std::time_t tm2 = std::time(NULL);
  cout << "\nRun duration: " << difftime(tm2,tm1) << " seconds." << endl;
  cout << " ... main() broadband version is done." << endl;

  return 0;
} // end of Broadband main();


AnyOption *parseInputOptions( int argc, char **argv ) {

  // parse input options
  AnyOption *opt = new AnyOption();

  opt->addUsage( "----------------------------------------------------------------------------" );
  opt->addUsage( "|                             NCPA Infrasound                              |" );  	
  opt->addUsage( "|                         Normal Modes Broadband                           |" );
  opt->addUsage( "|    Based on either: Effective Sound Speed Approximation - see ModESS     |" );
  opt->addUsage( "|                     Wide_Angle High-Mach code - see WMod                 |" );  
  opt->addUsage( "----------------------------------------------------------------------------" );	
  opt->addUsage( "Usage: " );
  opt->addUsage( "" );
  opt->addUsage( "The options below can be specified in a colon-separated file \"ModBB.options\"");
  opt->addUsage( "or at the command line.  Command-line options override file options." );
  opt->addUsage( "Be sure to precede all options with two minuses (--)." );
  opt->addUsage( " --help -h                Print this message and exit" );
  opt->addUsage( "" );
  opt->addUsage( "One of two algorithms can be used to perform pulse propagation." );
  opt->addUsage( "The first is based on the Effective Sound Speed Approximation (as in ModESS);" );
  opt->addUsage( "the second is based on the the Wide_Angle High-Mach solution ");
  opt->addUsage( "of the wave equation (see  implementation in WMod)." );
  opt->addUsage( "ModESS is faster but it is accurate for (launch) angles less than 30 deg and" );
  opt->addUsage( "low wind speeds. WMod extends the validity to higher angles" ); 
  opt->addUsage( "and high Mach numbers but it runs slower." );
  opt->addUsage( "Options --use_modess and --use_wmod allow the user to choose");
  opt->addUsage( "the desired algorithm when computing the dispersion data (see step 1 below)." );
  opt->addUsage( "" );
  opt->addUsage( "To propagate a pulse, 2 steps must be completed:");
  opt->addUsage( " 1. A dispersion file must be available or computed" );
  opt->addUsage( "     using the option --out_disp_src2rcv_file  ." );  
//  opt->addUsage( "     use either option --out_dispersion_files  or --out_disp_src2rcv_file" );
  opt->addUsage( " 2. Perform pulse propagation for one of several scenarios:");
  opt->addUsage( "    a. source-to-receiver at one range (option --pulse_prop_src2rcv)");	
  opt->addUsage( "    b. source-to-receiver at several equally spaced ranges " );
  opt->addUsage( "       (option --pulse_prop_src2rcv_grid)");
//  opt->addUsage( "    c. computing the whole 2D pressure field " );
//  opt->addUsage( "       (most expensive - option --pulse_prop_grid)" );
  opt->addUsage( "    The source type can be one of the following:" );
  opt->addUsage( "        delta function              -> see option --get_impulse_resp" );  
  opt->addUsage( "        built-in pulse choice 1     -> see option --use_builtin_pulse1" );
  opt->addUsage( "        built-in pulse choice 2     -> see option --use_builtin_pulse2" );
  opt->addUsage( "        user-provided spectrum file -> see option --src_spectrum_file" );
  opt->addUsage( "        user-provided waveform file -> see option --src_waveform_file" );
  opt->addUsage( "" );
//  opt->addUsage( "To compute a dispersion file: one of the following 2 options is REQUIRED:" );
  opt->addUsage( "To compute a dispersion file: the following option is REQUIRED:" );
  opt->addUsage( " --out_disp_src2rcv_file  <dispersion filename>");
  opt->addUsage( "                    Output dispersion curves and modal values for" );
  opt->addUsage( "                    source-to-receiver propagation to the specified file." );	
//  opt->addUsage( " --out_dispersion_files   <dispersion filename stub>");
//  opt->addUsage( "                    Output dispersion curves and modal values on a 2D grid" );
//  opt->addUsage( "                    to binary files at each frequency. The resulting filenames" );
//  opt->addUsage( "                    have the stub and frequency appended: " );
//  opt->addUsage( "                    e.g. <stub><freq>_nm.bin." );
//  opt->addUsage( "                    This option is computationally expensive." );
  opt->addUsage( "" );
//  opt->addUsage( " Examples (run in the 'samples' directory):" );
//  opt->addUsage( "" );
//  opt->addUsage( "   a. Compute dispersion file that will be used to compute the pressure pulse at 1 receiver. Assume that we want to end up with a pulse having a spectrum with a maximum frequency of f_max=0.5 Hz. Also assume that we want the pulse represented on a time record of T=512 seconds. The number of positive frequencies necessary for the calculation is T*f_max = 256 i.e.256 frequencies between 0 and 0.5 Hz. Thus we know f_max=0.5 Hz and f_step=f_max/256=0.001953125 Hz. The corresponding run command is:" );
  opt->addUsage( " Example (run in the 'samples' directory):" );
  opt->addUsage( "" );
  opt->addUsage( " Example to obtain a dispersion file that will be used to compute" );
  opt->addUsage( " the propagated pressure pulse at a distant receiver. " );
  opt->addUsage( " Assume that we want to end up with a pulse having a spectrum" );
  opt->addUsage( " with a maximum frequency of f_max=0.5 Hz. " );
  opt->addUsage( " Also assume that we want the pulse represented on a time record " );
  opt->addUsage( " of T=512 seconds. The number of positive frequencies necessary  " );
  opt->addUsage( " for the calculation is T*f_max = 256 i.e.256 frequencies " );
  opt->addUsage( " between f_min=f_step and 0.5 Hz. The step in frequency is  " );
  opt->addUsage( " f_step=f_max/256=0.001953125 Hz." );
  opt->addUsage( " The corresponding run command is:" );
  opt->addUsage( "" );
  opt->addUsage( "    >> ../bin/ModBB --out_disp_src2rcv_file myDispersionFile.dat " );
  opt->addUsage( "       --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp " );
  opt->addUsage( "       --skiplines 0 --azimuth 90 --f_step 0.001953125 --f_max 0.5 --use_modess" );
  opt->addUsage( "" );
  opt->addUsage( " Each line in this dispersion file has (4 + 4*n_modes) entries" );
  opt->addUsage( " in the following format:" );  
  opt->addUsage( " freq n_modes rho_src rho_rcv Re(k_pert_m) Im(k_pert_m) V_m(z_src) V_m(z_rcv)" );
  opt->addUsage( " where m varies from 1 to n_modes. k_pert_m, and V_m are" );
  opt->addUsage( " the m_th wavenumber and mode amplitude respectively." );  
  opt->addUsage( " z_src and z_rcv stand for source and receiver height respectively." );
  opt->addUsage( "" );
//  opt->addUsage( "   b. Compute dispersion files for propagation to all receivers on a 2D grid: for 256 frequencies from 0 to 0.5 Hz in steps of 0.5/256 Hz:" ); 
//  opt->addUsage( "" ); 
//  opt->addUsage( "    >> ../bin/ModBB --out_dispersion_files disprs --atmosfile NCPA_canonical_profile_zuvwtdp.dat --atmosfileorder zuvwtdp --skiplines 0 --azimuth 90 --f_step 0.001953125 --f_max 0.5 --use_modess" );
  opt->addUsage( "" );  
  opt->addUsage( "In addition the following options are REQUIRED:" );
  opt->addUsage( " --use_modess              Prompts the use of ModESS algorithm." );
  opt->addUsage( " --use_wmod                Prompts the use of WMod algorithm." );
  opt->addUsage( " Note that --use_modess and --use_wmod are mutually exclusive." );
  opt->addUsage( " --atmosfile  <filename>   Uses an ASCII atmosphere file" );
  opt->addUsage( "                           referenced to Mean Sea Level (MSL)." );  
  opt->addUsage( " --atmosfileorder          The order of the (z,t,u,v,w,p,d) fields in " );
  opt->addUsage( "                           the ASCII file (Ex: 'ztuvpd')" );
  opt->addUsage( " --skiplines               Lines at the beginning of the ASCII file to skip" );	
  opt->addUsage( " --azimuth                 Value in range [0,360), clockwise from North" );
  opt->addUsage( " --f_step                  The frequency step" );		
  opt->addUsage( " --f_max                   Maximum frequency to propagate" );
  opt->addUsage( "    Note that in this case the array of frequencies is [f_step:f_step:f_max]." );	  	
  opt->addUsage( "" );  
  opt->addUsage( "OPTIONAL [defaults]:" );
  opt->addUsage( " --f_min                   Minimum frequency [f_step Hz] " ); 
  opt->addUsage( " --maxheight_km            Calculation grid height in km above MSL [150 km]" );
  opt->addUsage( " --zground_km              Height of the ground level above MSL [0 km]" );  
  opt->addUsage( " --Nz_grid                 Number of points on the z-grid from the ground" );
  opt->addUsage( "                           to maxheight [20000]" );
  opt->addUsage( " --sourceheight_km         Source height in km Above Ground Level (AGL) [0]" );
  opt->addUsage( " --receiverheight_km       Receiver height in km AGL [0]" );   
  opt->addUsage( " --maxrange_km             Maximum horizontal distance from origin to propagate" );
  opt->addUsage( "                           [1000 km]" );
  opt->addUsage( " --ground_impedance_model  Name of the ground impedance models to be employed:" );
  opt->addUsage( "                           [rigid], TBD" );
  opt->addUsage( " --Lamb_wave_BC            For a rigid ground: if ==1 it sets" );
  opt->addUsage( "                           admittance= = -1/2*dln(rho)/dz; [ 0 ]" );
  opt->addUsage( " --wind_units              Use it to specify 'kmpersec' if the winds are given" );
  opt->addUsage( "                           in km/s [ mpersec ]" );
  opt->addUsage( " --use_attn_file           Option to specify a file name containing" );
  opt->addUsage( "                           user-provided attenuation coefficients " );
  opt->addUsage( "                           to be loaded instead of the default" );
  opt->addUsage( "                           Sutherland-Bass attenuation. " );
  opt->addUsage( "                           The text file should contain two columns: " );
  opt->addUsage( "                              height (km AGL) and " );
  opt->addUsage( "                              attenuation coefficients in np/m." );

  opt->addUsage( " --c_min                   Specify the minimum phase speed (in m/sec)." );
  opt->addUsage( "                           This is used in conjunction with the --wvnum_filter" );
  opt->addUsage( "                           flag which turns on wavenumber filtering by" );
  opt->addUsage( "                           phase speed. See also the --wvnum_filter flag" );
  opt->addUsage( "                           and the --c_max option." );
  opt->addUsage( " --c_max                   Specify the maximum phase speed (in m/sec)." );
	
  opt->addUsage( "" );
  opt->addUsage( "FLAGS (no value required after the flag itself):" );
  opt->addUsage( " --turnoff_WKB             Turn off the WKB least phase speed estimation," );
  opt->addUsage( "                           an approximation that speeds-up ground-to-ground" );
  opt->addUsage( "                           propagation. If either or both the source and " ); 
  opt->addUsage( "                           receiver are not on the ground [height=0] " );
  opt->addUsage( "                           the approximation is turned off by default" ); 
  opt->addUsage( "                           regardless whether the flag is mentioned" ); 
  opt->addUsage( "                           in the command line or not." ); 
  opt->addUsage( " --use_zero_attn           Set attenuation to zero.");
  opt->addUsage( "                           Reads the dispersion file and sets the");
  opt->addUsage( "                           imaginary part of the wavenumber to zero.");
  opt->addUsage( " --wvnum_filter            Applies wavenumber filtering by phase speed" );
  opt->addUsage( "                           and should be followed by specification of" ); 
  opt->addUsage( "                           the parameters:" );
  opt->addUsage( "                              --c_min   minimum phase speed (in m/sec)." );
  opt->addUsage( "                              --c_max   maximum phase speed (in m/sec)." );
  opt->addUsage( "" );
  opt->addUsage( "" );
  opt->addUsage( "Options for PULSE PROPAGATION:" );
  opt->addUsage( " --pulse_prop_src2rcv <dispersion filename> ");
  opt->addUsage( "                    Propagate pulse from source to 1 receiver");
  opt->addUsage( "                    at a distance specified by option --range_R_km; " );
  opt->addUsage( " --range_R_km       Propagate pulse to this range [km]" );
  opt->addUsage( " --waveform_out_file <waveform filename>   Name of the waveform output file" );
  opt->addUsage( "" );
  opt->addUsage( " --pulse_prop_src2rcv_grid <dispersion filename> ");
  opt->addUsage( "                    Propagate pulse from source to array of ");
  opt->addUsage( "                    horizontally equally-spaced receivers" );  
  opt->addUsage( "" );
  opt->addUsage( " REQUIRED additional options:" );
  opt->addUsage( " --f_center         The center frequency of the built-in pulse choices 1 and 2" );
  opt->addUsage( "                    (f_center<=f_max/5)" );
  opt->addUsage( " --R_start_km       Propagation from this range to R_end_km in DR_km steps." );
  opt->addUsage( " --R_end_km         Pulse is propagated from R_start_km to this range." );
  opt->addUsage( " --DR_km            Range step to propagate from R_start_km to R_end_km." );
  opt->addUsage( " --waveform_out_file <waveform filename> ");
  opt->addUsage( "                    Name of the waveform output file." );
  opt->addUsage( "" );
  opt->addUsage( " OPTIONAL [defaults]:" );
  opt->addUsage( " --max_celerity     Maximum celerity [340 m/s]." );
  opt->addUsage( " --nfft             Number of points used in the FFT computation. ");
  opt->addUsage( "                    Defaults to [4*f_max/f_step]." );		
  
  opt->addUsage( "" );	   
  opt->addUsage( "" );
  opt->addUsage( "SOURCE TYPE input options: Use one of the following 4 options" );
  opt->addUsage( "                           to specify the source type:" );  
  opt->addUsage( " --get_impulse_resp       Flag to use a delta function as source and" );
  opt->addUsage( "                          to output the impulse response." );
  opt->addUsage( " --use_builtin_pulse1     Flag to request the use of the built-in source pulse." );
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

  opt->addUsage( " --use_builtin_pulse2     Flag to request a second shape of built-in source pulse." );
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
  opt->addUsage( "                          Note that the frequencies provided in the file" );
  opt->addUsage( "                          must match the frequencies in the dispersion file." );
  opt->addUsage( " --src_waveform_file      Specify the file name of the user-provided " );
  opt->addUsage( "                          source waveform. The file must have 2 columns" );
  opt->addUsage( "                             |Time | Amplitude |" ); 
  opt->addUsage( "   If none of then source type options are specified the delta function source");
  opt->addUsage( "   is the default i.e. the output is the impulse response." );  
  opt->addUsage( "" );
  opt->addUsage( " Example: Pulse propagation to a point on the ground at range_R_km" ); 
  opt->addUsage( "          and output the impulse response:" );
  opt->addUsage( "" );
  opt->addUsage( "    ../bin/ModBB --pulse_prop_src2rcv myDispersionFile.dat --range_R_km 240 " );
  opt->addUsage( "                 --waveform_out_file mywavf.dat --get_impulse_resp" );
  opt->addUsage( "" );
  opt->addUsage( " Example: Pulse propagation to a point on the ground at range_R_km" );
  opt->addUsage( "          and employ the user-provided source spectrum:" );
  opt->addUsage( "" );
  opt->addUsage( "   ../bin/ModBB --pulse_prop_src2rcv myDispersionFile.dat --range_R_km 240 " );
  opt->addUsage( "                --waveform_out_file mywavf.dat --max_celerity 300 " );
  opt->addUsage( "                --src_spectrum_file source_spectrum_example.dat" );    
  opt->addUsage( "" );
  opt->addUsage( " Example: Pulse propagation to several points on the ground 20 km apart" );
  opt->addUsage( "          and employ the user-provided source waveform:" );
  opt->addUsage( "" );
  opt->addUsage( "   ../bin/ModBB --pulse_prop_src2rcv_grid myDispersionFile.dat " );
  opt->addUsage( "                --R_start_km 240 --DR_km 20 --R_end_km 300 --waveform_out_file mywavf.dat " );
  opt->addUsage( "                --src_waveform_file source_waveform_input_example.dat" );
  opt->addUsage( "" );
  opt->addUsage( "" );
  
/*
// This functionality works but was disabled because of too long computation times.
  opt->addUsage( " To compute a 2D field:" );
  opt->addUsage( " --pulse_prop_grid <dispersion directory name> ");
  opt->addUsage( "                    Compute/view pulse on the 2D spatial x-z grid of 'height_km'");
  opt->addUsage( "                    and 'width_km' starting at 'R_start_km" );
  opt->addUsage( "" );
  opt->addUsage( "   --------------------------------------------------------------" );
  opt->addUsage( " height_km      |                                     |" );
  opt->addUsage( "                |  Pressure field computed within     |" );
  opt->addUsage( "                |  a 2D (width_km x height_km) grid   |" );
  opt->addUsage( "                |          'ntsteps' times            |" );
  opt->addUsage( "                |      every 'tmstep' seconds         |" );
  opt->addUsage( "                |                                     |" );  
  opt->addUsage( "                |                                     |" );	
  opt->addUsage( "   -------------x------------------------------------------------" );	
  opt->addUsage( "              R_start_km" );
  opt->addUsage( "" );
  opt->addUsage(" Additional parameters:" );
  opt->addUsage( " --R_start_km       The grid (viewing window) starts at R_start_km" );
  opt->addUsage( " --width_km         Grid width" );
  opt->addUsage( " --max_celerity     Reference speed [m/s]; in conjunction with R_start_km");
  opt->addUsage( "                    it is determining where inside the grid the field is at");
  opt->addUsage( "                    a time step; a value smaller than the speed of sound");
  opt->addUsage( "                    at the ground is suggested." );
  opt->addUsage( " --tmstep           2D pressure field is calculated at this specified time step." );
  opt->addUsage( " --ntsteps          Number of times the 2D pressure field is calculated");
  opt->addUsage( "                    'tmstep' seconds apart." );
  opt->addUsage( "" );
  opt->addUsage( " OPTIONAL [defaults]:" );
  opt->addUsage( " --height_km        The height of the 2D grid. [maximum height]" );
  opt->addUsage( " --frame_file_stub  Each 2D grid is saved into a file with the name");
  opt->addUsage( "                    frame_file_stub_<time_of_start>; Default:[Pressure2D]." );
  opt->addUsage( "" );
  opt->addUsage( " Example: >> ../bin/ModBB --pulse_prop_grid mydispersionFolder --R_start_km 220 --width_km 50 --height_km 25 --max_celerity 300 --tmstep 30 --ntsteps 5 --frame_file_stub myPressure --use_builtin_pulse" );  
  opt->addUsage( "" );
*/
  

    


  // Set up the actual flags, etc.
  opt->setFlag( "help", 'h' );
  opt->setFlag( "use_modess" );
  opt->setFlag( "use_wmod" );
  opt->setFlag( "get_impulse_resp" );
  opt->setFlag( "use_builtin_pulse1" );
  opt->setFlag( "use_builtin_pulse2" );
  opt->setFlag( "turnoff_WKB");
  opt->setFlag( "plot");
  opt->setFlag( "use_zero_attn");
  opt->setFlag( "wvnum_filter");
  
  opt->setOption( "atmosfile" );
  opt->setOption( "atmosfileorder" );
  opt->setOption( "wind_units" );
  opt->setOption( "skiplines" );		
  opt->setOption( "azimuth" );
  opt->setOption( "maxrange_km" );
  opt->setOption( "sourceheight_km" );
  opt->setOption( "receiverheight_km" );
  opt->setOption( "maxheight_km" );
  opt->setOption( "zground_km" );  
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
  opt->setOption( "out_dispersion_files" );
  opt->setOption( "out_disp_src2rcv_file" );
  opt->setOption( "waveform_out_file" );
  opt->setOption( "width_km" );
  opt->setOption( "height_km" );
  opt->setOption( "tmstep" );
  opt->setOption( "ntsteps" );
  opt->setOption( "frame_file_stub" );
  opt->setOption( "disp_dirname" );
  opt->setOption( "src_spectrum_file" );
  opt->setOption( "src_waveform_file" );
  opt->setOption( "nfft" );
  opt->setOption( "use_attn_file" );
  opt->setOption( "c_min" );
  opt->setOption( "c_max" );

  // Process the command-line arguments
  opt->processFile( "./ModBB.options" );
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


 
