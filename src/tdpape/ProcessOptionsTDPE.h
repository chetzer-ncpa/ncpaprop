#ifndef _ProcessOptionsTDPE_H_
#define _ProcessOptionsTDPE_H_

#include "anyoption.h"


namespace NCPA {
  class ProcessOptionsTDPE {
    public:
      ProcessOptionsTDPE(AnyOption *opt);

      // get functions
      string getPape_out_dir();
      
      string getWaveform_out_file(); 
      string getSrcfile();         
            
      //bool   getPprop_grid_flg();
      bool   getPprop_s2r_flg();
      bool   getPprop_s2r_grid_flg();
      bool   getPlot_flg();
    
      int    getNrng_steps();
      int    getNfreq();
      int    getNtsteps();
      int    getSrc_flg();
      int    getNFFT();

      double getZ_min(); 
      double getMax_celerity();
      double getRR();
      double getR_start();
      double getR_end();
      double getDR();
      double getMaxrange();    
      double getF_min();
      double getF_step();
      double getF_max();
      double getF_center();
      double getR_start_km();
      double getTmstep(); 

    private:
      string  pape_out_dir;

      string waveform_out_file;
      string srcfile;         // file name of the user-provided source spectrum or source waveform       
      
      bool   pprop_s2r_flg;	
      bool   pprop_s2r_grid_flg;
      bool   plot_flg;

      int    Nrng_steps;      // number of range steps		
      int    Nfreq;           // number of positive frequencies 
      int    ntsteps;         // number of time steps
      int    src_flg;         // source flag; 0 for impulse response; 1 for built-in inpulse 
                              // 2 for source spectrum file; 3 for source waveform file provided
      int    NFFT;            // number of fft points 
      
      double RR;
      double R_start;
      double R_end;
      double DR;
      double R_start_km;
      double width_km;
      double height_km;
      double tmstep;
      double maxrange;        // meters

      double f_min;           // minimum frequency
      double f_step;          // frequency step
      double f_max;           // Nyquist frequency
      double max_cel;         // maximum celerity
      double f_center;        // specifies the center frequency of the pulse
		
  }; // mandatory semicolon here
}

#endif
