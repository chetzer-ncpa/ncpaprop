#ifndef _ProcessOptionsBB_H_
#define _ProcessOptionsBB_H_

#include "anyoption.h"


namespace NCPA {
  class ProcessOptionsBB {
    public:
      ProcessOptionsBB(AnyOption *opt);

      // get functions
      string getAtmosfile();
      string getAtmosfileorder();   
      string getGnd_imp_model();  
      string getWindUnits();
      string getPprop_s2r_disp_file();
      string getWaveform_out_file(); 
      string getOut_disp_file();
      string getPprop_grid_dirname();
      string getFrame_file_stub();
      string getSrcfile();
      string getUsrAttFile();      
            
      bool   getW_disp_src2rcv_flg();
      bool   getW_disp_flg();
      bool   getPprop_grid_flg();
      bool   getPprop_s2r_flg();
      bool   getPprop_s2r_grid_flg();
      bool   getUsemodess_flg();
      bool   getTurnoff_WKB();
      bool   getPlot_flg();
      bool   getZeroAttn_flg();
    
      int    getNrng_steps();
      int    getNz_grid();  
      int    getLamb_wave_BC();
      int    getNfreq();
      int    getSkiplines();
      int    getNtsteps();
      int    getSrc_flg();
      int    getNFFT();

      double getZ_min(); 
      double getMax_celerity();
      double getRR();
      double getR_start();
      double getR_end();
      double getDR();
      double getAzimuth();
      double getMaxrange();
      double getMaxheight();
      double getSourceheight();
      double getReceiverheight();      
      double getF_min();
      double getF_step();
      double getF_max();
      double getF_center();
      double getR_start_km();
      double getWidth_km();
      double getHeight_km();
      double getC_ref();
      double getTmstep(); 

    private:
      string out_disp_file;
      string pprop_grid_dirname;
      string pprop_s2r_disp_file;
      string frame_file_stub;
      string disp_dirname;  
      string waveform_out_file;
      string atmosfile;       // stores the atmospheric profile name 
      string atmosfileorder;  // order of column names in atmosfile
      string wind_units;      // default mpersec
      string gnd_imp_model;   // ("rigid");
      string srcfile;         // file name of the user-provided source spectrum or source waveform
      string usrattfile;          // user-provided attenuation filename  
      
      bool   w_disp_src2rcv_flg;
      bool   w_disp_flg;
      bool   pprop_grid_flg;
      bool   pprop_s2r_flg;	
      bool   pprop_s2r_grid_flg;
      bool   usemodess_flg;   // if ==1 it uses the Eff. Sound Speed approx.
      bool   usewmod_flg;     // if ==1 prompts the use of wmod (mutually exclusive with usemodess_flg)
      bool   turnoff_WKB;
      bool   plot_flg;
      bool   zero_attn_flg;   // if ==1 sets attenuation to zero
      
      int    Nz_grid;         // number of points on the z-grid
      int    Nrng_steps;      // number of range steps		
      int    Nfreq;           // number of positive frequencies 
      int    skiplines;       // number of lines to skip in atmosfile
      int    ntsteps;         // number of time steps
      int    Lamb_wave_BC;    // if ==1 then admittance = -1/2*dln(rho)/dz
      int    src_flg;         // source flag; 0 for impulse response; 1 for built-in impulse 
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
      double z_min;           // meters
      double azi;             // degrees
      double maxrange;        // meters
      double maxheight;       // meters
      double sourceheight;    // meters
      double receiverheight;  // meters
      double f_min;           // minimum frequency
      double f_step;          // frequency step
      double f_max;           // Nyquist frequency
      double max_cel;         // maximum celerity
      double f_center;        // specifies the center frequency of the pulse
		
  }; // mandatory semicolon here
}

#endif
