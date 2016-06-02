#ifndef _PROCESSOPTIONSNB_H_
#define _PROCESSOPTIONSNB_H_

#include "anyoption.h"

namespace NCPA {
  class ProcessOptionsNB {
    public:
      // constructor
      ProcessOptionsNB(AnyOption *opt);

      // get functions 
      string   getAtmosfile();
      string   getAtmosfileorder();
      string   getAtm_profile_dir();   
      string   getGnd_imp_model();
      string   getProfileRanges();  
      string   getWindUnits();
      string   getUsrAttFile();      

      int      getFiletype();
      int      getSkiplines();
      int      getNrng_steps();
      int      getNz_grid();
      int      getLamb_wave_BC();
            
      double   getFreq();
      double   getAzimuth();
      double   getMaxrange();
      double   getMaxheight();
      double   getSourceheight();
      double   getReceiverheight();
      double   getSlepcTolerance();
      double   getZ_min(); 
      double   getMax_celerity();
      double   getReq_profile_step();
                 
      bool     getWrite_2D_TLoss();
      bool     getWrite_phase_speeds();
      bool     getWrite_modes();
      bool     getWrite_dispersion(); 
      bool     getProfile_ranges_given_flag();
      bool     getTurnoff_WKB();
      bool     getPlot_flg();
      
      // print parameters
      void     printParams();

	
    private:
      string   atmosfile;           // stores the atmospheric profile name 
      string   atmosfileorder;      // order of column names in atmosfile
      string   atm_profile_dir;
      string   wind_units;          // m/s as default
      string   gnd_imp_model;       // ("rigid");
      string   prf_ranges_km;       // string specifying profile ranges
      string   usrattfile;          // user-provided attenuation filename

      int      filetype;            // the filetype: atmosfile, slicefile, etc.
      int      Nz_grid;             // number of points on the z-grid
      int      Nrng_steps;          // number of range steps		
      int      Nfreq;               // number of positive frequencies 
      int      skiplines;           // number of lines to skip in "atmosfile"
      int      Lamb_wave_BC;        // for rigid ground: if ==1 then admittance = -1/2*dln(rho)/dz
             
      double   freq;                // Hz	
      double   z_min;               // meters
      double   azi;                 // degrees
      double   maxrange;            // meters
      double   maxheight;           // meters
      double   sourceheight;        // meters
      double   receiverheight;      // meters
      double	 tol;                 // tolerance for Slepc calculations
      double   req_profile_step;    // the profiles are requested at equidistant intervals specified by this number      
           
      bool     write_2D_TLoss;
      bool     write_phase_speeds;
      bool     write_modes;
      bool     write_dispersion;
      bool     profile_ranges_given;// flag signaling that prf_ranges_km are given
      bool     turnoff_WKB;
      bool     plot_flg;
 
	}; // mandatory semicolon here
}

#endif
