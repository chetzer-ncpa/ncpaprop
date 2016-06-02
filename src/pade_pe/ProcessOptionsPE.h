#ifndef _ProcessOptionsPE_H_
#define _ProcessOptionsPE_H_

#include "anyoption.h"

namespace NCPA {
  class ProcessOptionsPE {
    public:
      // constructor
      ProcessOptionsPE(AnyOption *opt);
      
      // set functions
      void   setMaxheight(double newmax);

      // get functions
      void     printParams();
      
      string   getAtmosfile();
      string   getAtmosfileorder();
      string   getAtm_profile_dir();  
      string   getGrnd_imp_model();
      string   getWindUnits();
      string   getProfileRanges(); 
      string   getStarterType();
      string   getUsrAttFile();
      string   getModalStarterFile();
      //string   getWindUnits();
            
      int      getFiletype();
      int      getSkiplines();
      int      getNz_grid();
      int      getNpade();
            
      double   getFreq();
      double   getAzimuth();
      double   getMaxrange();
      double   getMaxheight();
      double   getSourceheight();
      double   getReceiverheight();
      double   getRngStep();
      double   getZ_min();
      double   getReq_profile_step();
      double   getMax_celerity();
               
      bool     getWrite_2D_TLoss();
      bool     getProfile_ranges_given_flag();
      bool     getNoabsorption();
      bool     getPlot_flg();

    private:
      string   atmosfile;           // stores the atmospheric profile name 
      string   atmosfileorder;      // order of column names in atmosfile
      string   atm_profile_dir;     // directory containing ASCII profiles
      string   prf_ranges_km;       // string specifying profile ranges
      string   grnd_imp_model;      // ("rigid");
      string   wind_units;          // default 'mpersec'
      string   starter_type;        // PE starter field type
      string   usrattfile;          // user-provided attenuation filename
      string   modstartfile;

      int      filetype;            // the filetype: atmosfile, slicefile, etc.
      int      Nz_grid;             // number of points on the z-grid		
      int      Nfreq;               // number of positive frequencies 
      int      skiplines;           // number of lines to skip in "atmosfile" 
      int      n_pade;              // number of Pade coefficients
                   
      double   freq;                // Hz	
      double   z_min;               // meters
      double   azi;                 // degrees
      double   maxrange;            // meters
      double   maxheight;           // meters
      double   sourceheight;        // meters
      double   receiverheight;      // meters
      double   rng_step;            // units of wavelength
      double   req_profile_step;    // the profiles are requested at equidistant intervals specified by this number      
      
      bool     ncpatoy;             // flag to use built-in NCPA canonical profiles   
      bool     write_2D_TLoss;      // compute 2D TL
      bool     profile_ranges_given;// flag signaling that prf_ranges_km are given    
      bool     do_lossless;         // flag to compute the with no atm. absorption
      bool     plot_flg;
	}; // mandatory semicolon here
}

#endif
