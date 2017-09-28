#ifndef _ProcessOptionsNB_H_
#define _ProcessOptionsNB_H_

#include "anyoption.h"

namespace NCPA {
  class ProcessOptionsNB {
    public:
      // constructor
      ProcessOptionsNB(AnyOption *opt);

      string   getAtmosfile();
      string   getAtmosfileorder();   
      string   getGnd_imp_model();  
      string   getWindUnits();
      string   getUsrAttFile();      
               
      int      getSkiplines();
      int      getNrng_steps();
      int      getNz_grid();
      int      getLamb_wave_BC();
           
      double   getFreq();
      double   getAzimuth();
      double   getAzimuthStart();
      double   getAzimuthEnd();
      double   getAzimuthStep();      
      double   getMaxrange();
      double   getMaxheight();
      double   getSourceheight();
      double   getReceiverheight();
      double   getSlepcTolerance();
      double   getZ_min();
      double   getMax_celerity();
      double   getC_min();
      double   getC_max();
             
      bool     getWrite_2D_TLoss();
      bool     getWrite_phase_speeds();
      bool     getWrite_modes();
      bool     getWrite_dispersion();
      bool     getNby2Dprop();
      bool     getWriteAtmProfile();
      bool     getTurnoff_WKB();
      bool     getPlot_flg();
      bool     getWvnum_filter_flg();

	
    private:
      string   atmosfile;           // stores the atmospheric profile name 
      string   atmosfileorder;      // order of column names in atmosfile
      string   gnd_imp_model;       // ("rigid");
      string   wind_units;          // default mpersec
      string   usrattfile;          // user-provided attenuation filename
     
      bool     write_2D_TLoss;
      bool     write_phase_speeds;
      bool     write_modes;
      bool     write_dispersion;
      bool     Nby2Dprop;
      bool     write_atm_profile;
      bool     turnoff_WKB;
      bool     plot_flg;
      bool     wvnum_filter_flg;    // wavenumber filtering flag
      
      int      Nz_grid;             // number of points on the z-grid
      int      Nrng_steps;          // number of range steps		
      int      Nfreq;               // number of positive frequencies 
      int      skiplines;           // number of lines to skip in "atmosfile"
      int      Lamb_wave_BC;        // for rigid ground: if ==1 then admittance = -1/2*dln(rho)/dz
           
      double   freq;                // Hz	
      double   z_min;               // meters
      double   azi;                 // degrees
      double   azi_min;             // degrees
      double   azi_max;             // degrees
      double   azi_step;            // degrees
      double   maxrange;            // meters
      double   maxheight;           // meters
      double   sourceheight;        // meters
      double   receiverheight;      // meters
      double   tol;                 // tolerance for Slepc calculations
      double   c_min;               // minimum sound speed requested by user to do wavenumber filtering
      double   c_max;               // maximum sound speed requested by user to do wavenumber filtering



	}; // mandatory semicolon here
}

#endif
