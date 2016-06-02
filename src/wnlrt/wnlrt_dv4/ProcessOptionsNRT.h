#ifndef _PROCESSOPTIONSNRT_H_
#define _PROCESSOPTIONSNRT_H_

#include "anyoption.h"

namespace NCPA {
  class ProcessOptionsNRT {
    public:
      // constructor
      ProcessOptionsNRT(AnyOption *opt);

      // get functions
      string   getAtmosfile();
      string   getAtmosfileorder(); 
        
      string   getGnd_imp_model();  
      string   getWindUnits();
      string   getUsrAttFile();
         
      int    getSkiplines();
      int    getNrng_steps();
      int    getNz_grid();
      int    getLamb_wave_BC();       
      double getFreq();
      
      double getAzimuthStart();
      double getAzimuthEnd();
      double getAzimuthStep();
      double getMaxrange();
      double getMaxheight();
      double getSourceheight();
      double getReceiverheight();
      double getSlepcTolerance();
      double getZ_min(); 
      double getMax_celerity();           
      bool   getWrite_2D_TLoss();
      bool   getWrite_phase_speeds();
      bool   getWrite_speeds();
      bool   getWrite_modes();
      bool   getWrite_dispersion();
      bool   getNby2Dprop();
      bool   getWriteAtmProfile();
      bool   getTurnoff_WKB();
      bool   getPlot_flg();
      
      
      
      
      string getEigenfile();
      string getWftype();
      double getInclination();
      double getAzimuth();
      double getSrcz();
      double getRcvx();
      double getRcvy();
      double getRcvz();
      double getRange();
      double getDth();
      double getDaz();
      double getTol();
      
      double getWfampl();
      double getWfduration();
      
      bool getFindEigenray();
	
    private:

      string   gnd_imp_model;       // ("rigid");
      string   wind_units;          // default mpersec
      string   usrattfile;          // user-provided attenuation filename
      int      Nz_grid;             // number of points on the z-grid
      int      Nrng_steps;          // number of range steps		
      int      Nfreq;               // number of positive frequencies 
      int      skiplines;           // number of lines to skip in "atmosfile"
      int      Lamb_wave_BC;        // for rigid ground: if ==1 then admittance = -1/2*dln(rho)/dz        
      double   freq;                // Hz	
      double   z_min;               // meters
      
      
 
      string   atmosfile;           // stores the atmospheric profile name 
      string   atmosfileorder;      // order of column names in atmosfile
      string   eigenfile;
      string   wftype;
           
      double theta; //degrees
      double   azi;                 // degrees
      double src_z;
      double rcv_x;
      double rcv_y;
      double rcv_z;
      double range;
      double dth;
      double daz;
      double tol;
      double wfAmpl;
      double wfDuration;
      
      
      bool findEigenray;
      //bool useEigenfile;
      bool shootRay;
      
      
      
      
      
      
      
      double   azi_min;             // degrees
      double   azi_max;             // degrees
      double   azi_step;            // degrees
      double   maxrange;            // meters
      double   maxheight;           // meters
      double   sourceheight;        // meters
      double   receiverheight;      // meters

      bool     write_2D_TLoss;
      bool     write_phase_speeds;
      bool     write_speeds;
      bool     write_modes;
      bool     write_dispersion;
      bool     Nby2Dprop;
      bool     write_atm_profile;
      bool     turnoff_WKB;
      bool     plot_flg;
	}; // mandatory semicolon here
}

#endif
