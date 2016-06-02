#include <stdlib.h>
#include <fstream>
#include <cmath>
#include "WMod_lib.h"
#include "Atmosphere.h"
#include "ProcessOptionsNB.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif

// utility functions

double **dmatrix(long nr, long nc) {
  // allocate a double matrix
 	double **v;
	v = new double* [nr];
	for (long i=0; i<nr; i++) {
			v[i] = new double [nc];
	}
	return v;
}

int free_dmatrix(double**v, long nr, long nc) {
  // free a double matrix
	for (long i=0; i<nr; i++) {
			delete v[i];
	}
	delete v;
	return 0;
}


/*
// plot from a GNUPlot script
int plotwGNUplot_from_script() {
  FILE *gnuplotPipe, *scrpt;
  const char *scriptFile = "gnuplot_Modess_output.p";
  
  scrpt = fopen(scriptFile, "r");
  if (scrpt==NULL) {
      perror("in plotwGNUplot(): Error: no gnuplot script found.");
      exit(1);
  }
  gnuplotPipe = popen("gnuplot -persist","w"); //dv
  if (gnuplotPipe) {
      printf("Plotting with gnuplot...\n");
      fprintf(gnuplotPipe,"load \"%s\"\n",scriptFile);
      fflush(gnuplotPipe);
      fprintf(gnuplotPipe,"exit \n");  
      pclose(gnuplotPipe);  
      printf("...the plots should be persistent on the screen.\n");
  } else {        
      printf("gnuplot not found!\n");    
  }      
  return 0;
}


int plotwGNUplot(double freq, bool write_2D_TLoss, bool Nby2Dprop, bool write_phase_speeds) {
  FILE *pipe; 
  if (!Nby2Dprop) {  
      // open a pipe to gnuplot and plot the results
      pipe = popen("gnuplot -persist","w");
      printf("Plotting with gnuplot...\n");              
      fprintf(pipe, "set title 'WMod (Wide Angle - High Mach) - Transmission Loss; Freq = %g Hz'\n", freq);
      fprintf(pipe, "set xlabel 'Range [km]'\n");
      fprintf(pipe, "set ylabel 'TL [dB]'\n");
      fprintf(pipe, "set grid; show grid;\n");
      fprintf(pipe, "plot './wtloss_1d.lossless.nm' using 1:(10*log10($2**2 + $3**2)) with lines lt 3 title 'lossless' , \\\n");
      fprintf(pipe, "    './wtloss_1d.nm' using 1:(10*log10($2**2 + $3**2)) with lines lt 1 title 'lossy'\n");      
      pclose(pipe);
  }
  
  if (write_2D_TLoss & (!Nby2Dprop)) {
      //To plot a surface plot:
      pipe = popen("gnuplot -persist","w");
      //fprintf(pipe,"set term wxt 1\n");
      fprintf(pipe,"set pm3d map\n");
      fprintf(pipe,"set cbrange [-140:-100]\n");
      fprintf(pipe,"set xlabel 'Range [km]'\n");
      fprintf(pipe,"set ylabel 'Height [km]'\n");
      fprintf(pipe, "set title 'WMod (Wide Angle - High Mach) - Transmission loss; Frequency %g Hz'\n", freq);
      fprintf(pipe,"splot './wtloss_2d.nm' using 1:2:(20*log10(sqrt($3**2 + $4**2)))\n");
      pclose(pipe);
  }
  
  
  if (Nby2Dprop) {
      pipe = popen("gnuplot -persist","w");
      //fprintf(pipe,"set term wxt 1\n");
      fprintf(pipe,"set pm3d map\n");
      fprintf(pipe,"set palette rgb  7, 5,15;\n");
      fprintf(pipe,"set cbrange [-140:-100]\n");
      fprintf(pipe,"set xlabel 'km'\n");
      fprintf(pipe,"set ylabel 'km'\n");
      fprintf(pipe, "set title 'WMod (Wide Angle - High Mach) - Transmission loss; Frequency %g Hz'\n", freq);
      fprintf(pipe,"splot './Nby2D_wtloss_1d.nm' using ($1*(sin($2*pi/180))):($1*(cos($2*pi/180))):(10*log10($3**2 + $4**2)) title ''\n");    
      pclose(pipe);        
  
  }      
  
  
  if (write_phase_speeds) {
      // open a pipe to gnuplot and plot the results
      pipe = popen("gnuplot -persist","w");
      printf("Plotting phase speeds...\n");
      fprintf(pipe, "set title 'Phase speeds; Frequency %g Hz'\n", freq);
      fprintf(pipe, "set xlabel 'Mode #'\n");
      fprintf(pipe, "set ylabel 'm/s'\n");
      fprintf(pipe, "set grid; show grid;\n");
      fprintf(pipe, "plot './wphasespeeds.nm' using 1:2 with lines; \\\n");      
      pclose(pipe);
  }     
  return 0;
}


int plotwGNUplot(NCPA::ProcessOptionsNB *oNB) {
  FILE *pipe; 
  if (!oNB->getNby2Dprop()) {  
      // open a pipe to gnuplot and plot the results
      pipe = popen("gnuplot -persist","w");
      printf("Plotting with gnuplot...\n");
      fprintf(pipe, "set title 'WMod (Wide Angle - High Mach) - Transmission loss; Freq. %g Hz'\n", oNB->getFreq());
      fprintf(pipe, "set xlabel 'Range [km]'\n");
      fprintf(pipe, "set ylabel 'TL [dB]'\n");
      fprintf(pipe, "set grid; show grid;\n");              
      fprintf(pipe, "plot './wtloss_1d.lossless.nm' using 1:(10*log10($2**2 + $3**2)) with lines lt 3 title 'lossless' , \\\n");
      fprintf(pipe, "    './wtloss_1d.nm' using 1:(10*log10($2**2 + $3**2)) with lines lt 1 title 'lossy'\n");        

      pclose(pipe);
  }
  
  if (oNB->getWrite_2D_TLoss() & (!oNB->getNby2Dprop())) {
      //To plot a surface plot:
      pipe = popen("gnuplot -persist","w");
      //fprintf(pipe,"set term wxt 1\n");
      fprintf(pipe,"set pm3d map\n");
      fprintf(pipe,"set cbrange [-140:-100]\n");
      fprintf(pipe,"set xlabel 'Range [km]'\n");
      fprintf(pipe,"set ylabel 'Height [km]'\n");
      fprintf(pipe, "set title 'NWMod (Wide Angle - High Mach) - Transmission loss; Freq. %g Hz'\n", oNB->getFreq());
      fprintf(pipe,"splot './wtloss_2d.nm' using 1:2:(20*log10(sqrt($3**2 + $4**2)))\n");
      pclose(pipe);
  }
  
  
  if (oNB->getNby2Dprop()) {
      pipe = popen("gnuplot -persist","w");
      //fprintf(pipe,"set term wxt 1\n");
      fprintf(pipe,"set pm3d map\n");
      fprintf(pipe,"set palette rgb  7, 5,15;\n");
      fprintf(pipe,"set cbrange [-140:-100]\n");
      fprintf(pipe,"set size square\n");
      fprintf(pipe,"set xlabel 'km'\n");
      fprintf(pipe,"set ylabel 'km'\n");
      fprintf(pipe, "set title 'WMod (Wide Angle - High Mach) - Transmission loss; Freq. %g Hz'\n", oNB->getFreq());
      fprintf(pipe,"splot './Nby2D_wtloss_1d.nm' using ($1*(sin($2*pi/180))):($1*(cos($2*pi/180))):(10*log10($3**2 + $4**2)) title ''\n");
      pclose(pipe);        
  
  }
 
  if (oNB->getWrite_phase_speeds()) {
      // open a pipe to gnuplot and plot the results
      pipe = popen("gnuplot -persist","w");
      printf("Plotting phase speeds...\n");
      //fprintf(pipe, "set data style lines\n");
      //fprintf(pipe, "set yrange [%f:%f]\n",min(alt,nAlt),max(alt,nAlt));       
      fprintf(pipe, "set title 'WMod (Wide Angle - High Mach) - Phase speeds; Freq. %g Hz'\n", oNB->getFreq());
      fprintf(pipe, "set xlabel 'Mode #'\n");
      fprintf(pipe, "set ylabel 'm/s'\n");
      fprintf(pipe, "set grid; show grid;\n");
      fprintf(pipe, "plot './wphasespeeds.nm' using 1:2 with lines; \\\n");
      pclose(pipe);
  }     
  return 0;
}
*/


/*
int saveAtm_profile(NCPA::SampledProfile *p) {
  int i, nz;
  double z, azi; //dz_km;
  nz  = p->nz();
  azi = p->getPropagationAzimuth();
  FILE *fp = fopen("atm_profile.nm", "w");
  
  for (i=0; i<nz; i++) {
      z = p->z(i);
      fprintf(fp, "%f %f %f %f %f %f %f %f %f\n",\
          z, p->u(z), p->v(z), p->w(z), p->t(z), p->rho(z), p->p(z), \
          p->c0(z)*1000.0, p->ceff(z, azi)*1000.0);
  }
  fclose(fp);
  printf("Interpolated atmospheric profiles (z,u,v,w,t,d,p,c,c_eff) saved in: atm_profile.nm\n");
  return 0;
}
*/


int saveAtm_profile(NCPA::SampledProfile *p, std::string wind_units) {
  int i, nz;
  double z, u, v, c, ceff, azi_rad; //dz_km; 
  double kmps2mps = 1.0;
  if (!wind_units.compare("kmpersec")) {
      kmps2mps = 1000.0;
  }
  
  nz  = p->nz();
  azi_rad = p->getPropagationAzimuth()*Pi/180.0;
  
  FILE *fp = fopen("atm_profile.nm", "w");
  for (i=0; i<nz; i++) {
      z    = p->z(i);
      u    = p->u(z)*kmps2mps;
      v    = p->v(z)*kmps2mps;
      c    = p->c0(z)*1000.0;
      ceff = c + u*sin(azi_rad) + v*cos(azi_rad);
      fprintf(fp, "%9.3f %10.3e %10.3e %10.3e %9.3f %10.3e %9.3f %8.3f %8.3f\n",\
          z, u, v, p->w(z), p->t(z), p->rho(z), p->p(z), p->c0(z)*1000.0, ceff);
  }
  fclose(fp);
  printf("Interpolated atmospheric profiles (z,u,v,w,t,d,p,c,c_eff) saved in: atm_profile.nm\n");
  return 0;
}


