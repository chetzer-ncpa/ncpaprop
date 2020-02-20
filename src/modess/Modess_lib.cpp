#include <stdlib.h>
#include <fstream>
#include <cmath>
#include "Modess_lib.h"
#include "Atmosphere.h"

#ifndef Pi
#define Pi 3.141592653589793
#endif

//utility functions

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

