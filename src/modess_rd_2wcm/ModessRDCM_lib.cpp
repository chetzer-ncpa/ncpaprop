#include <iostream>
#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <cmath>
#include <complex>
#include <dirent.h>
#include <vector>
#include <string>
#include "Atmosphere.h"
#include "ModessRDCM_lib.h"
#include "binaryreader.h"

#include <petscksp.h>
#include <sys/stat.h>
#include <ctime>
//#include <vector>

#ifndef Pi
#define Pi 3.141592653589793
#endif

using namespace NCPA;
using namespace std;

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

complex<double> **cmatrix(long nr, long nc) {
  // allocate a complex<double> matrix
  complex<double> **v;
  v = new complex<double>* [nr];
  for (long i=0; i<nr; i++) {
      v[i] = new complex<double> [nc];
  }
  return v;
}

int free_cmatrix(complex<double>**v, long nr, long nc) {
  // free a complex<double> matrix
  for (long i=0; i<nr; i++) {
      delete v[i];
  }
  delete v;
  return 0;
}


/*
int plotwGNUplot(double freq, bool write_2D_TLoss) {
    if (1) {      
        string str1 = "./tloss_rd2wcm_1d.lossless.nm";
        string str2 = "./tloss_rd2wcm_1d.nm";
        // open a pipe to gnuplot and plot the results
        FILE *pipe = popen("gnuplot -persist","w");
        printf("Plotting with gnuplot...\n");
        //fprintf(pipe, "set data style lines\n");
        //fprintf(pipe, "set yrange [%f:%f]\n",min(alt,nAlt),max(alt,nAlt));        
        fprintf(pipe, "set title 'Range-Dep. Two-Way Coupled Modes: Transmission loss; Frequency %g Hz'\n", freq);
        fprintf(pipe, "set xlabel 'Range [km]'\n");
        fprintf(pipe, "set ylabel 'TL [dB]'\n");
        fprintf(pipe, "set grid; show grid;\n");
        fprintf(pipe, "plot   '%s' using 1:(10*log10($2**2 + $3**2)) with lines lt 3 title 'lossless'\n", str1.c_str());         
        fprintf(pipe, "replot   '%s' using 1:(10*log10($2**2 + $3**2)) with lines lt 1 title 'lossy'\n", str2.c_str());        
        pclose(pipe);
        
        if (write_2D_TLoss) {
            //To plot a surface plot:
            pipe = popen("gnuplot -persist","w");
            //fprintf(pipe,"set term wxt 1\n");
            fprintf(pipe,"set pm3d map\n");
            fprintf(pipe,"set cbrange [-140:-100]\n");
            fprintf(pipe,"set xlabel 'Range [km]'\n");
            fprintf(pipe,"set ylabel 'Height [km]'\n");
            fprintf(pipe, "set title 'Range-Dep. Two-Way Coupled Modes: Transmission loss (no attenuation case); Frequency %g Hz'\n", freq);
            fprintf(pipe,"splot './tloss_rd2wcm_2d.lossless.nm' using 1:2:(20*log10(sqrt($3**2 + $4**2))) title ''\n");
            pclose(pipe);
        }
        //printf("Plotting done.\n");
    }
  return 0;
}
*/


NCPA::SampledProfile * get_RngDepnd_profile(string env_file, double R_meters) {

  double R_km = R_meters/1000.0;

  ifstream *prof_stream = new ifstream( env_file.c_str(), ios_base::in );
  if (!prof_stream->good()) {
    cerr << "Problem with input stream from file " << env_file << ".  Exiting..." << endl;
  }

  BinaryReader *binread = new BinaryReader();
  int *envinfo = new int[ 2 ];
  binread->readLittleIntArray( prof_stream, 2, envinfo );
  int nProfiles = envinfo[0];
  int nAlt = envinfo[1];
  delete [] envinfo;

  double *lats      = new double[ nProfiles ];
  double *lons      = new double[ nProfiles ];
  double *bazs      = new double[ nProfiles ];
  double *ranges_km = new double[ nProfiles ];
  double *ground_m  = new double[ nProfiles ];
  double *alts_km   = new double[ nAlt ];

  binread->readLittleDoubleArray( prof_stream, nProfiles, lats );
  binread->readLittleDoubleArray( prof_stream, nProfiles, lons );
  binread->readLittleDoubleArray( prof_stream, nProfiles, bazs );
  binread->readLittleDoubleArray( prof_stream, nProfiles, ranges_km );
  binread->readLittleDoubleArray( prof_stream, nProfiles, ground_m );
  binread->readLittleDoubleArray( prof_stream, nAlt, alts_km );

  //cout << "Profiles: " << nProfiles << endl << "Altitude Points: "
  //     << nAlt << endl;
  //cout << "Range step = "    << ranges_km[1] - ranges_km[0] << " km" << endl;
  //cout << "Altitude step = " << alts_km[1] - alts_km[0]     << " km" << endl;

  // find the profile right before R_km; if R_km> any ranges_km then use the last profile
  int i;
  double baz = 0.0; // backazimuth
  int nthprof = nProfiles-1;
  for (i=0; i<nProfiles; i++) {
    if (ranges_km[i]>R_km) {
        nthprof = i-1;  
        baz = bazs[nthprof];
        cout << "Using profile # " << nthprof
             << " given at range " << ranges_km[i-1] << " km;" 
             << " backazimuth = " << baz << " degrees" << endl;          
        break;
    }
  }
  //cout << "R_km = " << R_km << "; using profile nthprof= " << nthprof <<  endl;

  delete [] lats;
  delete [] lons;
  delete [] bazs;
  delete [] ranges_km;
  delete [] ground_m;

  //double *alt            = new double[ nAlt ];
  double *temp_K         = new double[ nAlt ];
  double *density_gpcm3  = new double[ nAlt ];
  double *pressure_hPa   = new double[ nAlt ];
  double *wind_along_mps = new double[ nAlt ];
  double *wind_cross_mps = new double[ nAlt ];
  double *wind_vert_mps  = new double[ nAlt ];

  // Now we read in each of the applicable data vectors.  There is only
  // one altitude vector but there are nProfiles vectors for each of the
  // other quantities.  So, we have to skip around a bit.
  //binread->readLittleDoubleArray( prof_stream, nAlt, alt );

  int bytesToSkipBefore = 8 * nthprof * nAlt;
  int bytesToSkipAfter  = 8 * nAlt * (nProfiles - nthprof - 1);

  prof_stream->seekg( bytesToSkipBefore, ios_base::cur );
  binread->readLittleDoubleArray( prof_stream, nAlt, temp_K );
  prof_stream->seekg( bytesToSkipAfter, ios_base::cur );

  prof_stream->seekg( bytesToSkipBefore, ios_base::cur );
  binread->readLittleDoubleArray( prof_stream, nAlt, density_gpcm3 );
  prof_stream->seekg( bytesToSkipAfter, ios_base::cur );

  prof_stream->seekg( bytesToSkipBefore, ios_base::cur );
  binread->readLittleDoubleArray( prof_stream, nAlt, pressure_hPa );
  prof_stream->seekg( bytesToSkipAfter, ios_base::cur );

  prof_stream->seekg( bytesToSkipBefore, ios_base::cur );
  binread->readLittleDoubleArray( prof_stream, nAlt, wind_along_mps );
  prof_stream->seekg( bytesToSkipAfter, ios_base::cur );

  prof_stream->seekg( bytesToSkipBefore, ios_base::cur );
  binread->readLittleDoubleArray( prof_stream, nAlt, wind_cross_mps );
  prof_stream->seekg( bytesToSkipAfter, ios_base::cur );

  prof_stream->seekg( bytesToSkipBefore, ios_base::cur );
  binread->readLittleDoubleArray( prof_stream, nAlt, wind_vert_mps );
  prof_stream->seekg( bytesToSkipAfter, ios_base::cur );

  prof_stream->close();

  // convert winds to units of km/s  as required in SampledProfile; this should be changed!!!
  double *uu, *vv, *ww;

  uu = new double [nAlt];
  vv = new double [nAlt];
  ww = new double [nAlt];
  double az = (baz+180.0)*Pi/180.0;
  
  /*
  for (i=0; i<nAlt; i++) {
      uu[i] = wind_along_mps[i]/1000.0*sin(az) - wind_cross_mps[i]/1000.0*cos(az);
      vv[i] = wind_along_mps[i]/1000.0*cos(az) + wind_cross_mps[i]/1000.0*sin(az);
      ww[i] = wind_vert_mps[i]/1000.0;
      alts_km[i] = alts_km[i]-alts_km[0]; // adjust altitudes to start from zero (no terrain considered yet)
      //printf("i=%d; alt=%g \n", i, alts_km[i]);
  }
  */
  
  // if we want to extract the along and cross winds
  //printf("in get_RngDepnd_profile(): saving the along- and cross-winds in units of km/sec\n");
  for (i=0; i<nAlt; i++) {
      uu[i] = wind_along_mps[i]; ///1000.0;
      vv[i] = wind_cross_mps[i]; ///1000.0;
      ww[i] = wind_vert_mps[i]; ///1000.0;
      alts_km[i] = alts_km[i]-alts_km[0]; // adjust altitudes to start from zero (no terrain considered yet)
      //printf("i=%d; alt=%g \n", i, alts_km[i]);
  }  

  // assemble the Sampled profile and set its azimuth
  SampledProfile *atm_prof;
  double z0 = 0.0; // for flat ground use z0=0
  atm_prof = new SampledProfile(nAlt, alts_km, temp_K, uu, vv, ww, density_gpcm3, pressure_hPa, z0);
  atm_prof->setPropagationAzimuth(az*180.0/Pi);

  if (0) {
      saveSampledProfile("./profile.int", atm_prof);
      printf("in get_RngDepnd_profile(): saving the interpolated profile to file profile.int\n");
  }

  delete prof_stream;;
  delete binread;
  delete [] alts_km;
  delete [] temp_K;
  delete [] density_gpcm3;
  delete [] pressure_hPa;
  delete [] wind_along_mps;
  delete [] wind_cross_mps;
  delete [] wind_vert_mps; 

  delete[] uu;
  delete[] vv;
  delete[] ww;

  return atm_prof;
}


int getFile_list(std::string dir, std::list<string> &files, std::string pattern) {
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
      a   = string(dirp->d_name);
      pos = a.find(pattern);
      //cout << "pos = " << pos <<  endl;
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


NCPA::SampledProfile * get_RngDepndProfiles_ascii( \
              int N, string atmosfileorder, uint skiplines, \
              string dirname, string pattern, bool print_flg) {
              
// get the Nth profile stored in a file named "profile000M" residing in a specified
// directory. If N is greater than the number of files found then the last file
// is used repeatedly as necessary             

  list<string> files;
  list<string>::iterator it;
  
  // get and sort the files (they should have names such as profile0001.dat, etc)
  getFile_list(dirname, files, pattern);
  files.sort();
  if (print_flg) { //print the sorted file list
      cout << "sorted atm. file list contains:" << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }
  
  int j = 1;
  string s;
  if (N<(int)(files.size()+1)) {
      for (it=files.begin(); it!=files.end(); ++it) {
          if (j>=N) {
              //cout << "j= " << j << ";  *it= " << *it << endl;
              s = (*it);
              break;
          }
          j++;
      }
  }
  else {
      s = files.back();
  } 

  string atmosfile = dirname + "/" + s;
  cout << "Making SampledProfile object from " << atmosfile << endl;

  SampledProfile *atm_prof;
  atm_prof = new SampledProfile(atmosfile, atmosfileorder.c_str(), skiplines);

  return atm_prof;
}


int saveSampledProfile(string filename, SampledProfile *p) {

  int i, Nz;
  double *zz, *uu, *vv, *ww, *tt, *dd, *pp;
  
  printf("Nz = %d\n", p->nz());
  Nz = p->nz();
  
  // read/store the Nth profile
  zz = new double [Nz];
  uu = new double [Nz];
  vv = new double [Nz];
  ww = new double [Nz];
  tt = new double [Nz];
  dd = new double [Nz];
  pp = new double [Nz];
  
  p->get_z(zz, Nz);
  p->get_u(uu, Nz);
  p->get_v(vv, Nz);
  p->get_w(ww, Nz);
  p->get_t(tt, Nz);
  p->get_rho(dd, Nz);
  p->get_p(pp, Nz);
  
  FILE *fp;
  
  fp = fopen(filename.c_str(), "w");
  for (i=0; i<Nz; i++) {
    fprintf(fp, "%8.3f  %8.3e %8.3e %8.3e %8.3f  %8.4e  %8.4e\n", zz[i], uu[i], vv[i], ww[i], tt[i], dd[i], pp[i]);
  }
  fclose(fp);
  
  delete[] zz;
  delete[] uu;
  delete[] vv;
  delete[] ww;
  delete[] tt;
  delete[] dd;
  delete[] pp;

  return 0;
} 


/*
  int getR1234(SampledProfile *atm_profile1, SampledProfile *atm_profile2, SolveModNBRDCM *a1, SolveModNBRDCM *a2, int Nz_grid, int Nm, double dz, double Rj1, double Rj2);
{
  int i,l,m;
  double *rho_curr, v1v2dz, d1ovd2, sRj2ovRj1;
  complex<double> C_lm, k1ovk2, I_Rj1_Rj2, *H1;
  complex<double> I (0.0, 1.0);

  rho_curr = new double [Nz_grid];
  for (i=0; i<Nz_grid; i++) {
      rho_curr[i] = atm_profile1->rho(i*dz/1000.0)*1000.0; // store density in rho_prev
  }
  
  rho_next = new double [Nz_grid];
  for (i=0; i<Nz_grid; i++) {
      rho_next[i] = atm_profile2->rho(i*dz/1000.0)*1000.0; // store density in rho_prev
  }  

  sRj2ovRj1 = sqrt(Rj2/Rj1);
  I_Rj1_Rj2 = I*(Rj1-Rj2);
  
  // the trace of H1 (H1 is diagonal)
  H1 = new complex<double> [Nm];
  for (m=0; m<Nm; m++) {
      H1[m] = sRj2ovRj1*exp(I_Rj1_Rj2*k_prev[m]); 
      //printf("H1[%d]=%g %g\n", m, real(H1[m]), imag(H1[m]));
  }
  
  
  // compute A_curr = R1*A_prev; the matrix R1 = 1/2(C_tilde+ C_hat)*H1
  for (l=0; l<Nm; l++) {    
      A_curr[l] = 0.0; 
      for (m=0; m<Nm; m++) {
          k1ovk2 = k_prev[m]/k_curr[l];
          C_lm = 0.0; // 1/2(C_tilde + C_hat);
          for (i=0; i<Nz_grid; i++) { // compute integrals
              d1ovd2 = sqrt(rho_prev[i]/rho_curr[i]);
              v1v2dz = v_curr[i][l]*v_prev[i][m]*dz;
              C_lm   = C_lm + (d1ovd2+ k1ovk2/d1ovd2)*v1v2dz/2.0; // eqs. 31, 32 DV's notes
          }
          //Rmx[l][m] = C_lm*H1[m]; 
          A_curr[l] = A_curr[l] + C_lm*H1[m]*A_prev[m]; // eq. 38 in DV notes
      }
  }
  
  
  delete[] rho_curr;
  delete[] H1;
  return 0;
}
*/




int getAcurr(int Nz_grid, int Nm, double dz, complex<double> *A_prev, double **v_prev, complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, complex<double> *k_curr, SampledProfile *atm_profile, complex<double> *A_curr)
{
  int i,l,m;
  double *rho_curr, v1v2dz, d1ovd2, sRj2ovRj1;
  complex<double> C_lm, k1ovk2, I_Rj1_Rj2, *H1;
  complex<double> I (0.0, 1.0);

  rho_curr = new double [Nz_grid];
  for (i=0; i<Nz_grid; i++) {
      rho_curr[i] = atm_profile->rho(i*dz/1000.0)*1000.0; // store density in rho_prev
  }

  sRj2ovRj1 = sqrt(Rj2/Rj1);
  I_Rj1_Rj2 = I*(Rj1-Rj2);
  
  // the trace of H1 (H1 is diagonal)
  H1 = new complex<double> [Nm];
  for (m=0; m<Nm; m++) {
      H1[m] = sRj2ovRj1*exp(I_Rj1_Rj2*k_prev[m]); 
      //printf("H1[%d]=%g %g\n", m, real(H1[m]), imag(H1[m]));
  }
  
  
  // compute A_curr = R1*A_prev; the matrix R1 = 1/2(C_tilde+ C_hat)*H1
  for (l=0; l<Nm; l++) {    
      A_curr[l] = 0.0; 
      for (m=0; m<Nm; m++) {
          k1ovk2 = k_prev[m]/k_curr[l];
          C_lm = 0.0; // 1/2(C_tilde + C_hat);
          for (i=0; i<Nz_grid; i++) { // compute integrals
              d1ovd2 = sqrt(rho_prev[i]/rho_curr[i]);
              v1v2dz = v_curr[i][l]*v_prev[i][m]*dz;
              C_lm   = C_lm + (d1ovd2+ k1ovk2/d1ovd2)*v1v2dz/2.0; // eqs. 31, 32 DV's notes
          }
          //Rmx[l][m] = C_lm*H1[m]; 
          A_curr[l] = A_curr[l] + C_lm*H1[m]*A_prev[m]; // eq. 38 in DV notes
      }
  }
  
  
  delete[] rho_curr;
  delete[] H1;
  return 0;
}

// lossless version
int getAcurr_ll(int Nz_grid, int Nm, double dz, complex<double> *A_prev_ll, double **v_prev, complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, complex<double> *k_curr, SampledProfile *atm_profile, complex<double> *A_curr_ll)
{
  int i,l,m;
  double *rho_curr, v1v2dz, d1ovd2, sRj2ovRj1;
  complex<double> C_lm, k1ovk2, I_Rj1_Rj2, *H1;
  complex<double> I (0.0, 1.0);

  rho_curr = new double [Nz_grid];
  for (i=0; i<Nz_grid; i++) {
      rho_curr[i] = atm_profile->rho(i*dz/1000.0)*1000.0; // store density in rho_prev
  }

  sRj2ovRj1 = sqrt(Rj2/Rj1);
  I_Rj1_Rj2 = I*(Rj1-Rj2);
  
  // the trace of H1 (H1 is diagonal)
  H1 = new complex<double> [Nm];
  for (m=0; m<Nm; m++) {
      H1[m] = sRj2ovRj1*exp(I_Rj1_Rj2*real(k_prev[m])); 
      //printf("H1[%d]=%g %g\n", m, real(H1[m]), imag(H1[m]));
  }
  
  // compute A_curr = R1*A_prev; the matrix R1 = 1/2(C_tilde+ C_hat)*H1
  for (l=0; l<Nm; l++) {    
      A_curr_ll[l] = 0.0; 
      for (m=0; m<Nm; m++) {
          k1ovk2 = real(k_prev[m])/real(k_curr[l]);
          C_lm = 0.0; // 1/2(C_tilde + C_hat);
          for (i=0; i<Nz_grid; i++) { // compute integrals
              d1ovd2 = sqrt(rho_prev[i]/rho_curr[i]);
              v1v2dz = v_curr[i][l]*v_prev[i][m]*dz;
              C_lm   = C_lm + (d1ovd2+ k1ovk2/d1ovd2)*v1v2dz/2.0; // eqs. 31, 32 DV's notes
          }
          //Rmx[l][m] = C_lm*H1[m]; 
          A_curr_ll[l] = A_curr_ll[l] + C_lm*H1[m]*A_prev_ll[m]; // eq. 38 in DV notes
      }
  }
  
  delete[] rho_curr;
  delete[] H1;
  return 0;
}


void parseReqRanges(std::string str, std::vector<double>& retVal) {
  // parse a string of the form "2_20_34_566_34.35" to extract the numbers into an array
  //cout << "in parseReqRanges\n";

  size_t  ix, it, i, p1; 

  // first pass to find out the number of underscores
  ix = -1;
  it =  0;
  ix = str.find("_");
  //cout << "ix= " << ix << endl;
  if (ix==string::npos) { // no underscore found
      vector<double> x (2,-1);
      x[1] = atof(str.c_str())*1000.0; // in meters
      x[0] = 0;
      retVal.swap(x);
      return; //NULL;
  }
  else { //  if (ix!=string::npos) i.e. underscore found
  it = 1;
      while (ix!=string::npos) {
          ix = str.find("_", ix+1);
          it++;
          //cout << "ix= " << ix << endl;
          //cout << "it= " << it << endl;
      }
  }
  //cout << "it= " << it << endl;

  // second pass to populate x
  vector<double> x (it,100);
  ix = str.find("_");
  //cout << "ix= " << ix << endl;
  if (ix==0) {
      throw invalid_argument("profile_ranges should not start with an underscore.");
      exit(1);
  }
      
  if (ix!=string::npos) {     
    //printf("%s\n", str.substr(0,ix).c_str());
     
      x[0] = atof(str.substr(0,ix).c_str())*1000.0;
      //printf("x[0]=%g\n", x[0]); 

      for (i=1; i<it; i++) {
          p1 = ix+1;
          ix = str.find("_", p1);
          x[i] = atof(str.substr(p1,ix-p1).c_str())*1000.0; // in meters
          //printf("x[%d]=%g\n", i, x[i]);  
      }
      //
  }
  retVal.swap(x);
  return;
}


int saveTLoss1D (int Nm, int n_zsrc, double sqrtrho_s, int it, vector<double> Rv, double rng_step, double **v_prev, complex<double> *k_prev, complex<double> *A_prev, complex<double> *A_prev_ll, const char *wa, double *prng) {

  int m;
  double rng;
  complex<double> PP, PP_ll;
  complex<double> I (0.0, 1.0);
  FILE *fp1dll, *fp1dloss;

  rng = (*prng); //rng_step; // current (marching) range is initialized here

  fp1dloss  = fopen("tloss_rd_1d.nm", wa);
  fp1dll    = fopen("tloss_rd_1d.lossless.nm", wa);
  while (rng<=Rv[it]) {
      PP    = 0.0 + I*0.0;
      PP_ll = 0.0 + I*0.0;
      if (it==1) {
          for (m=0; m<Nm; m++) {
              // pressure on the ground   - this version of PP doesn't work for big rng because NaN appear in the calculation
              //PP    = PP + 4.0*Pi*A_prev[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[0]/rng)*exp(I*k_prev[m]*(rng-Rv[0]));
              //PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[0]/rng)*exp(I*real(k_prev[m])*(rng-Rv[0]));
              
              // see eq. 11b in DV notes NMRD-OWCM: this should coincide with formula 5.1.14 in Ocean Acoustics page 274
              PP    = PP    + v_prev[0][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              PP_ll = PP_ll + v_prev[0][m]*v_prev[n_zsrc][m]*exp(I*real(k_prev[m])*rng)/sqrt(real(k_prev[m]));
          }  
          PP    = 4.0*Pi*exp(I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP;
          PP_ll = 4.0*Pi*exp(I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP_ll;
      }
      else {
          for (m=0; m<Nm; m++) {         
              // pressure on the ground (or TLoss if a 4*PI factor appears here) - eq 33 in DV's notes 
              PP    = PP + 4.0*Pi*A_prev[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1])); 
              PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[it-1]/rng)*exp(I*real(k_prev[m])*(rng-Rv[it-1])); 
          }
      }
      fprintf(fp1dloss, "%10.3f %18.8e  %18.8e\n", rng/1000.0, real(PP), imag(PP));
      fprintf(fp1dll  , "%10.3f %18.8e  %18.8e\n", rng/1000.0, real(PP_ll), imag(PP_ll));
      // update current range
      rng = rng + rng_step;   
  }
  fclose(fp1dloss);
  fclose(fp1dll);
  
  (*prng) = rng;
  return 0;
}


int saveTLoss2D (int Nm, int Nz_grid, int n_zsrc, double dz_km, int stepj, \
                  double sqrtrho_s, int it, vector<double> Rv, double rng_step, \
                  double **v_prev, complex<double> *k_prev, \
                  complex<double> *A_prev, complex<double> *A_prev_ll, \
                  const char *wa, double *prng) {

  int m,j;
  double rng, rng_km;
  complex<double> PP, eIpir;
  complex<double> I (0.0, 1.0);
  FILE *fp2d;

  rng    = (*prng); // re-initialize current (marching) range
  rng_km = rng/1000.0;
  
  // stepj controls the vertical sampling of 2D data saved
  // ensure it's never less than 1; necessary for the loop below
  if (stepj==0) {
      cout << "stepj = " << stepj << endl;
      throw invalid_argument( "stepj must be > 0:" );	
  }

  fp2d    = fopen("tloss_rd_2d.nm", wa);
  while (rng<=Rv[it]) {     
      if (it==1) { // region 1 only
          eIpir = 4.0*Pi*exp(I*Pi/4.0)/sqrt(8.0*Pi*rng);
          for (j=0; j<Nz_grid; j=j+stepj) { // note: stepj must be >0
              PP    = 0.0 + I*0.0;
              for (m=0; m<Nm; m++) {
                  // see eq. 11b in DV notes NMRD-OWCM: this should coincide with formula 5.1.14 in Ocean Acoustics page 274
                  PP    = PP + v_prev[j][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              }
              PP = eIpir*PP; 
              fprintf(fp2d,"%10.3f %10.3f %18.8e %18.8e\n", rng_km, j*dz_km, real(PP), imag(PP));     
          }
          fprintf(fp2d  , "\n");
      }
      else { // regions 2,3 ...
          for (j=0; j<Nz_grid; j=j+stepj) {
              PP    = 0.0 + I*0.0;
              for (m=0; m<Nm; m++) {
                  // see eq. 25 in DV notes NMRD-OWCM
                  PP = PP + 4.0*Pi*A_prev[m]*sqrtrho_s*v_prev[j][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1]));     
              }
              fprintf(fp2d,"%10.3f %10.3f %18.8e %18.8e\n", rng_km, j*dz_km, real(PP), imag(PP));       
          }
          fprintf(fp2d  , "\n"); 
      }
      rng = rng + rng_step; // update current range
      rng_km = rng/1000.0;
  }
  fclose(fp2d);
  
  (*prng) = rng; // returned value
  return 0;
}


int getNmodesNz(string fn, int *n_modes, int *nz, double *dz)
{
  int i; double freq;
  // read freq n_modes, nz, dz, k_pert, rho, v.
  FILE *f = fopen(fn.c_str(), "r");
  i = fscanf(f, "%lf %d %d %lf\n", &freq, n_modes, nz, dz);
  fclose(f); 
  return 0;
} 


int readEigenValVecs(string fn, complex<double> *k_pert, double *rho, double **v, int how_many)
{
  // reads from file all k_pert but only how_many eigenvectors
  //char mode_output[40];
  int i,j, n_modes, nz, jj, m;
  double tmp;
  double freq;
  double dz_km;
  double Re_k, Im_k;
  complex<double> I (0.0, 1.0);
  
  // read freq n_modes, nz, dz, k_pert, rho, v.
  
  FILE *f = fopen(fn.c_str(), "r");
  m = fscanf(f, "%lf %d %d %lf\n", &freq, &n_modes, &nz, &dz_km);
  if (how_many>n_modes) {
          std::ostringstream es;
      es << "Number of requested eigenvectors (" << how_many << ")"
         << " is greater than what is stored in the file: " << n_modes << endl;
      throw invalid_argument(es.str());
      return 0;
  }
  
  // read the wavenumbers
  for (j=0; j< how_many; j++) {
      m = fscanf(f, "%d %le %le\n", &jj, &Re_k, &Im_k);
      k_pert[j] = Re_k + I*Im_k;
  }
  
  // skip lines
  for (j=how_many; j< n_modes; j++) {
      m = fscanf(f, "%d %le %le\n", &jj, &Re_k, &Im_k);
  }  
  
  // read rho
  for (i=0; i<nz; i++) {
      m = fscanf(f,"%lf %le\n", &dz_km, &rho[i]);
  }

  // read the modes
  for (j=0; j<how_many; j++) {
      for (i=0; i<nz; i++) {
	        m = fscanf(f,"%lf %le\n", &dz_km, &tmp);
	        v[i][j] = tmp;
      } 
  }
  
  // zero the rest of elements from n_modes to Nmax
  //for (j=n_modes; j<Nmax; j++) {
  //   k_pert[j] = 0.0;
  //   for (i=0; i<nz; i++) {
	//        v[i][j] = 0.0;
  //    } 
  //}
  
  fclose(f);
  //printf(" Reading eigenvals and modes from file: %s (it has %d modes in total)\n", fn.c_str(), n_modes);  
  return 0;
}


// ----------------------------------------------

int getRRmats(string fn, double r1, double r2, int n_modes, int nz, double dz,  \
              complex<double> *k_curr, double *rho_curr, double **v_curr, \
              complex<double> *k_next, double *rho_next, double **v_next, \
              complex<double> **RRR)
{

  int i,l,m;
  double v1v2dz, d1ovd2, sr1ovr2, reknextlx2;
  complex<double> R1, R2, R3, R4, cknextl;
  complex<double> Ch_lm, Ct_lm, I_r2_r1, *H1, *H2;
  complex<double> I (0.0, 1.0);

  sr1ovr2 = sqrt(r1/r2);
  I_r2_r1 = I*(r2-r1);
  
  // the traces of H1, H2 (diagonal) at r=r2
  H1 = new complex<double> [n_modes];
  H2 = new complex<double> [n_modes];
  //printf("Warning: in getRRmats(): H1, H2 use only the real part of k\n");
  for (m=0; m<n_modes; m++) {
      //H1[m] = sr1ovr2*exp(I_r2_r1*real(k_curr[m])); 
      //H2[m] = sr1ovr2*exp(-I_r2_r1*real(k_curr[m]));

      H1[m] = sr1ovr2*exp(I_r2_r1*(real(k_curr[m])+I*imag(k_curr[m])));  // right-going
      H2[m] = sr1ovr2*exp(-I_r2_r1*(real(k_curr[m])-I*imag(k_curr[m]))); // still decaying for left going wave

      //H1[m] = sr1ovr2*exp(I_r2_r1*k_curr[m]); 
      //H2[m] = sr1ovr2*exp(-I_r2_r1*k_curr[m]);
      //printf("r1=%g; r2=%g; sr1ovr2 = %g; H1[%d]=%g %g\n", r1, r2, sr1ovr2, m, real(H1[m]), imag(H1[m]));
      //printf("k_curr[%d] = %g+I*%g\n", m, real(k_curr[m]), imag(k_curr[m]));
  }

 //complex<double> k1ovk2, Cp_lm, tmp1;
  // 
  for (l=0; l<n_modes; l++) {
      cknextl  = conj(k_next[l]);
      reknextlx2 = real(k_next[l])*2.0;  
      for (m=0; m<n_modes; m++) {
      
          //k1ovk2 = k_curr[m]/k_next[l];
          //Cp_lm = 0.0;
      
      
          Ct_lm = 0.0; // C_tilde;
          Ch_lm = 0.0; // C_hat;
          

          for (i=0; i<nz; i++) { // compute integrals
              d1ovd2 = sqrt(rho_curr[i]/rho_next[i]);
              v1v2dz = v_curr[i][m]*v_next[i][l]*dz;
              Ct_lm   = Ct_lm + d1ovd2*v1v2dz; // eqs. 87 DV's notes
              Ch_lm   = Ch_lm + v1v2dz/d1ovd2; // eqs. 88 DV's notes    
              
              
              
              //Cp_lm   = Cp_lm + (d1ovd2 + k1ovk2/d1ovd2)*v1v2dz/2.0; // eqs. 66 DV's notes
              //Cm_lm   = Cm_lm + (d1ovd2 - k1ovk2/d1ovd2)*v1v2dz/2.0; // eqs. 67 DV's notes          
          }
          
          //tmp1 = Cp_lm*H1[m]; // e.g. eq. 59 in DV notes: NMRD-CM pg 9
          
          
          
          
          

          R1 = (k_curr[m]*Ch_lm + cknextl*Ct_lm)*H1[m]/reknextlx2;   // eq. 83       
          R2 = (-conj(k_curr[m])*Ch_lm + cknextl*Ct_lm)*H2[m]/reknextlx2; // eq. 83      
          R3 = (-k_curr[m]*Ch_lm + k_next[l]*Ct_lm)*H1[m]/reknextlx2; // eq. 84  
          R4 = (conj(k_curr[m])*Ch_lm + k_next[l]*Ct_lm)*H2[m]/reknextlx2;  // eq. 84
          
          
          //if (R1!=tmp1) {
          //  printf("l=%d; m=%d; R1 = %g + I*%g\n", l, m, real(R1), imag(R1));
           // printf("l=%d; m=%d; tmp1 = %g + I*%g\n", l, m, real(tmp1), imag(tmp1));
          //}
          
          
          
          
          RRR[l][m]                 = R1;
          RRR[l][m+n_modes]         = R2;
          RRR[l+n_modes][m]         = R3;
          RRR[l+n_modes][m+n_modes] = R4;
      }
  } 
  
  //
  // save to file
  //
  FILE *f = fopen(fn.c_str(), "w");
  for (l=0; l<2*n_modes; l++) {    
      for (m=0; m<2*n_modes; m++) {
          fprintf(f, "%le %le\n", real(RRR[l][m]), imag(RRR[l][m]));
      }
  }
  fclose(f);     

  delete[] H1;
  delete[] H2;

  return 0;
}


// ----------------------------------------------
// lossless version
int getRRmats_ll(string fn, double r1, double r2, int n_modes, int nz, double dz,  \
              complex<double> *k_curr, double *rho_curr, double **v_curr, \
              complex<double> *k_next, double *rho_next, double **v_next, \
              complex<double> **RRR)
{

  int i,l,m;
  double v1v2dz, d1ovd2, sr1ovr2, reknextlx2;
  complex<double> R1, R2, R3, R4, cknextl;
  complex<double> Ch_lm, Ct_lm, I_r2_r1, *H1, *H2;
  complex<double> I (0.0, 1.0);

  sr1ovr2 = sqrt(r1/r2);
  I_r2_r1 = I*(r2-r1);
  
  // the traces of H1, H2 (diagonal) at r=r2
  H1 = new complex<double> [n_modes];
  H2 = new complex<double> [n_modes];
  //printf("Warning: in getRRmats(): H1, H2 use only the real part of k\n");
  for (m=0; m<n_modes; m++) {
      H1[m] = sr1ovr2*exp(I_r2_r1*(real(k_curr[m])));  // right-going
      H2[m] = sr1ovr2*exp(-I_r2_r1*(real(k_curr[m]))); // left going wave
  }

 //complex<double> k1ovk2, Cp_lm, tmp1;
  // 
  for (l=0; l<n_modes; l++) {
      cknextl  = conj(k_next[l]);
      reknextlx2 = real(k_next[l])*2.0;  
      for (m=0; m<n_modes; m++) {
          Ct_lm = 0.0; // C_tilde;
          Ch_lm = 0.0; // C_hat;

          for (i=0; i<nz; i++) { // compute integrals
              d1ovd2 = sqrt(rho_curr[i]/rho_next[i]);
              v1v2dz = v_curr[i][m]*v_next[i][l]*dz;
              Ct_lm   = Ct_lm + d1ovd2*v1v2dz; // eqs. 87 DV's notes
              Ch_lm   = Ch_lm + v1v2dz/d1ovd2; // eqs. 88 DV's notes           
          }
          R1 = (k_curr[m]*Ch_lm + cknextl*Ct_lm)*H1[m]/reknextlx2;   // eq. 83       
          R2 = (-conj(k_curr[m])*Ch_lm + cknextl*Ct_lm)*H2[m]/reknextlx2; // eq. 83      
          R3 = (-k_curr[m]*Ch_lm + k_next[l]*Ct_lm)*H1[m]/reknextlx2; // eq. 84  
          R4 = (conj(k_curr[m])*Ch_lm + k_next[l]*Ct_lm)*H2[m]/reknextlx2;  // eq. 84

          RRR[l][m]                 = R1;
          RRR[l][m+n_modes]         = R2;
          RRR[l+n_modes][m]         = R3;
          RRR[l+n_modes][m+n_modes] = R4;
      }
  } 
  
  //
  // save to file
  //
  FILE *f = fopen(fn.c_str(), "w");
  for (l=0; l<2*n_modes; l++) {    
      for (m=0; m<2*n_modes; m++) {
          fprintf(f, "%le %le\n", real(RRR[l][m]), imag(RRR[l][m]));
      }
  }
  fclose(f);     

  delete[] H1;
  delete[] H2;

  return 0;
}


int readRRmat(string fn, complex<double> **RRR, int n) {
  int l,m, tmp;
  double re, im;
  complex<double> I (0.0, 1.0);
  // read from file
  FILE *f = fopen(fn.c_str(), "r");
  for (l=0; l<n; l++) {    
      for (m=0; m<n; m++) {
          tmp = fscanf(f, "%lf %lf\n", &re, &im);
          RRR[l][m] = re + I*im;
      }
  }
  fclose(f);
  
  return 0;
}


int MatCMultiply(complex<double> **A, complex<double> **B, complex<double> **C, int M1, int N1, int M2, int N2)
{
  // matrix multiplication: A is a (M1 by N1), B is (M2 by N2) C will be (M1 by N2)
  if (N1!=M2) {
      std::ostringstream es;
      es << "Cannot do matrix multiplication: "
         << "Number of columns of A does not match number of lines of B" << endl;
      throw invalid_argument(es.str());
      return 0;
   }
  
  int i,j,k;
  complex<double> s;
  complex<double> I (0.0, 1.0);
  
  for (i=0; i<M1; i++) {
    for (j=0; j<N2; j++) {
      s = 0.0 + I*0.0;
      for (k=0; k<N1; k++) {
         s = s + A[i][k]*B[k][j];
      }
      C[i][j] = s;
    }
  }
  
  return 0;
}


int MatVecCMultiply(complex<double> **A, complex<double> *B, complex<double> *C, int M1, int N1, int M2)
{
  // mulptiply matrix by column vector: A is a (M1 by N1), B is (M2 by 1) C will be (M1 by 1)
  if (N1!=M2) {
      std::ostringstream es;
      es << "Cannot multiply A*B: "
         << "Number of columns of A does not match number of lines of B" << endl;
      throw invalid_argument(es.str());
      return 0;
   }
  
  int i,k;
  complex<double> s;
  complex<double> I (0.0, 1.0);
  
  for (i=0; i<M1; i++) {
      s = 0.0 + I*0.0;
      for (k=0; k<N1; k++) {
         s = s + A[i][k]*B[k];
      }
      C[i] = s;
  }
  return 0;
}


int SolveLinSys2(PetscScalar **AA, PetscScalar *bb, PetscScalar *yy, PetscInt n)
{ 
  // Solve: Ax=b; the output is the array yy (also contained in x - a Vec type container)
 
  Vec            x, b, u;          // approx solution, RHS, exact solution
  Mat            A;                // linear system matrix
  KSP            ksp;              // linear solver context
  PC             pc;               // preconditioner context
  //PetscReal      tol=1.e-14;  // norm of solution error
  PetscErrorCode ierr;
  PetscInt       i,col[n];
  PetscMPIInt    size;
  PetscScalar    neg_one = -1.0,one = 1.0;
  PetscBool      nonzeroguess = PETSC_FALSE;
 
  //static char help[] = "Solves a linear system with KSP.\n\n";
 
  //PetscInitialize(NULL,NULL,(char *)0,(char *)0);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  //printf("after MPI_Comm_size\n");
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
  //printf("after PetscOptionsGetInt\n");
  ierr = PetscOptionsGetBool(PETSC_NULL,"-nonzero_guess",&nonzeroguess,PETSC_NULL);CHKERRQ(ierr);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //       Compute the matrix and right-hand-side vector that define
  //       the linear system, Ax = b.
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // 
  //   Create vectors.  Note that we form 1 vector from scratch and
  //   then duplicate as needed.
  //
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&u);CHKERRQ(ierr);
  

  // 
  //   Create matrix.  When using MatCreate(), the matrix format can
  //   be specified at runtime.
  //
  //   Performance tuning note:  For problems of substantial size,
  //   preallocation of matrix memory is crucial for attaining good 
  //   performance. See the matrix chapter of the users manual for details.
  //
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  // 
  //   Assemble matrix A
  //
  //printf("assembling matrix A in SolveLinSys2()\n");
  for (i=0; i<n; i++) {
    col[i] = i;
  }
  
  // set values row by row    
  for (i=0; i<n; i++) {
    ierr = MatSetValues(A,1,&i,n,col,AA[i],INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
  
  //PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO);
  //PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);
  //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

  VecCreateSeqWithArray(MPI_COMM_SELF, 1,n,bb,&b);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  //VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // 
  //   Set exact solution; then compute right-hand-side vector.
  //
  ierr = VecSet(u,one);CHKERRQ(ierr);
  //ierr = MatMult(A,u,b);CHKERRQ(ierr);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //              Create the linear solver and set various options
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // 
  //   Create linear solver context
  //
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  // 
  //   Set operators. Here the matrix that defines the linear system
  //   also serves as the preconditioning matrix.
  //
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  // 
  //   Set linear solver defaults for this problem (optional).
  //   - By extracting the KSP and PC contexts from the KSP context,
  //     we can then directly call any KSP and PC routines to set
  //     various options.
  //   - The following four statements are optional; all of these
  //     parameters could alternatively be specified at runtime via
  //     KSPSetFromOptions();
  //
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  // 
  //  Set runtime options, e.g.,
  //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  //  These options will override those specified above as long as
  //  KSPSetFromOptions() is called _after_ any other customization
  //  routines.
  //
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  if (nonzeroguess) {
    PetscScalar p = .5;
    ierr = VecSet(x,p);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }
 
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //                    Solve the linear system
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // 
  //   Solve linear system
  //
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr); 
  // VecView(x, PETSC_VIEWER_STDOUT_WORLD); // view solution (as type Vec)
  
  // copy array from Vec x to yy  - the output of this routine
  PetscScalar *y;
  VecGetArray(x,&y);
  for (i=0; i<n; i++) {
    yy[i] = y[i];
    //PetscPrintf(PETSC_COMM_WORLD, "avec[%d] = %g\n",i,avec[i]);
  }

  //
  //   View solver info; we could instead use the option -ksp_view to
  //   print this info to the screen at the conclusion of KSPSolve().
  //
  //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //                    Check solution and clean up
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // 
  //   Check the error
  //
  //ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  //ierr  = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  //ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  //if (norm > tol){
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %G, Iterations %D\n",
  //                   norm,its);CHKERRQ(ierr);
  //}

  // 
  //   Free work space.  All PETSc objects should be destroyed when they
  //   are no longer needed.
  //
  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  //
  //   Always call PetscFinalize() before exiting a program.  This routine
  //     - finalizes the PETSc libraries as well as MPI
  //     - provides summary and diagnostic information if certain runtime
  //       options are chosen (e.g., -log_summary).
  //
  //ierr = PetscFinalize();
  return 0;
}


int save_printCMatrix(complex<double> **A, int nr, int nc, std::string filename, bool flag) {

  int i,j;
  if (flag) { // print to screen 
      for (i=0; i<nr; i++) {    
          for (j=0; j<nc; j++) {
              printf("(%e %e)  ", real(A[i][j]), imag(A[i][j]));
          }
          printf("\n");
      }
  }
  else { // save to file
      FILE *f = fopen(filename.c_str(), "w");
      for (i=0; i<nr; i++) {    
          for (j=0; j<nc; j++) {
              fprintf(f, "%18.15e %18.15e\n", real(A[i][j]), imag(A[i][j]));
          }
          //fprintf(f,"\n");
      }
      fclose(f);
  }
      
  return 0;      
}

int save_printCVector(complex<double> *A, int n, std::string filename, bool flag) {

  int i;
  if (flag) { // print to screen 
      for (i=0; i<n; i++) {    
              printf("%e %e\n", real(A[i]), imag(A[i]));
      }
  }
  else { // save to file
      FILE *f = fopen(filename.c_str(), "w");
      for (i=0; i<n; i++) {    
              fprintf(f,"%18.15e %18.15e\n", real(A[i]), imag(A[i]));
      }
      fclose(f);
  }
      
  return 0;      
}


int computePressure_1D( int it, int Nm, int n_z, vector<double> Rv, double rng_step, double sqrtrho_z, complex<double> *k_curr,  double **v_curr, complex<double> *ab_curr, complex<double> *ab_curr_ll, const char *wa, double *prng_curr)
{
  int m;
  double rng, rj_1;
  complex<double> H1, H2, PP, PP_ll;
  complex<double> I (0.0, 1.0);
  FILE *fp1dloss, *fp1dll;

  rng = (*prng_curr); // current (marching) range is initialized here      
  rj_1 = Rv[it-1]; // 'it' is the Region we are in
  //printf("rng = %g; Rv[it-1] = %g\n", rng, Rv[it-1]);
  
  fp1dloss  = fopen("tloss_rd2wcm_1d.nm", wa);
  fp1dll    = fopen("tloss_rd2wcm_1d.lossless.nm", wa);
  
  //printf("!!! Warning: in computePressure_1D(): H1, H2 use only the real part of k\n");

  while (rng<=Rv[it]) { 
    PP    = 0.0 + I*0.0;
    PP_ll = 0.0 + I*0.0;
    for (m=0; m<Nm; m++) {
        H1= sqrt(rj_1/rng)*exp(I*k_curr[m]*(rng - rj_1));
        H2= sqrt(rj_1/rng)*exp(-I*k_curr[m]*(rng - rj_1));
        
        // pressure (or TLoss if a 4*PI factor appears here) - eq 33 in DV's notes 
        PP = PP + (ab_curr[m]*H1 + ab_curr[m+Nm]*H2)*v_curr[n_z][m];

        // lossless case
        H1= sqrt(rj_1/rng)*exp(I*real(k_curr[m])*(rng - rj_1));
        H2= sqrt(rj_1/rng)*exp(-I*real(k_curr[m])*(rng - rj_1));         
        PP_ll = PP_ll + (ab_curr_ll[m]*H1 + ab_curr_ll[m+Nm]*H2)*v_curr[n_z][m];
        
        //// use the next 2 lines only if the modes need scaling by sqrt(rho)
        //PP = PP + (ab_curr[m]*H1 + ab_curr[m+Nm]*H2)*sqrtrho_z*v_curr[n_z][m];
        //PP_ll = PP_ll + (ab_curr_ll[m]*H1 + ab_curr_ll[m+Nm]*H2)*sqrtrho_z*v_curr[n_z][m];
    } 
    fprintf(fp1dloss, "%10.3f %18.8e  %18.8e\n", rng/1000.0, real(4.0*Pi*PP), imag(4.0*Pi*PP));
    fprintf(fp1dll  , "%10.3f %18.8e  %18.8e\n", rng/1000.0, real(4.0*Pi*PP_ll), imag(4.0*Pi*PP_ll));
    // update current range
    rng = rng + rng_step;   
  }
  fclose(fp1dloss);
  fclose(fp1dll);
   (*prng_curr) = rng; // this will be returned to the caller
  return 0; 
}


int computePressure_2D( int it, int Nm, int Nz, int stepn, vector<double> Rv, double rng_step, double dz, double *rho_curr, complex<double> *k_curr,  double **v_curr, complex<double> *ab_curr, complex<double> *ab_curr_ll, const char *wa, double *prng_curr)
{
  int m, n;
  double rng, rj_1, dz_km;
  double *sqrtrho;
  complex<double> PP, PP_ll;
  complex<double> *H1, *H2, *H1_ll, *H2_ll;
  complex<double> I (0.0, 1.0);
  FILE *fp2dloss, *fp2dll;
  
  sqrtrho = new double [Nm];
  H1      = new complex<double> [Nm];
  H2      = new complex<double> [Nm];
  H1_ll   = new complex<double> [Nm];
  H2_ll   = new complex<double> [Nm];
  
  // stepn controls the vertical sampling of 2D data saved
  // ensure it's never less than 1; necessary for the loop below
  if (stepn==0) {
      cout << "stepn = " << stepn << endl;
      throw invalid_argument( "stepn must be > 0:" );	
  }      
  
  dz_km = dz/1000.0;
  rng = (*prng_curr); // current (marching) range is initialized here
  rj_1 = Rv[it-1]; // 'it' is the Region we are in
  //printf("rng = %g; Rv[it-1] = %g\n", rng, Rv[it-1]);  

  fp2dloss  = fopen("tloss_rd2wcm_2d.nm", wa);
  fp2dll    = fopen("tloss_rd2wcm_2d.lossless.nm", wa);
  //printf("!!! Warning: in computePressure_2D(): H1, H2 use only the real part of k\n");

  for (m=0; m<Nm; m++) { //pre-compute sqrt(rho)
      sqrtrho[m] = sqrt(rho_curr[m]);
  }

  while (rng<=Rv[it]) {
    for (m=0; m<Nm; m++) {
      H1[m]   = sqrt(rj_1/rng)*exp(I*k_curr[m]*(rng - rj_1));
      H2[m]   = sqrt(rj_1/rng)*exp(-I*k_curr[m]*(rng - rj_1));
      H1_ll[m]= sqrt(rj_1/rng)*exp(I*real(k_curr[m])*(rng - rj_1));
      H2_ll[m]= sqrt(rj_1/rng)*exp(-I*real(k_curr[m])*(rng - rj_1)); 
    }
  
    for (n=0; n<Nz; n = n+stepn) {
        PP    = 0.0 + I*0.0;
        PP_ll = 0.0 + I*0.0;
        for (m=0; m<Nm; m++) {
            // pressure (or TLoss if a 4*PI factor appears here) - eq 33 in DV's notes 
            PP    = PP + (ab_curr[m]*H1[m] + ab_curr[m+Nm]*H2[m])*v_curr[n][m];
            PP_ll = PP_ll + (ab_curr_ll[m]*H1_ll[m] + ab_curr_ll[m+Nm]*H2_ll[m])*v_curr[n][m];
            
            //// use the next 2 lines only if the modes need scaling by sqrt(rho)
            //PP    = PP + (ab_curr[m]*H1[m] + ab_curr[m+Nm]*H2[m])*sqrtrho[m]*v_curr[n][m];
            //PP_ll = PP_ll + (ab_curr_ll[m]*H1_ll[m] + ab_curr_ll[m+Nm]*H2_ll[m])*sqrtrho[m]*v_curr[n][m];
        } 
        fprintf(fp2dloss, "%10.3f %10.3f %18.8e  %18.8e\n", rng/1000.0, n*dz_km, real(4.0*Pi*PP), imag(4.0*Pi*PP));
        fprintf(fp2dll  , "%10.3f %10.3f %18.8e  %18.8e\n", rng/1000.0, n*dz_km, real(4.0*Pi*PP_ll), imag(4.0*Pi*PP_ll));  
    }
    // update current range
    rng = rng + rng_step;
    fprintf(fp2dloss, "\n");
    fprintf(fp2dll, "\n");
  }
  fclose(fp2dloss);
  fclose(fp2dll);
    
  delete [] sqrtrho;  
  delete [] H1;
  delete [] H2;
  delete [] H1_ll;
  delete [] H2_ll;
  (*prng_curr) = rng; // this will be returned to the caller
  return 0; 
}


int makeYYYY_MM_DD_subdir(std::string *subdir) {
  std::time_t tm1 = std::time(NULL);
  char sdir[25];
  if (std::strftime(sdir, 25, "%F_%H_%M_%S", std::localtime(&tm1))) {
	    //cout << sdir << endl;
  }
  else {
      throw invalid_argument("time string YYYY_MM_DD_HH_MM_SS needed for dir name could not be formed.");
  }
  (*subdir)  = string(sdir);             // subdirectory name
  // make a dir (i.e. subdirectory subdir in this case)
  // if the wrong directory permissions are set for this dir look up
  // man 2 umask and man 2 stat
  if (mkdir((*subdir).c_str(), S_IRWXU | S_IRWXG | S_IRWXO)==-1) {
      std::ostringstream es;
      es << "attempt to make dir << " << (*subdir) << " failed.";
      throw invalid_argument(es.str());
  }
  return 0;
}


int getRegionBoundaries(bool flg, double maxrange, double req_profile_step, string prf_ranges_km, int *Nprofiles, vector<double> *R) {

// this function populates vector R with the Region boundaries.
// if in the main() the following option is used e.g.:
//  --use_profile_ranges_km 30_100_500_700.5 
// then Rv will be [30, 30, 100, 500, 700.5]. Note that the first entry is doubled - see algorithm
// otherwise
// Rv is constructed from req_profile_step: [step:step:maxrange]

  int i;
  vector<double> Rv;
  
  if (flg) { // if option --use_profile_ranges_km such as 30_100_500_700.5 is provided
      //oNB->getProfileRanges(prf_ranges_km);
      //printf("prf_ranges_km=%s\n", prf_ranges_km.c_str());
      parseReqRanges(prf_ranges_km, Rv);
      
      if (Rv[0]!=0.0) { // if the first element is not zero insert 0 at the beginning
          vector<double>::iterator it;
          it = Rv.begin();
          Rv.insert(it, Rv[0]);
      }
      
      if (Rv[Rv.size()-1]>=maxrange) {//if Rv(end) is >= maxrange force the last elem. =maxrange
          Rv[Rv.size()-1] = maxrange;
      }
      else { //add maxrange at end
          Rv.push_back(maxrange);
      }
      (*Nprofiles) = Rv.size()-1;
      //for (i=0; i<Rv.size(); i++) cout << " " << Rv[i] << endl;
  }
  else {
      (*Nprofiles) = (int) ceil(maxrange/req_profile_step); // at least one profile is ensured
      vector<double> Rv1 ((*Nprofiles)+1, 0.0);
      for (i=1; i<(int) Rv1.size(); i++) {
          Rv1[i] = i*req_profile_step;  
      } 
      Rv1[Rv1.size()-1] = maxrange; // last element adjusted to maxrange
      //Rv1[0] = Rv1[1]; // the convention is the first 2 elements are equal - see algorithm
      Rv.swap(Rv1);
  }
  Rv[0] = Rv[1]; // by convention the first 2 elements are equal - see algorithm
  //cout << "Nprofiles = " << Nprofiles << endl;
  //cout << "Rv size   = " << Rv.size() << endl;
  //for (i=0; i<Rv.size(); i++) {
  //  printf("Rv[%i] = %g\n", i, Rv[i]);
  //}
  
  (*R).swap(Rv);
  return 0;     
}


int getCeffMinMax(NCPA::SampledProfile *p, double *ceffmin, double *ceffmax) {
  // obtain ceffmin, ceffmax of a profile

  int    i, nz;
  double azi;
  double *u, *v, *rho, *P, *ceff;
  nz  = p->nz();
  azi = p->getPropagationAzimuth()*Pi/180.0;

  double gamma = 1.4;

  u    = new double [nz];
  v    = new double [nz];
  rho  = new double [nz];
  P    = new double [nz];
  ceff = new double [nz];
  
  p->get_u(u, nz);
  p->get_v(v, nz);
  p->get_rho(rho, nz);
  p->get_p(P, nz);
  
  ceff[0] = sqrt(gamma*P[0]*100.0/(rho[0]*1000.0)) + (u[0]*sin(azi)+v[0]*cos(azi));
  (*ceffmax) = ceff[0];
  (*ceffmin) = ceff[0];

  for (i=1; i<nz; i++) {
      ceff[i] = sqrt(gamma*P[i]*100.0/(rho[i]*1000.0)) + (u[i]*sin(azi)+v[i]*cos(azi));
      if (ceff[i]>(*ceffmax)) (*ceffmax) = ceff[i];
      if (ceff[i]<(*ceffmin)) (*ceffmin) = ceff[i];
  }
  
  delete [] u;
  delete [] v;
  delete [] rho;
  delete [] P;
  delete [] ceff;
  return 0;
}



int getGlobalceffMinMax(int Nprofiles, vector<double> Rv, double azi, string atm_profile_dir, string atmosfile, string atmosfileorder, int skiplines, double *ceffMin, double *ceffMax) {

  int j;
  double cmin, cmax;
  SampledProfile *atm_profile; 
  
  // the case j=1
  j = 1;
  //printf("\nRegion %d (%g to %g km)\n", j, 0.0, Rv[j]/1000.0);
  if (!atm_profile_dir.empty()) { // get ascii 1D profiles from files in directory 'atm_profile_dir'
          atm_profile = get_RngDepndProfiles_ascii(j, atmosfileorder, skiplines, atm_profile_dir, "profile", 0);
  }
  else { // get profiles from the .env file
      atm_profile = get_RngDepnd_profile(atmosfile, 0.0); // if the very first profile is required in the first region
  }
  atm_profile->setPropagationAzimuth(azi);
  getCeffMinMax(atm_profile, &cmin, &cmax);
  (*ceffMin) = cmin;
  (*ceffMax) = cmax;
  
  delete atm_profile;

  // the rest of profiles j>1
  for (j=2; j<=Nprofiles; j++) {
      //printf("\nRegion %d (%g to %g km)\n", j, Rv[j-1]/1000.0, Rv[j]/1000.0);
      if (!atm_profile_dir.empty()) { // get ascii 1D profiles from files in directory 'atm_profile_dir'
          atm_profile = get_RngDepndProfiles_ascii(j, atmosfileorder, skiplines, atm_profile_dir, "profile", 0);
      }
      else { // get profiles from the .env file
          atm_profile = get_RngDepnd_profile(atmosfile, Rv[j-1]);
      }
      atm_profile->setPropagationAzimuth(azi);
      getCeffMinMax(atm_profile, &cmin, &cmax);
      if (cmin<(*ceffMin)) (*ceffMin) = cmin;
      if (cmax>(*ceffMax)) (*ceffMax) = cmax;

      delete atm_profile;
  }

  return 0;     
}

