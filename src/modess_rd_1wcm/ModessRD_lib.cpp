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
#include "ModessRD_lib.h"

#include "binaryreader.h"

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


/*
int plotwGNUplot(double freq, bool write_2D_TLoss) {
    if (1) {   
        // open a pipe to gnuplot and plot the results
        FILE *pipe = popen("gnuplot -persist","w");
        printf("Plotting with gnuplot...\n");
        //fprintf(pipe, "set data style lines\n");
        //fprintf(pipe, "set yrange [%f:%f]\n",min(alt,nAlt),max(alt,nAlt));       
       
        fprintf(pipe, "set title 'One-way Coupled Modes - Transmission loss; Frequency %g Hz'\n", freq);
        fprintf(pipe, "set xlabel 'Range [km]'\n");
        fprintf(pipe, "set ylabel 'TL [dB]'\n");
        fprintf(pipe, "set grid; show grid;\n");
        fprintf(pipe, "plot './tloss_rd_1d.lossless.nm' using 1:(10*log10($2**2 + $3**2)) with lines lt 3 title 'lossless' , \\\n");
        fprintf(pipe, "    './tloss_rd_1d.nm' using 1:(10*log10($2**2 + $3**2)) with lines lt 1 title 'lossy'\n");        
        pclose(pipe);
        
        if (write_2D_TLoss) {
            //To plot a surface plot:
            pipe = popen("gnuplot -persist","w");
            //fprintf(pipe,"set term wxt 1\n");
            fprintf(pipe, "set pm3d map\n");
            fprintf(pipe, "set cbrange [-140:-100]\n");
            fprintf(pipe, "set xlabel 'Range [km]'\n");
            fprintf(pipe, "set ylabel 'Height [km]'\n");
            fprintf(pipe, "set title 'One-way Coupled Modes - Transmission loss; Frequency %g Hz'\n", freq);            
            fprintf(pipe,"splot './tloss_rd_2d.nm' using 1:2:(20*log10(sqrt($3**2 + $4**2))) title ''\n");
            pclose(pipe);
        }
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

  // convert winds to units of km/s  as required in SampledProfile; this could change eventually!!!
  double *uu, *vv, *ww;

  uu = new double [nAlt];
  vv = new double [nAlt];
  ww = new double [nAlt];
  double az = (baz+180.0)*Pi/180.0; // forward azimuth in radians
  
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
  // printf("in get_RngDepnd_profile(): saving the along and cross winds in units of km/sec\n");
  for (i=0; i<nAlt; i++) {
      uu[i] = wind_along_mps[i]; ///1000.0;
      vv[i] = wind_cross_mps[i]; ///1000.0;
      ww[i] = wind_vert_mps[i]; ///1000.0;
      alts_km[i] = alts_km[i]-alts_km[0]; // adjust altitudes to start from zero (no terrain considered yet)
  }

  // assemble the Sampled profile and set its azimuth
  SampledProfile *atm_prof;
  double z0 = 0.0; // for flat ground use z0=0
  atm_prof = new SampledProfile(nAlt, alts_km, temp_K, uu, vv, ww, density_gpcm3, pressure_hPa, z0);
  atm_prof->setPropagationAzimuth(az*180.0/Pi);

  if (1) {
      saveSampledProfile("./profile.int", atm_prof);
      //printf("in get_RngDepnd_profile(): saving the interpolated profile\n");
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
              int N, string atmosfileorder, int skiplines, \
              string dirname, string pattern) {
              
// get the Nth profile stored in a file named "profile000M" residing in a specified
// directory. If N is greater than the number of files found then the last file
// is used repeatedly as necessary             

  list<string> files;
  list<string>::iterator it;
  //string pattern ("prof");
  
  // get and sort the files (they should have names such as profile0001.dat, etc)
  getFile_list(dirname, files, pattern);
  files.sort();
  if (1) { //print the sorted file list
      cout << "Sorted file list from directory: " << dirname << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }
  // cout << "filessize = " << files.size() << endl;
  
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
  
  //printf("Nz = %d\n", p->nz());
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
  
  //for (i=0; i<Nz; i++) {
  //  printf("zz=%g\n", zz[i]);
  //}
  
  FILE *fp;
  
  fp = fopen(filename.c_str(), "w");
  for (i=0; i<Nz; i++) {
    fprintf(fp, "%8.3f %12.3e %12.3e %12.3e %8.3f  %8.4e  %8.4e\n", zz[i], uu[i], vv[i], ww[i], tt[i], dd[i], pp[i]);
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


int computeTLoss1D(int Nmodes, double r, double RR, int n_zsrc, complex<double> *k_pert, double **v_s, int *signv2, complex<double> *Kintg_atR)
{
  int m;
  double modal_sum_i, modal_sum_i_ll;
  complex<double> modal_sum_c, modal_sum_c_ll;
  complex<double> kint, expov8pi;
  complex<double> I (0.0, 1.0);
  
  expov8pi = exp(I*Pi*0.25)*sqrt(1./8./Pi); // I*exp(-I*Pi/4) = +exp(I*Pi/4)

  FILE *tloss_1d    = fopen("tloss_rd_1d.nm","a");
  FILE *tloss_ll_1d = fopen("tloss_rd_1d.lossless.nm","a");

  //for (i=0; i<=n_r; i++) {
  //    r = r1 + i*dr;
      modal_sum_c    = 0.;  // coherent sum
      modal_sum_i    = 0.;  // incoherent sum
      modal_sum_c_ll = 0;
      modal_sum_i_ll = 0;
      
      for (m=0; m<Nmodes; m++) {
          kint = Kintg_atR[m] + k_pert[m]*(r-RR);
          modal_sum_c    = modal_sum_c + signv2[m]*v_s[n_zsrc][m]*v_s[0][m]*exp(I*kint)/sqrt(k_pert[m]);
          modal_sum_i    = modal_sum_i + signv2[m]*pow(v_s[n_zsrc][m]*v_s[0][m],2)*exp(-2*imag(kint))/abs(k_pert[m]);
          modal_sum_c_ll = modal_sum_c_ll + signv2[m]*v_s[n_zsrc][m]*v_s[0][m]*exp(I*real(kint))/sqrt(real(k_pert[m]));
          modal_sum_i_ll = modal_sum_i_ll + signv2[m]*pow(v_s[n_zsrc][m]*v_s[0][m],2)/real(k_pert[m]);
      }
      
      modal_sum_c    = 4*Pi*modal_sum_c*expov8pi/sqrt(r); // I*exp(-I*Pi/4) = +exp(I*Pi/4)
      modal_sum_i    = 4*Pi*sqrt(modal_sum_i)*sqrt(1./8./Pi/r);
      modal_sum_c_ll = 4*Pi*modal_sum_c_ll*expov8pi/sqrt(r);
      modal_sum_i_ll = 4*Pi*sqrt(modal_sum_i_ll)*sqrt(1./8./Pi/r);
        
     /* 
     for (m=0; m<Nmodes; m++) {
          kint = Kintg_atR[m] + k_pert[m]*(r-RR);
          modal_sum_c    = modal_sum_c + v_s[n_zsrc][m]*v_s[0][m];
          modal_sum_i    = modal_sum_i + pow(v_s[n_zsrc][m]*v_s[0][m],2)*exp(-2*imag(kint))/abs(k_pert[m]);
          modal_sum_c_ll = modal_sum_c_ll + v_s[n_zsrc][m]*v_s[0][m];
          modal_sum_i_ll = modal_sum_i_ll + signv2[m]*pow(v_s[n_zsrc][m]*v_s[0][m],2)/real(k_pert[m]);
      }
    */
 
      fprintf(tloss_1d,"%f %20.12e %20.12e %20.12e\n", r/1000, real(modal_sum_c), imag(modal_sum_c), modal_sum_i);
      fprintf(tloss_ll_1d,"%f %20.12e %20.12e %20.12e\n", r/1000, real(modal_sum_c_ll), imag(modal_sum_c_ll), modal_sum_i_ll);
  //}
  fclose(tloss_1d);
  fclose(tloss_ll_1d);
  printf("           appended to file tloss_rd_1d.nm\n");
  printf("           appended to file tloss_rd_1d.lossless.nm\n");
  return 0;
}


int writeDispersion(int select_modes, int n_zsrc, double freq, complex<double> *k_pert, double **v_s, string file_stub)
{
  int i;
  //int n_zsrc = ceil(z_src/dz)+1;
  char dispersion_file[256];

  sprintf(dispersion_file,"%s_%e.nm",file_stub.c_str(), freq);
  FILE *dispersion = fopen(dispersion_file,"w");
  fprintf(dispersion,"%e   %d",freq,select_modes);
  for(i=0;i<select_modes;i++) {
      fprintf(dispersion,"   %.12e   %.12e",real(k_pert[i]),imag(k_pert[i]));
      //fprintf(dispersion,"   %.12e   %.12e",real(v_s[n_zsrc][i]),real(v_s[0][i]));
      fprintf(dispersion,"   %.12e   %.12e",(v_s[n_zsrc][i]),(v_s[0][i]));
  }
  fprintf(dispersion,"\n");
  fclose(dispersion);
  printf("           file %s created\n", dispersion_file);
  return 0;
}


int writeEigenVec(int nz, int select_modes, double dz, double **v_s, string file_stub)
{
  char mode_output[256];
  int i,j;

  for (j=0; j<select_modes; j++) {
      sprintf(mode_output,"%s_mode_%d.dat", file_stub.c_str(),j);
      FILE *eigenfunction= fopen(mode_output, "w");
      double chk = 0;
      for (i=0; i<nz; i++) {
	        fprintf(eigenfunction,"%f %15.8e\n", (i+1)*dz, v_s[i][j]);
	        chk = chk + v_s[i][j]*v_s[i][j]*dz;
      }
      if (fabs(1.-chk) > 0.1) { printf("Check if eigenfunction %d is normalized!\n", j); }
      fclose(eigenfunction);
  }
  printf("           files mode_<mode_number> created (%d in total)\n", select_modes);  
  return 0;
} 


int getAcurr(int Nz_grid, int Nm_prev, int Nm, double dz, complex<double> *A_prev, double **v_prev, complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, complex<double> *k_curr, SampledProfile *atm_profile, complex<double> *A_curr) {

  int i,l,m;
  double *rho_curr, v1v2dz, d1ovd2, sRj2ovRj1;
  complex<double> C_lm, k1ovk2, I_Rj1_Rj2, *H1;
  complex<double> I (0.0, 1.0);

  rho_curr = new double [Nz_grid];
  for (i=0; i<Nz_grid; i++) {
      rho_curr[i] = atm_profile->rho(i*dz/1000.0)*1000.0;
  }

  sRj2ovRj1 = sqrt(Rj2/Rj1);
  I_Rj1_Rj2 = I*(Rj1-Rj2);
  
  // the trace of H1 (H1 is diagonal)
  H1 = new complex<double> [Nm_prev];
  for (m=0; m<Nm_prev; m++) {
      H1[m] = sRj2ovRj1*exp(I_Rj1_Rj2*k_prev[m]); 
      //printf("H1[%d]=%g %g\n", m, real(H1[m]), imag(H1[m]));
  }

  // compute A_curr = R1*A_prev; the matrix R1 = 1/2(C_tilde+ C_hat)*H1
  for (l=0; l<Nm; l++) {
      A_curr[l] = 0.0; 
      for (m=0; m<Nm_prev; m++) {
        
          k1ovk2 = k_prev[m]/k_curr[l];
          C_lm = 0.0; // 1/2(C_tilde + C_hat);
          for (i=0; i<Nz_grid; i++) { // compute integrals
              d1ovd2 = sqrt(rho_prev[i]/rho_curr[i]);
              v1v2dz = v_curr[i][l]*v_prev[i][m]*dz;
              C_lm   = C_lm + (d1ovd2 + k1ovk2/d1ovd2)*v1v2dz/2.0; // eqs. 31, 32 DV's NMRD-OWCM notes
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
int getAcurr_ll(int Nz_grid, int Nm_prev, int Nm, double dz, complex<double> *A_prev_ll, double **v_prev, complex<double> *k_prev, double *rho_prev, double Rj1, double Rj2, double **v_curr, complex<double> *k_curr, SampledProfile *atm_profile, complex<double> *A_curr_ll)
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
  H1 = new complex<double> [Nm_prev];
  for (m=0; m<Nm_prev; m++) {
      H1[m] = sRj2ovRj1*exp(I_Rj1_Rj2*real(k_prev[m])); 
      //printf("H1[%d]=%g %g\n", m, real(H1[m]), imag(H1[m]));
  }
  
  // compute A_curr = R1*A_prev; the matrix R1 = 1/2(C_tilde+ C_hat)*H1
  for (l=0; l<Nm; l++) {    
      A_curr_ll[l] = 0.0; 
      for (m=0; m<Nm_prev; m++) {
          k1ovk2 = real(k_prev[m])/real(k_curr[l]);
          C_lm = 0.0; // 1/2(C_tilde + C_hat);
          for (i=0; i<Nz_grid; i++) { // compute integrals
              d1ovd2 = sqrt(rho_prev[i]/rho_curr[i]);
              v1v2dz = v_curr[i][l]*v_prev[i][m]*dz;
              C_lm   = C_lm + (d1ovd2+ k1ovk2/d1ovd2)*v1v2dz/2.0; // eqs. 31, 32 DV's NMRD-OWCM notes
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
  if (ix==string::npos) { // no underscore found but assumes atof() returns a valid double
      //vector<double> x (1,-1);
      //retVal.swap(x);      
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


/*
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
*/


/*
int saveTLoss1D (int Nm, int n_zsrc, double sqrtrho_s, int n_zrcv, double sqrtrho_r, int it, vector<double> Rv, double rng_step, double **v_prev, complex<double> *k_prev, complex<double> *A_prev, complex<double> *A_prev_ll, const char *wa, double *prng) {

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
      if (it==1) { // for region 1 only
          for (m=0; m<Nm; m++) {
              // pressure on the ground   - this version of PP doesn't work for big rng because NaN appear in the calculation
              //PP    = PP + 4.0*Pi*A_prev[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[0]/rng)*exp(I*k_prev[m]*(rng-Rv[0]));
              //PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[0]/rng)*exp(I*real(k_prev[m])*(rng-Rv[0]));
              
              // see eq. 11b in DV notes NMRD-OWCM: this should coincide with formula 5.1.14 in Ocean Acoustics page 274
              
              //// use the next 2 lines if modes are scaled by 1/sqrt(rho)
              //PP    = PP    + sqrtrho_r/sqrtrho_s*v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              //PP_ll = PP_ll + sqrtrho_r/sqrtrho_s*v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*real(k_prev[m])*rng)/sqrt(real(k_prev[m]));
              // just the modal sum - to correspond with the results of range-independent Modess             
              PP    = PP    + v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              PP_ll = PP_ll + v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*real(k_prev[m])*rng)/sqrt(real(k_prev[m]));              
              
          }  
          PP    = 4.0*Pi*I*exp(-I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP;
          PP_ll = 4.0*Pi*I*exp(-I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP_ll;
      }
      else { // for regions 2,3,...
          for (m=0; m<Nm; m++) {         
              // pressure on the ground (or TLoss if a 4*PI factor appears here) - eq 33 in DV's notes
              //// use the next 2 lines if modes are scaled by 1/sqrt(rho) 
              //PP    = PP + 4.0*Pi*A_prev[m]*sqrtrho_r*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1])); 
              //PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*sqrtrho_r*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*real(k_prev[m])*(rng-Rv[it-1]));
              
              // here we save the modal sum = (4*Pi*reduced pressure*sqrt(rho(zs))); 
              // p_red(r,z) = p(r,z)/sqrt(rho(z)); see eq. 25 on page 3 in in DV's Modess notes
              // i.e. to obtain the actual physical pressure multiply the modal sum below by sqrt(rho(z)/rho(zs))
              PP    = PP    + 4.0*Pi*A_prev[m]*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1])); 
              PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*real(k_prev[m])*(rng-Rv[it-1])); 
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

*/



// DV 20150930 - has the option to have sqrt(rho(zrcv)/rho(zsrc)) factor 
// in the Ocean Acoustics book the modes are normalized as 
// integral( [v_m/sqrt(rho(z))]^2*dz ) = 1
// but in this code the modes are normalized as integral(V_m^2*dz) = 1 
// then we need the factor sqrt(rho(zrcv)/rho(zsrc))

int saveTLoss1D (int Nm, int n_zsrc, double sqrtrho_s, int n_zrcv, double sqrtrho_r, int it, vector<double> Rv, double rng_step, double **v_prev, complex<double> *k_prev, complex<double> *A_prev, complex<double> *A_prev_ll, const char *wa, double *prng) {

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
      if (it==1) { // for region 1 only
          for (m=0; m<Nm; m++) {
              // pressure on the ground   - this version of PP doesn't work for big rng 
              // because NaN appear in the calculation
              //PP    = PP + 4.0*Pi*A_prev[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[0]/rng)*exp(I*k_prev[m]*(rng-Rv[0]));
              //PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*sqrtrho_s*v_prev[0][m]*sqrt(Rv[0]/rng)*exp(I*real(k_prev[m])*(rng-Rv[0]));
              
              // see eq. 11b in DV notes NMRD-OWCM: this should coincide with formula 5.1.14 in Ocean Acoustics page 274
              
              //// use the next 2 lines if modes are scaled by 1/sqrt(rho)
              //PP    = PP    + sqrtrho_r/sqrtrho_s*v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              //PP_ll = PP_ll + sqrtrho_r/sqrtrho_s*v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*real(k_prev[m])*rng)/sqrt(real(k_prev[m]));
              
              // just the modal sum -            
              PP    = PP    + v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              PP_ll = PP_ll + v_prev[n_zrcv][m]*v_prev[n_zsrc][m]*exp(I*real(k_prev[m])*rng)/sqrt(real(k_prev[m]));              
              
          }         
          
          // no sqrtrho_r/sqrtrho_s factor
          if (1) {
          PP    = 4.0*Pi*I*exp(-I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP;
          PP_ll = 4.0*Pi*I*exp(-I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP_ll;
          }
          
          if (0) {
          PP    = 4.0*Pi*(sqrtrho_r/sqrtrho_s)*I*exp(-I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP;
          PP_ll = 4.0*Pi*(sqrtrho_r/sqrtrho_s)*I*exp(-I*Pi/4.0)/sqrt(8.0*Pi*rng)*PP_ll;
          }
          
      }
      else { // for regions 2,3,...
          for (m=0; m<Nm; m++) {         
              // pressure on the ground (or TLoss if a 4*PI factor appears here) - eq 33 in DV's notes
              //// use the next 2 lines if modes are scaled by 1/sqrt(rho) 
              //PP    = PP + 4.0*Pi*A_prev[m]*sqrtrho_r*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1])); 
              //PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*sqrtrho_r*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*real(k_prev[m])*(rng-Rv[it-1]));
              
              // here we save the modal sum = (4*Pi*reduced pressure*sqrt(rho(zs))); 
              // p_red(r,z) = p(r,z)/sqrt(rho(z)); see eq. 25 on page 3 in in DV's Modess notes
              // i.e. to obtain the actual physical pressure multiply the modal sum below by sqrt(rho(z)/rho(zs))
              PP    = PP    + 4.0*Pi*A_prev[m]*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1])); 
              PP_ll = PP_ll + 4.0*Pi*A_prev_ll[m]*v_prev[n_zrcv][m]*sqrt(Rv[it-1]/rng)*exp(I*real(k_prev[m])*(rng-Rv[it-1])); 
          }
          
          // scaling by (sqrtrho_r/sqrtrho_s) to get actual pressure
          if (0) {
          PP    = (sqrtrho_r/sqrtrho_s)*PP;
          PP_ll = (sqrtrho_r/sqrtrho_s)*PP_ll;
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







/*
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
*/


/*
int saveTLoss2D (int Nm, int Nz_grid, int n_zsrc, double dz_km, int stepj, \
                  double sqrtrho_s, int it, vector<double> Rv, double rng_step, \
                  double *rho_prev, double **v_prev, complex<double> *k_prev, \
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
                  // use next line if modes are scaled by 1/sqrt(rho)
                  //PP    = PP + sqrt(rho_prev[j])/sqrtrho_s*v_prev[j][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);

                  PP    = PP + v_prev[j][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              }
              PP = eIpir*PP; // saving the modal sum * 4*pi - see eq. 25 in DV's Modess notes
              fprintf(fp2d,"%10.3f %10.3f %18.8e %18.8e\n", rng_km, j*dz_km, real(PP), imag(PP));     
          }
          fprintf(fp2d  , "\n");
      }
      else { // regions 2,3 ...
          for (j=0; j<Nz_grid; j=j+stepj) {
              PP    = 0.0 + I*0.0;
              for (m=0; m<Nm; m++) {
                  // see eq. 25 in DV notes NMRD-OWCM
                  // use next line if modes are scaled by 1/sqrt(rho)
                  //PP = PP + 4.0*Pi*A_prev[m]*sqrt(rho_prev[j])*v_prev[j][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1]));
                  PP = PP + 4.0*Pi*A_prev[m]*v_prev[j][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1]));
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
*/


// DV 20150930 
// in the Ocean Acoustics book the modes are normalized as 
// integral( [v_m/sqrt(rho(z))]^2*dz ) = 1
// but in this code the modes are normalized as integral(V_m^2*dz) = 1 
// then we need the factor sqrt(rho_prev[j])/sqrtrho_s
int saveTLoss2D (int Nm, int Nz_grid, int n_zsrc, double dz_km, int stepj, \
                  double sqrtrho_s, int it, vector<double> Rv, double rng_step, \
                  double *rho_prev, double **v_prev, complex<double> *k_prev, \
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
                  // use next line if modes are scaled by 1/sqrt(rho)
                  //PP    = PP + sqrt(rho_prev[j])/sqrtrho_s*v_prev[j][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);

                  PP    = PP + v_prev[j][m]*v_prev[n_zsrc][m]*exp(I*k_prev[m]*rng)/sqrt(k_prev[m]);      
              }
              PP = eIpir*PP; // saving the modal sum * 4*pi 

              // in the Ocean Acoustics book the modes are normalized as 
              // integral( [v_m/sqrt(rho(z))]^2*dz ) = 1
              // but in this code the modes are normalized as integral(V_m^2*dz) = 1 
              // then we need the factor sqrt(rho_prev[j])/sqrtrho_s
              if (0) {
              PP = PP*sqrt(rho_prev[j])/sqrtrho_s;  // see eq. 25 in DV's Modess notes
              }
       
              fprintf(fp2d,"%10.3f %10.3f %18.8e %18.8e\n", rng_km, j*dz_km, real(PP), imag(PP));     
          }
          fprintf(fp2d  , "\n");
      }
      else { // regions 2,3 ...
          for (j=0; j<Nz_grid; j=j+stepj) {
              PP    = 0.0 + I*0.0;
              for (m=0; m<Nm; m++) {
                  // see eq. 25 in DV notes NMRD-OWCM
                  PP = PP + 4.0*Pi*A_prev[m]*v_prev[j][m]*sqrt(Rv[it-1]/rng)*exp(I*k_prev[m]*(rng-Rv[it-1]));
              }
              
              // in the Ocean Acoustics book the modes are normalized as 
              // integral( [v_m/sqrt(rho(z))]^2*dz ) = 1
              // but in this code the modes are normalized as integral(V_m^2*dz) = 1 
              // then we need the factor sqrt(rho_prev[j])/sqrtrho_s
              if (0) {
              PP = PP*sqrt(rho_prev[j])/sqrtrho_s;
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


// new prototype: DV 20150929
int getRegionBoundaries(bool flg, double maxrange, double req_profile_step, string prf_ranges_km, int & Nprofiles, vector<double> *R) {

// this function populates vector R with the Region boundaries.
// if in the main() the following option is used e.g.:
//  --use_profile_ranges_km 30_100_500_700.5 
// then Rv will be [30, 30, 100, 500, 700.5]. Note that the first entry is doubled - see algorithm
// otherwise
// Rv is constructed from req_profile_step: [step:step:maxrange]

  unsigned int i;
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
      Nprofiles = Rv.size()-1;
      //for (i=0; i<Rv.size(); i++) cout << " " << Rv[i] << endl;
  }
  else {

      // this next line works on PCs but apparently not on all Macs?
      Nprofiles = (int) ceil(maxrange/req_profile_step); // at least one profile is ensured

      //// an alternative to the previous line
      //if ( fabs( (maxrange/req_profile_step) - ((int) (maxrange/req_profile_step) ) ) < 1e-5 ) { 
      //// if the ratio maxrange/req_profile_step is an integer
      //  Nprofiles = (int)  (maxrange/req_profile_step) ; // Nprofiles should be at least 1
      //}
      //else {
      //  Nprofiles = (int) ( (maxrange/req_profile_step) + 1.0 ); // Nprofiles should be at least 1
      //}

      //cout << "In getRegionBoundaries(...): " << endl;
      //cout << "maxrange = " << maxrange << endl;
      //cout << "req_profile_step = " << req_profile_step << endl;
      //cout << "Nprofiles = " << Nprofiles << endl;

      vector<double> Rv1 (Nprofiles+1, 0.0);
      for (i=1; i< Rv1.size(); i++) {
          Rv1[i] = i*req_profile_step;  
      } 
      Rv1[Rv1.size()-1] = maxrange; // last element adjusted to maxrange
      //Rv1[0] = Rv1[1]; // the convention is the first 2 elements are equal - see algorithm
      Rv.swap(Rv1);
  }
  Rv[0] = Rv[1]; // by convention the first 2 elements are equal - see algorithm
  //cout << "Nprofiles = " << (*Nprofiles) << endl;
  cout << "Rv size   = " << Rv.size() << endl;
  for (i=0; i<Rv.size(); i++) {
    printf("Rv[%i] = %g\n", i, Rv[i]);
  }
  
  (*R).swap(Rv);
  
  //cout << "End of getRegionBoundaries(...): " << endl;
  
  return 0;     
}
