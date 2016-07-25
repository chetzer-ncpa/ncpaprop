#include <stdexcept>
#include <ctime>
#include <dirent.h>
#include <fftw3.h>
#include "anyoption.h"
#include "ModBB_lib.h"

using namespace NCPA;
using namespace std;

#ifndef Pi
#define Pi 3.141592653589793
#endif

#define FFTN 32*1024
#define MAX_MODES 4000
//#define MIN(X,Y) ((X) < (Y) ? : (X) : (Y)) 



int writeDispersion(int select_modes, double dz, double z_src, double freq, complex<double> *k_pert, double **v_s)
{
  int i;
  int n_zsrc = ceil(z_src/dz)+1;
  char dispersion_file[40];

  sprintf(dispersion_file,"dispersion_%e.nm", freq);

  FILE *dispersion = fopen(dispersion_file,"w");
  fprintf(dispersion,"%e   %d",freq,select_modes);
  for(i=0;i<select_modes;i++) {
      fprintf(dispersion,"   %.12e   %.12e",real(k_pert[i]),imag(k_pert[i]));
      fprintf(dispersion,"   %.12e   %.12e",(v_s[n_zsrc][i]),(v_s[0][i]));
  }
  fprintf(dispersion,"\n");
  fclose(dispersion);

  printf("           file %s created\n", dispersion_file);
  return 0;
}


int readDispersion_bb_ascii(string filen, double *f_vec, int *mode_count, \
								double *prho_zsrc, double *prho_zrcv, double **re_k, double **im_k, \
								double **mode_S, double **mode_R)
{
  // this function reads the file created with writeDispersion_bb_ascii()
  // reads freq #modes rho(@z=src), rho(@z=rcv), k_horiz, modes(@z=src), modes(@z=receiverheight)
  int i,m,n;
  double x1,x2,y1,y2;
  FILE *fp = fopen(filen.c_str(),"r");
  if(fp==NULL){
      std::ostringstream es;
      es << "file << " << filen.c_str() << " could not be opened.\n";
      throw invalid_argument(es.str());
  }
  //read_header(f); // use it if file has header
  printf("--> Opening dispersion file %s ...if the program chokes here the file may be corrupt.\n",\
                             filen.c_str());
  // read line by line
  n=0;
  while(fscanf(fp,"%lf",&f_vec[n])==1){
      i = fscanf(fp,"%d",&mode_count[n]);
      i = fscanf(fp,"%lf",prho_zsrc);
      i = fscanf(fp,"%lf",prho_zrcv);
      for(m=0;m<mode_count[n];m++){
          if(fscanf(fp,"%lf %lf %lf %lf",&x1,&x2,&y1,&y2)!=4){
              fprintf(stderr,"%s: read_dispersion_data format error\n",filen.c_str());
              exit(1);
          }
          re_k[n][m]  =x1;
          im_k[n][m]  =x2;
          mode_S[n][m]=y1;
          mode_R[n][m]=y2;
      }
      ++n;
  }
  fclose(fp);
  printf("--> Found modes at %d positive frequencies in file %s\n",n, filen.c_str());
  return 0; 
}

int count_rows_arbcol(const char *filename) {
  int answer,c;
  FILE *f=fopen(filename,"r");

  if(f==NULL){
      std::ostringstream es;
      es << "file << " << filename << " could not be opened.\n";
      throw invalid_argument(es.str());
  }
  answer=0;
  //read_header(f);
  while((c=getc(f))!=EOF){
      if(c=='\n') answer=answer+1;
  }
  fclose(f);
  return answer;
}
 
 
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


double half_hann(int begin,int end,int i) {
  double answer;
  if(i<begin) answer=1.0;
  else if((begin<=i) && (i<=end)){
      answer=0.5*(cos(Pi*((i-begin))/((int)(end-begin)))+1.0);
  }
  else answer=0.0;
  return answer;
}


// DV 20151021: Added some code to save fft and ifftw results to compare with Matlab
// we have established the relationship between FFTW_BACKWARD and Matlab ifft 
//      FFTW_BACKWARD(x,NFFT) = NFFT*ifft_matlab(x, NFFT)
// and  FFTW_FORWARD(x,NFFT)  = fft_matlab(x, NFFT)
// Thus spectrum = FFTW_BACKWARD(pulse,NFFT)*dt
//      pulse    = FFTW_FORWARD(spectrum,NFFT)*df
// where df = f_max/f_step and dt = 1/(NFFT*df)
//
// Added NFFT as argument; normalizes the builtin pulse
// the source spectrum is loaded/computed in this function
int get_source_spectrum( \
								int n_freqs, int NFFT, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile) 
{
  int i;
  double dt, fmx, scale;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  fftw_plan p;

  fmx = ((double)NFFT)*f_step; // max frequency; the resulting time interval is dt = 1/fmx
  dt  = 1.0/fmx;
  
  if (src_flg==0) { // use spectrum of a delta function => to get the impulse response;
      for(i=0; i<n_freqs; i++) {
          dft_vec[i] = 1.0 + I*0.0; // 'ones' for all positive frequencies
          dft_vec[i] = dft_vec[i]*exp(I*2.0*Pi*f_vec[i]*(double)n_freqs/2.0); // add delay - DV 20150720
      }
  }
  else if (src_flg==1) { //use the built-in Roger pulse synthesised from a given spectrum          
          
      // make sure the dispersion file contains what we need to generate a good pulse
      // pulse center frequency f_center should be set <= f_max/5
      
      // f_center is initialized with negative value if the option --f_center was not used
      // if that's the case redefine f_center to be = f_max/5;
      if (f_center < 0) {
          f_center = f_vec[n_freqs-1]/5.0;
      }
      
      // get the scale for the built-in pulse
      if ((f_center > 0) & (f_center <= f_vec[n_freqs-1]/5.0+1.0E-16)) {	
          scale = 1/f_center;
      }
      else {
          std::ostringstream es;
          es << endl << "f_center = " << f_center << " Hz is too large." << endl
             << "For the built-in pulse f_center should be set smaller than f_max/5 = " 
             << f_vec[n_freqs-1]/5.0;
          throw invalid_argument(es.str());
      }
      
      f=fopen("builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          //dft_vec[i]=f_step*pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale);
          //dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale); // source spectrum
          dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale)*exp(I*2.0*Pi*f_vec[i]*5.0*scale/4.0); // DV 20150720 - time delay of '5.0*scale/4.0' added; source spectrum
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "--> The built-in source spectrum is saved in file 'builtin_source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;


      // The built-in source spectrum 'dft_vec' gives a time domain source pulse that 
      // is not normalized. Next we will normalize it to have a max positive
      // amplitude of one. So 1. we go into time domain; 2. normalize the pulse
      // then 3. come back to frequency domain with a 'normalized' dft_vec[].

      // Step 1 (of 3): from dft_vec to time domain pulse_vec
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; no need for the left zero pad  if fmin=f_step
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec has src spectrum in it
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // temporarily save the non-normalized pulse
      if (1) {
        f = fopen("unnormalized_pulse_input.dat","w");
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "--> Unnormalized Initial source waveform saved in file " 
             << "'unnormalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }

      // Step 2 (of 3) - normalize the builtin pulse to an amplitude of 1
      double pmax = 0.0; // will store the max(abs) value
      int    imax = 0;   // index of max
      for (i=0; i<NFFT; i++) {
        if ( pmax<abs(real(pulse_vec[i])) ) {
          pmax = abs(real(pulse_vec[i]));
          imax = i;
        };
      }
      
      // the pulse maximum should not be zero at this point; if it is, something went wrong
      if (pmax==0.0) {
        std::ostringstream es;
        es << endl << "The built-in pulse appears to be all zeros." << endl
           << "Check the parameters defining the built-in source spectrum "
           << " such as f_max, f_center." << endl;
        throw std::invalid_argument(es.str());
      }
      
      // normalize the time domain pulse; the positive max amplitude will be =1.
      // first redefine pmax to have the sign of real(pulse_vec[i])
      pmax = real(pulse_vec[imax]);
      for (i=0; i<NFFT; i++) {
        pulse_vec[i] = real(pulse_vec[i])/pmax; 
      }

      // temporarily save the normalized pulse
      if (0) {
        f = fopen("normalized_pulse_input.dat","w");
          //for(i=0;i<(NFFT*f_step*5*scale);i++) {
          for(i=0;i<NFFT;i++) {
            	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "Normalized Initial source waveform saved in file "
             << "'normalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }
          
      // Step 3 (of 3) - go from the normalized pulse back to the frequency domain
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      //
      // we have established the relationship between FFTW_BACKWARD and Matlab ifft 
      // FFTW_BACKWARD(x,NFFT) = NFFT*ifft_matlab(x, NFFT)
      //

      // temporary save of pure FFTW_BACKWARD(normalized pulse)
      if (0) {
              f = fopen("FFTW_BACKWARD_of_norm_pulse.dat","w");
          //for(i=0;i<(NFFT*f_step*5*scale);i++) {
          for(i=0;i<NFFT;i++) {
            	fprintf(f,"%14.9f %15.6e %15.6e\n", i*f_step, real(arg_vec[i]), imag(arg_vec[i]));
        }      
        fclose(f);
        cout << "FFTW_BACKWARD (normalized pulse) time domain waveform saved in file "
             << "'FFTW_BACKWARD_of_norm_pulse.dat' " << endl
             << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;  
             
      }
      
      
      // temporary test showing FFTW_FORWARD( FFTW_BACKWARD(normalized_pulse_vec)).
      // Our choice here:
      // FFTW_BACKWARD(pulse) takes us to frequency domain
      // FFTW_FORWARD() take us to time domain
      if (0) {
      // here arg_vec = FFTW_BACKWARD(normalized_pulse_vec)
      //
      // perform FFTW_FORWARD to obtain the pulse at the source (time domain) 
      // of NFFT points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p); 
      
      // now pulse_vec stores the result of 
      // FFTW_FORWARD( FFTW_BACKWARD(normalized_pulse_vec)). Call this P2.
      // P2 is in fact equal to P2 = NFFT*normalized_pulse_vec i.e.
      // FFTW_FORWARD( FFTW_BACKWARD(normalized_pulse_vec))= NFFT*normalized_pulse_vec 

        f = fopen("FFTW_FORWARD_arg_vec.dat","w");
          //for(i=0;i<(NFFT*f_step*5*scale);i++) {
          for(i=0;i<NFFT;i++) {
            	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "FFTW_FORWARD_arg_vec.dat time domain waveform saved in file "
             << "'FFTW_FORWARD_arg_vec.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;      
      
      }    

      
      // We just computed arg_vec = FFTW_BACKWARD(normalized pulse_vec))
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt; // multiply by the dt factor
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
                 
      //f=fopen("norm_builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      f=fopen("source_spectrum_input.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "--> The built-in source spectrum of the built-in normalized source pulse"
           << "is saved in file 'source_spectrum_input.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;           
       
  }
  else if (src_flg==2) { //use custom source spectrum from a file
      cout << "Loading source spectrum from file " << srcfile << endl;
      double *freqv;
      freqv = new double [n_freqs];
      load_source_spectrum(srcfile, freqv, dft_vec, n_freqs);
      if (fabs(freqv[0]-f_vec[0])>1.0e-10) { // compare frequencies
          cerr << "The frequencies from the source spectrum file " << srcfile
               << " do not appear to match the ones from the provided dispersion file" 
               << endl << " ... aborting." << endl;
               exit(1);
      }
      delete [] freqv;
  }

  // obtain the pulse at the source: either from builtin or user provided
  if ((src_flg<3)) { // common block only if the source spectrum (at positive frequencies) 
                     // is available: either the builtin spectrum or loaded from a file
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; the left zero pad
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec contains src spectrum times df
    
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      // The multiplication my df to complete the Fouruier integral should have
      // been done above.
      // Note though that since arg_vec contains only the positive freq spectrum
      // we'll have to double the result when we convert to the time domain
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      
  }
  else if (src_flg==3) { // use custom source pulse (time domain) from a file
      double dt;
      vector<double> t, tdp;
      load_source_pulse_td(srcfile, t, tdp);
      dt = t[1]-t[0];
      for (i=0; i<(int)t.size(); i++) {
          pulse_vec[i] = tdp[i]; // automatic conversion to complex<double>
          //cout << i << "  " << t[i] << " " << tdp[i] << endl;
      }
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt;
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
  }
  else { // the code should never get here
      cerr << "Unknown option " << src_flg << endl;
      exit(1);
  }
  
  // save initial pulse: pulse_vec
  if (1) {
      f = fopen("source_waveform_input.dat","w");
      if (src_flg==1) { // the builtin pulse
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 because we only had the spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, 2*real(pulse_vec[i])); 
        }
      fclose(f);
      }
      else { // if not the builtin pulse
        double fct = 2;
        if (src_flg==3) fct = 1;  // use custom source pulse (time domain) from a file
        
        for(i=0;i<NFFT/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 if we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, fct*real(pulse_vec[i]));
        }    
      }
      cout << "--> Initial source waveform saved in file 'source_waveform_input.dat' " << endl
           << " with format: | Time (s) | Amplitude |" << endl;
           //<< " with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
  }

  //// save arg_vec to file
  //printf("in pulse_prop: saving arg_vec to file arg_vec.dat; n_freqs=%d\n", n_freqs);
  //fmx = ((double)FFTN)*f_step;
  //f   = fopen("arg_vec.dat","w");
  //for(i=0;i<FFTN; i++) {
  //    fprintf(f,"%e   %e   %e\n", 1.0*i/fmx, real(arg_vec[i]), imag(arg_vec[i]));
  //}
  //fclose(f);
  //

  return 0;
}



/*
// DV 20151019: Added NFFT as argument; normalizes the builtin pulse
// the source spectrum is loaded/computed in this function
int get_source_spectrum( \
								int n_freqs, int NFFT, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile) 
{
  int i;
  double fmx, scale;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  fftw_plan p;

  fmx = ((double)NFFT)*f_step; // max frequency
  
  if (src_flg==0) { // use spectrum of a delta function => to get the impulse response;
      for(i=0; i<n_freqs; i++) {
          dft_vec[i] = 1.0 + I*0.0; // 'ones' for all positive frequencies
          dft_vec[i] = dft_vec[i]*exp(I*2.0*Pi*f_vec[i]*(double)n_freqs/2.0); // add delay - DV 20150720
      }
  }
  else if (src_flg==1) { //use the built-in Roger pulse synthesised from a given spectrum          
          
      // make sure the dispersion file contains what we need to generate a good pulse
      // pulse center frequency f_center should be set <= f_max/5
      
      // f_center is initialized with negative value if the option --f_center was not used
      // if that's the case redefine f_center to be = f_max/5;
      if (f_center < 0) {
          f_center = f_vec[n_freqs-1]/5.0;
      }
      
      // get the scale for the built-in pulse
      if ((f_center > 0) & (f_center <= f_vec[n_freqs-1]/5.0+1.0E-16)) {	
          scale = 1/f_center;
      }
      else {
          std::ostringstream es;
          es << endl << "f_center = " << f_center << " Hz is too large." << endl
             << "For the built-in pulse f_center should be set smaller than f_max/5 = " 
             << f_vec[n_freqs-1]/5.0;
          throw invalid_argument(es.str());
      }
      
      f=fopen("builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          //dft_vec[i]=f_step*pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale);
          //dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale); // source spectrum
          dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale)*exp(I*2.0*Pi*f_vec[i]*5.0*scale/4.0); // DV 20150720 - time delay of '5.0*scale/4.0' added; source spectrum
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The built-in source spectrum is saved in file 'builtin_source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;


      // The built-in source spectrum 'dft_vec' gives a time domain source pulse that 
      // is not normalized. Next we will normalize it to have a max positive
      // amplitude of one. So 1. we go into time domain; 2. normalize the pulse
      // then 3. come back to frequency domain with a 'normalized' dft_vec[].

      // Step 1 (of 3): from dft_vec to time domain pulse_vec
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; no need for the left zero pad  if fmin=f_step
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec has src spectrum in it
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // temporarily save the non-normalized pulse
      if (1) {
        f = fopen("unnormalized_pulse_input.dat","w");
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "Unnormalized Initial source waveform saved in file " 
             << "'unnormalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }

      // Step 2 (of 3) - normalize the builtin pulse to an amplitude of 1
      double pmax = 0.0; // will store the max(abs) value
      int    imax = 0;   // index of max
      for (i=0; i<NFFT; i++) {
        if ( pmax<abs(real(pulse_vec[i])) ) {
          pmax = abs(real(pulse_vec[i]));
          imax = i;
        };
      }
      
      // the pulse maximum should not be zero at this point; if it is, something went wrong
      if (pmax==0.0) {
        std::ostringstream es;
        es << endl << "The built-in pulse appears to be all zeros." << endl
           << "Check the parameters defining the built-in source spectrum "
           << " such as f_max, f_center." << endl;
        throw std::invalid_argument(es.str());
      }
      
      // normalize the time domain pulse; the positive max amplitude will be =1.
      // first redefine pmax to have the sign of real(pulse_vec[i])
      pmax = real(pulse_vec[imax]);
      for (i=0; i<NFFT; i++) {
        pulse_vec[i] = real(pulse_vec[i])/pmax; 
      }

      // temporarily save the normalized pulse
      if (1) {
        f = fopen("normalized_pulse_input.dat","w");
          for(i=0;i<(NFFT*f_step*5*scale);i++) {
            	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, real(pulse_vec[i])); 
        }      
        fclose(f);
        cout << "Normalized Initial source waveform saved in file "
             << "'normalized_pulse_input.dat' " << endl
             << " with format: | Time (s) | Amplitude |" << endl;
      }
          
      // Step 3 (of 3) - go from the normalized pulse back to the frequency domain
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      //double dt = f_vec[n_freqs-1]/n_freqs;
      //double dt = double (n_freqs)/((double)NFFT);
      
      // We just computed arg_vec = FFTW_BACKWARD(normalized pulse_vec))
      // Since FFTW has the property: FFTW_FORWARD(FFTW_BACKWARD(x)) = NFFT*x then
      // FFTW_FORWARD(arg_vec) = NFFT*normalized_pulse. 
      // Thus if we divide by arg_vec by NFFT below we'll reobtain the
      // normalized_pulse when we compute FFTW_FORWARD(arg_vec)
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]/((double)NFFT); // divide by NFFT in FFTW
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
                 
      f=fopen("norm_builtin_source_spectrum.dat","w"); // save the source spectrum in this file
      for(i=0; i<n_freqs; i++) {
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The built-in source spectrum is saved in file 'norm_builtin_source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;           
           
           
           
           
           
  }
  else if (src_flg==2) { //use custom source spectrum from a file
      cout << "Loading source spectrum from file " << srcfile << endl;
      double *freqv;
      freqv = new double [n_freqs];
      load_source_spectrum(srcfile, freqv, dft_vec, n_freqs);
      if (fabs(freqv[0]-f_vec[0])>1.0e-10) { // compare frequencies
          cerr << "The frequencies from the source spectrum file " << srcfile
               << " do not appear to match the ones from the provided dispersion file" 
               << endl << " ... aborting." << endl;
               exit(1);
      }
      delete [] freqv;
  }

  if ((src_flg<3)) { // common block only if the source spectrum (at positive frequencies) 
                     // is available: either the builtin spectrum or loaded from a file
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; the left zero pad
      int i0 = i;
      //for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec contains src spectrum
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]; // arg_vec contains src spectrum     
      for( ; i<NFFT; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of NFFT points: 'pulse_vec'
      // Note though that since arg_vec contains only the positive freq spectrum
      // we'll have to double the result when we convert to the time domain
      //
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      
  }
  else if (src_flg==3) { // use custom source pulse (time domain) from a file
      double dt;
      vector<double> t, tdp;
      load_source_pulse_td(srcfile, t, tdp);
      dt = t[1]-t[0];
      for (i=0; i<(int)t.size(); i++) {
          pulse_vec[i] = tdp[i]; // automatic conversion to complex<double>
          //cout << i << "  " << t[i] << " " << tdp[i] << endl;
      }
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with NFFT points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( NFFT, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt;
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
  }
  else { // the code should never get here
      cerr << "Unknown option " << src_flg << endl;
      exit(1);
  }
  
  // save initial pulse: pulse_vec
  if (1) {
      f = fopen("source_waveform_input.dat","w");
      if (src_flg==1) {
        for(i=0;i<(NFFT*f_step*5*scale);i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 because we only had the spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, 2*real(pulse_vec[i])); 
        }
      fclose(f);
      }
      else { // if not the builtin pulse
        double fct = 2;
        if (src_flg==3) fct = 1;  // use custom source pulse (time domain) from a file
        
        for(i=0;i<NFFT/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 if we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, fct*real(pulse_vec[i]));
        }    
      }
      cout << "Initial source waveform saved in file 'source_waveform_input.dat' " << endl
           << " with format: | Time (s) | Amplitude |" << endl;
           //<< " with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
  }

  //// save arg_vec to file
  //printf("in pulse_prop: saving arg_vec to file arg_vec.dat; n_freqs=%d\n", n_freqs);
  //fmx = ((double)FFTN)*f_step;
  //f   = fopen("arg_vec.dat","w");
  //for(i=0;i<FFTN; i++) {
  //    fprintf(f,"%e   %e   %e\n", 1.0*i/fmx, real(arg_vec[i]), imag(arg_vec[i]));
  //}
  //fclose(f);
  //

  return 0;
}
*/


// the source spectrum is loaded/computed in this function
// this version uses FFTN from #define FFTN; 
// the newer version above has NFFT as argument
int get_source_spectrum( \
								int n_freqs, double f_step, double *f_vec, double f_center, \
								complex<double> *dft_vec, complex<double> *pulse_vec, \
								complex<double> *arg_vec, int src_flg, string srcfile) 
{
  int i;
  double fmx, scale;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  fftw_plan p;

  fmx = ((double)FFTN)*f_step; // max frequency
  
  if (src_flg==0) { // use spectrum of a delta function => to get impulse response;
      for(i=0; i<n_freqs; i++) {
          dft_vec[i] = 1.0 + I*0.0; // 'ones' for all positive freqencies
          dft_vec[i] = dft_vec[i]*exp(I*2.0*Pi*f_vec[i]*(double)n_freqs/2.0); // add delay - DV 20150720
      }
  }
  else if (src_flg==1) { //use the built-in Roger pulse synthesised from a given spectrum          
          
      // make sure the dispersion file contains what we need to generate a good pulse
      // pulse center ferquency f_center should be set <= f_max/5
      
      // f_center is initialized with negative value if the option --f_center was not used
      // if that's the case redefine f_center to be =f_max/5;
      if (f_center < 0) {
          f_center = f_vec[n_freqs-1]/5.0;
      }
      
      // get the scale for the built-in pulse
      if ((f_center > 0) & (f_center <= f_vec[n_freqs-1]/5.0+1.0E-16)) {	
          scale = 1/f_center;
      }
      else {
          std::ostringstream es;
          es << endl << "f_center = " << f_center << " Hz is too large." << endl
             << "For the built-in pulse f_center should be set smaller than f_max/5 = " 
             << f_vec[n_freqs-1]/5.0;
          throw invalid_argument(es.str());
      }
      
      f=fopen("source_spectrum.dat","w"); // save the source spectrum
      for(i=0; i<n_freqs; i++) {
          //dft_vec[i]=f_step*pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale);
          //dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale); // source spectrum
          dft_vec[i] = pulse_spec_fit(scale,f_vec[i])*exp(I*2.0*Pi*f_vec[i]*scale)*exp(I*2.0*Pi*f_vec[i]*5.0*scale/4.0); // DV 20150720 - time delay of '5.0*scale/4.0' added; source spectrum
          fprintf(f,"%14.9f %15.6e %15.6e\n", f_vec[i], real(dft_vec[i]), imag(dft_vec[i]));    
      }
      fclose(f);
      cout << "The source spectrum is saved in file 'source_spectrum.dat' " << endl
           << "with format: | Freq (Hz) | Re(S) | Imag(S) |" << endl;
  }
  else if (src_flg==2) { //use custom source spectrum from a file
      cout << "Loading source spectrum from file " << srcfile << endl;
      double *freqv;
      freqv = new double [n_freqs];
      load_source_spectrum(srcfile, freqv, dft_vec, n_freqs);
      if (fabs(freqv[0]-f_vec[0])>1.0e-10) { // compare frequencies
          cerr << "The frequencies from the source spectrum file " << srcfile
               << " do not appear to match the ones from the provided dispersion file" 
               << endl << " ... aborting." << endl;
               exit(1);
      }
      delete [] freqv;
  }

  if ((src_flg<3)) { // common block only if the source spectrum (at positive frequencies) is available
      
      //for(i=0;i*f_step<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum; //dv20131016 - commented out this line and inserted the next line: i=0;
      i = 0; //dv20131016 - this is a new line; the left zero pad above was unnecessary
      int i0 = i;
      for(i=i0; i<(i0+n_freqs); i++) arg_vec[i] = dft_vec[i-i0]*f_step; // arg_vec has src spectrum in it
      for( ; i<FFTN; i++) arg_vec[i]=0.0; //arg_vec[] zero pad to the right of the src spectrum
      
      //
      // perform fft to obtain the pulse at the source (time domain) of FFTN points: 'pulse_vec'
      //
      p=fftw_plan_dft_1d( FFTN, reinterpret_cast<fftw_complex*> (arg_vec), \
		        										reinterpret_cast<fftw_complex*> (pulse_vec), \
			        									FFTW_FORWARD,FFTW_ESTIMATE );
      fftw_execute(p);
      fftw_destroy_plan(p);
  }
  else if (src_flg==3) { // use custom source pulse (time domain) from a file
      double dt;
      vector<double> t, tdp;
      load_source_pulse_td(srcfile, t, tdp);
      dt = t[1]-t[0];
      for (i=0; i<(int)t.size(); i++) {
          pulse_vec[i] = tdp[i]; // automatic conversion to complex<double>
          //cout << i << "  " << t[i] << " " << tdp[i] << endl;
      }
      //
      // perform integral of (pulse(t))*exp(+iwt)*dt via fft
      // to obtain the source spectrum (with FFTN points): 'arg_vec'
      // FFTW_BACKWARD performs (pulse(t))*exp(+iwt); will multiply by dt factor afterwards
      // Note that there is no need to multiply by N as in Matlab;
      // note: this FFTW_BACKWARD lacks the 1/N factor that Matlab ifft() has
      // in other words: FFTW_FORWARD(FFTW_BACKWARD(x)) = N*x
      // 
      p=fftw_plan_dft_1d( FFTN, reinterpret_cast<fftw_complex*> (pulse_vec), \
            										reinterpret_cast<fftw_complex*> (arg_vec), \
	            									FFTW_BACKWARD,FFTW_ESTIMATE );    
      fftw_execute(p);
      fftw_destroy_plan(p);
      
      // multiply by the dt factor to complete the Fourier integral over time
      // note: it is expected that the energy in the pulse is concentrated well below f_max
      for (i=0; i<n_freqs; i++) {
          dft_vec[i] = arg_vec[i]*dt;
          //printf("dft[%d] = %g + %g i\n", i, real(dft_vec[i]), imag(dft_vec[i]));
      }
  }
  else {
      cerr << "Unknown option " << src_flg << endl;
      exit(1);
  }
  
  // save initial pulse: pulse_vec
  if (1) {
      f = fopen("source_waveform.dat","w");
      if (src_flg==1) {
        for(i=0;i<(FFTN*f_step*5*scale);i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 because we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, 2*real(pulse_vec[i])); 
        }
      fclose(f);
      }
      else { // if not the builtin pulse
        double fct = 2;
        if (src_flg==3) fct = 1;
        
        for(i=0;i<FFTN/2;i++) {
          	//fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx, real(pulse_vec[i]), imag(pulse_vec[i]));
          	// save 2*real(pulse_vec[i]); factor of 2 if we only had spectrum for positive freqs
          	fprintf(f,"%12.6f %15.6e\n", 1.0*i/fmx, fct*real(pulse_vec[i]));
        }    
      }
      cout << "Initial source waveform saved in file 'source_waveform.dat' " << endl
           << " with format: | Time (s) | Amplitude |" << endl;
           //<< " with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;
  }

  //// save arg_vec to file
  //printf("in pulse_prop: saving arg_vec to file arg_vec.dat; n_freqs=%d\n", n_freqs);
  //fmx = ((double)FFTN)*f_step;
  //f   = fopen("arg_vec.dat","w");
  //for(i=0;i<FFTN; i++) {
  //    fprintf(f,"%e   %e   %e\n", 1.0*i/fmx, real(arg_vec[i]), imag(arg_vec[i]));
  //}
  //fclose(f);
  //

  return 0;
}



// 20130710 - 'fft_pulse_prop' with comments
// dft_vec has the source spectrum (for positive frequencies)
// arg_vec will hold the (delayed Fourier pressure component)*source_spectrum*df
// pulse_vec is the fft(arg_vec) i.e. the time domain propagated waveform (pulse)
void fft_pulse_prop(double t0, int n_freqs, double df, double *f_vec, \
                    double range, complex<double> *dft_vec, \
                    complex<double> *pulse_vec, \
                    int *mode_count, double rho_zsrc, double rho_zrcv, \
                    double **re_k, double **im_k, \
                    double **mode_S, double **mode_R)
{
  int i,i0,j,smooth_space;
  double sqrt_rho_ratio = sqrt(rho_zrcv/rho_zsrc);
  complex<double> cup,pot,k_H,mode_prod,t_phase,*arg_vec;
  complex<double> I (0.0, 1.0);
  complex<double> expov8pir = I*exp(-I*Pi*0.25)/sqrt(8.0*Pi*range);
  
  fftw_plan p;

  if(FFTN < n_freqs){
      throw invalid_argument("fft too short (i.e. FFTN < n_freqs), exiting.");
  }
  
  arg_vec = new complex<double> [FFTN];
  
  
  // old left zero-pad; sometimes with a C compiler it gives i=2 even if fvec[0]=df
  // so I rewrote it below; DV
  //for(i=0;i*df<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum;
  
  // new left zero pad up to f_min present in the spectrum;
  for(i=0;(fabs(i*df-f_vec[0])>1.0e-10) & (i*df<f_vec[0]); i++) arg_vec[i]=0.0;
  //printf("in fft_pulse_prop:f_vec[0] = %f; df = %f; i=%d\n", f_vec[0], df, i);

  i0 = i;
  for(i=0;i<n_freqs;i++){
      t_phase=exp(-I*2.0*Pi*f_vec[i]*t0); // corresponding to reduced time t0
      cup=0.0;
      for(j=0;j<mode_count[i];j++){
          k_H = re_k[i][j]+I*im_k[i][j];
          pot = exp(I*range*k_H);
          mode_prod = mode_S[i][j]*mode_R[i][j];
          pot = mode_prod*pot;
          pot = pot/sqrt(k_H);
          cup = cup+pot; // modal sum: sum( exp(ikr)/sqrt(k)*Vr*Vs )
      }
      cup=cup*t_phase;
      
      // up to a factor ((delayed Fourier pressure component)*source_spectrum*df) 
      // (see eq. 5.14 pg 274 and eq. 8.9 pg 480 in Comp. Oc. Acoust. 1994 ed.)
      cup=cup*dft_vec[i]*df;

      // the reduced pressure is defined as: 
      // (actual pressure)/sqrt(rho(z)): p_reduced(r,z) = p(r,z)/sqrt(rho(z))

      // Note on mode scaling in this code vs. the modes in Computational Oc. Acoust. 1994: 
      // V_book =  sqrt(rho)*V_in_this_code; 
      // Thus the formula for the Fourier component of pressure using the modes in this code
      // is given in DV Modess notes eq. 25 pg. 3 and contains the factor sqrt_rho_ratio
      // as opposed to the factor 1/rho(z_s) in the book eq. 5.14 pg 274.
      arg_vec[i0+i] = expov8pir*cup*sqrt_rho_ratio; // note sqrt_rho_ratio
  }

  smooth_space=(int)floor(0.1*n_freqs); // smoothly zero out on right; as in RW (July 2012)

  for(i=n_freqs-smooth_space;i<n_freqs;i++){
      //arg_vec[i]=arg_vec[i]*half_hann(n_freqs-smooth_space,n_freqs-df,i); // old code : has df in it
      arg_vec[i]=arg_vec[i]*half_hann(n_freqs-smooth_space,n_freqs-1,i); // changed df to 1 as df doesn't make sense here
  }
  for(;i<FFTN;i++) arg_vec[i]=0.0; // right zero pad
  
  // save arg_vec for debugging purposes
  if (0) {
      FILE *fp = fopen("arg_vec.dat", "w");
      for (i=0; i<FFTN; i++) {
          fprintf(fp, "%12.6f %15.6e %15.6e\n", i*df, real(arg_vec[i]), imag(arg_vec[i]));
      }
      fclose(fp);
      printf("arg_vec saved in file 'arg_vec.dat' with format: i*df | Re(arg_vec) | Im(arg_vec)\n");
  }
   
  // 
  //perform fft of arg_vec to obtain the propagated time domain pulse
  //
  p=fftw_plan_dft_1d(FFTN, reinterpret_cast<fftw_complex*>(arg_vec), \
												   reinterpret_cast<fftw_complex*>(pulse_vec), \
												   FFTW_FORWARD,FFTW_ESTIMATE); 												 
  fftw_execute(p); 
  fftw_destroy_plan(p);
  delete [] arg_vec;
}  // end of debug version of 'fft_pulse_prop'



// -----------------------------------------------------------------------
// DV 20151017: Added NFFT as argument 
// 20130710 - 'fft_pulse_prop' with comments
// dft_vec has the source spectrum (for positive frequencies)
// arg_vec will hold the (delayed Fourier pressure component)*source_spectrum*df
// pulse_vec is the fft(arg_vec) i.e. the time domain propagated waveform (pulse)
void fft_pulse_prop(double t0, int n_freqs, int NFFT, double df, double *f_vec, \
                    double range, complex<double> *dft_vec, \
                    complex<double> *pulse_vec, \
                    int *mode_count, double rho_zsrc, double rho_zrcv, \
                    double **re_k, double **im_k, \
                    double **mode_S, double **mode_R)
{
  int i,i0,j,smooth_space;
  double sqrt_rho_ratio = sqrt(rho_zrcv/rho_zsrc);
  complex<double> cup,pot,k_H,mode_prod,t_phase,*arg_vec;
  complex<double> I (0.0, 1.0);
  complex<double> expov8pir = I*exp(-I*Pi*0.25)/sqrt(8.0*Pi*range);
  
  fftw_plan p;

  if(NFFT < n_freqs){
      throw invalid_argument("fft too short (i.e. NFFT < n_freqs), exiting.");
  }
  
  arg_vec = new complex<double> [NFFT];
  
  
  // old left zero-pad; sometimes with a C compiler it gives i=2 even if fvec[0]=df
  // so I rewrote it below; DV
  //for(i=0;i*df<f_vec[0];i++) arg_vec[i]=0.0; // left zero pad up to f_min present in the spectrum;
  
  // new left zero pad up to f_min present in the spectrum;
  for(i=0;(fabs(i*df-f_vec[0])>1.0e-10) & (i*df<f_vec[0]); i++) arg_vec[i]=0.0;
  //printf("in fft_pulse_prop:f_vec[0] = %f; df = %f; i=%d\n", f_vec[0], df, i);

  i0 = i;
  for(i=0;i<n_freqs;i++){
      t_phase=exp(-I*2.0*Pi*f_vec[i]*t0); // corresponding to reduced time t0
      cup=0.0;
      for(j=0;j<mode_count[i];j++){
          k_H = re_k[i][j]+I*im_k[i][j];
          pot = exp(I*range*k_H);
          mode_prod = mode_S[i][j]*mode_R[i][j];
          pot = mode_prod*pot;
          pot = pot/sqrt(k_H);
          cup = cup+pot; // modal sum: sum( exp(ikr)/sqrt(k)*Vr*Vs )
      }
      cup=cup*t_phase;
      
      // up to a factor ((delayed Fourier pressure component)*source_spectrum*df) 
      // (see eq. 5.14 pg 274 and eq. 8.9 pg 480 in Comp. Oc. Acoust. 1994 ed.)
      cup=cup*dft_vec[i]*df;

      // the reduced pressure is defined as: 
      // (actual pressure)/sqrt(rho(z)): p_reduced(r,z) = p(r,z)/sqrt(rho(z))

      // Note on mode scaling in this code vs. the modes in Computational Oc. Acoust. 1994: 
      // V_book =  sqrt(rho)*V_in_this_code; 
      // Thus the formula for the Fourier component of pressure using the modes in this code
      // is given in DV Modess notes eq. 25 pg. 3 and contains the factor sqrt_rho_ratio
      // as opposed to the factor 1/rho(z_s) in the book eq. 5.14 pg 274.
      arg_vec[i0+i] = expov8pir*cup*sqrt_rho_ratio; // note sqrt_rho_ratio
  }

  smooth_space=(int)floor(0.1*n_freqs); // smoothly zero out on right; as in RW (July 2012)

  for(i=n_freqs-smooth_space;i<n_freqs;i++){
      //arg_vec[i]=arg_vec[i]*half_hann(n_freqs-smooth_space,n_freqs-df,i); // old code : has df in it
      arg_vec[i]=arg_vec[i]*half_hann(n_freqs-smooth_space,n_freqs-1,i); // changed df to 1 as df doesn't make sense here
  }
  for(;i<NFFT;i++) arg_vec[i]=0.0; // right zero pad
  
  // save arg_vec for debugging purposes
  if (0) {
      FILE *fp = fopen("arg_vec.dat", "w");
      for (i=0; i<NFFT; i++) {
          fprintf(fp, "%12.6f %15.6e %15.6e\n", i*df, real(arg_vec[i]), imag(arg_vec[i]));
      }
      fclose(fp);
      printf("arg_vec saved in file 'arg_vec.dat' with format: i*df | Re(arg_vec) | Im(arg_vec)\n");
  }
   
  // 
  //perform fft of arg_vec to obtain the propagated time domain pulse
  //
  p=fftw_plan_dft_1d(NFFT, reinterpret_cast<fftw_complex*>(arg_vec), \
												   reinterpret_cast<fftw_complex*>(pulse_vec), \
												   FFTW_FORWARD,FFTW_ESTIMATE); 												 
  fftw_execute(p); 
  fftw_destroy_plan(p);
  delete [] arg_vec;
}  // end of debug version of 'fft_pulse_prop'
//-----------------------------------------------------------------------

// DV 20151017: Added NFFT as argument
// DV 20160717: Added no_attenuation flag
int pulse_prop_src2rcv_grid2(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, int NFFT,  double f_step, double *f_vec, \
					double f_center, int *mode_count, double rho_zsrc, double rho_zrcv, \
					double **re_k, double **im_k, double **mode_S, double **mode_R, \
					int src_flg, string srcfile, int pprop_src2rcv_flg, bool zero_attn_flg) 
{
  int i,n;
  double rr, tskip, fmx, t0;	
  complex<double> cup,*dft_vec,*pulse_vec,*arg_vec;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  
  // if the option --nfft was not passed to the main program then the 
  // NFFT default was -1. Otherwise it has the requested positive value. 
  // Here we make sure that whatever the current value of NFFT is 
  // it's not less than 4*n_freqs
  if (NFFT<n_freqs) {
    NFFT = 4*n_freqs;
  }
  
  
  // if zero attenuation requested set the imaginary part of the wavenumber = 0
  if (zero_attn_flg) {
    for (n=0; n<n_freqs; n++) {
      for(int m=0;m<mode_count[n];m++){
          im_k[n][m]  = 0.0;
      }
    }
  }

  dft_vec   = new complex<double> [n_freqs];
  pulse_vec = new complex<double> [NFFT];
  arg_vec   = new complex<double> [NFFT];

  fmx = ((double)NFFT)*f_step; // max frequency
								      
  get_source_spectrum(n_freqs, NFFT, f_step, f_vec,  f_center, \
								      dft_vec, pulse_vec, arg_vec, src_flg,  srcfile);								      
								      
	
	if (pprop_src2rcv_flg) { // propagation to one receiver at distance RR from source
	    cout << "--> Doing pulse propagation source-to-receiver at one range: " 
	         << R_start/1000.0 << " km" << endl;
	         
	    t0 = R_start/max_cel;
	    //
	    // fft propagation: note that here 'pulse_vec' is overwritten with the propagated pulse
	    //							
//      fft_pulse_prop( t0, n_freqs, f_step, f_vec, R_start, \
								      dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, \
								      re_k, im_k, mode_S, mode_R);
								      
      fft_pulse_prop( t0, n_freqs, NFFT, f_step, f_vec, R_start, \
								      dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, \
								      re_k, im_k, mode_S, mode_R);								      
      
      // save propagated pulse to file
	    // DV 20150514 - factor of two necessary because we only used 
	    // the positive freq. spectrum in the fft														
      f = fopen(filename,"w");
      for(i=0; i<NFFT; i++) {
          //fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, real(pulse_vec[i]), imag(pulse_vec[i]));
          fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, 2.0*real(pulse_vec[i]), 2.0*imag(pulse_vec[i]));
      }
      fclose(f);

      cout << "--> Propagated pulse saved in file: '" << filename << "'" << endl
           << "' with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;      

	}
	else { // propagation to several receivers
	    cout << "--> Propagating pulse from source-to-receivers on grid ..." << endl;					      
      printf("----------------------------------------------\n");
      printf("max_celerity     t0             R\n");
      printf("    m/s          sec            km\n");
      printf("----------------------------------------------\n");

      f=fopen(filename,"w");
      tskip = 0.0;
      for(n=0; n<=(int)(floor((R_end-R_start)/DR)); n++) {
          rr= R_start + DR*n;
          t0=tskip+rr/max_cel;
          printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

//          fft_pulse_prop( t0, n_freqs, f_step, f_vec, rr, \
					            dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, re_k, im_k, \
					            mode_S, mode_R);
					            
          fft_pulse_prop( t0, n_freqs, NFFT, f_step, f_vec, rr, \
					            dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, re_k, im_k, \
					            mode_S, mode_R);					            
							              
          for(i=0;i<NFFT;i++){
              //fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]));
              fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, 2.0*real(pulse_vec[i])); // factor of 2; DV20150930
          }
          fprintf(f,"\n");
      }
      fclose(f);
      printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
      printf("Time array length = %d; delta_t = %g s\n", NFFT, 1.0/fmx);
      printf("Propagation results saved in file: %s\n", filename);
      printf("with columns: R (km) | time (s) | pulse(R,t) |\n");
  }

  delete [] dft_vec;
  delete [] arg_vec;
  delete [] pulse_vec;

  return 0;
}

/*
// DV 20151017: Added NFFT as argument
int pulse_prop_src2rcv_grid2(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, int NFFT,  double f_step, double *f_vec, \
					double f_center, int *mode_count, double rho_zsrc, double rho_zrcv, \
					double **re_k, double **im_k, double **mode_S, double **mode_R, \
					int src_flg, string srcfile, int pprop_src2rcv_flg) 
{
  int i,n;
  double rr, tskip, fmx, t0;	
  complex<double> cup,*dft_vec,*pulse_vec,*arg_vec;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;
  
  // if the option --nfft was not passed to the main program then the 
  // NFFT default was -1. Otherwise it has the requested positive value. 
  // Here we make sure that whatever the current value of NFFT is 
  // it's not less than 4*n_freqs
  if (NFFT<n_freqs) {
    NFFT = 4*n_freqs;
  }

  dft_vec   = new complex<double> [n_freqs];
  pulse_vec = new complex<double> [NFFT];
  arg_vec   = new complex<double> [NFFT];

  fmx = ((double)NFFT)*f_step; // max frequency
								      
  get_source_spectrum(n_freqs, NFFT, f_step, f_vec,  f_center, \
								      dft_vec, pulse_vec, arg_vec, src_flg,  srcfile);								      
								      
	
	if (pprop_src2rcv_flg) { // propagation to one receiver at distance RR from source
	    cout << "--> Doing pulse propagation source-to-receiver at one range: " 
	         << R_start/1000.0 << " km" << endl;
	         
	    t0 = R_start/max_cel;
	    //
	    // fft propagation: note that here 'pulse_vec' is overwritten with the propagated pulse
	    //							
//      fft_pulse_prop( t0, n_freqs, f_step, f_vec, R_start, \
								      dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, \
								      re_k, im_k, mode_S, mode_R);
								      
      fft_pulse_prop( t0, n_freqs, NFFT, f_step, f_vec, R_start, \
								      dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, \
								      re_k, im_k, mode_S, mode_R);								      
      
      // save propagated pulse to file
	    // DV 20150514 - factor of two necessary because we only used 
	    // the positive freq. spectrum in the fft														
      f = fopen(filename,"w");
      for(i=0; i<NFFT; i++) {
          //fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, real(pulse_vec[i]), imag(pulse_vec[i]));
          fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, 2.0*real(pulse_vec[i]), 2.0*imag(pulse_vec[i]));
      }
      fclose(f);

      cout << "--> Propagated pulse saved in file: '" << filename << "'" << endl
           << "' with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;      

	}
	else { // propagation to several receivers
	    cout << "--> Propagating pulse from source-to-receivers on grid ..." << endl;					      
      printf("----------------------------------------------\n");
      printf("max_celerity     t0             R\n");
      printf("    m/s          sec            km\n");
      printf("----------------------------------------------\n");

      f=fopen(filename,"w");
      tskip = 0.0;
      for(n=0; n<=(int)(floor((R_end-R_start)/DR)); n++) {
          rr= R_start + DR*n;
          t0=tskip+rr/max_cel;
          printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

//          fft_pulse_prop( t0, n_freqs, f_step, f_vec, rr, \
					            dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, re_k, im_k, \
					            mode_S, mode_R);
					            
          fft_pulse_prop( t0, n_freqs, NFFT, f_step, f_vec, rr, \
					            dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, re_k, im_k, \
					            mode_S, mode_R);					            
							              
          for(i=0;i<NFFT;i++){
              //fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]));
              fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, 2.0*real(pulse_vec[i])); // factor of 2; DV20150930
          }
          fprintf(f,"\n");
      }
      fclose(f);
      printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
      printf("Time array length = %d; delta_t = %g s\n", NFFT, 1.0/fmx);
      printf("Propagation results saved in file: %s\n", filename);
      printf("with columns: R (km) | time (s) | pulse(R,t) |\n");
  }

  delete [] dft_vec;
  delete [] arg_vec;
  delete [] pulse_vec;

  return 0;
}
*/


// This version uses the #define FFTN; obsolete as of 20151017 when
// the --nfft option was added to the code.  ... But will keep it for now
int pulse_prop_src2rcv_grid2(\
          const char *filename,double max_cel, \
          double R_start,double DR,double R_end, \
					int n_freqs, double f_step, double *f_vec, \
					double f_center, int *mode_count, double rho_zsrc, double rho_zrcv, \
					double **re_k, double **im_k, double **mode_S, double **mode_R, \
					int src_flg, string srcfile, int pprop_src2rcv_flg) 
{
  int i,n;
  double rr, tskip, fmx, t0;	
  complex<double> cup,*dft_vec,*pulse_vec,*arg_vec;
  complex<double> I = complex<double> (0.0, 1.0);
  FILE *f;

  dft_vec   = new complex<double> [n_freqs];
  pulse_vec = new complex<double> [FFTN];
  arg_vec   = new complex<double> [FFTN];

  fmx = ((double)FFTN)*f_step; // max frequency
  
  
  get_source_spectrum(n_freqs,  f_step, f_vec,  f_center, \
								      dft_vec, pulse_vec, arg_vec, src_flg,  srcfile);
								      
	
	if (pprop_src2rcv_flg) { // propagation to one receiver at distance RR from source
	    cout << "--> Doing pulse propagation source-to-receiver at one range: " 
	         << R_start/1000.0 << " km" << endl;
	         
	    t0 = R_start/max_cel;
	    //
	    // fft propagation: note that here 'pulse_vec' is overwritten with the propagated pulse
	    //							
      fft_pulse_prop( t0, n_freqs, f_step, f_vec, R_start, \
								      dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, \
								      re_k, im_k, mode_S, mode_R);
      
      // save propagated pulse to file
	    // DV 20150514 - factor of two necessary because we only used 
	    // the positive freq. spectrum in the fft														
      f = fopen(filename,"w");
      for(i=0; i<FFTN; i++) {
          //fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, real(pulse_vec[i]), imag(pulse_vec[i]));
          fprintf(f,"%12.6f %15.6e %15.6e\n", 1.0*i/fmx+t0, 2.0*real(pulse_vec[i]), 2.0*imag(pulse_vec[i]));
      }
      fclose(f);

      cout << "--> Propagated pulse saved in file: '" << filename << "'" << endl
           << "' with format: | Time (s) | Re(pulse) | Imag(pulse) |" << endl;      

      /*
      // plot with gnuplot
      if (1) {
          double *xData, *yData;
          xData = new double [FFTN];
          yData = new double [FFTN];
          for(i=0; i<FFTN; i++) {
            xData[i] = 1.0*i/fmx + t0; 
            yData[i] = real(pulse_vec[i]);
          }
          plotXY(xData, yData, FFTN, R_start/1000.0);
          delete [] xData;
          delete [] yData;
      }
      */
	}
	else { // propagation to several receivers
	    cout << "--> Propagating pulse from source-to-receivers on grid ..." << endl;			      
      printf("----------------------------------------------\n");
      printf("max_celerity     t0             R\n");
      printf("    m/s          sec            km\n");
      printf("----------------------------------------------\n");

      f=fopen(filename,"w");
      tskip = 0.0;
      for(n=0; n<=(int)(floor((R_end-R_start)/DR)); n++) {
          rr= R_start + DR*n;
          t0=tskip+rr/max_cel;
          printf("%8.3f     %9.3f      %9.3f\n", max_cel, t0, rr/1000.0);

          fft_pulse_prop( t0, n_freqs, f_step, f_vec, rr, \
					            dft_vec, pulse_vec, mode_count, rho_zsrc, rho_zrcv, re_k, im_k, \
					            mode_S, mode_R);
							              
          for(i=0;i<FFTN;i++){
              //fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, real(pulse_vec[i]));
              fprintf(f,"%10.3f %12.6f %15.6e\n", rr/1000.0, 1.0*i/fmx, 2.0*real(pulse_vec[i])); // factor of 2; DV20150930
          }
          fprintf(f,"\n");
      }
      fclose(f);
      printf("f_step = %f   1/f_step = %f\n", f_step, 1.0/f_step);
      printf("Time array length = %d; delta_t = %g s\n", FFTN, 1.0/fmx);
      printf("Propagation results saved in file: %s\n", filename);
      printf("with columns: R (km) | time (s) | pulse(R,t) |\n");
  }

  delete [] dft_vec;
  delete [] arg_vec;
  delete [] pulse_vec;

  return 0;
}


complex<double> pulse_spec_fit(double scale, double x) { 
  double fnorm = ((sqrt(4.0+0.25)-0.5)*0.5);
  //double fnorm = 1.0;
  complex<double> answer, I;
  I = complex<double> (0.0, 1.0);
  // fit to effective source spectrum for the propane cannon
  // scale down the frequency to get infrasound

  x=scale*fnorm*x;

  // fit to fourier transform of signal at 10 m
  answer=0.2421405*x*exp((-x-x*x)/2)*exp(I*(-0.5*Pi));
  /* calibrate (times 10 divided by 2 for pressure doubling) and time shift by 0.01 sec */
  /* multiply by 4 pi to get a delta func source strength */
  answer=4.0*Pi*5.0*answer*exp(I*2.0*Pi*x*(0.01));

  return answer;
}


/*
int plotXY(double* xData, double* yData, int dataSize, double RR) {
  // helper function to plot with gnuplot from within C code
  FILE *gnuplotPipe,*tempDataFile;
  double x,y;
  int i;
  const char *tempDataFileName = "tempData.dat";
  //gnuplotPipe = popen("c:\\gnuplot\\bin\\pgnuplot -persist","w");
  gnuplotPipe = popen("gnuplot -persist","w"); //dv
  if (gnuplotPipe) {
      printf("Plotting with gnuplot...\n");
      fprintf(gnuplotPipe, "plot \"%s\" with lines lt 1 title ''\n",tempDataFileName);
      fprintf(gnuplotPipe, "set xlabel \"Time [s]\"\n");
      fprintf(gnuplotPipe, "set ylabel \"Amplitude\"\n");
      fprintf(gnuplotPipe, "set title 'Propagated pulse to range = %g km'; set grid; show grid \nreplot\n", RR);
      fflush(gnuplotPipe);
      tempDataFile = fopen(tempDataFileName,"w");
      for (i=0; i < dataSize; i++) {
          x = xData[i];
          y = yData[i];            
          fprintf(tempDataFile,"%e %e\n",x,y);     
      }        
      fclose(tempDataFile);               
      //remove(tempDataFileName);        
      fprintf(gnuplotPipe,"exit \n");  
      pclose(gnuplotPipe);  
      printf("...the plot should be persistent on the screen\n");
  } else {        
      printf("gnuplot not found...\n");    
  }
  return 0;
}


// developer's function
int plotInitialPulse(int src_flg, double freq, double fmax) {
  FILE *pipe; 
      // open a pipe to gnuplot and plot the results
      pipe = popen("gnuplot -persist","w");
  if (pipe) {
      char title[128];
      //save to a png file
      //fprintf(pipe, "set terminal png\n");
      //fprintf(pipe, "set output 'initial_pulse.png'\n", freq);
      
      if (src_flg==0) {
          sprintf(title, "Initial pulse (band-limited delta function)");
      }
      else if (src_flg==1) {
          if (freq<0) freq=fmax/5.0; // the built-in pulse has fcenter=fmax/5 (if fcenter not specified)
          sprintf(title, "Initial pulse; center frequency: %g Hz", freq);
      }
      else if (src_flg>1) {
          sprintf(title, "Source waveform");
      }
          
      //fprintf(pipe, "set title 'Initial pulse; center frequency: %g Hz'\n", freq);
      fprintf(pipe, "set title '%s'\n", title);
      fprintf(pipe, "set xlabel 'Time [s]'\n");
      fprintf(pipe, "set ylabel 'Amplitude'\n");
      fprintf(pipe, "set grid; show grid;\n");
      fprintf(pipe, "plot './source_waveform.dat' using 1:2 with lines lt 3 title 'initial pulse';\n");

      pclose(pipe);
  } else {        
      printf("gnuplot not found...\n");    
  }
  return 0;
}
*/


int getFile_list(string dir, list<string> &files, string pattern)
{
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
      a = string(dirp->d_name);
      pos = a.find(pattern);
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


complex<double> ***getPz1z2 (int I1, int I2, double r1, double r2, double dr, string dirn, list<string> files, bool zero_attn_flg) 
{
  int lenI, ier;
  int ir, j, jj, iz;
  int Nfreq, Nz_grid, Nz_subgrid, NN, nmods;
  int Nrng_steps = (int) (r2-r1)/dr+1;  

  double rho_zsrc, *rho_rcv, sqrtrho;
  double *kr, *ki;
  double *vs, **vr;
  double rr, df, dz, delZ, z_min, freq;
  complex<double> kH, epio8, epio8ovr, sumM;
  complex<double> II (0.0, 1.0);
  complex<double> ***Pfzr;  
  FILE *fp;
  list<string>::iterator it;  

  int    szi = sizeof(int);
  double szd = sizeof(double);
  int Nfiles = files.size();
  std::string fullname;

  lenI = I2-I1+1;
  epio8 = II*exp(-II*(Pi/4.0))/sqrt(8.0*Pi);

  // allocate 3-rank tensor
  Pfzr = c3Darray(Nfiles, Nrng_steps, lenI);

  // iterate over sorted file list
  jj = 0;
  for (it=files.begin(); it!=files.end(); ++it) {
  
      cout << "Iteration " << jj << "; processing file: " << (*it) << endl;

      fullname = dirn + "/" + (*it);
      fp = fopen(fullname.c_str(), "r");
      if (fp==NULL) {
          std::ostringstream es;
          es << "file << " << fullname << " could not be opened.\n";
          throw invalid_argument(es.str());
      }

      ier = fread(&freq,       szd, 1, fp);
      ier = fread(&nmods,      szi, 1, fp);		
      ier = fread(&Nfreq,      szi, 1, fp);
      ier = fread(&df,         szd, 1, fp);
      ier = fread(&Nz_grid,    szi, 1, fp);
      ier = fread(&dz,         szd, 1, fp);      
      ier = fread(&Nz_subgrid, szi, 1, fp);
      ier = fread(&delZ,       szd, 1, fp);
      ier = fread(&z_min,      szd, 1, fp);
      ier = fread(&rho_zsrc,   szd, 1, fp);
      
      NN = Nz_grid/Nz_subgrid; // integer division

      if (nmods>0) {
          kr      = new double [nmods];
          ki      = new double [nmods];
          vs      = new double [nmods];
          rho_rcv = new double [lenI];
          vr      = dmatrix(lenI, nmods);

          // read 4(ints) + 6(doubles) = 4*4+6*8 = 64 bytes so far
          ier = fread(kr, szd, nmods, fp);
          ier = fread(ki, szd, nmods, fp);
          ier = fread(vs, szd, nmods, fp);

          // at each z on the z-subgrid we saved (nmods+1) doubles i.e.
          // we offset to the requested altitude specified by I1: I1*(nmods+1)*szd bytes
          fseek(fp, (I1)*(nmods+1)*szd, SEEK_CUR);
            
          for (iz=0; iz<lenI; iz++) { // at each z level reading (nmods+1)*szd bytes
              ier = fread(&rho_rcv[iz], szd, 1, fp);  // read rho at all subgrid heights
              ier = fread(vr[iz], szd, nmods, fp);
          }

          rr = r1;
          for (ir=0; ir<Nrng_steps; ir++) {
              epio8ovr = epio8/sqrt(rr*rho_zsrc);
              for (iz=0; iz<lenI; iz++) {
                  sumM = complex<double> (0.0, 0.0);
                  sqrtrho = sqrt(rho_rcv[iz]);
                  for (j=0; j<nmods; j++) {
                      //kH =  kr[j] + II*ki[j];
                      // 20160717 taking zero_attn_flg into consideration 
                      kH =  kr[j] + (II*ki[j])*(1.0-(double)zero_attn_flg); // 20160717					
                      sumM += vs[j]*vr[iz][j]*sqrtrho*exp(II*kH*rr)/sqrt(kH);     
                  }
                  Pfzr[jj][ir][iz] = epio8ovr*sumM;
              }
              rr = rr + dr;
          }

          delete[] kr;
          delete[] ki;
          delete[] vs;
          delete[] rho_rcv;
          free_dmatrix(vr,lenI,nmods);
      }
      fclose(fp);
      jj++;
  } // end for loop by files
  //printf(" ... Pfzr done.\n");
  return Pfzr;
}


complex<double> **cmatrix(long nr, long nc) {
  // allocate a complex<double> matrix
  complex<double> **v;
  v = new complex<double> *[nr];
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


complex<double> ***c3Darray(size_t xlen, size_t ylen, size_t zlen)
{
  complex<double> ***p;
  size_t i, j;

  if ((p = (complex<double> ***) calloc(xlen, (sizeof(*p)))) == NULL) {
      perror("calloc 1");
      return NULL;
  }

  for (i=0; i < xlen; ++i)
      p[i] = NULL;

  for (i=0; i < xlen; ++i)
      if ((p[i] = (complex<double> **) calloc(ylen, (sizeof (*p[i])))) == NULL) {
          perror("calloc 2");
          free_c3Darray(p, xlen, ylen);
          return NULL;
      }

  for (i=0; i < xlen; ++i)
      for (j=0; j < ylen; ++j)
          p[i][j] = NULL;

  for (i=0; i < xlen; ++i)
      for (j=0; j < ylen; ++j)
          if ((p[i][j] = (complex<double> *) calloc(zlen, (sizeof (*p[i][j])))) == NULL) {
              perror("calloc 3");
              free_c3Darray(p, xlen, ylen);
              return NULL;
          }
  return p;
}


void free_c3Darray(complex<double> ***data, size_t xlen, size_t ylen)
{
  size_t i, j;
  for (i=0; i < xlen; ++i) {
      if (data[i] != NULL) {
          for (j=0; j < ylen; ++j)
              free(data[i][j]);
          free(data[i]);
      }
  }
  free(data);
}

// comparison, freq in filename.
bool compare_freq (string first, string second)
{
  string a, b;
  string patt1 ("_");
  string patt2 ("_nm.bin");
  size_t pos1, pos2;
  double x, y;
    
  pos2 = first.find(patt2);
  pos1 = first.rfind(patt1, pos2-1) + patt1.length();
  
  a = first.substr(pos1, pos2-pos1);
  b = second.substr(pos1, pos2-pos1);
  
  x = strtod(a.c_str(), NULL);
  y = strtod(b.c_str(), NULL);

  if (x<y) {return true; }
  else     {return false;}
}


// DV 20160717 - added zero_atten flag
int process2DPressure(double R_start_km, double width_km, double height_km, \
                      double c_ref, double tmstep, int ntsteps, double f_center, \
                      string framefn, string dirn, int src_flg, string srcfile, \
                      bool zero_attn_flg) 
{
  // compute the 2D spatial pressure field in (width_km by height_km) window 
  // starting at R_start_km
  // cout << "in process2DPressure():\n";
  
  char   filename [256]; // holds the frame file name(s)
  int    i, j, ier, m, n, k, Nx, Nz, NN, zs;
  int    I1, I2, lenI;
  int    Nfreq, Nz_grid, Nz_subgrid, nmods, Nfiles;
  double tmin, tmax, tt;
  double freq, fr, df, f_max;	
  double dz, delZ, z_min;
  double r1, r2, dr, sumP;
  double *f_vec;
  double **pp;

  complex<double> *S, *pulse_vec, *arg_vec; //source spectrum
  complex<double> II = complex<double> (0.0, 1.0);	
  complex<double> ***Pfzr;
  string fullname;
  FILE *fp;

  list<string> files;
  list<string>::iterator it;
  string pattern (".bin");

  // get and sort the files (they have the frequency in the name)
  getFile_list(dirn, files, pattern);
  files.sort(compare_freq);

  if (0) { //print the sorted file list
      cout << "sorted list contains:" << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }
  
  // get necessary info from 1 file
  fullname = dirn + "/" + (*files.begin());
  fp = fopen(fullname.c_str(), "r");
  if (fp==NULL) {
        std::ostringstream es;
        es << "file << " << fullname << " could not be opened.";
        throw invalid_argument(es.str());
  }
  ier = fread(&freq,       sizeof(double), 1, fp);
  ier = fread(&nmods,      sizeof(int),    1, fp);		
  ier = fread(&Nfreq,      sizeof(int),    1, fp);
  ier = fread(&df,         sizeof(double), 1, fp);
  ier = fread(&Nz_grid,    sizeof(int),    1, fp);
  ier = fread(&dz,         sizeof(double), 1, fp);  
  ier = fread(&Nz_subgrid, sizeof(int),    1, fp);
  ier = fread(&delZ,       sizeof(double), 1, fp);
  ier = fread(&z_min,      sizeof(double), 1, fp);
  fclose(fp);

  f_max = Nfreq*df; // fmax
  Nfiles = files.size();
  
  cout << "Nfiles    = " << Nfiles << endl;
  cout << "Nfreq     = " << Nfreq  << endl;
  cout << "max freq  = " << f_max  << endl;
  cout << "delta_Z   = " << delZ   << endl;
  cout << "freq step = " << df     << endl;

  if (Nfreq != Nfiles) {
      throw invalid_argument("Error: number of frequencies is not equal to the number of files" );
  }

  // get the source spectrum with function get_source_spectrum
  // for this we need some intermediary arrays
  S         = new complex<double> [Nfreq];
  pulse_vec = new complex<double> [FFTN];
  arg_vec   = new complex<double> [FFTN];
  
  f_vec     = new double [Nfreq];
  for(i=0; i<Nfreq; i++) {
      f_vec[i] = (i+1)*df;
  }
  
  get_source_spectrum(Nfreq, df, f_vec, f_center, \
								      S, pulse_vec, arg_vec, src_flg,  srcfile);
	
	// only the source spectrum S is necessary from now on; delete the other arrays							      
	delete [] pulse_vec;
	delete [] arg_vec;
	delete [] f_vec;						      

  // the subgrid (window) parameters
  r1   = R_start_km*1000;
  r2   = r1 + width_km*1000;
  dr   = delZ;          // make horizontal spacing equal to the vertical spacing
  Nx   = (int)((r2-r1)/dr)+1;
  Nz   = (int) (height_km*1000.0/delZ)+1;
  tmin = r1/c_ref;
  tmax = r2/c_ref;

  // check if enough memory is available; if not, will iterate more than once
  long maxel = 4000000; // maximum number of elements in the window (matrix)
  zs = Nz;              // zs = how many z-levels to be processed at a time; 
  NN = 0;               // NN is number of iterations-1  (see loop below)
  if (Nx*zs>maxel) {    // zs is adjusted down when not enough memory is available 
      NN  =(Nx*zs)/maxel + 1;  // integer division (floor)
      zs = zs/NN;              // integer division    
  }

  pp = dmatrix(zs,Nx); // allocate 2D pressure field
  tt = tmin;
  //cout <<"Computing the 2D pressure field..." << endl;
  for (j=0; j<NN+1; j++) { // loop to get pressure spectrum at points on the jth subgrid of size [zs, Nx]
      I1   = j*zs;         // zero-based index to pass to getPz1z2()
      I2   = min(Nz-1, I1+zs-1);
      lenI = I2-I1+1;
      cout << "iteration = " << j+1 << " of " << NN+1 << endl;
      
      // get pressure field from height z(I1) to z(I2)
      // Pfzr = getPz1z2 (I1, I2, r1, r2, dr, dirn, files); //Pfzr = c3Darray(Nfiles, Nrng_steps, lenI);
      Pfzr = getPz1z2 (I1, I2, r1, r2, dr, dirn, files, zero_attn_flg); //Pfzr = c3Darray(Nfiles, Nrng_steps, lenI);
			
      for (k=0; k < ntsteps; k++) { // iterate over time steps
          tt = tmin + k*tmstep;
          sprintf(filename, "%s_%d.bin", framefn.c_str(), (int)floor(tt)); // filename

          if (j==0)	{ // info in the file header - appears only once
              fp = fopen(filename, "a+b");			
              fwrite(&Nz,   1, sizeof(int),    fp);
              fwrite(&Nx,   1, sizeof(int),    fp);
              fwrite(&delZ, 1, sizeof(double), fp);
              fwrite(&r1,   1, sizeof(double), fp);							
              fwrite(&tt,   1, sizeof(double), fp);
              fclose(fp);
          }
          for (n=0; n<Nx; n++) {       // for ranges on the x-grid
              for (m=0; m<lenI; m++) { // for altitudes on z-grid							
                  fr       = 0.0;
                  sumP     = 0.0;
                  pp[m][n] = 0.0;
                  for (i=0; i<Nfreq; i++) { //integral over all positive frequencies to get pressure at (x,y,tt)
		                  fr   = fr + df;
		                  sumP = sumP + real(2.0*df*S[i]*Pfzr[i][n][m]*exp(-II*2.0*Pi*fr*tt));						
                  }
                  pp[m][n] = sumP;
              }
          }
          cout << " ... saving 2D field at time " << tt << " secs to file " << filename << endl;
          saveMatrix2bin(filename, pp, lenI, Nx);
      }
      free_c3Darray(Pfzr, Nfiles, Nx);		
  }

  delete [] S;
  free_dmatrix(pp, zs, Nx);
  return 0;
}

/*
int process2DPressure(double R_start_km, double width_km, double height_km, \
                      double c_ref, double tmstep, int ntsteps, double f_center, \
                      string framefn, string dirn, int src_flg, string srcfile) 
{
  // compute the 2D spatial pressure field in (width_km by height_km) window 
  // starting at R_start_km
  // cout << "in process2DPressure():\n";
  
  char   filename [256]; // holds the frame file name(s)
  int    i, j, ier, m, n, k, Nx, Nz, NN, zs;
  int    I1, I2, lenI;
  int    Nfreq, Nz_grid, Nz_subgrid, nmods, Nfiles;
  double tmin, tmax, tt;
  double freq, fr, df, f_max;	
  double dz, delZ, z_min;
  double r1, r2, dr, sumP;
  double *f_vec;
  double **pp;

  complex<double> *S, *pulse_vec, *arg_vec; //source spectrum
  complex<double> II = complex<double> (0.0, 1.0);	
  complex<double> ***Pfzr;
  string fullname;
  FILE *fp;

  list<string> files;
  list<string>::iterator it;
  string pattern (".bin");

  // get and sort the files (they have the frequency in the name)
  getFile_list(dirn, files, pattern);
  files.sort(compare_freq);

  if (0) { //print the sorted file list
      cout << "sorted list contains:" << endl;
      for (it=files.begin(); it!=files.end(); ++it) {
          cout << *it << endl;
      }
      cout << endl;
  }
  
  // get necessary info from 1 file
  fullname = dirn + "/" + (*files.begin());
  fp = fopen(fullname.c_str(), "r");
  if (fp==NULL) {
        std::ostringstream es;
        es << "file << " << fullname << " could not be opened.";
        throw invalid_argument(es.str());
  }
  ier = fread(&freq,       sizeof(double), 1, fp);
  ier = fread(&nmods,      sizeof(int),    1, fp);		
  ier = fread(&Nfreq,      sizeof(int),    1, fp);
  ier = fread(&df,         sizeof(double), 1, fp);
  ier = fread(&Nz_grid,    sizeof(int),    1, fp);
  ier = fread(&dz,         sizeof(double), 1, fp);  
  ier = fread(&Nz_subgrid, sizeof(int),    1, fp);
  ier = fread(&delZ,       sizeof(double), 1, fp);
  ier = fread(&z_min,      sizeof(double), 1, fp);
  fclose(fp);

  f_max = Nfreq*df; // fmax
  Nfiles = files.size();
  
  cout << "Nfiles    = " << Nfiles << endl;
  cout << "Nfreq     = " << Nfreq  << endl;
  cout << "max freq  = " << f_max  << endl;
  cout << "delta_Z   = " << delZ   << endl;
  cout << "freq step = " << df     << endl;

  if (Nfreq != Nfiles) {
      throw invalid_argument("Error: number of frequencies is not equal to the number of files" );
  }

  // get the source spectrum with function get_source_spectrum
  // for this we need some intermediary arrays
  S         = new complex<double> [Nfreq];
  pulse_vec = new complex<double> [FFTN];
  arg_vec   = new complex<double> [FFTN];
  
  f_vec     = new double [Nfreq];
  for(i=0; i<Nfreq; i++) {
      f_vec[i] = (i+1)*df;
  }
  
  get_source_spectrum(Nfreq, df, f_vec, f_center, \
								      S, pulse_vec, arg_vec, src_flg,  srcfile);
	
	// only the source spectrum S is necessary from now on; delete the other arrays							      
	delete [] pulse_vec;
	delete [] arg_vec;
	delete [] f_vec;						      

  // the subgrid (window) parameters
  r1   = R_start_km*1000;
  r2   = r1 + width_km*1000;
  dr   = delZ;          // make horizontal spacing equal to the vertical spacing
  Nx   = (int)((r2-r1)/dr)+1;
  Nz   = (int) (height_km*1000.0/delZ)+1;
  tmin = r1/c_ref;
  tmax = r2/c_ref;

  // check if enough memory is available; if not, will iterate more than once
  long maxel = 4000000; // maximum number of elements in the window (matrix)
  zs = Nz;              // zs = how many z-levels to be processed at a time; 
  NN = 0;               // NN is number of iterations-1  (see loop below)
  if (Nx*zs>maxel) {    // zs is adjusted down when not enough memory is available 
      NN  =(Nx*zs)/maxel + 1;  // integer division (floor)
      zs = zs/NN;              // integer division    
  }

  pp = dmatrix(zs,Nx); // allocate 2D pressure field
  tt = tmin;
  //cout <<"Computing the 2D pressure field..." << endl;
  for (j=0; j<NN+1; j++) { // loop to get pressure spectrum at points on the jth subgrid of size [zs, Nx]
      I1   = j*zs;         // zero-based index to pass to getPz1z2()
      I2   = min(Nz-1, I1+zs-1);
      lenI = I2-I1+1;
      cout << "iteration = " << j+1 << " of " << NN+1 << endl;
      
      // get pressure field from height z(I1) to z(I2)
      Pfzr = getPz1z2 (I1, I2, r1, r2, dr, dirn, files); //Pfzr = c3Darray(Nfiles, Nrng_steps, lenI);
			
      for (k=0; k < ntsteps; k++) { // iterate over time steps
          tt = tmin + k*tmstep;
          sprintf(filename, "%s_%d.bin", framefn.c_str(), (int)floor(tt)); // filename

          if (j==0)	{ // info in the file header - appears only once
              fp = fopen(filename, "a+b");			
              fwrite(&Nz,   1, sizeof(int),    fp);
              fwrite(&Nx,   1, sizeof(int),    fp);
              fwrite(&delZ, 1, sizeof(double), fp);
              fwrite(&r1,   1, sizeof(double), fp);							
              fwrite(&tt,   1, sizeof(double), fp);
              fclose(fp);
          }
          for (n=0; n<Nx; n++) {       // for ranges on the x-grid
              for (m=0; m<lenI; m++) { // for altitudes on z-grid							
                  fr       = 0.0;
                  sumP     = 0.0;
                  pp[m][n] = 0.0;
                  for (i=0; i<Nfreq; i++) { //integral over all positive frequencies to get pressure at (x,y,tt)
		                  fr   = fr + df;
		                  sumP = sumP + real(2.0*df*S[i]*Pfzr[i][n][m]*exp(-II*2.0*Pi*fr*tt));						
                  }
                  pp[m][n] = sumP;
              }
          }
          cout << " ... saving 2D field at time " << tt << " secs to file " << filename << endl;
          saveMatrix2bin(filename, pp, lenI, Nx);
      }
      free_c3Darray(Pfzr, Nfiles, Nx);		
  }

  delete [] S;
  free_dmatrix(pp, zs, Nx);
  return 0;
}
*/

int saveMatrix2bin(const char *filename, double **M, int nr, int nc) {
  int i,j, szd;
  szd = sizeof(double);
  FILE *fp = fopen(filename, "a+b");
  for (i=0; i<nr; i++) {
      for (j=0; j<nc; j++) {
          fwrite(&M[i][j], szd, 1, fp);
      }
  }
  fclose(fp);
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
  
  FILE *fp = fopen("atm_profile.bbnm", "w");
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
  printf("Interpolated atmospheric profiles (z,u,v,w,t,d,p,c,c_eff) saved in: atm_profile.bbnm\n");
  return 0;
}


int load_source_pulse_td(string srcpulsetdfn, vector<double> &t, vector<double> &tdp ) {
  // load time domain source pulse: time (s) | amplitude
  cout << "Loading time domain source pulse from file " << srcpulsetdfn << endl;

  double d1, d2; //, d3;
  ifstream indata;

  indata.open(srcpulsetdfn.c_str());	
  if (!indata) {
      cerr << "file " << srcpulsetdfn << " could not be opened.\n";
      exit(1);
  }

  //indata >> d1 >> d2 >> d3; // read first line in file
  indata >> d1 >> d2; // read first line in file
  while (!indata.eof()) {
      t.push_back(d1);
      tdp.push_back(d2);
      //indata >> d1 >> d2 >> d3;
      indata >> d1 >> d2;
  }
  indata.close();

  //// print data
  //for (unsigned i=0; i<t.size(); i++) {
  //    cout << i << "  " << t[i] << endl;
  //}

  return 0;
}


int load_source_spectrum(string srcspfn, double *freqv, complex<double> *dft_vec, int fftn) {
  // load source spectrum
  //cout << "Loading source spectrum from file " << srcspfn << endl;
  
  int i;
  double d1, d2, d3;
  complex<double> I (0.0, 1.0);
  ifstream indata;

  indata.open(srcspfn.c_str());	
  if (!indata) {
      cerr << "file " << srcspfn << " could not be opened.\n";
      exit(1);
  }

  i = 0;
  indata >> d1 >> d2 >> d3; // read first line in file
  while (!indata.eof()) {
      freqv[i] = d1;
      dft_vec[i] = d2 + I*d3;
      indata >> d1 >> d2 >> d3;
      if (i>(fftn-1)) {
          cerr << "file " << srcspfn << "appears to have more points than maximum allowed of " << endl
               << "fftn = " << fftn << endl
               << " ... aborting."  << endl;
          exit(1);
      }
      i++;
  }
  indata.close();

  return 0;
}	  





