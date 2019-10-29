#include "atmlib.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "binaryreader.h"

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace std;
typedef unsigned char byte;

  AtmLibrary::AtmLibrary() { 
}

void AtmLibrary::readG2SBinary( ifstream* g2s_env, double *alt, double *range, double **zw_2D, double **mw_2D, double **T_2D, double **rho_2D, double **pr_2D)
 {
  g2s_env->seekg( 0 , ios_base::beg );             // go back to beginning of file
  NCPA::BinaryReader *binread = new NCPA::BinaryReader();
  int *filevars = new int[ 2 ];

  // get # profiles and # layers per profile from .env file
  binread->readLittleIntArray( g2s_env, 2, filevars );
  int atm_r = filevars[ 0 ];
  int atm_z = filevars[ 1 ];
  delete [] filevars;

  // get physical information per profile. this gives us range and altitude arrays
  double *lats     = new double[ atm_r ];
  double *lons     = new double[ atm_r ];
  double *bazs     = new double[ atm_r ];
  double *ground_m = new double[ atm_r ];
  binread->readLittleDoubleArray( g2s_env, atm_r, lats );
  binread->readLittleDoubleArray( g2s_env, atm_r, lons );
  binread->readLittleDoubleArray( g2s_env, atm_r, bazs );
  binread->readLittleDoubleArray( g2s_env, atm_r, range );
  binread->readLittleDoubleArray( g2s_env, atm_r, ground_m );
  binread->readLittleDoubleArray( g2s_env, atm_z, alt );
  delete [] lats;
  delete [] lons;
  delete [] bazs;
  delete [] ground_m;

  // get atmospheric per array. for now, just store c_eff values
  double *temp_K        = new double[ atm_z ];
  double *density_gpcm3 = new double[ atm_z ];
  double *pressure_hPa  = new double[ atm_z ];
  double *wind_zw_mps   = new double[ atm_z ];
  double *wind_mw_mps   = new double[ atm_z ];
  double *wind_vert_mps = new double[ atm_z ];
  int bytesToSkipFromBeg = 8 * 5 * atm_r + 8 * atm_z + 8;
  for (int i = 0; i < atm_r; i++) {
   int bytesToSkipBefore  = 8 * i * atm_z;
   int bytesToSkipAfter   = 8 * atm_z * (atm_r - i - 1);

   // for subsequent iterations, skip to correct beginning
   if (i > 0) {
     g2s_env->seekg( bytesToSkipFromBeg, ios_base::beg );
    }

   g2s_env->seekg( bytesToSkipBefore, ios_base::cur );
   binread->readLittleDoubleArray( g2s_env, atm_z, temp_K );
   g2s_env->seekg( bytesToSkipAfter, ios_base::cur );
    
   g2s_env->seekg( bytesToSkipBefore, ios_base::cur );
   binread->readLittleDoubleArray( g2s_env, atm_z, density_gpcm3 );
   g2s_env->seekg( bytesToSkipAfter, ios_base::cur );

   g2s_env->seekg( bytesToSkipBefore, ios_base::cur );
   binread->readLittleDoubleArray( g2s_env, atm_z, pressure_hPa );
   g2s_env->seekg( bytesToSkipAfter, ios_base::cur );

   g2s_env->seekg( bytesToSkipBefore, ios_base::cur );
   binread->readLittleDoubleArray( g2s_env, atm_z, wind_zw_mps );
   g2s_env->seekg( bytesToSkipAfter, ios_base::cur );

   g2s_env->seekg( bytesToSkipBefore, ios_base::cur );
   binread->readLittleDoubleArray( g2s_env, atm_z, wind_mw_mps );
   g2s_env->seekg( bytesToSkipAfter, ios_base::cur );

   g2s_env->seekg( bytesToSkipBefore, ios_base::cur );
   binread->readLittleDoubleArray( g2s_env, atm_z, wind_vert_mps );
   g2s_env->seekg( bytesToSkipAfter, ios_base::cur );

   for(int j=0; j < atm_z; j++) {
    if (i == atm_r-1) {
     alt[j]    = alt[j]*1000;
    }
    zw_2D[j][i]  = wind_zw_mps[j];
    mw_2D[j][i]  = wind_mw_mps[j];
    T_2D[j][i]   = temp_K[j];
    rho_2D[j][i] = density_gpcm3[j]*1000;     
    pr_2D[j][i]  = pressure_hPa[j]*100;     
   }
   range[i] = range[i]*1000;
  }

  delete [] temp_K;
  delete [] density_gpcm3;
  delete [] pressure_hPa;
  delete [] wind_zw_mps;
  delete [] wind_mw_mps;
  delete [] wind_vert_mps;
  
  delete binread;
 }

  void AtmLibrary::readG2SAscii( ifstream* infile, double *alt, double *zw_, double *mw_, double *T_, double *rho_, double *pr_)
{
   int i = 0;
   float dummy;
   while(*infile && (!infile->eof()) ) {
     *infile >> alt[i] >> zw_[i] >> mw_[i] >> dummy >> T_[i] >> rho_[i] >> pr_[i];
    alt[i]  = alt[i]*1000;
    rho_[i] = rho_[i]*1000;
    pr_[i]  = pr_[i]*100;
    i++;
   }
}

  void AtmLibrary::getNumberOfLines( ifstream* infile, int *atm_rows)
{
   int i = -1;
   float dummy;
   while(*infile && (!infile->eof()) ) {
     *infile >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
    i++;
   }
 *atm_rows = i;
 infile->clear();
 infile->seekg(0, ios::beg);
}

void AtmLibrary::getBinaryG2SDimensions( ifstream* g2s_env, int *atm_nz, int *atm_nr)
 {
  g2s_env->seekg( 0 , ios_base::beg );             // go back to beginning of file
  NCPA::BinaryReader *binread = new NCPA::BinaryReader();
  int *filevars = new int[ 2 ];

  // get # profiles and # layers per profile from .env file
  binread->readLittleIntArray( g2s_env, 2, filevars );
  *atm_nr = filevars[ 0 ];
  *atm_nz = filevars[ 1 ];
  delete [] filevars;
  delete binread;
 }

 void AtmLibrary::readEOFs( ifstream* eoffile, int n_eof, int n_levels,  double **eof_data)
  {
   for (int i=0; i < n_levels; i++) 
    {
     for (int j=0; j < n_eof+1; j++)   // #n_eof eof's + 1 average field
      {
       *eoffile >> eof_data[i][j];
      }
    }
  } 

 void AtmLibrary::makeProfile( int n_levels, int n_eof, double *beta, double **eof_data, double *profile)
  {
   for (int i=0; i < n_levels; i++) 
    {
     profile[i] = eof_data[i][0];
     for (int j=1; j < n_eof+1; j++) 
      {
       profile[i] = profile[i] + beta[j-1]*eof_data[i][j];
      }
    }
  }

 void AtmLibrary::makeZeroProfile( int n_levels, double *profile)
  {
   for (int i=0; i < n_levels; i++) 
    {
     profile[i] = 0.0;
    }
  }

 void AtmLibrary::makeGaussianProfile( int n_levels, double aj, double hj, double wj, double *alt, double *profile)
  {
   double gaussarg;
 
   for (int i=0; i < n_levels; i++) 
    {
     gaussarg   = -(alt[i]-hj)*(alt[i]-hj)/(2*wj*wj);
     profile[i] = aj*exp(gaussarg);
    }
  }

/*
 void AtmLibrary::writeProfile( char *fid, int n_levels, double *alt, double *u, double *v, double *T, double *rho, double *P )
  {
   FILE *profile = fopen(fid,"w");
   for (int i=0; i< n_levels; i++)
    {
     fprintf(profile,"%9.3f %8.3f %8.3f %8.3f %10.3f %11.4e %11.4e\n", alt[i]/1000, u[i], v[i], 0.0, T[i], rho[i]/1000, P[i]/100); 
    } 
   fclose(profile);
  }
*/


// DV 20130417: added zmin argument to accomodate profiles referred to Mean Sea Level
void AtmLibrary::writeProfile( char *fid, int n_levels, double zmin, double *alt, double *u, double *v, double *T, double *rho, double *P )  {
   FILE *profile = fopen(fid,"w");
   for (int i=0; i< n_levels; i++)
    {
     fprintf(profile,"%9.3f %8.3f %8.3f %8.3f %10.3f %11.4e %11.4e\n", (alt[i]+zmin)/1000, u[i], v[i], 0.0, T[i], rho[i]/1000, P[i]/100); 
    } 
   fclose(profile);
}


/*
  void AtmLibrary::getAbsorptionCoefficients(int n, double freq, double *alt, double *T, double *pr, double *c, double *alpha)
 {
  // Expressions based on Bass and Sutherland, JASA 2004
  // Computes alpha(freq) for given G2S output
  // Subroutine can easily be modified to include dispersion effects
  // In that case, make alpha a complex array and
  // use the real part as absorption and the imaginary part for dispersion
  int m, ii;
  double T_o, P_o, S, z;
  double X[7], X_ON, Z_rot[2], Z_rot_;
  double sigma, nn, chi, cchi, mu, nu, mu_o;
  double beta_0, beta_1, beta_2, alpha_1, alpha_2;
  double a_cl, a_rot, a_diff, a_vib;
  double T_z, P_z, c_snd_z, gamma;
  double A1, A2, B, C, D, E, F, G, H, I, J, K, L, ZZ, hu;
  double f_vib[4], a_vib_c[4], Cp_R[4], Cv_R[4], theta[4], C_R, A_max, Tr;
  
  double tweak_abs = 1; //0.3; // tweak absorption alpha by this factor

  // Atmospheric composition constants
  mu_o  = 18.192E-6;       // Reference viscosity [kg/(m*s)]
  T_o   = T[0];            // Reference temperature [K]
  P_o   = pr[0];           // Reference pressure [Pa]
  S     = 117;             // Sutherland constant [K]        

  Cv_R[0] = 5/2;                                    // Heat capacity|volume (O2)
  Cv_R[1] = 5/2;                                    // Heat capacity|volume (N2)
  Cv_R[2] = 3;                                      // Heat capacity|volume (CO2)
  Cv_R[3] = 3;                                      // Heat capacity|volume (O3)
  Cp_R[0] = 7/2;                                    // Heat capacity|pressure (O2)
  Cp_R[1] = 7/2;                                    // Heat capacity|pressure (N2)
  Cp_R[2] = 4;                                      // Heat capacity|pressure (CO2)
  Cp_R[3] = 4;                                      // Heat capacity|pressure (O3)
  theta[0]= 2239.1;                                 // Charact. temperature (O2)
  theta[1]= 3352;                                   // Charact. temperature (N2)
  theta[2]= 915;                                    // Charact. temperature (CO2)
  theta[3]= 1037;                                   // Charact. temperature (O3)

  gamma   = 1.4;
  for (ii=0; ii<n; ii++) {
   z       = alt[ii]/1000; // in km
   T_z     = T[ii];
   P_z     = pr[ii];
   c_snd_z = c[ii]; 
   mu      = mu_o*sqrt(T_z/T_o)*((1+S/T_o)/(1+S/T_z)); // Viscosity [kg/(m*s)]
   nu      = (8*PI*freq*mu)/(3*P_z);                   // Nondimensional frequency
   
  //-------- Gas fraction polynomial fits -----------------------------------
  if (z > 90.)                                         // O2 profile
   X[0] = pow(10,49.296-(1.5524*z)+(1.8714E-2*pow(z,2))-(1.1069E-4*pow(z,3))+(3.199E-7*pow(z,4))-(3.6211E-10*pow(z,5)));
  else
   X[0] = pow(10,-0.67887);

  if (z > 76.)                                         // N2 profile
   X[1] = pow(10,(1.3972E-1)-(5.6269E-3*z)+(3.9407E-5*pow(z,2))-(1.0737E-7*pow(z,3)));
  else
   X[1] = pow(10,-0.10744);

  X[2]  = pow(10,-3.3979);                             // CO2 profile

  if (z > 80. )                                        // O3 profile
   X[3] = pow(10,-4.234-(3.0975E-2*z));
  else
   X[3] = pow(10,-19.027+(1.3093*z)-(4.6496E-2*pow(z,2))+(7.8543E-4*pow(z,3))-(6.5169E-6*pow(z,4))+(2.1343E-8*pow(z,5)));

  if (z > 95. )                                        // O profile
   X[4] = pow(10,-3.2456+(4.6642E-2*z)-(2.6894E-4*pow(z,2))+(5.264E-7*pow(z,3)));
  else
   X[4] = pow(10,-11.195+(1.5408E-1*z)-(1.4348E-3*pow(z,2))+(1.0166E-5*pow(z,3)));

                                                       // N profile 
  X[5]  = pow(10,-53.746+(1.5439*z)-(1.8824E-2*pow(z,2))+(1.1587E-4*pow(z,3))-(3.5399E-7*pow(z,4))+(4.2609E-10*pow(z,5)));

  if (z > 30. )                                         // H2O profile
   X[6] = pow(10,-4.2563+(7.6245E-2*z)-(2.1824E-3*pow(z,2))-(2.3010E-6*pow(z,3))+(2.4265E-7*pow(z,4))-(1.2500E-09*pow(z,5)));
  else 
   {
    if (z > 100.)
     X[6] = pow(10,-0.62534-(8.3665E-2*z));
    else
     X[6] = pow(10,-1.7491+(4.4986E-2*z)-(6.8549E-2*pow(z,2))+(5.4639E-3*pow(z,3))-(1.5539E-4*pow(z,4))+(1.5063E-06*pow(z,5)));
   }
  X_ON = (X[0] + X[1])/0.9903;
  
  //-------- Rotational collision number-------------------------------------
  Z_rot[0] = 54.1*exp(-17.3*(pow(T_z,-1./3.)));   // O2
  Z_rot[1] = 63.3*exp(-16.7*(pow(T_z,-1./3.)));   // N2
  Z_rot_   = 1./((X[1]/Z_rot[1])+(X[0]/Z_rot[0]));
  
  //-------- Nondimensional atmospheric quantities---------------------------
  sigma = 5./sqrt(21.);
  nn = (4./5.)*sqrt(3./7.)*Z_rot_;
  chi=3.*nn*nu/4.;
  cchi=2.36*chi;
 
  //---------Classical + rotational loss/dispersion--------------------------
  beta_0  = 2*PI*freq/c_snd_z; 
  beta_1  = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/(1+pow(nu,2)));
  beta_2  = beta_0*sqrt((1+pow(chi,2))/(1+pow((sigma*chi),2))); 
  alpha_1 = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))-1)/(1+pow(nu,2)));
  alpha_2 = beta_0*(((sigma/2-1/(2*sigma))*chi)/(sqrt((1+pow(chi,2))*(1+pow(sigma*chi,2)))));
  //a_cl    = alpha_1*(beta_2/beta_0);
  //a_rot = alpha_2*(beta_1/beta_0)*X_ON;

  a_cl    = (2*PI*freq/c_snd_z)*sqrt(0.5*(sqrt(1+pow(nu,2))-1)*(1+pow(cchi,2))/((1+pow(nu,2))*(1+pow(sigma*cchi,2))));
  a_rot   = (2*PI*freq/c_snd_z)*X_ON*((pow(sigma,2)-1)*chi/(2*sigma))*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/((1+pow(nu,2))*(1+pow(cchi,2))));
  a_diff  = 0.003*a_cl;

  //---------Vibrational relaxation-------------------------------------------
  Tr = pow(T_z/T_o,-1./3.)-1;
  A1 = (X[0]+X[1])*24*exp(-9.16*Tr);
  A2 = (X[4]+X[5])*2400;
  B  = 40400*exp(10*Tr);
  C  = 0.02*exp(-11.2*Tr);
  D  = 0.391*exp(8.41*Tr);
  E  = 9*exp(-19.9*Tr);
  F  = 60000;
  G  = 28000*exp(-4.17*Tr);
  H  = 22000*exp(-7.68*Tr);
  I  = 15100*exp(-10.4*Tr);
  J  = 11500*exp(-9.17*Tr);
  K  = (8.48E08)*exp(9.17*Tr);
  L  = exp(-7.72*Tr);
  ZZ = H*X[2]+I*(X[0]+0.5*X[4])+J*(X[1]+0.5*X[5])+K*(X[6]+X[3]);
  hu = 100*(X[3]+X[6]);
  f_vib[0] = (P_z/P_o)*(mu_o/mu)*(A1+A2+B*hu*(C+hu)*(D+hu));
  f_vib[1] = (P_z/P_o)*(mu_o/mu)*(E+F*X[3]+G*X[6]);
  f_vib[2] = (P_z/P_o)*(mu_o/mu)*ZZ;
  f_vib[3] = (P_z/P_o)*(mu_o/mu)*(1.2E5)*L;

  a_vib = 0.;
  for (m=0; m<4; m++)
   {
    C_R          = ((pow(theta[m]/T_z,2))*exp(-theta[m]/T_z))/(pow(1-exp(-theta[m]/T_z),2));
    A_max        = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
    a_vib_c[m]   = (A_max/c_snd_z)*((2*(pow(freq,2))/f_vib[m])/(1+pow(freq/f_vib[m],2)));
    a_vib        = a_vib + a_vib_c[m];
   }

   alpha[ii] = a_cl + a_rot + a_diff + a_vib;
   
   // alter alpha for some studies
	 alpha[ii] = tweak_abs*alpha[ii];
  }
  printf(" -> Using Sutherland-Bass absorption times a factor of %g\n", tweak_abs);
 }
 */
 

// updated: will accept attenuation coeff. loaded from a file
 void AtmLibrary::getAbsorptionCoefficients(int n, double freq, double *alt, double *T, double *pr, double *c, string usrattfile, double *alpha)
{

  if (usrattfile.empty()) {
  // Expressions based on Bass and Sutherland, JASA 2004
  // Computes alpha(freq) for given G2S output
  // Subroutine can easily be modified to include dispersion effects
  // In that case, make alpha a complex array and
  // use the real part as absorption and the imaginary part for dispersion
  int m, ii;
  double T_o, P_o, S, z;
  double X[7], X_ON, Z_rot[2], Z_rot_;
  double sigma, nn, chi, cchi, mu, nu, mu_o;
  double beta_0, beta_1, beta_2, alpha_1, alpha_2;
  double a_cl, a_rot, a_diff, a_vib;
  double T_z, P_z, c_snd_z, gamma;
  double A1, A2, B, C, D, E, F, G, H, I, J, K, L, ZZ, hu;
  double f_vib[4], a_vib_c[4], Cp_R[4], Cv_R[4], theta[4], C_R, A_max, Tr;
  
  double tweak_abs = 1; //0.3; // tweak absorption alpha by this factor

  // Atmospheric composition constants
  mu_o  = 18.192E-6;       // Reference viscosity [kg/(m*s)]
  T_o   = T[0];            // Reference temperature [K]
  P_o   = pr[0];           // Reference pressure [Pa]
  S     = 117;             // Sutherland constant [K]        

  Cv_R[0] = 5/2;                                    // Heat capacity|volume (O2)
  Cv_R[1] = 5/2;                                    // Heat capacity|volume (N2)
  Cv_R[2] = 3;                                      // Heat capacity|volume (CO2)
  Cv_R[3] = 3;                                      // Heat capacity|volume (O3)
  Cp_R[0] = 7/2;                                    // Heat capacity|pressure (O2)
  Cp_R[1] = 7/2;                                    // Heat capacity|pressure (N2)
  Cp_R[2] = 4;                                      // Heat capacity|pressure (CO2)
  Cp_R[3] = 4;                                      // Heat capacity|pressure (O3)
  theta[0]= 2239.1;                                 // Charact. temperature (O2)
  theta[1]= 3352;                                   // Charact. temperature (N2)
  theta[2]= 915;                                    // Charact. temperature (CO2)
  theta[3]= 1037;                                   // Charact. temperature (O3)

  gamma   = 1.4;
  for (ii=0; ii<n; ii++) {
   z       = alt[ii]/1000; // in km
   T_z     = T[ii];
   P_z     = pr[ii];
   c_snd_z = c[ii]; 
   mu      = mu_o*sqrt(T_z/T_o)*((1+S/T_o)/(1+S/T_z)); // Viscosity [kg/(m*s)]
   nu      = (8*PI*freq*mu)/(3*P_z);                   // Nondimensional frequency
   
  //-------- Gas fraction polynomial fits -----------------------------------
  if (z > 90.)                                         // O2 profile
   X[0] = pow(10,49.296-(1.5524*z)+(1.8714E-2*pow(z,2))-(1.1069E-4*pow(z,3))+(3.199E-7*pow(z,4))-(3.6211E-10*pow(z,5)));
  else
   X[0] = pow(10,-0.67887);

  if (z > 76.)                                         // N2 profile
   X[1] = pow(10,(1.3972E-1)-(5.6269E-3*z)+(3.9407E-5*pow(z,2))-(1.0737E-7*pow(z,3)));
  else
   X[1] = pow(10,-0.10744);

  X[2]  = pow(10,-3.3979);                             // CO2 profile

  if (z > 80. )                                        // O3 profile
   X[3] = pow(10,-4.234-(3.0975E-2*z));
  else
   X[3] = pow(10,-19.027+(1.3093*z)-(4.6496E-2*pow(z,2))+(7.8543E-4*pow(z,3))-(6.5169E-6*pow(z,4))+(2.1343E-8*pow(z,5)));

  if (z > 95. )                                        // O profile
   X[4] = pow(10,-3.2456+(4.6642E-2*z)-(2.6894E-4*pow(z,2))+(5.264E-7*pow(z,3)));
  else
   X[4] = pow(10,-11.195+(1.5408E-1*z)-(1.4348E-3*pow(z,2))+(1.0166E-5*pow(z,3)));

                                                       // N profile 
  X[5]  = pow(10,-53.746+(1.5439*z)-(1.8824E-2*pow(z,2))+(1.1587E-4*pow(z,3))-(3.5399E-7*pow(z,4))+(4.2609E-10*pow(z,5)));

  if (z > 30. )                                         // H2O profile
   X[6] = pow(10,-4.2563+(7.6245E-2*z)-(2.1824E-3*pow(z,2))-(2.3010E-6*pow(z,3))+(2.4265E-7*pow(z,4))-(1.2500E-09*pow(z,5)));
  else 
   {
    if (z > 100.)
     X[6] = pow(10,-0.62534-(8.3665E-2*z));
    else
     X[6] = pow(10,-1.7491+(4.4986E-2*z)-(6.8549E-2*pow(z,2))+(5.4639E-3*pow(z,3))-(1.5539E-4*pow(z,4))+(1.5063E-06*pow(z,5)));
   }
  X_ON = (X[0] + X[1])/0.9903;
  
  //-------- Rotational collision number-------------------------------------
  Z_rot[0] = 54.1*exp(-17.3*(pow(T_z,-1./3.)));   // O2
  Z_rot[1] = 63.3*exp(-16.7*(pow(T_z,-1./3.)));   // N2
  Z_rot_   = 1./((X[1]/Z_rot[1])+(X[0]/Z_rot[0]));
  
  //-------- Nondimensional atmospheric quantities---------------------------
  sigma = 5./sqrt(21.);
  nn = (4./5.)*sqrt(3./7.)*Z_rot_;
  chi=3.*nn*nu/4.;
  cchi=2.36*chi;
 
  //---------Classical + rotational loss/dispersion--------------------------
  beta_0  = 2*PI*freq/c_snd_z; 
  beta_1  = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/(1+pow(nu,2)));
  beta_2  = beta_0*sqrt((1+pow(chi,2))/(1+pow((sigma*chi),2))); 
  alpha_1 = beta_0*sqrt(0.5*(sqrt(1+pow(nu,2))-1)/(1+pow(nu,2)));
  alpha_2 = beta_0*(((sigma/2-1/(2*sigma))*chi)/(sqrt((1+pow(chi,2))*(1+pow(sigma*chi,2)))));
  //a_cl    = alpha_1*(beta_2/beta_0);
  //a_rot = alpha_2*(beta_1/beta_0)*X_ON;

  a_cl    = (2*PI*freq/c_snd_z)*sqrt(0.5*(sqrt(1+pow(nu,2))-1)*(1+pow(cchi,2))/((1+pow(nu,2))*(1+pow(sigma*cchi,2))));
  a_rot   = (2*PI*freq/c_snd_z)*X_ON*((pow(sigma,2)-1)*chi/(2*sigma))*sqrt(0.5*(sqrt(1+pow(nu,2))+1)/((1+pow(nu,2))*(1+pow(cchi,2))));
  a_diff  = 0.003*a_cl;

  //---------Vibrational relaxation-------------------------------------------
  Tr = pow(T_z/T_o,-1./3.)-1;
  A1 = (X[0]+X[1])*24*exp(-9.16*Tr);
  A2 = (X[4]+X[5])*2400;
  B  = 40400*exp(10*Tr);
  C  = 0.02*exp(-11.2*Tr);
  D  = 0.391*exp(8.41*Tr);
  E  = 9*exp(-19.9*Tr);
  F  = 60000;
  G  = 28000*exp(-4.17*Tr);
  H  = 22000*exp(-7.68*Tr);
  I  = 15100*exp(-10.4*Tr);
  J  = 11500*exp(-9.17*Tr);
  K  = (8.48E08)*exp(9.17*Tr);
  L  = exp(-7.72*Tr);
  ZZ = H*X[2]+I*(X[0]+0.5*X[4])+J*(X[1]+0.5*X[5])+K*(X[6]+X[3]);
  hu = 100*(X[3]+X[6]);
  f_vib[0] = (P_z/P_o)*(mu_o/mu)*(A1+A2+B*hu*(C+hu)*(D+hu));
  f_vib[1] = (P_z/P_o)*(mu_o/mu)*(E+F*X[3]+G*X[6]);
  f_vib[2] = (P_z/P_o)*(mu_o/mu)*ZZ;
  f_vib[3] = (P_z/P_o)*(mu_o/mu)*(1.2E5)*L;

  a_vib = 0.;
  for (m=0; m<4; m++)
   {
    C_R          = ((pow(theta[m]/T_z,2))*exp(-theta[m]/T_z))/(pow(1-exp(-theta[m]/T_z),2));
    A_max        = (X[m]*(PI/2)*C_R)/(Cp_R[m]*(Cv_R[m]+C_R));
    a_vib_c[m]   = (A_max/c_snd_z)*((2*(pow(freq,2))/f_vib[m])/(1+pow(freq/f_vib[m],2)));
    a_vib        = a_vib + a_vib_c[m];
   }

   alpha[ii] = a_cl + a_rot + a_diff + a_vib;
   
   // alter alpha for some studies
	 alpha[ii] = tweak_abs*alpha[ii];
  }
  printf(" -> Using Sutherland-Bass absorption times a factor of %g\n", tweak_abs);
  
	}
	else {
	
		    // load atten. coeff from a text file with columns: z (km AGL) | attn
	    cout << " -> Loading attenuation coefficients from file " << usrattfile << endl;
	
	    double d1, d2;
	    vector<double> zatt, att;
	    ifstream indata;

	    indata.open(usrattfile.c_str());	
	    if (!indata) {
	        cerr << "file " << usrattfile << " could not be opened.\n";
	        exit(1);
	    }
	
	    indata >> d1 >> d2; // read first line in file
	    while (!indata.eof()) {
	        zatt.push_back(d1*1000.0);
	        att.push_back(d2);
	        indata >> d1 >> d2;
	    }
	    indata.close();
	    
	    //// print data
	    //for (unsigned i=0; i<zatt.size(); i++) {
	    //    cout << i << "  " << zatt[i] << endl;
	    //}

	    // do spline interpolation of attenuation coefficient
	    gsl_interp_accel *acc_;
	    acc_ = gsl_interp_accel_alloc();
	    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, zatt.size());
	    gsl_spline_init(spline, zatt.data(), att.data(), zatt.size());
	    
	    // do spline up to the max height hatt found in the provided atten. file
	    int it = 0;
	    while ((alt[it]< zatt[zatt.size()-1]) && (it<n)) {
	        alpha[it] = gsl_spline_eval(spline, alt[it], acc_ );
	        it++;
	    }
	    // if the grid extends above hatt => extend alpha with the last value up to grid height 
      for (int i=it; i< n; i++) {
          alpha[i] = alpha[it-1];
      }
            
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc_);
  }
     
  // save alpha?
  if (1) {
      FILE *fp = fopen("attn.pe", "w");
      for (int i=0; i<n; i++) {
          fprintf(fp, "%8.3f  %14.6e\n", alt[i], alpha[i]);
      }
      printf(" -> Attenuation coefficients saved in 'attn.pe'\n");
      fclose(fp);
  }
  
 }
  
