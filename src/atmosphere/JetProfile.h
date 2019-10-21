#ifndef _NCPA_JETPROFILE_H__
#define _NCPA_JETPROFILE_H__

#include "AtmosphericProfile.h"
#include "util.h"
#include <iostream>
#include <vector>

namespace NCPA {
	
	struct _jetStruct {
		double azimuth, speed, width, height;
	};

	/**
	  * Atmospheric profile consisting of three wind jets atop the thermodynamic sound speed.
	  */
	class JetProfile : public NCPA::AtmosphericProfile {

		private:
			const static double A1, A2, A3, A4, A5, A6, A7, A8;
			const static double B1, B2, B3, B4, B5, B6, B7, B8;
			/*
			const static double A1 =     -3.9082017e-2;
			const static double A2 =     -1.1526465e-3;
			const static double A3 =     3.2891937e-5;
			const static double A4 =     -2.0494958e-7;
			const static double A5 =     -4.7087295e-2;
			const static double A6 =     1.2506387e-3;
			const static double A7 =     -1.5194498e-5;
			const static double A8 =     6.518877e-8;

			const static double B1 =     -4.9244637e-3;
			const static double B2 =     -1.2984142e-6;
			const static double B3 =     -1.5701595e-6;
			const static double B4 =     1.5535974e-8;
			const static double B5 =     -2.7221769e-2;
			const static double B6 =     4.247473e-4;
			const static double B7 =     -3.958318e-6;
			const static double B8 =     1.7295795e-8;
			*/

			/*
			// Prototyping and testing variables, currently unused
			double  v_strat, v_jet, v_noct,
				z_strat, z_jet, z_noct,
				strat_jet_width, jet_width, noct_width;
			double azimuth, jet_azimuth, noct_azimuth;
			*/
			
			std::vector< _jetStruct > _jets;
			
			double T0, rho0;

		public:

			/**
			  * Default constructor.  Initializes the profile values with defaults suitable
			  * for testing.
			  */
			JetProfile();

			/**
			  * Initializes profile with values from param file.  Not yet implemented.
			  * @param filename The parameter file to use for values.
			  * @todo Implement this method.
			  */
			JetProfile( std::string filename );

			double t( double z );
                        double u( double z );
                        double v( double z );
                        double w( double z );
                        double ceff( double z, double phi );
			// @CHH: Use parent object's method
			//double c0( double z );    
			double p( double z );
                        double rho( double z );
			double z0();
			
			void addJet( double height, double speed, double width, double azimuth );


	};  // class JetProfile

}  // namespace NCPA


#endif  // #ifndef _NCPA_JETPROFILE_H__
