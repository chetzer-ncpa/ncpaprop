#include "AtmosphericSpecification.h"
#include "geographic.h"
#include "AtmosphericProfile.h"
#include <vector>

namespace NCPA {
	
	class Slice : public AtmosphericSpecification {
		
		protected:
		//public:
			std::vector< AtmosphericProfile * > profiles_;
			double pathAz_;
			std::vector< double > ranges_;
			bool strict_;
			double azTolerance_;
			
			// find the indices on either side of the requested point
			void getBrackets_( double r, int &nearIndex, int &farIndex );
			void init_();
			void clearOut();
			
			
		public:
			// constructors
			Slice();
			virtual ~Slice();
			
			virtual void readSummaryFile( std::string summaryFile );
			virtual AtmosphericProfile *getProfile( double lat, double lon, bool exact = false );
			virtual bool stratified() const;
			virtual double z0( double x, double y );
			void strict( bool st );
			void verticalTolerance( double eps_z );
			double pathAzimuth() const;
			
			virtual double t( double x, double y, double z );
			virtual double u( double x, double y, double z );
			virtual double v( double x, double y, double z );
			virtual double w( double x, double y, double z );
			virtual double p( double x, double y, double z );
			virtual double rho( double x, double y, double z );
		
			// spatial derivatives: temperature
			virtual double dtdz( double x, double y, double z );
			virtual double ddtdzdz( double x, double y, double z );
			virtual inline double dtdx( double x, double y, double z ) { return 0.0; }; 
			virtual inline double dtdy( double x, double y, double z ) { return 0.0; };  
			virtual inline double ddtdxdz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddtdydz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddtdxdx( double x, double y, double z ) { return 0.0; };
			virtual inline double ddtdydy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddtdxdy( double x, double y, double z ) { return 0.0; };

			// Spatial derivatives: effective sound speed
			virtual double dceffdz( double x, double y, double z, double phi );
			virtual double ddceffdzdz( double x, double y, double z, double phi );
			virtual inline double dceffdx( double x, double y, double z, double phi ) { return 0.0; };
			virtual inline double dceffdy( double x, double y, double z, double phi ) { return 0.0; };
			virtual inline double ddceffdxdz( double x, double y, double z, double phi ) { return 0.0; };
			virtual inline double ddceffdydz( double x, double y, double z, double phi ) { return 0.0; };
			virtual inline double ddceffdxdx( double x, double y, double z, double phi ) { return 0.0; };
			virtual inline double ddceffdydy( double x, double y, double z, double phi ) { return 0.0; };
			virtual inline double ddceffdxdy( double x, double y, double z, double phi ) { return 0.0; };

			// Spatial derivatives: thermodynamic sound speed
			virtual double dc0dz( double x, double y, double z );
			virtual double ddc0dzdz( double x, double y, double z );
			virtual inline double dc0dx( double x, double y, double z ) { return 0.0; };
			virtual inline double dc0dy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddc0dydy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddc0dxdx( double x, double y, double z ) { return 0.0; };
			virtual inline double ddc0dxdy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddc0dxdz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddc0dydz( double x, double y, double z ) { return 0.0; };

			// Spatial derivatives: zonal (E-W) winds
			virtual double dudz( double x, double y, double z );
			virtual double ddudzdz( double x, double y, double z );
			virtual inline double dudx( double x, double y, double z ) { return 0.0; };
			virtual inline double dudy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddudxdz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddudydz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddudxdx( double x, double y, double z ) { return 0.0; };
			virtual inline double ddudydy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddudxdy( double x, double y, double z ) { return 0.0; };

			// Spatial derivatives: meridional (N-S) winds
			virtual double dvdz( double x, double y, double z );
			virtual double ddvdzdz( double x, double y, double z );
			virtual inline double dvdx( double x, double y, double z ) { return 0.0; };
			virtual inline double dvdy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddvdxdz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddvdydz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddvdxdx( double x, double y, double z ) { return 0.0; };
			virtual inline double ddvdydy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddvdxdy( double x, double y, double z ) { return 0.0; };

			// Spatial derivatives: vertical winds
			virtual double dwdz( double x, double y, double z );
			virtual double ddwdzdz( double x, double y, double z );
			virtual inline double dwdx( double x, double y, double z ) { return 0.0; };
			virtual inline double dwdy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddwdxdz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddwdydz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddwdxdx( double x, double y, double z ) { return 0.0; };
			virtual inline double ddwdydy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddwdxdy( double x, double y, double z ) { return 0.0; };

			// Spatial derivatives: pressure
			virtual double dpdz( double x, double y, double z );
			virtual double ddpdzdz( double x, double y, double z );
			virtual inline double dpdx( double x, double y, double z ) { return 0.0; };
			virtual inline double dpdy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddpdxdz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddpdydz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddpdxdx( double x, double y, double z ) { return 0.0; };
			virtual inline double ddpdydy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddpdxdy( double x, double y, double z ) { return 0.0; };

			// Spatial derivatives: density
			virtual double drhodz( double x, double y, double z );
			virtual double ddrhodzdz( double x, double y, double z );
			virtual inline double drhodx( double x, double y, double z ) { return 0.0; };
			virtual inline double drhody( double x, double y, double z ) { return 0.0; };
			virtual inline double ddrhodxdz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddrhodydz( double x, double y, double z ) { return 0.0; };
			virtual inline double ddrhodxdx( double x, double y, double z ) { return 0.0; };
			virtual inline double ddrhodydy( double x, double y, double z ) { return 0.0; };
			virtual inline double ddrhodxdy( double x, double y, double z ) { return 0.0; };

		
	};  // end class Slice
	
}  // end namespace NCPA