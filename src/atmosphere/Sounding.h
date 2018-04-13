#ifndef __SOUNDING_H__
#define __SOUNDING_H__

#include "AtmosphericProfile.h"
#include "AtmosphericSpecification.h"
#include "geographic.h"
#include <string>

namespace NCPA {

	class Sounding : public AtmosphericSpecification {

		protected:
			
			AtmosphericProfile *profile_;
			virtual void init_();
			virtual void clearOut();

		public:
			Sounding();
			virtual ~Sounding();
			Sounding( double lat, double lon, int nz, double *z, double *t,
                                double *u, double *v, double *w = 0, double *rho = 0,
				double *p = 0, double z0 = -99999.0 );
			Sounding( int nz, double *z, double *t, double *u, double *v, double *w = 0,
                                double *rho = 0, double *p = 0, double z0 = -99999.0 );
			Sounding( AtmosphericProfile *p );
			
			NCPA::AtmosphericProfile *getProfile( double lat, double lon, bool exact = false );
			

			bool stratified() const;
			double lat() const;
			double lon() const;

			virtual double t( double x, double y, double z );
			virtual double u( double x, double y, double z );
			virtual double v( double x, double y, double z );
			virtual double w( double x, double y, double z );
			virtual double p( double x, double y, double z );
			virtual double rho( double x, double y, double z );
			virtual double c0( double x, double y, double z );
			virtual double z0( double x, double y );
			
                        virtual double dtdz( double x, double y, double z );
                        virtual double dceffdz( double x, double y, double z, double phi );
                        virtual double dc0dz( double x, double y, double z );
                        virtual double dudz( double x, double y, double z );
                        virtual double dvdz( double x, double y, double z );
                        virtual double dwdz( double x, double y, double z );
                        virtual double dpdz( double x, double y, double z );
                        virtual double drhodz( double x, double y, double z );
                        virtual double ddtdzdz( double x, double y, double z );
                        virtual double ddceffdzdz( double x, double y, double z, double phi );
                        virtual double ddc0dzdz( double x, double y, double z );
                        virtual double ddudzdz( double x, double y, double z );
                        virtual double ddvdzdz( double x, double y, double z );
                        virtual double ddwdzdz( double x, double y, double z );
                        virtual double ddpdzdz( double x, double y, double z );
                        virtual double ddrhodzdz( double x, double y, double z );


                        // Spatial derivatives, all == 0 because Sounding is 1-D
                        double dtdx( double x, double y, double z ); /**< @return 0.0 */
                        double dtdy( double x, double y, double z ); /**< @return 0.0 */
                        double ddtdxdz( double x, double y, double z ); /**< @return 0.0 */
                        double ddtdydz( double x, double y, double z ); /**< @return 0.0 */
                        double ddtdxdx( double x, double y, double z ); /**< @return 0.0 */
                        double ddtdydy( double x, double y, double z ); /**< @return 0.0 */
                        double ddtdxdy( double x, double y, double z ); /**< @return 0.0 */
                        double dceffdx( double x, double y, double z, double phi ); /**< @return 0.0 */
                        double dceffdy( double x, double y, double z, double phi ); /**< @return 0.0 */
                        double ddceffdxdz( double x, double y, double z, double phi ); /**< @return 0.0 */
                        double ddceffdydz( double x, double y, double z, double phi ); /**< @return 0.0 */
                        double ddceffdxdx( double x, double y, double z, double phi ); /**< @return 0.0 */
                        double ddceffdydy( double x, double y, double z, double phi ); /**< @return 0.0 */
                        double ddceffdxdy( double x, double y, double z, double phi ); /**< @return 0.0 */
                        double dc0dx( double x, double y, double z ); /**< @return 0.0 */
                        double dc0dy( double x, double y, double z ); /**< @return 0.0 */
                        double ddc0dydy( double x, double y, double z ); /**< @return 0.0 */
                        double ddc0dxdx( double x, double y, double z ); /**< @return 0.0 */
                        double ddc0dxdy( double x, double y, double z ); /**< @return 0.0 */
                        double ddc0dxdz( double x, double y, double z ); /**< @return 0.0 */
                        double ddc0dydz( double x, double y, double z ); /**< @return 0.0 */
                        double dudx( double x, double y, double z ); /**< @return 0.0 */
                        double dudy( double x, double y, double z ); /**< @return 0.0 */
                        double ddudxdz( double x, double y, double z ); /**< @return 0.0 */
                        double ddudydz( double x, double y, double z ); /**< @return 0.0 */
                        double ddudxdx( double x, double y, double z ); /**< @return 0.0 */
                        double ddudydy( double x, double y, double z ); /**< @return 0.0 */
                        double ddudxdy( double x, double y, double z ); /**< @return 0.0 */
                        double dvdx( double x, double y, double z ); /**< @return 0.0 */
                        double dvdy( double x, double y, double z ); /**< @return 0.0 */
                        double ddvdxdz( double x, double y, double z ); /**< @return 0.0 */
                        double ddvdydz( double x, double y, double z ); /**< @return 0.0 */
                        double ddvdxdx( double x, double y, double z ); /**< @return 0.0 */
                        double ddvdydy( double x, double y, double z ); /**< @return 0.0 */
                        double ddvdxdy( double x, double y, double z ); /**< @return 0.0 */
                        double dwdx( double x, double y, double z ); /**< @return 0.0 */
                        double dwdy( double x, double y, double z ); /**< @return 0.0 */
                        double ddwdxdz( double x, double y, double z ); /**< @return 0.0 */
                        double ddwdydz( double x, double y, double z ); /**< @return 0.0 */
                        double ddwdxdx( double x, double y, double z ); /**< @return 0.0 */
                        double ddwdydy( double x, double y, double z ); /**< @return 0.0 */
                        double ddwdxdy( double x, double y, double z ); /**< @return 0.0 */
                        double dpdx( double x, double y, double z ); /**< @return 0.0 */
                        double dpdy( double x, double y, double z ); /**< @return 0.0 */
                        double ddpdxdz( double x, double y, double z ); /**< @return 0.0 */
                        double ddpdydz( double x, double y, double z ); /**< @return 0.0 */
                        double ddpdxdx( double x, double y, double z ); /**< @return 0.0 */
                        double ddpdydy( double x, double y, double z ); /**< @return 0.0 */
                        double ddpdxdy( double x, double y, double z ); /**< @return 0.0 */
                        double drhodx( double x, double y, double z ); /**< @return 0.0 */
                        double drhody( double x, double y, double z ); /**< @return 0.0 */
                        double ddrhodxdz( double x, double y, double z ); /**< @return 0.0 */
                        double ddrhodydz( double x, double y, double z ); /**< @return 0.0 */
                        double ddrhodxdx( double x, double y, double z ); /**< @return 0.0 */
                        double ddrhodydy( double x, double y, double z ); /**< @return 0.0 */
                        double ddrhodxdy( double x, double y, double z ); /**< @return 0.0 */

			
			
	};
}

#endif
