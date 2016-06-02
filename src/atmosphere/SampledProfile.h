#ifndef __SAMPLEDPROFILE_H__
#define __SAMPLEDPROFILE_H__

#include "AtmosphericProfile.h"
#include "geographic.h"
#include <string>
#include <iostream>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

namespace NCPA {
	
	class SampledProfile : public AtmosphericProfile {
		
	protected:
		unsigned int nz_;
		double z0_, propAz_, minZ_, maxZ_;
		double *z_, *t_, *u_, *v_, *w_, *rho_, *p_, *c0_, *ceff_;

		// Spline objects and accelerators
		gsl_interp_accel *t_acc, *u_acc, *v_acc, *w_acc, *p_acc, *rho_acc, *c0_acc, *ceff_acc;
		gsl_spline *t_spline, *u_spline, *v_spline, *w_spline, *p_spline, *rho_spline, *c0_spline, *ceff_spline;
		
		virtual void init_();
		virtual void clearOut();
		virtual void initSplines();
		
		virtual double d_dz( double z, double *q );
		virtual double dd_dzdz( double z, double *q );
		
	public:
		SampledProfile();
		virtual ~SampledProfile();
		SampledProfile( double lat, double lon, int nz, double *z, double *t,
			double *u, double *v, double *w = 0, double *rho = 0, 
			double *p = 0, double z0 = -99999.0 );
		SampledProfile( int nz, double *z, double *t, double *u, double *v, double *w = 0,
			double *rho = 0, double *p = 0, double z0 = -99999.0 );
		SampledProfile( SampledProfile &s );
		
		/**
		 * File constructor.  Takes a filename and an optional number of lines to skip, and
		 * reads the values in in CSV format.  The order in which the values are presented
		 * is represented by a string of characters such as "ztuvwpd", where:
		 *      z: height (km)
		 *      t: temperature (K)
		 *      u: zonal winds (km/s)
		 *      v: meridional winds (km/s)
		 *      w: vertical winds (km/s)
		 *      p: pressure (hPa)
		 *      d: density (g/cm3)
		 *      Any other character will cause the field to be skipped.  If there are not enough
		 *      fields in the line, the extra order characters will be ignored.
		 * @param filename The name of the file to read
		 * @param order The order in which fields are present
		 * @param propAz The azimuth of propagation, for calculating the effective sound speed
		 * @param skiplines The number of lines to skip before reading
		 * @param inMPS Flag to indicate that data is in mps, not kmps
		 */
		SampledProfile( std::string filename, const char *order, int skiplines = 0 );
		SampledProfile( std::string filename, const char *order, int skiplines, bool inMPS );
		//SampledProfile( std::istream instream, const char *order, int skiplines = 0 );
		
		
		virtual int nz() const;
		//virtual double dz_min() const;
		double lat() const;
		double lon() const;
		//void resample( double newDZ );
		void setPropagationAzimuth( double newAzimuth );
		double getPropagationAzimuth() const;
		
		virtual double t( double z );
		virtual double u( double z );
		virtual double v( double z );
		virtual double w( double z );
		virtual double p( double z );
		virtual double rho( double z );
		virtual double c0( double z );
		virtual double ceff( double z, double phi );
		
		virtual double z0();
		//virtual unsigned int firstValidIndex();
		virtual double z( unsigned int zind ) const;
		
		virtual double dtdz( double z );
		virtual double dceffdz( double z, double phi );
		virtual double dc0dz( double z );
		virtual double dudz( double z );
		virtual double dvdz( double z );
		virtual double dwdz( double z );
		virtual double dpdz( double z );
		virtual double drhodz( double z );
		virtual double ddtdzdz( double z );
		virtual double ddceffdzdz( double z, double phi );
		virtual double ddc0dzdz( double z );
		virtual double ddudzdz( double z );
		virtual double ddvdzdz( double z );
		virtual double ddwdzdz( double z );
		virtual double ddpdzdz( double z );
		virtual double ddrhodzdz( double z );
		
		virtual void get_z( double *target, int nz ) const;
		virtual void get_t( double *target, int nz ) const;
		virtual void get_u( double *target, int nz ) const;
		virtual void get_v( double *target, int nz ) const;
		virtual void get_w( double *target, int nz ) const;
		virtual void get_p( double *target, int nz ) const;
		virtual void get_rho( double *target, int nz ) const;
		virtual void get_c0( double *target, int nz ) const;
		virtual void get_ceff( double *target, int nz ) const;
		/*
		virtual void set_t( double *source, int nz ) const;
		virtual void set_u( double *source, int nz ) const;
		virtual void set_v( double *source, int nz ) const;
		virtual void set_w( double *source, int nz ) const;
		virtual void set_p( double *source, int nz ) const;
		virtual void set_rho( double *source, int nz ) const;
		*/
		
		virtual unsigned int z2ind_( double z_in ) const;
		virtual unsigned int z2ind_floor_( double z_in ) const;
	};
}

#endif
