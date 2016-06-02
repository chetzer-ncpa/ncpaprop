#ifndef __PROFILEGROUP_H__
#define __PROFILEGROUP_H__

#include "AtmosphericSpecification.h"
#include "AtmosphericProfile.h"
#include <vector>

namespace NCPA {

	class ProfileGroup : public AtmosphericSpecification {

		protected:
			std::vector< NCPA::AtmosphericProfile * > profiles_, buffer_;
			double lat0_, lon0_;
			
			

		public:
			/**
			  * Virtual destructor.
			  */
			virtual ~ProfileGroup();

			virtual void setOrigin( double lat, double lon );
			virtual NCPA::AtmosphericProfile *getProfile( double lat, double lon, bool exact = false );

			/**
			  * Minimum valid altitude.  Returns the minimum altitude for which the profile can be queried for a
			  * valid value.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @return The minimum valid altitude, in km.
			  */
			virtual double z0( double x, double y );

			/**
			  * Temperature.  Returns the Kelvin temperature at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the temperature, in km relative to MSL.
			  * @return The temperature at point (x,y,z), in Kelvin.
			  */
			virtual double t( double x, double y, double z );
			
			/**
			  * Zonal wind.  Returns the zonal (positive to East) wind speed at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The zonal wind speed at point (x,y,z), in km/s.
			  */
			virtual double u( double x, double y, double z );

			/**
			  * Meridional wind.  Returns the meridional (positive to North) wind speed at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The meridional wind speed at point (x,y,z), in km/s.
			  */
			virtual double v( double x, double y, double z );

			/**
			  * Vertical wind.  Returns the vertical (positive up) wind speed at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The vertical wind speed at point (x,y,z), in km/s.
			  */
			virtual double w( double x, double y, double z );

			/**
			  * Density.  Returns the air density at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the density, in km relative to MSL.
			  * @return The density at point (x,y,z), in g/cm3 (?).
			  * @todo Confirm the density units.
			  */
			virtual double rho( double x, double y, double z );
			
			/**
			  * Pressure.  Returns the air pressure at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the density, in km relative to MSL.
			  * @return The pressure at point (x,y,z), in hPa (?).
			  * @todo Confirm the pressure units.
			  */
			virtual double p( double x, double y, double z );

			// Calculated properties of atmospheric, may be overridden

			// Spatial derivatives: temperature
			virtual double dtdz( double x, double y, double z );
			virtual double ddtdzdz( double x, double y, double z );
			virtual double ddtdxdz( double x, double y, double z );
			virtual double ddtdydz( double x, double y, double z );

			// Spatial derivatives: effective sound speed
			virtual double dceffdz( double x, double y, double z, double phi );
			virtual double ddceffdzdz( double x, double y, double z, double phi );
			virtual double ddceffdxdz( double x, double y, double z, double phi );
			virtual double ddceffdydz( double x, double y, double z, double phi );

			// Spatial derivatives: thermodynamic sound speed
			virtual double dc0dz( double x, double y, double z );
			virtual double ddc0dzdz( double x, double y, double z );
			virtual double ddc0dxdz( double x, double y, double z );
			virtual double ddc0dydz( double x, double y, double z );

			// Spatial derivatives: zonal (E-W) winds
			virtual double dudz( double x, double y, double z );
			virtual double ddudzdz( double x, double y, double z );
			virtual double ddudxdz( double x, double y, double z );
			virtual double ddudydz( double x, double y, double z );

			// Spatial derivatives: meridional (N-S) winds
			virtual double dvdz( double x, double y, double z );
			virtual double ddvdzdz( double x, double y, double z );
			virtual double ddvdxdz( double x, double y, double z );
			virtual double ddvdydz( double x, double y, double z );

			// Spatial derivatives: vertical winds
			virtual double dwdz( double x, double y, double z );
			virtual double ddwdzdz( double x, double y, double z );
			virtual double ddwdxdz( double x, double y, double z );
			virtual double ddwdydz( double x, double y, double z );
			
			// Spatial derivatives: pressure
			virtual double dpdz( double x, double y, double z );
			virtual double ddpdzdz( double x, double y, double z );
			virtual double ddpdxdz( double x, double y, double z );
			virtual double ddpdydz( double x, double y, double z );
			
			// Spatial derivatives: density
			virtual double drhodz( double x, double y, double z );
			virtual double ddrhodzdz( double x, double y, double z );
			virtual double ddrhodxdz( double x, double y, double z );
			virtual double ddrhodydz( double x, double y, double z );

			// horizontal derivatives of quantities not required by base class
			virtual double dpdx( double x, double y, double z );
			virtual double dpdy( double x, double y, double z );
			virtual double drhodx( double x, double y, double z );
			virtual double drhody( double x, double y, double z );




	};   // class ProfileGroup

}   // namespace NCPA




#endif    // #ifndef __SOUNDINGGROUP_H__
