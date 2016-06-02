#ifndef __ATMOSPHERIC_SPECIFICATION_H__
#define __ATMOSPHERIC_SPECIFICATION_H__

#include "AtmosphericProfile.h"

/**
 * The NCPA namespace.  Used to indicate NCPA-written code.
 */
namespace NCPA {
	
	struct AtmosphericPoint {
		double  u, v, w, p, t, rho, c0,
			dudx, dudy, dudz,
			dvdx, dvdy, dvdz,
			dwdx, dwdy, dwdz,
			dpdx, dpdy, dpdz,
			dtdx, dtdy, dtdz,
			drhodx, drhody, drhodz,
			dc0dx, dc0dy, dc0dz,
			ddudxdx, ddudxdy, ddudydy, ddudxdz, ddudydz, ddudzdz,
			ddvdxdx, ddvdxdy, ddvdydy, ddvdxdz, ddvdydz, ddvdzdz,
			ddwdxdx, ddwdxdy, ddwdydy, ddwdxdz, ddwdydz, ddwdzdz,
			ddpdxdx, ddpdxdy, ddpdydy, ddpdxdz, ddpdydz, ddpdzdz,
			ddtdxdx, ddtdxdy, ddtdydy, ddtdxdz, ddtdydz, ddtdzdz,
			ddrhodxdx, ddrhodxdy, ddrhodydy, ddrhodxdz, ddrhodydz, ddrhodzdz,
			ddc0dxdx, ddc0dxdy, ddc0dydy, ddc0dxdz, ddc0dydz, ddc0dzdz,
			c00;
	};

	/**
	 * The AtmosphericSpecification abstract base class.  All classes that are to be used as atmospheric specifications
	 * for propagation modeling and the like should be children of this class.
	 */
	class AtmosphericSpecification {

		protected:
			bool good_; /**< Status flag.  Indicates whether or not the profile contains valid data. */
			double  eps_x, /**< Horizontal tolerance.  Used for computing derivatives in the horizontal directions. */
				eps_z; /**< Vertical tolerance.  Used for computing vertical derivatives. Must be >= dz for stratified specifications.*/

		public:
			/**
			  * Virtual destructor.
			  */
			virtual ~AtmosphericSpecification();

			/**
			  * Status test.  Used to test whether or not the specification contains valid data and is ready for use.
			  * @return TRUE if the specification is ready for use, FALSE otherwise.
			  */
			virtual bool good();
			
			/**
			 * Pull a 1-D atmospheric profile from the Specification.
			 * @param lat The latitude of the desired profile
			 * @param lon The longitude of the desired profile
			 * @param exact Flag to indicate whether the exact location (if available) must be returned (true), or 
			 * 		any Profile within eps_x km may be returned (false).
			 */
			virtual NCPA::AtmosphericProfile *getProfile( double lat, double lon, bool exact = false ) = 0;

			/**
			  * Stratification test.  Signals whether or not the specification is stratified (i.e. vertically non-continuous).
			  * @return TRUE if the specification is stratified, FALSE otherwise.
			  */
			virtual bool stratified() const = 0;

			// Properties of the atmospheric specification, must be overridden by subclass
			/**
			  * Minimum valid altitude.  Returns the minimum altitude for which the specification can be queried for a
			  * valid value.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @return The minimum valid altitude, in km.
			  */
			virtual double z0( double x, double y ) = 0;

			/**
			  * Temperature.  Returns the Kelvin temperature at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the temperature, in km relative to MSL.
			  * @return The temperature at point (x,y,z), in Kelvin.
			  */
			virtual double t( double x, double y, double z ) = 0;
			
			/**
			  * Zonal wind.  Returns the zonal (positive to East) wind speed at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The zonal wind speed at point (x,y,z), in km/s.
			  */
			virtual double u( double x, double y, double z ) = 0;

			/**
			  * Meridional wind.  Returns the meridional (positive to North) wind speed at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The meridional wind speed at point (x,y,z), in km/s.
			  */
			virtual double v( double x, double y, double z ) = 0;

			/**
			  * Vertical wind.  Returns the vertical (positive up) wind speed at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The vertical wind speed at point (x,y,z), in km/s.
			  */
			virtual double w( double x, double y, double z ) = 0;

			/**
			  * Density.  Returns the air density at the indicated point.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the density, in km relative to MSL.
			  * @return The density at point (x,y,z), in g/cm3 (?).
			  * @todo Confirm the density units.
			  */
			virtual double rho( double x, double y, double z ) = 0;

			/**
			 * Pressure.  Returns the atmospheric pressure at the indicated point.
			 * @param x The distance from the origin in the X direction, in km.
			 * @param y The distance from the origin in the Y direction, in km.
			 * @param z The altitude for which to return the pressure, in km relative to MSL.
			 * @return The pressure at point (x,y,z), in hPa.
			 */
			virtual double p( double x, double y, double z ) = 0;
			
			// Calculated properties of atmospheric, may be overridden

			/**
			  * Effective sound speed.  Returns the effective (i.e. advected) sound speed in the indicated direction.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the sound speed, in km relative to MSL.
			  * @param phi The direction in which to determine the effective soundspeed, in degrees clockwise from North
			  * @return The effective sound speed along azimuth phi at point (x,y,z), in km/s.
			  */
			virtual double ceff( double x, double y, double z, double phi );

			/**
			  * Static sound speed.  Returns the thermodynamic sound speed, not including wind effects.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the sound speed, in km relative to MSL.
			  * @return The static (thermodynamic) sound speed at point (x,y,z), in km/s.
			  */
			virtual double c0( double x, double y, double z );

			/**
			  * Wind Speed.  Returns the magnitude of the wind vector.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The scalar wind speed at point (x,y,z), in km/s.
			  */
			virtual double wspeed( double x, double y, double z );

			/**
			  * Wind Direction.  Returns the 2-D geographic direction of the wind vector.
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind direction, in km relative to MSL.
			  * @return The direction of the projection of the wind vector into the X-Y plane, in degrees clockwise from North.
			  */
			virtual double wdirection( double x, double y, double z );

			/**
			  * Wind Component.  Returns the magnitude of the wind vector projected first into the X-Y plane, then in the specified direction
			  * @param x The distance from the origin in the X direction, in km.
			  * @param y The distance from the origin in the Y direction, in km.
			  * @param z The altitude for which to return the wind direction, in km relative to MSL.
			  * @param phi The direction, in degrees clockwise from North, into which to project the wind vector
			  * @return The magnitude of the projection of the horizonal wind vector into the phi direction, in km/s.
			  */
			virtual double wcomponent( double x, double y, double z, double phi );

			// Spatial derivatives: temperature

			virtual double dtdx( double x, double y, double z ); 
			virtual double dtdy( double x, double y, double z );  
			virtual double dtdz( double x, double y, double z );
			virtual double ddtdzdz( double x, double y, double z );
			virtual double ddtdxdz( double x, double y, double z );
			virtual double ddtdydz( double x, double y, double z );
			virtual double ddtdxdx( double x, double y, double z );
			virtual double ddtdydy( double x, double y, double z );
			virtual double ddtdxdy( double x, double y, double z );

			// Spatial derivatives: effective sound speed
			virtual double dceffdx( double x, double y, double z, double phi );
			virtual double dceffdy( double x, double y, double z, double phi );
			virtual double dceffdz( double x, double y, double z, double phi );
			virtual double ddceffdzdz( double x, double y, double z, double phi );
			virtual double ddceffdxdz( double x, double y, double z, double phi );
			virtual double ddceffdydz( double x, double y, double z, double phi );
			virtual double ddceffdxdx( double x, double y, double z, double phi );
			virtual double ddceffdydy( double x, double y, double z, double phi );
			virtual double ddceffdxdy( double x, double y, double z, double phi );

			// Spatial derivatives: thermodynamic sound speed
			virtual double dc0dx( double x, double y, double z );
			virtual double dc0dy( double x, double y, double z );
			virtual double dc0dz( double x, double y, double z );
			virtual double ddc0dzdz( double x, double y, double z );
			virtual double ddc0dydy( double x, double y, double z );
			virtual double ddc0dxdx( double x, double y, double z );
			virtual double ddc0dxdy( double x, double y, double z );
			virtual double ddc0dxdz( double x, double y, double z );
			virtual double ddc0dydz( double x, double y, double z );

			// Spatial derivatives: zonal (E-W) winds
			virtual double dudx( double x, double y, double z );
			virtual double dudy( double x, double y, double z );
			virtual double dudz( double x, double y, double z );
			virtual double ddudzdz( double x, double y, double z );
			virtual double ddudxdz( double x, double y, double z );
			virtual double ddudydz( double x, double y, double z );
			virtual double ddudxdx( double x, double y, double z );
			virtual double ddudydy( double x, double y, double z );
			virtual double ddudxdy( double x, double y, double z );

			// Spatial derivatives: meridional (N-S) winds
			virtual double dvdx( double x, double y, double z );
			virtual double dvdy( double x, double y, double z );
			virtual double dvdz( double x, double y, double z );
			virtual double ddvdzdz( double x, double y, double z );
			virtual double ddvdxdz( double x, double y, double z );
			virtual double ddvdydz( double x, double y, double z );
			virtual double ddvdxdx( double x, double y, double z );
			virtual double ddvdydy( double x, double y, double z );
			virtual double ddvdxdy( double x, double y, double z );

			// Spatial derivatives: vertical winds
			virtual double dwdx( double x, double y, double z );
			virtual double dwdy( double x, double y, double z );
			virtual double dwdz( double x, double y, double z );
			virtual double ddwdzdz( double x, double y, double z );
			virtual double ddwdxdz( double x, double y, double z );
			virtual double ddwdydz( double x, double y, double z );
			virtual double ddwdxdx( double x, double y, double z );
			virtual double ddwdydy( double x, double y, double z );
			virtual double ddwdxdy( double x, double y, double z );

			// Spatial derivatives: pressure
			virtual double dpdx( double x, double y, double z );
			virtual double dpdy( double x, double y, double z );
			virtual double dpdz( double x, double y, double z );
			virtual double ddpdzdz( double x, double y, double z );
			virtual double ddpdxdz( double x, double y, double z );
			virtual double ddpdydz( double x, double y, double z );
			virtual double ddpdxdx( double x, double y, double z );
			virtual double ddpdydy( double x, double y, double z );
			virtual double ddpdxdy( double x, double y, double z );

			// Spatial derivatives: density
			virtual double drhodx( double x, double y, double z );
			virtual double drhody( double x, double y, double z );
			virtual double drhodz( double x, double y, double z );
			virtual double ddrhodzdz( double x, double y, double z );
			virtual double ddrhodxdz( double x, double y, double z );
			virtual double ddrhodydz( double x, double y, double z );
			virtual double ddrhodxdx( double x, double y, double z );
			virtual double ddrhodydy( double x, double y, double z );
			virtual double ddrhodxdy( double x, double y, double z );

			virtual NCPA::AtmosphericPoint take_snapshot(double x, double y, double z);
			
	};   // class AtmosphericSpecification

}   // namespace NCPA




#endif    // #ifndef __ATMOSPHERIC_SPECIFICATION_H__
