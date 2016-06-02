#ifndef __ATMOSPHERIC_PROFILE_H__
#define __ATMOSPHERIC_PROFILE_H__

#include "geographic.h"

/**
 * The NCPA namespace.  Used to indicate NCPA-written code.
 */
namespace NCPA {

	/**
	 * The AtmosphericSpecification abstract base class.  All classes that are to be used as atmospheric specifications
	 * for propagation modeling and the like should be children of this class.
	 */
	class AtmosphericProfile {

		protected:
			bool good_; /**< Status flag.  Indicates whether or not the profile contains valid data. */
			double eps_z; /**< Vertical tolerance.  Used for computing vertical derivatives. Must be >= dz for stratified profiles.*/
			Location *origin_;

		public:
			/**
			  * Virtual destructor.
			  */
			virtual ~AtmosphericProfile();

			
			virtual void setOrigin( double lat, double lon );
			virtual double lat();
			virtual double lon();
			
			/**
			  * Status test.  Used to test whether or not the specification contains valid data and is ready for use.
			  * @return TRUE if the specification is ready for use, FALSE otherwise.
			  */
			virtual bool good();

			// Properties of the atmospheric specification, must be overridden by subclass
			/**
			  * Minimum valid altitude.  Returns the minimum altitude for which the specification can be queried for a
			  * valid value.
			  * @return The minimum valid altitude, in km.
			  */
			virtual double z0() = 0;

			/**
			  * Temperature.  Returns the Kelvin temperature at the indicated point.
			  * @param z The altitude for which to return the temperature, in km relative to MSL.
			  * @return The temperature at point (x,y,z), in Kelvin.
			  */
			virtual double t( double z ) = 0;
			
			/**
			  * Zonal wind.  Returns the zonal (positive to East) wind speed at the indicated point.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The zonal wind speed at point (x,y,z), in km/s.
			  */
			virtual double u( double z ) = 0;

			/**
			  * Meridional wind.  Returns the meridional (positive to North) wind speed at the indicated point.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The meridional wind speed at point (x,y,z), in km/s.
			  */
			virtual double v( double z ) = 0;

			/**
			  * Vertical wind.  Returns the vertical (positive up) wind speed at the indicated point.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The vertical wind speed at point (x,y,z), in km/s.
			  */
			virtual double w( double z ) = 0;
			
			/**
			 * Pressure.  Returns the atmospheric pressure at the indicated point.
			 * @param z The altitude for which to return the pressure, in km relative to MSL.
			 * @return The pressure in hPa
			 */
			virtual double p( double z ) = 0;

			/**
			  * Density.  Returns the air density at the indicated point.
			  * @param z The altitude for which to return the density, in km relative to MSL.
			  * @return The density at point (x,y,z), in g/cm3.
			  */
			virtual double rho( double z ) = 0;

			// Calculated properties of atmospheric, may be overridden

			/**
			  * Effective sound speed.  Returns the effective (i.e. advected) sound speed in the indicated direction.
			  * @param z The altitude for which to return the sound speed, in km relative to MSL.
			  * @param phi The direction in which to determine the effective soundspeed, in degrees clockwise from North
			  * @return The effective sound speed along azimuth phi at point (x,y,z), in km/s.
			  */
			virtual double ceff( double z, double phi );

			/**
			  * Static sound speed.  Returns the thermodynamic sound speed, not including wind effects.
			  * @param z The altitude for which to return the sound speed, in km relative to MSL.
			  * @return The static (thermodynamic) sound speed at point (x,y,z), in km/s.
			  */
			virtual double c0( double z );

			/**
			  * Wind Speed.  Returns the magnitude of the wind vector.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The scalar wind speed at point (x,y,z), in km/s.
			  */
			virtual double wspeed( double z );

			/**
			  * Wind Direction.  Returns the 2-D geographic direction of the wind vector.
			  * @param z The altitude for which to return the wind direction, in km relative to MSL.
			  * @return The direction of the projection of the wind vector into the X-Y plane, in degrees clockwise from North.
			  */
			virtual double wdirection( double z );

			/**
			  * Wind Component.  Returns the magnitude of the wind vector projected first into the X-Y plane, then in the specified direction
			  * @param z The altitude for which to return the wind direction, in km relative to MSL.
			  * @param phi The direction, in degrees clockwise from North, into which to project the wind vector
			  * @return The magnitude of the projection of the horizonal wind vector into the phi direction, in km/s.
			  */
			virtual double wcomponent( double z, double phi );

			// Spatial derivatives: temperature

			virtual double dtdz(  double z );
			virtual double ddtdzdz(  double z );

			// Spatial derivatives: effective sound speed
			virtual double dceffdz(  double z, double phi );
			virtual double ddceffdzdz(  double z, double phi );

			// Spatial derivatives: thermodynamic sound speed
			virtual double dc0dz(  double z );
			virtual double ddc0dzdz(  double z );

			// Spatial derivatives: zonal (E-W) winds
			virtual double dudz(  double z );
			virtual double ddudzdz(  double z );

			// Spatial derivatives: meridional (N-S) winds
			virtual double dvdz(  double z );
			virtual double ddvdzdz(  double z );

			// Spatial derivatives: vertical winds
			virtual double dwdz(  double z );
			virtual double ddwdzdz(  double z );

			// Spatial derivatives: pressure
			virtual double dpdz(  double z );
			virtual double ddpdzdz(  double z );
			
			// Spatial derivatives: pressure
			virtual double drhodz(  double z );
			virtual double ddrhodzdz(  double z );
			

	};   // class AtmosphericProfile

}   // namespace NCPA




#endif    // #ifndef __ATMOSPHERIC_PROFILE_H__
