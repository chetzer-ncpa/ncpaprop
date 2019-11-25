/*
Abstract base class AtmosphericProfile - subclasses must override at least:
double minimumAltitude() - return minimum valid altitude
double maximumAltitude() - return maximum valid altitude
void convertUnits( ATMOSPHERIC_QUANTITY, UNITS_TYPE ) - change internal units
double get( ATMOSPHERIC_QUANTITY quantity, double z ) - return a single value from the profile
double get( ATMOSPHERIC_QUANTITY quantity, double z, UNITS_TYPE new_units ) - return a single
	value from the profile, converted into specific units
void get( ATMOSPHERIC_QUANTITY quantity, unsigned int n_z, const double *z, double *target )
	- extract a vector of values from the profile
void get( ATMOSPHERIC_QUANTITY quantity, unsigned int n_z, const double *z, UNITS_TYPE new_units, 
	double *target ) - extract a vector of values from the profile, converted into specific units
bool has( ATMOSPHERIC_QUANTITY quantity ) - returns whether the profile contains the quantity
*/

#ifndef ATMOSPHERIC_PROFILE_H__
#define ATMOSPHERIC_PROFILE_H__

#include <map>
#include "geographic.h"
#include "units.h"


/**
 * The NCPA namespace.  Used to indicate NCPA-written code.
 */
namespace NCPA {
	
	enum ATMOSPHERIC_QUANTITY : unsigned int {
		QUANTITY_NONE = 0,
		QUANTITY_ALTITUDE,
		QUANTITY_TEMPERATURE,
		QUANTITY_WIND_SPEED_WEST_TO_EAST,
		QUANTITY_WIND_SPEED_SOUTH_TO_NORTH,
		QUANTITY_WIND_SPEED_VERTICAL,
		QUANTITY_AIR_PRESSURE,
		QUANTITY_AIR_DENSITY,
		QUANTITY_STATIC_SOUND_SPEED,
		QUANTITY_EFFECTIVE_SOUND_SPEED,
		QUANTITY_WIND_SPEED,
		QUANTITY_WIND_SPEED_HORIZONTAL,
		QUANTITY_WIND_DIRECTION
	};

	/**
	 * The AtmosphericProfile abstract base class.  All classes that are to be used as atmospheric profiles
	 * for propagation modeling and the like should be children of this class.  This class specifically
	 * represents 1-D atmospheric profiles, whereas the AtmosphericSpecification ABC can be 2-D or 3-D.
	 */
	class AtmosphericProfile {

		protected:
			
			// Methods
			
			/** calculate c0 using internal logic to determine how */
			virtual double calculate_c0_( double z );
			
			/** calculate c0 using sqrt( gamma * R * T ) */
			virtual double calculate_c0_using_t_( double z );
			
			/** calculate c0 using sqrt( gamma * p / rho ) */
			virtual double calculate_c0_using_p_( double z );	
			
			/** calculate total (including vertical) wind speed */
			virtual double calculate_total_wind_speed_( double z );
			
			/** calculate total horiztonal wind speed */
			virtual double calculate_horizontal_wind_speed_( double z );
			
			/** calculate the direction of the horizontal wind vector */
			virtual double calculate_wind_direction_( double z );
			
			
			// Members
			
			/** Status flag.  Indicates whether or not the profile contains valid data. */
			bool good_;
			
			/** 
			Vertical tolerance.  
			Used for computing vertical derivatives. Must be >= dz for stratified profiles.
			*/
			double eps_z;
			
			/**
			Geographic location of profile.
			*/
			Location *origin_;
			
			//bool 	hasW_,		/**< Indicates whether the specification includes vertical winds. */
			//	hasP_,		/**< Indicates whether the specification includes pressure. */
			//	hasRho_;	/**< Indicates whether the specification includes density. */
				
			/**
			Units tracker.
			Keeps track of the internal units of all atmospheric quantities.
			*/
			std::map< ATMOSPHERIC_QUANTITY, UNITS_TYPE > units_map_;


		public:
			
			/**
			  * Destructor.
			  */
			virtual ~AtmosphericProfile();

			virtual void setOrigin( double lat, double lon );
			virtual double lat() const;
			virtual double lon() const;
			
			/**
			  * Specify the internal units for atmospheric quantities.
			  * Does not do any conversion, just tells the Profile object "any numbers
			  * you currently have or that are provided to you are in these units."
			  * No validity checking is performed to ensure that e.g. wind speeds
			  * are not expressed in temperature units.  It is on the user not to do
			  * anything stupid.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @param units	The units that quantity is considered to be in
			  */
			virtual void setUnits( ATMOSPHERIC_QUANTITY quantity, UNITS_TYPE units );
			
			/**
			  * Returns the current internal units for a given quantity.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @return		The units that quantity is considered to be in
			  * @throws out_of_bounds if no units have been specified for the quantity
			  */
			virtual UNITS_TYPE getUnits( ATMOSPHERIC_QUANTITY quantity );
			
			
			/**
			  * Status test.  Used to test whether or not the specification contains valid data and is ready for use.
			  * @return TRUE if the specification is ready for use, FALSE otherwise.
			  */
			virtual bool good();
			
			
			// Vertical derivatives
			virtual double get_dZ( ATMOSPHERIC_QUANTITY quantity, double z );
			virtual double get_ddZ( ATMOSPHERIC_QUANTITY quantity, double z );
			
			
			
			// Pure virtual methods, must be overridden
			
			/**
			  * Returns the minimum valid altitude for the profile.
			  * @return	The minimum valid altitude
			  */
			virtual double minimumAltitude() = 0;
			
			/**
			  * Returns the maximum valid altitude for the profile.
			  * @return	The maximum valid altitude
			  */
			virtual double maximumAltitude() = 0;
		
			/**
			  * Converts units internally.
			  * Converts the internal values of the Profile to a new unit.  Specifics are
			  * implementation-dependent.  Should throw an out_of_bounds exception if an
			  * invalid conversion is specified.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @param new_units	The units to convert into (old units should have already
			  *			been specified using setUnits().)
			  * @throws std::out_of_bounds if an undefined conversion is requested.
			  */
			virtual void convertUnits( ATMOSPHERIC_QUANTITY quantity, UNITS_TYPE new_units ) = 0;
		
			/**
			  * Retrieves an atmospheric quantity at a single altitude.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @param z		The altitude, in internal units
			  * @return		The requested quantity
			  */
			virtual double get( ATMOSPHERIC_QUANTITY quantity, double z ) = 0;
		
			/**
			  * Retrieves an atmospheric quantity at a single altitude, with a temporary
			  * units conversion.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @param z		The altitude, in internal units
			  * @param units	The units to return in
			  * @return		The requested quantity
			  * @throws out_of_bounds if the unit conversion is undefined.
			  */
			virtual double get( ATMOSPHERIC_QUANTITY quantity, double z, UNITS_TYPE new_units ) = 0;
		
			/**
			  * Retrieves an atmospheric quantity at a vector of altitudes.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @param n_z		The number of points to retrieve
			  * @param z		The altitudes, in internal units
			  * @param target	An allocated array for the values to be returned
			  */
			virtual void get( ATMOSPHERIC_QUANTITY quantity, unsigned int n_z, const double *z,
					  double *target ) = 0;
					  
			/**
			  * Retrieves an atmospheric quantity at a vector of altitudes, with a temporary
			  * unit conversion.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @param n_z		The number of points to retrieve
			  * @param z		The altitudes, in internal units
			  * @param units	The units to return in
			  * @param target	An allocated array for the values to be returned
			  */
			virtual void get( ATMOSPHERIC_QUANTITY quantity, unsigned int n_z, const double *z, 
					  UNITS_TYPE new_units, double *target ) = 0;
		
			/**
			  * Notifies the user if a quantity exists in the profile.
			  *
			  * @param quantity	The atmospheric quantity (e.g. temperature)
			  * @return		true if the quantity exists in the profile, false otherwise
			  */
			virtual bool has( ATMOSPHERIC_QUANTITY quantity ) = 0;
		
			
			
			
			
			/*
			FOLLOWING FUNCTIONS ARE DEPRECATED
			*/
			
			/**
			  * returns whether the specification includes vertical winds. 
			  */
			[[deprecated("use has(QUANTITY_WIND_SPEED_VERTICAL) instead")]]
			virtual bool hasW();
		
			/**
			  * returns whether the specification includes pressure. 
			  */
			[[deprecated("use has(QUANTITY_AIR_PRESSURE) instead")]]
			virtual bool hasP();
		
			/**
			  * returns whether the specification includes density. 
			  */
			[[deprecated("use has(QUANTITY_AIR_DENSITY) instead")]]
			virtual bool hasRho();
			
			/**
			  * Effective sound speed.  Returns the effective (i.e. advected) sound speed 
			  * in the indicated direction.
			  * @param z The altitude for which to return the sound speed, in km relative to MSL.
			  * @param phi The direction in which to determine the effective soundspeed, in degrees clockwise from North
			  * @return The effective sound speed along azimuth phi at point (x,y,z), in km/s.
			  */
			[[deprecated("use setPropagationDirection(phi);get(QUANTITY_EFFECTIVE_SOUND_SPEED,z) instead")]]
			virtual double ceff( double z, double phi );

			/**
			  * Static sound speed.  Returns the  sound speed, not including wind effects.  Sound speed is calculated as
			  * sqrt( gamma * P / rho ) if P and rho are available, sqrt( gamma * R * T ) otherwise.
			  * @param z The altitude for which to return the sound speed, in km relative to MSL.
			  * @return The static sound speed at point (x,y,z), in km/s.
			  */
			[[deprecated("use get(QUANTITY_STATIC_SOUND_SPEED,z) instead")]]
			virtual double c0( double z );

			/**
			  * Wind Speed.  Returns the magnitude of the wind vector.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The scalar wind speed at point (x,y,z), in km/s.
			  */
			[[deprecated("use get(QUANTITY_WIND_SPEED,z) instead")]]
			virtual double wspeed( double z );

			/**
			  * Wind Direction.  Returns the 2-D geographic direction of the wind vector.
			  * @param z The altitude for which to return the wind direction, in km relative to MSL.
			  * @return The direction of the projection of the wind vector into the X-Y plane, in degrees clockwise from North.
			  */
			[[deprecated("use get(QUANTITY_WIND_DIRECTION,z) instead")]]
			virtual double wdirection( double z );

			/**
			  * Wind Component.  Returns the magnitude of the wind vector projected first into the X-Y plane, then in the specified direction
			  * @param z The altitude for which to return the wind direction, in km relative to MSL.
			  * @param phi The direction, in degrees clockwise from North, into which to project the wind vector
			  * @return The magnitude of the projection of the horizonal wind vector into the phi direction, in km/s.
			  */
			[[deprecated("use setPropagationDirection(phi);get(QUANTITY_WIND_COMPONENT,z) instead")]]
			virtual double wcomponent( double z, double phi );
			
			

			// Spatial derivatives: temperature
			[[deprecated("use get_dZ( QUANTITY_TEMPERATURE, z ) instead")]]
			virtual double dtdz(  double z );
			[[deprecated("use get_ddZ( QUANTITY_TEMPERATURE, z ) instead")]]
			virtual double ddtdzdz(  double z );

			// Spatial derivatives: effective sound speed
			[[deprecated("use setPropagationDirection(phi);get_dZ(QUANTITY_EFFECTIVE_SOUND_SPEED,z) instead")]]
			virtual double dceffdz(  double z, double phi );
			[[deprecated("use setPropagationDirection(phi);get_ddZ(QUANTITY_EFFECTIVE_SOUND_SPEED,z) instead")]]
			virtual double ddceffdzdz(  double z, double phi );

			// Spatial derivatives: thermodynamic sound speed
			[[deprecated("use get_dZ( QUANTITY_STATIC_SOUND_SPEED, z ) instead")]]
			virtual double dc0dz(  double z );
			[[deprecated("use get_ddZ( QUANTITY_STATIC_SOUND_SPEED, z ) instead")]]
			virtual double ddc0dzdz(  double z );

			// Spatial derivatives: zonal (E-W) winds
			[[deprecated("use get_dZ( QUANTITY_WIND_SPEED_WEST_TO_EAST, z ) instead")]]
			virtual double dudz(  double z );
			[[deprecated("use get_ddZ( QUANTITY_WIND_SPEED_WEST_TO_EAST, z ) instead")]]
			virtual double ddudzdz(  double z );

			// Spatial derivatives: meridional (N-S) winds
			[[deprecated("use get_dZ( QUANTITY_WIND_SPEED_SOUTH_TO_NORTH, z ) instead")]]
			virtual double dvdz(  double z );
			[[deprecated("use get_ddZ( QUANTITY_WIND_SPEED_SOUTH_TO_NORTH, z ) instead")]]
			virtual double ddvdzdz(  double z );

			// Spatial derivatives: vertical winds
			[[deprecated("use get_dZ( QUANTITY_WIND_SPEED_VERTICAL, z ) instead")]]
			virtual double dwdz(  double z );
			[[deprecated("use get_ddZ( QUANTITY_WIND_SPEED_VERTICAL, z ) instead")]]
			virtual double ddwdzdz(  double z );

			// Spatial derivatives: pressure
			[[deprecated("use get_dZ( QUANTITY_AIR_PRESSURE, z ) instead")]]
			virtual double dpdz(  double z );
			[[deprecated("use get_ddZ( QUANTITY_AIR_PRESSURE, z ) instead")]]
			virtual double ddpdzdz(  double z );
		
			// Spatial derivatives: pressure
			[[deprecated("use get_dZ( QUANTITY_AIR_DENSITY, z ) instead")]]
			virtual double drhodz(  double z );
			[[deprecated("use get_ddZ( QUANTITY_AIR_DENSITY, z ) instead")]]
			virtual double ddrhodzdz(  double z );
			
			/**
			  * Minimum valid altitude.  Returns the minimum altitude for which the specification can be queried for a
			  * valid value.
			  * @return The minimum valid altitude, in km.
			  */
			[[deprecated("use getMinimumAltitude() instead")]]
			virtual double z0();
			

			/**
			  * Temperature.  Returns the Kelvin temperature at the indicated point.
			  * @param z The altitude for which to return the temperature, in km relative to MSL.
			  * @return The temperature at point (x,y,z), in Kelvin.
			  */
			[[deprecated("use get(QUANTITY_TEMPERATURE,z) instead")]]
			virtual double t( double z );
			
			/**
			  * Zonal wind.  Returns the zonal (positive to East) wind speed at the indicated point.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The zonal wind speed at point (x,y,z), in km/s.
			  */
			[[deprecated("use get(QUANTITY_WIND_SPEED_WEST_TO_EAST,z) instead")]]
			virtual double u( double z );

			/**
			  * Meridional wind.  Returns the meridional (positive to North) wind speed at the indicated point.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The meridional wind speed at point (x,y,z), in km/s.
			  */
			[[deprecated("use get(QUANTITY_WIND_SPEED_SOUTH_TO_NORTH,z) instead")]]
			virtual double v( double z );

			/**
			  * Vertical wind.  Returns the vertical (positive up) wind speed at the indicated point.
			  * @param z The altitude for which to return the wind speed, in km relative to MSL.
			  * @return The vertical wind speed at point (x,y,z), in km/s.
			  */
			[[deprecated("use get(QUANTITY_WIND_SPEED_VERTICAL,z) instead")]]
			virtual double w( double z );
			
			/**
			 * Pressure.  Returns the atmospheric pressure at the indicated point.
			 * @param z The altitude for which to return the pressure, in km relative to MSL.
			 * @return The pressure in hPa
			 */
			[[deprecated("use get(QUANTITY_AIR_PRESSURE,z) instead")]]
			virtual double p( double z );

			/**
			  * Density.  Returns the air density at the indicated point.
			  * @param z The altitude for which to return the density, in km relative to MSL.
			  * @return The density at point (x,y,z), in g/cm3.
			  */
			[[deprecated("use get(QUANTITY_AIR_DENSITY,z) instead")]]
			virtual double rho( double z );

	};   // class AtmosphericProfile

}   // namespace NCPA




#endif    // #ifndef __ATMOSPHERIC_PROFILE_H__
