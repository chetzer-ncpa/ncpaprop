#include "units.h"
#include <iostream>

using namespace std;
using namespace NCPA;

int main( int argc, char **argv ) {
	
	// test individual conversions
	double tempC, tempK, tempF, tempSAME, distM, distKM, speedMPS, speedKMPS, speedSAME, 
		pressPA, pressMBAR, pressSAME, densityKGPM3, densityGPCM3, densitySAME, 
		directionCWFN, directionCCWFE, directionSAME, angleDEG, angleRAD, angleSAME;

	// celsius conversions
	tempC = 20.0;
	tempK = Units::convert( tempC, UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_KELVIN );
	Units::convert( &tempC, 1, UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_FAHRENHEIT, &tempF );
	tempSAME = Units::convert( tempC, UNITS_TEMPERATURE_CELSIUS, UNITS_TEMPERATURE_CELSIUS );
	cout << tempC << " " << toStr( UNITS_TEMPERATURE_CELSIUS ) << " = " << tempK << " " << toStr( UNITS_TEMPERATURE_KELVIN ) << endl;
	cout << tempC << " " << toString( UNITS_TEMPERATURE_CELSIUS ) << " = " << tempF << " " << toString( UNITS_TEMPERATURE_FAHRENHEIT ) << endl;
	cout << tempC << " " << toString( UNITS_TEMPERATURE_CELSIUS ) << " = " << tempSAME << " " << toString( UNITS_TEMPERATURE_CELSIUS ) << endl;

	// Kelvin conversions
	tempK = 307.4;
	tempF = Units::convert( tempK, UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_FAHRENHEIT );
	Units::convert( &tempK, 1, UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_CELSIUS, &tempC );
	tempSAME = Units::convert( tempK, UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_KELVIN );
	cout << tempK << " " << toStr( UNITS_TEMPERATURE_KELVIN ) << " = " << tempF << " " << toStr( UNITS_TEMPERATURE_FAHRENHEIT ) << endl;
	cout << tempK << " " << toString( UNITS_TEMPERATURE_KELVIN ) << " = " << tempC << " " << toString( UNITS_TEMPERATURE_CELSIUS ) << endl;
	cout << tempK << " " << toString( UNITS_TEMPERATURE_KELVIN ) << " = " << tempSAME << " " << toString( UNITS_TEMPERATURE_KELVIN ) << endl;

	// Fahrenheit conversions
	tempF = 32.0;
	tempC = Units::convert( tempF, UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_CELSIUS );
	Units::convert( &tempF, 1, UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_KELVIN, &tempK );
	tempSAME = Units::convert( tempF, UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_FAHRENHEIT );
	cout << tempF << " " << toStr( UNITS_TEMPERATURE_FAHRENHEIT ) << " = " << tempC << " " << toStr( UNITS_TEMPERATURE_CELSIUS ) << endl;
	cout << tempF << " " << toString( UNITS_TEMPERATURE_FAHRENHEIT ) << " = " << tempK << " " << toString( UNITS_TEMPERATURE_KELVIN ) << endl;
	cout << tempF << " " << toString( UNITS_TEMPERATURE_FAHRENHEIT ) << " = " << tempSAME << " " << toString( UNITS_TEMPERATURE_FAHRENHEIT ) << endl;

	// Distance conversions
	distM = 15.0;
	distKM = Units::convert( distM, UNITS_DISTANCE_METERS, UNITS_DISTANCE_KILOMETERS );
	cout << distM << " " << toStr( UNITS_DISTANCE_METERS ) << " = " << distKM << " " << toStr( UNITS_DISTANCE_KILOMETERS ) << endl;
	distKM = 42.0;
	distM = Units::convert( distKM, UNITS_DISTANCE_KILOMETERS, UNITS_DISTANCE_METERS );
	cout << distKM << " " << toString( UNITS_DISTANCE_KILOMETERS ) << " = " << distM << " " << toString( UNITS_DISTANCE_METERS ) << endl;

	// Speed conversions
	speedMPS = 0.1;
	speedKMPS = Units::convert( speedMPS, UNITS_SPEED_METERS_PER_SECOND, UNITS_SPEED_KILOMETERS_PER_SECOND );
	cout << speedMPS << " " << toStr( UNITS_SPEED_METERS_PER_SECOND ) << " = " << speedKMPS << " " << toStr( UNITS_SPEED_KILOMETERS_PER_SECOND ) << endl;
	speedKMPS = 12.4;
	speedMPS = Units::convert( speedKMPS, UNITS_SPEED_KILOMETERS_PER_SECOND, UNITS_SPEED_METERS_PER_SECOND );
	cout << speedKMPS << " " << toStr( UNITS_SPEED_KILOMETERS_PER_SECOND ) << " = " << speedMPS << " " << toStr( UNITS_SPEED_METERS_PER_SECOND ) << endl;
	
	

	return 1;
}