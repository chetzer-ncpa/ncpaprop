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
	cout << tempC << " " << Units::toStr( UNITS_TEMPERATURE_CELSIUS ) << " = " << tempK << " " << Units::toStr( UNITS_TEMPERATURE_KELVIN ) << endl;
	cout << tempC << " " << Units::toString( UNITS_TEMPERATURE_CELSIUS ) << " = " << tempF << " " << Units::toString( UNITS_TEMPERATURE_FAHRENHEIT ) << endl;
	cout << tempC << " " << Units::toString( UNITS_TEMPERATURE_CELSIUS ) << " = " << tempSAME << " " << Units::toString( UNITS_TEMPERATURE_CELSIUS ) << endl;

	// Kelvin conversions
	tempK = 307.4;
	tempF = Units::convert( tempK, UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_FAHRENHEIT );
	Units::convert( &tempK, 1, UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_CELSIUS, &tempC );
	tempSAME = Units::convert( tempK, UNITS_TEMPERATURE_KELVIN, UNITS_TEMPERATURE_KELVIN );
	cout << tempK << " " << Units::toStr( UNITS_TEMPERATURE_KELVIN ) << " = " << tempF << " " << Units::toStr( UNITS_TEMPERATURE_FAHRENHEIT ) << endl;
	cout << tempK << " " << Units::toString( UNITS_TEMPERATURE_KELVIN ) << " = " << tempC << " " << Units::toString( UNITS_TEMPERATURE_CELSIUS ) << endl;
	cout << tempK << " " << Units::toString( UNITS_TEMPERATURE_KELVIN ) << " = " << tempSAME << " " << Units::toString( UNITS_TEMPERATURE_KELVIN ) << endl;

	// Fahrenheit conversions
	tempF = 32.0;
	tempC = Units::convert( tempF, UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_CELSIUS );
	Units::convert( &tempF, 1, UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_KELVIN, &tempK );
	tempSAME = Units::convert( tempF, UNITS_TEMPERATURE_FAHRENHEIT, UNITS_TEMPERATURE_FAHRENHEIT );
	cout << tempF << " " << Units::toStr( UNITS_TEMPERATURE_FAHRENHEIT ) << " = " << tempC << " " << Units::toStr( UNITS_TEMPERATURE_CELSIUS ) << endl;
	cout << tempF << " " << Units::toString( UNITS_TEMPERATURE_FAHRENHEIT ) << " = " << tempK << " " << Units::toString( UNITS_TEMPERATURE_KELVIN ) << endl;
	cout << tempF << " " << Units::toString( UNITS_TEMPERATURE_FAHRENHEIT ) << " = " << tempSAME << " " << Units::toString( UNITS_TEMPERATURE_FAHRENHEIT ) << endl;

	// Distance conversions
	distM = 15.0;
	distKM = Units::convert( distM, UNITS_DISTANCE_METERS, UNITS_DISTANCE_KILOMETERS );
	cout << distM << " " << Units::toStr( UNITS_DISTANCE_METERS ) << " = " << distKM << " " << Units::toStr( UNITS_DISTANCE_KILOMETERS ) << endl;
	distKM = 42.0;
	distM = Units::convert( distKM, UNITS_DISTANCE_KILOMETERS, UNITS_DISTANCE_METERS );
	cout << distKM << " " << Units::toString( UNITS_DISTANCE_KILOMETERS ) << " = " << distM << " " << Units::toString( UNITS_DISTANCE_METERS ) << endl;

	// Speed conversions
	speedMPS = 0.1;
	speedKMPS = Units::convert( speedMPS, UNITS_SPEED_METERS_PER_SECOND, UNITS_SPEED_KILOMETERS_PER_SECOND );
	cout << speedMPS << " " << Units::toStr( UNITS_SPEED_METERS_PER_SECOND ) << " = " << speedKMPS << " " << Units::toStr( UNITS_SPEED_KILOMETERS_PER_SECOND ) << endl;
	speedKMPS = 12.4;
	speedMPS = Units::convert( speedKMPS, UNITS_SPEED_KILOMETERS_PER_SECOND, UNITS_SPEED_METERS_PER_SECOND );
	cout << speedKMPS << " " << Units::toStr( UNITS_SPEED_KILOMETERS_PER_SECOND ) << " = " << speedMPS << " " << Units::toStr( UNITS_SPEED_METERS_PER_SECOND ) << endl;
	
	cout << endl << endl;

	Units::list_recognized_strings( cout );

	return 1;
}