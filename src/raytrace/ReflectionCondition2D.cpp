#include "ReflectionCondition2D.h"
#include <stdexcept>
#include <cmath>


NCPA::ReflectionCondition2D::ReflectionCondition2D( NCPA::AtmosphericSpecification *s, double az, unsigned int maxb ) {
	propAzimuth = az;
	spec = s;
	maxbounces = maxb;
	bounces = 0;
	message = "Maximum number of skips reached.";
	triggered_ = false;
}

void NCPA::ReflectionCondition2D::reset() {
	bounces = 0;
	triggered_ = false;
}

bool NCPA::ReflectionCondition2D::triggered() const {
	return triggered_;
}

bool NCPA::ReflectionCondition2D::shouldBreak(double** solution, int k) {
	double cart_angle = NCPA::deg2rad( 90 - this->propAzimuth );
	double x = solution[k][0] * std::cos( cart_angle );
	double y = solution[k][0] * std::sin( cart_angle );
	double z0 = spec->z0( x, y );
	
	if (solution[k][1] < z0) {
		if (k < 2) {
			std::runtime_error e( "Not enough steps yet taken to compute reflection conditions.  This usually means something has gone wrong, because you shouldn't be reflecting yet." );
			throw e;
		}
		double conditions[ 6 ];
		
		conditions[0] = solution[k-1][0] + (solution[k][0] - solution[k-1][0])/(solution[k-1][1] - solution[k][1])*(solution[k-1][1] - z0)
				+ (solution[k][0] + solution[k-2][0] - 2*solution[k-1][0])/std::pow(solution[k-1][1] - solution[k][1],2.0)*std::pow(solution[k-1][1]-z0,2.0);
		
		// recalculate the incidence angle
		x = conditions[0] * std::cos( cart_angle );
		y = conditions[0] * std::sin( cart_angle );
		z0 = spec->z0( x, y );
		double ceff0 = spec->ceff( x,y,z0,propAzimuth );
		double theta_new = spec->stratified() ? launchAngle : std::atan2( solution[k-1][1] - z0, conditions[0] - solution[k-1][0] );
		
		conditions[1] = z0;
		conditions[2] = std::sin( theta_new ) / ceff0;
		
		conditions[3] = solution[k-1][3] + (solution[k][3] - solution[k-1][3])/(solution[k-1][1] - solution[k][1])*(solution[k-1][1] - z0)
				+ (solution[k][3] + solution[k-2][3] - 2*solution[k-1][3])/std::pow(solution[k-1][1] - solution[k][1],2.0)*std::pow(solution[k-1][1]-z0,2.0);
		
		conditions[4] = -(solution[k-1][4] + (solution[k][4] - solution[k-1][4])/(solution[k-1][1] - solution[k][1])*(solution[k-1][1] - z0)
				+ (solution[k][4] + solution[k-2][4] - 2*solution[k-1][4])/std::pow(solution[k-1][1] - solution[k][1],2.0)*std::pow(solution[k-1][1]-z0,2.0));
		
		conditions[5] = std::cos( theta_new ) / ceff0 - conditions[4] * spec->dceffdz(x,y,z0,propAzimuth) / ceff0 / ceff0 / std::sin( theta_new );
		
		for (unsigned int i = 0; i < 6; i++) {
			solution[k][i] = conditions[i];
		}
		
		if (maxbounces > 0 && ++bounces == maxbounces) {
			//bounces = 0;
			triggered_ = true;
			return true;
		}
		
	}
	// shouldn't break at all
	return false;
}

void NCPA::ReflectionCondition2D::setLaunchAzimuthDegrees(double az) {
	while (az < 0) {
		az += 360;
	}
	while (az >= 360) {
		az -= 360;
	}
	propAzimuth = az;
	bounces = 0;
}

void NCPA::ReflectionCondition2D::setLaunchAngleRadians(double radians)
{
	launchAngle = radians;
	bounces = 0;
}

unsigned int NCPA::ReflectionCondition2D::countBounces() const {
	return bounces;
}