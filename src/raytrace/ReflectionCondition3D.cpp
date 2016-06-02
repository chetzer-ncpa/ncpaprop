/**
 * Handles the 3D reflection condition
 * @version 1.0.0
 * @author Claus Hetzer, claus@olemiss.edu
 * @date 2012-04-02
 * 
 * Changelog:
 */

#include "ReflectionCondition3D.h"
#include <stdexcept>
#include <cmath>


NCPA::ReflectionCondition3D::ReflectionCondition3D( NCPA::AtmosphericSpecification *s, double az, unsigned int maxb ) {
	propAzimuth = az;
	spec = s;
	maxbounces = maxb;
	bounces = 0;
	triggered_ = false;
	message = "Maximum number of skips reached.";
}

void NCPA::ReflectionCondition3D::reset() {
	bounces = 0;
	triggered_ = false;
}

bool NCPA::ReflectionCondition3D::triggered() const {
	return triggered_;
}


bool NCPA::ReflectionCondition3D::shouldBreak(double** soln, int k) {
	//double cart_angle = NCPA::deg2rad( 90 - this->propAzimuth );
	double x = soln[k][0];
	double y = soln[k][1];
	double z0 = spec->z0( x, y );
	
	if (soln[k][2] < z0) {
		if (k < 2) {
			std::runtime_error e( "Not enough steps yet taken to compute reflection conditions.  This usually means something has gone wrong, because you shouldn't be reflecting yet." );
			throw e;
		}
		
		double conditions[ 18 ];
		double delta_z = soln[k-1][2] - z0;
		double nu1_approx, nu2_approx, nu3_approx, theta_ref, phi_ref;
		if (spec->stratified()) {
			nu1_approx = soln[0][3];
			nu2_approx = soln[0][4];
			nu3_approx = -soln[0][5];
		} else {
			nu1_approx = soln[k-1][3] 
				- (soln[k][3] - soln[k-1][3])/(soln[k][2] - soln[k-1][2])*delta_z
				+ (soln[k][3] + soln[k-2][3] - 2.0*soln[k-1][3])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
			nu2_approx = soln[k-1][4] 
				- (soln[k][4] - soln[k-1][4])/(soln[k][2] - soln[k-1][2])*delta_z
				+ (soln[k][4] + soln[k-2][4] - 2.0*soln[k-1][4])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
			nu3_approx = soln[k-1][5] 
				- (soln[k][5] - soln[k-1][5])/(soln[k][2] - soln[k-1][2])*delta_z
				+ (soln[k][5] + soln[k-2][5] - 2.0*soln[k-1][5])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
		}
		theta_ref = -std::asin(nu3_approx);
		phi_ref = std::atan2(nu2_approx,nu1_approx);

		// determine x and y at the reflection point using a Taylor approximation
		conditions[0] = soln[k-1][0] 
			- (soln[k][0] - soln[k-1][0])/(soln[k][2] - soln[k-1][2])*delta_z 
			+ (soln[k][0] + soln[k-2][0] - 2.0*soln[k-1][0])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
		conditions[1] = soln[k-1][1] 
			- (soln[k][1] - soln[k-1][1])/(soln[k][2] - soln[k-1][2])*delta_z 
			+ (soln[k][1] + soln[k-2][1] - 2.0*soln[k-1][1])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
		conditions[2] = z0;
		
		// dx/dtheta and dy/dtheta are continuous across the reflection, while dz/dtheta is continuous with a sign change
		conditions[6] = soln[k-1][6] 
			- (soln[k][6] - soln[k-1][6])/(soln[k][2] - soln[k-1][2])*delta_z
			+ (soln[k][6] + soln[k-2][6] - 2.0*soln[k-1][6])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
		conditions[7] = soln[k-1][7] 
			- (soln[k][7] - soln[k-1][7])/(soln[k][2] - soln[k-1][2])*delta_z
			+ (soln[k][7] + soln[k-2][7] - 2.0*soln[k-1][7])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
		conditions[8] = -(
			soln[k-1][8] 
			- (soln[k][8] - soln[k-1][8])/(soln[k][2] - soln[k-1][2])*delta_z
			+ (soln[k][8] + soln[k-2][8] - 2.0*soln[k-1][8])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z
			);

		// dx/dphi and dy/dphi are continuous across the reflection, while dz/dphi is continuous with a sign change
		conditions[12] = soln[k-1][12] 
			- (soln[k][12] - soln[k-1][12])/(soln[k][2] - soln[k-1][2])*delta_z
			+ (soln[k][12] + soln[k-2][12] - 2.0*soln[k-1][12])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
		conditions[13] = soln[k-1][13] 
			- (soln[k][13] - soln[k-1][13])/(soln[k][2] - soln[k-1][2])*delta_z
			+ (soln[k][13] + soln[k-2][13] - 2.0*soln[k-1][13])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z;
		conditions[14] = -(
			soln[k-1][14] 
			- (soln[k][14] - soln[k-1][14])/(soln[k][2] - soln[k-1][2])*delta_z
			+ (soln[k][14] + soln[k-2][14] - 2.0*soln[k-1][14])/pow(soln[k][2] - soln[k-1][2],2.0)*delta_z*delta_z
			);

		// calculate new nu vectors
		conditions[3] = std::cos(theta_ref)*std::cos(phi_ref);
		conditions[4] = std::cos(theta_ref)*std::sin(phi_ref);
		conditions[5] = std::sin(theta_ref);
		
		// derivatives of nu components wrt theta.  First, precalculate some atmospheric characteristics
		double c = spec->c0(conditions[0],conditions[1],conditions[2]);
		double dcdx = spec->dc0dx(conditions[0],conditions[1],conditions[2]);
		double dcdy = spec->dc0dy(conditions[0],conditions[1],conditions[2]);
		double dcdz = spec->dc0dz(conditions[0],conditions[1],conditions[2]);
		conditions[9] = -std::sin(theta_ref)*std::cos(phi_ref) + dcdx / c * conditions[8]/conditions[5];
		conditions[10] = -std::sin(theta_ref)*std::sin(phi_ref) + dcdy / c * conditions[8]/conditions[5];
		conditions[11] = std::cos(theta_ref) - dcdz / c * conditions[8]/conditions[5];
		
		// derivatives of nu components wrt phi
		conditions[15] = -std::cos(theta_ref)*std::sin(phi_ref) + dcdx/c * conditions[14]/conditions[5];
		conditions[16] = std::cos(theta_ref)*std::cos(phi_ref) + dcdy/c * conditions[14]/conditions[5];
		conditions[17] = -dcdz/c * conditions[14]/conditions[5];	
		
		
		
		
		
		for (unsigned int i = 0; i < 18; i++) {
			soln[k][i] = conditions[i];
		}
		
		if (maxbounces > 0 && ++bounces == maxbounces) {
			triggered_ = true;
			return true;
		}
		
	}
	// shouldn't break at all
	return false;
}

void NCPA::ReflectionCondition3D::setAzimuth(double az) {
	while (az < 0) {
		az += 360;
	}
	while (az >= 360) {
		az -= 360;
	}
	propAzimuth = az;
}
