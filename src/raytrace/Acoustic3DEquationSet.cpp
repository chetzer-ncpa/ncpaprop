#include "Acoustic3DEquationSet.h"
#include "AtmosphericSpecification.h"
#include <cmath>
#include "util.h"
#include <iostream>

#ifndef Pi
#define Pi 3.14159
#endif

NCPA::Acoustic3DEquationSet::Acoustic3DEquationSet( NCPA::AtmosphericSpecification *p, 
		double takeoffangle, double propagation_azimuth,
		bool rangeDependent ) {
	profile = p;
	azimuth = propagation_azimuth;
	takeoff = takeoffangle;
	rangeDependent_ = rangeDependent;
}

int NCPA::Acoustic3DEquationSet::numberOfEquations() const { return 18; }

void NCPA::Acoustic3DEquationSet::setupInitialConditions( double *initialConditions, double zmin ) {
	double az = deg2rad( 90 - azimuth );
	initialConditions[ 0 ] = 0;
	initialConditions[ 1 ] = 0;
	initialConditions[ 2 ] = zmin;
	initialConditions[ 3 ] = std::cos( takeoff ) * std::cos( az );
	initialConditions[ 4 ] = std::cos( takeoff ) * std::sin( az );
	initialConditions[ 5 ] = std::sin( takeoff );
	initialConditions[ 6 ] = 0;
	initialConditions[ 7 ] = 0;
	initialConditions[ 8 ] = 0;
	initialConditions[ 9 ] = -std::sin( takeoff ) * std::cos( az );
	initialConditions[ 10 ] = -std::sin( takeoff ) * std::sin( az );
	initialConditions[ 11 ] = std::cos( takeoff );
	initialConditions[ 12 ] = 0;
	initialConditions[ 13 ] = 0;
	initialConditions[ 14 ] = 0;
	initialConditions[ 15 ] = -std::cos( takeoff ) * std::sin( az );
	initialConditions[ 16 ] = std::cos( takeoff ) * std::cos( az );
	initialConditions[ 17 ] = 0;
}
	
// Magnitude of the nu vector
double NCPA::Acoustic3DEquationSet::mag_nu(double r1, double r2, double r3, double nu1, double nu2, double nu3) const {
	return profile->c0(0.0,0.0,profile->z0(0.0,0.0))/profile->c0(r1,r2,r3)
		* fabs( 1.0 - (nu1*profile->u(r1,r2,r3) + nu2*profile->v(r1,r2,r3) + nu3*profile->w(r1,r2,r3)) / profile->c0(0.0,0.0,profile->z0(0.0,0.0)) );
}

double NCPA::Acoustic3DEquationSet::mag_nu(double r1, double r2, double r3, double nu1, double nu2, double nu3,
					   double c0, double c0_0, double u, double v, double w ) const {
	return c0_0 / c0 * std::fabs( 1.0 - (nu1 * u + nu2 * v + nu3 * w ) / c0_0 );
}

// Group velocity vector component, x direction
double NCPA::Acoustic3DEquationSet::cg1(double r1, double r2, double r3, double nu1, double nu2, double nu3) const {
        return profile->u(r1,r2,r3) + profile->c0(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3)*nu1;
}

double NCPA::Acoustic3DEquationSet::cg1(double r1, double r2, double r3, double nu1, double nu2, double nu3, 
					double c0, double c0_0, double u, double v, double w ) const {
        return u + c0 / mag_nu(r1,r2,r3,nu1,nu2,nu3, c0, c0_0, u, v, w )*nu1;
}


// Group velocity vector component, y direction
double NCPA::Acoustic3DEquationSet::cg2(double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double c0, double c0_0, double u, double v, double w)  const {
        return v + c0/mag_nu(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w)*nu2;
}

double NCPA::Acoustic3DEquationSet::cg2(double r1, double r2, double r3, double nu1, double nu2, double nu3)  const {
        return profile->v(r1,r2,r3) + profile->c0(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3)*nu2;
}


// Group velocity vector component, z direction
double NCPA::Acoustic3DEquationSet::cg3(double r1, double r2, double r3, double nu1, double nu2, double nu3) const {
        return profile->w(r1,r2,r3) + profile->c0(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3)*nu3;
}

double NCPA::Acoustic3DEquationSet::cg3(double r1, double r2, double r3, double nu1, double nu2, double nu3,
					double c0, double c0_0, double u, double v, double w) const {
        return w + c0/mag_nu(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w)*nu3;
}



// Group velocity vector magnitude
double NCPA::Acoustic3DEquationSet::cg_mag(double r1, double r2, double r3, double nu1, double nu2, double nu3) const {
        double cg1_sqr = pow(fabs(profile->u(r1,r2,r3) + profile->c0(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3)*nu1), 2.0);
        double cg2_sqr = pow(fabs(profile->v(r1,r2,r3) + profile->c0(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3)*nu2), 2.0);
        double cg3_sqr = pow(fabs(profile->w(r1,r2,r3) + profile->c0(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3)*nu3), 2.0);

        return sqrt(cg1_sqr + cg2_sqr + cg3_sqr);
}

double NCPA::Acoustic3DEquationSet::cg_mag(double r1, double r2, double r3, double nu1, double nu2, double nu3,
					   double c0, double c0_0, double u, double v, double w) const {
        double cg1_sqr = pow(fabs(u + c0/mag_nu(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w)*nu1), 2.0);
        double cg2_sqr = pow(fabs(v + c0/mag_nu(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w)*nu2), 2.0);
        double cg3_sqr = pow(fabs(w + c0/mag_nu(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w)*nu3), 2.0);

        return sqrt(cg1_sqr + cg2_sqr + cg3_sqr);
}



// Angular derivative: c
double NCPA::Acoustic3DEquationSet::dc_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle ) const {
        return dr1_angle*profile->dc0dx(r1,r2,r3) 
		+ dr2_angle*profile->dc0dy(r1,r2,r3)
		+ dr3_angle*profile->dc0dz(r1,r2,r3);
}

double NCPA::Acoustic3DEquationSet::dc_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					     double dc0dx, double dc0dy, double dc0dz ) const {
        return dr1_angle*dc0dx + dr2_angle*dc0dy + dr3_angle*dc0dz;
}

// Angular derivative: u
double NCPA::Acoustic3DEquationSet::du_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle) const {
        return dr1_angle*profile->dudx(r1,r2,r3) 
		+ dr2_angle*profile->dudy(r1,r2,r3)
		+ dr3_angle*profile->dudz(r1,r2,r3);
}

double NCPA::Acoustic3DEquationSet::du_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					     double dudx, double dudy, double dudz) const {
        return dr1_angle*dudx + dr2_angle*dudy + dr3_angle*dudz;
}


// Angular derivative: v
double NCPA::Acoustic3DEquationSet::dv_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle) const {
        return dr1_angle*profile->dvdx(r1,r2,r3) 
		+ dr2_angle*profile->dvdy(r1,r2,r3)
		+ dr3_angle*profile->dvdz(r1,r2,r3);
}

double NCPA::Acoustic3DEquationSet::dv_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					     double dvdx, double dvdy, double dvdz) const {
        return dr1_angle*dvdx + dr2_angle*dvdy + dr3_angle*dvdz;
}



// Angular derivative: w
double NCPA::Acoustic3DEquationSet::dw_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle) const {
        return dr1_angle*profile->dwdx(r1,r2,r3) 
		+ dr2_angle*profile->dwdy(r1,r2,r3)
		+ dr3_angle*profile->dwdz(r1,r2,r3);
}

double NCPA::Acoustic3DEquationSet::dw_angle(double r1, double r2, double r3, double dr1_angle, double dr2_angle, double dr3_angle,
					     double dwdx, double dwdy, double dwdz) const {
        return dr1_angle*dwdx + dr2_angle*dwdy + dr3_angle*dwdz;
}


// Angular derivative: group velocity
double NCPA::Acoustic3DEquationSet::dCg_angle(double r1, double r2, double r3, double nu1, double nu2, double nu3, 
					      double dr1_angle, double dr2_angle, double dr3_angle, 
					      double dnu1_angle, double dnu2_angle, double dnu3_angle) const {
						      
        double nu_dot_dnu = nu1*dnu1_angle + nu2*dnu2_angle + nu3*dnu3_angle;

        double sub11 = profile->u(r1,r2,r3)/profile->c0(r1,r2,r3) + nu1/mag_nu(r1,r2,r3,nu1,nu2,nu3);
        double sub12 = 1/profile->c0(r1,r2,r3)*du_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle) - profile->u(r1,r2,r3)/pow(profile->c0(r1,r2,r3),2)*dc_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle) + dnu1_angle/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu1/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*nu_dot_dnu;

        double sub21 = profile->v(r1,r2,r3)/profile->c0(r1,r2,r3) + nu2/mag_nu(r1,r2,r3,nu1,nu2,nu3);
        double sub22 = 1/profile->c0(r1,r2,r3)*dv_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle) - profile->v(r1,r2,r3)/pow(profile->c0(r1,r2,r3),2)*dc_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle) + dnu2_angle/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu2/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*nu_dot_dnu;

        double sub31 = profile->w(r1,r2,r3)/profile->c0(r1,r2,r3) + nu3/mag_nu(r1,r2,r3,nu1,nu2,nu3);
        double sub32 = 1/profile->c0(r1,r2,r3)*dw_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle) - profile->w(r1,r2,r3)/pow(profile->c0(r1,r2,r3),2)*dc_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle) + dnu3_angle/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu3/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*nu_dot_dnu;

        return cg_mag(r1,r2,r3,nu1,nu2,nu3)*(1/profile->c0(r1,r2,r3)*dc_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle) + pow(profile->c0(r1,r2,r3)/cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*(sub11*sub12 + sub21*sub22 + sub31*sub32) );
}

double NCPA::Acoustic3DEquationSet::dCg_angle(double r1, double r2, double r3, double nu1, double nu2, double nu3, 
					      double dr1_angle, double dr2_angle, double dr3_angle, 
					      double dnu1_angle, double dnu2_angle, double dnu3_angle,
					      double c0, double c0_0, double u, double v, double w,
					      double dc0dx, double dc0dy, double dc0dz,
					      double dudx, double dudy, double dudz,
					      double dvdx, double dvdy, double dvdz,
					      double dwdx, double dwdy, double dwdz
     					) const {
						      
        double nu_dot_dnu = nu1*dnu1_angle + nu2*dnu2_angle + nu3*dnu3_angle;
	double mag_nu_pc = mag_nu(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w);
	double dc_angle_pc = dc_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle,dc0dx,dc0dy,dc0dz);
	double cg_mag_pc = cg_mag(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w);

        double sub11 = u/c0 + nu1/mag_nu_pc;
        double sub12 = 1/c0*du_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle,dudx,dudy,dudz) 
			- u/pow(c0,2)*dc_angle_pc 
			+ dnu1_angle/mag_nu_pc - nu1/pow(mag_nu_pc,3.0)*nu_dot_dnu;

        double sub21 = v/c0 + nu2/mag_nu_pc;
        double sub22 = 1/c0*dv_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle,dvdx,dvdy,dvdz) 
			- v/pow(c0,2)*dc_angle_pc
			+ dnu2_angle/mag_nu_pc - nu2/pow(mag_nu_pc,3.0)*nu_dot_dnu;

        double sub31 = w/c0 + nu3/mag_nu_pc;
        double sub32 = 1/c0*dw_angle(r1,r2,r3,dr1_angle,dr2_angle,dr3_angle,dwdx,dwdy,dwdz) 
			- w/pow(c0,2)*dc_angle_pc
			+ dnu3_angle/mag_nu_pc - nu3/pow(mag_nu_pc,3.0)*nu_dot_dnu;

        return cg_mag_pc*(1/c0 * dc_angle_pc + pow(c0/cg_mag_pc,2.0)*(sub11*sub12 + sub21*sub22 + sub31*sub32) );
}





/*
double NCPA::Acoustic3DEquationSet::result( double s, double *current_values, int which ) {

	// Identify the current values by name, for ease of reading/debugging
	double r1 = current_values[ 0 ];
	double r2 = current_values[ 1 ];
	double r3 = current_values[ 2 ];
	double nu1 = current_values[ 3 ];
	double nu2 = current_values[ 4 ];
	double nu3 = current_values[ 5 ];
	double dr1t = current_values[ 6 ];
	double dr2t = current_values[ 7 ];
	double dr3t = current_values[ 8 ];
	double dnu1t = current_values[ 9 ];
	double dnu2t = current_values[ 10 ];
	double dnu3t = current_values[ 11 ];
	double dr1p = current_values[ 12 ];
	double dr2p = current_values[ 13 ];
	double dr3p = current_values[ 14 ];
	double dnu1p = current_values[ 15 ];
	double dnu2p = current_values[ 16 ];
	double dnu3p = current_values[ 17 ];

	double sub_term1, sub_term2, sub_term3;

	switch (which) {
		case 0:
			// define source function dr1/ds
			return cg1(r1,r2,r3,nu1,nu2,nu3) / cg_mag(r1,r2,r3,nu1,nu2,nu3);

		case 1:
			// define source function dr2/ds
			return cg2(r1,r2,r3,nu1,nu2,nu3) / cg_mag(r1,r2,r3,nu1,nu2,nu3);

		case 2:
			// define source function dr3/ds
			return cg3(r1,r2,r3,nu1,nu2,nu3) / cg_mag(r1,r2,r3,nu1,nu2,nu3);

		case 3:
			//define source function dnu1/ds
			return -1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)
				*(mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dx(r1,r2,r3) + nu1*profile->dudx(r1,r2,r3) + nu2*profile->dvdx(r1,r2,r3) + nu3*profile->dwdx(r1,r2,r3));

		case 4:
			//define source function dnu2/ds
			return -1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)
				*(mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dy(r1,r2,r3) + nu1*profile->dudy(r1,r2,r3) + nu2*profile->dvdy(r1,r2,r3) + nu3*profile->dwdy(r1,r2,r3) );

		case 5:
			//define source function dnu3/ds
			return -1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)
				*(mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dz(r1,r2,r3) + nu1*profile->dudz(r1,r2,r3) + nu2*profile->dvdz(r1,r2,r3) + nu3*profile->dwdz(r1,r2,r3) );

		case 6:
			// define source function d/ds[ dr1/dtheta ]
			sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg1(r1,r2,r3,nu1,nu2,nu3);
			sub_term2 = du_angle(r1,r2,r3,dr1t,dr2t,dr3t) + nu1/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1t,dr2t,dr3t);
			sub_term3 = profile->c0(r1,r2,r3)*(dnu1t/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu1/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
			return sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 7:
			// define source function d/ds[ dr2/dtheta ]
			sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg2(r1,r2,r3,nu1,nu2,nu3);
			sub_term2 = dv_angle(r1,r2,r3,dr1t,dr2t,dr3t) + nu2/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1t,dr2t,dr3t);
			sub_term3 = profile->c0(r1,r2,r3)*(dnu2t/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu2/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
			return  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 8:
			// define source function d/ds[ dr3/dtheta ]
			sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg3(r1,r2,r3,nu1,nu2,nu3);
			sub_term2 = dw_angle(r1,r2,r3,dr1t,dr2t,dr3t) + nu3/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1t,dr2t,dr3t);
			sub_term3 = profile->c0(r1,r2,r3)*(dnu3t/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu3/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
			return sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 9:
			// define source function d/ds[ dnu1/dtheta ]
			sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dx(r1,r2,r3) 
				+ nu1*profile->dudx(r1,r2,r3) + nu2*profile->dvdx(r1,r2,r3)
				+ nu3*profile->dwdx(r1,r2,r3);
			sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*profile->dc0dx(r1,r2,r3)
				/ mag_nu(r1,r2,r3,nu1,nu2,nu3) 
				+ dnu1t*profile->dudx(r1,r2,r3) + dnu2t*profile->dvdx(r1,r2,r3)
				+ dnu3t*profile->dwdx(r1,r2,r3);
			sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
				* (dr1t*profile->ddc0dxdx(r1,r2,r3) + dr2t*profile->ddc0dxdy(r1,r2,r3) + dr3t*profile->ddc0dxdz(r1,r2,r3)) 
				+ nu1*(dr1t*profile->ddudxdx(r1,r2,r3) + dr2t*profile->ddudxdy(r1,r2,r3) + dr3t*profile->ddudxdz(r1,r2,r3)) 
				+ nu2*(dr1t*profile->ddvdxdx(r1,r2,r3) + dr2t*profile->ddvdxdy(r1,r2,r3) + dr3t*profile->ddvdxdz(r1,r2,r3))
				+ nu3*(dr1t*profile->ddwdxdx(r1,r2,r3) + dr2t*profile->ddwdxdy(r1,r2,r3) + dr3t*profile->ddwdxdz(r1,r2,r3));
			return  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 10:
			// define source function d/ds[ dnu2/dtheta ]
			sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dy(r1,r2,r3) 
				+ nu1*profile->dudy(r1,r2,r3) + nu2*profile->dvdy(r1,r2,r3)
				+ nu3*profile->dwdy(r1,r2,r3);
			sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*profile->dc0dy(r1,r2,r3)
				/ mag_nu(r1,r2,r3,nu1,nu2,nu3) 
				+ dnu1t*profile->dudy(r1,r2,r3) + dnu2t*profile->dvdy(r1,r2,r3)
				+ dnu3t*profile->dwdy(r1,r2,r3);
			sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
				* (dr1t*profile->ddc0dxdy(r1,r2,r3) + dr2t*profile->ddc0dydy(r1,r2,r3) + dr3t*profile->ddc0dydz(r1,r2,r3)) 
				+ nu1*(dr1t*profile->ddudxdy(r1,r2,r3) + dr2t*profile->ddudydy(r1,r2,r3) + dr3t*profile->ddudydz(r1,r2,r3)) 
				+ nu2*(dr1t*profile->ddvdxdy(r1,r2,r3) + dr2t*profile->ddvdydy(r1,r2,r3) + dr3t*profile->ddvdydz(r1,r2,r3))
				+ nu3*(dr1t*profile->ddwdxdy(r1,r2,r3) + dr2t*profile->ddwdydy(r1,r2,r3) + dr3t*profile->ddwdydz(r1,r2,r3));
			return  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 11:
			// define source function d/ds[ dnu3/dtheta ]
			sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dz(r1,r2,r3) 
				+ nu1*profile->dudz(r1,r2,r3) + nu2*profile->dvdz(r1,r2,r3)
				+ nu3*profile->dwdz(r1,r2,r3);
			sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*profile->dc0dz(r1,r2,r3)
				/ mag_nu(r1,r2,r3,nu1,nu2,nu3) 
				+ dnu1t*profile->dudz(r1,r2,r3) + dnu2t*profile->dvdz(r1,r2,r3)
				+ dnu3t*profile->dwdz(r1,r2,r3);
			sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
				* (dr1t*profile->ddc0dxdz(r1,r2,r3) + dr2t*profile->ddc0dydz(r1,r2,r3) + dr3t*profile->ddc0dzdz(r1,r2,r3)) 
				+ nu1*(dr1t*profile->ddudxdz(r1,r2,r3) + dr2t*profile->ddudydz(r1,r2,r3) + dr3t*profile->ddudzdz(r1,r2,r3)) 
				+ nu2*(dr1t*profile->ddvdxdz(r1,r2,r3) + dr2t*profile->ddvdydz(r1,r2,r3) + dr3t*profile->ddvdzdz(r1,r2,r3))
				+ nu3*(dr1t*profile->ddwdxdz(r1,r2,r3) + dr2t*profile->ddwdydz(r1,r2,r3) + dr3t*profile->ddwdzdz(r1,r2,r3));
			return  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 12:
			// define source function d/ds[ dr1/dphi ]
			sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg1(r1,r2,r3,nu1,nu2,nu3);
			sub_term2 = du_angle(r1,r2,r3,dr1p,dr2p,dr3p) + nu1/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1p,dr2p,dr3p);
			sub_term3 = profile->c0(r1,r2,r3)*(dnu1p/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu1/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
			return  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 13:
			// define source function d/ds[ dr2/dphi ]
			sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg2(r1,r2,r3,nu1,nu2,nu3);
			sub_term2 = dv_angle(r1,r2,r3,dr1p,dr2p,dr3p) + nu2/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1p,dr2p,dr3p);
			sub_term3 = profile->c0(r1,r2,r3)*(dnu2p/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu2/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
			return  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 14:
			// define source function d/ds[ dr3/dphi ]
			sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg3(r1,r2,r3,nu1,nu2,nu3);
			sub_term2 = dw_angle(r1,r2,r3,dr1p,dr2p,dr3p) + nu3/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1p,dr2p,dr3p);
			sub_term3 = profile->c0(r1,r2,r3)*(dnu3p/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu3/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
			return  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 15:
			// define source function d/ds[ dnu1/dphi ]
			sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dx(r1,r2,r3) 
				+ nu1*profile->dudx(r1,r2,r3) 
				+ nu2*profile->dvdx(r1,r2,r3)
				+ nu3*profile->dwdx(r1,r2,r3);
			sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)
				* profile->dc0dx(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3) 
				+ dnu1p*profile->dudx(r1,r2,r3) 
				+ dnu2p*profile->dvdx(r1,r2,r3)
				+ dnu3p*profile->dwdx(r1,r2,r3);
			sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
				* (
					dr1p*profile->ddc0dxdx(r1,r2,r3) 
					+ dr2p*profile->ddc0dxdy(r1,r2,r3) 
					+ dr3p*profile->ddc0dxdz(r1,r2,r3)
				  ) 
				+ nu1 * (
					dr1p*profile->ddudxdx(r1,r2,r3) 
					+ dr2p*profile->ddudxdy(r1,r2,r3) 
					+ dr3p*profile->ddudxdz(r1,r2,r3)
					) 
				+ nu2 * (
					dr1p*profile->ddvdxdx(r1,r2,r3) 
					+ dr2p*profile->ddvdxdy(r1,r2,r3) 
					+ dr3p*profile->ddvdxdz(r1,r2,r3)
					)
				+ nu3 * (
					dr1p*profile->ddwdxdx(r1,r2,r3) 
					+ dr2p*profile->ddwdxdy(r1,r2,r3) 
					+ dr3p*profile->ddwdxdz(r1,r2,r3)
					);
			return  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 16:
			// define source function d/ds[ dnu2/dphi ]
			sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dy(r1,r2,r3) 
				+ nu1*profile->dudy(r1,r2,r3) 
				+ nu2*profile->dvdy(r1,r2,r3)
				+ nu3*profile->dwdy(r1,r2,r3);
			sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)
				* profile->dc0dy(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3) 
				+ dnu1p*profile->dudy(r1,r2,r3) 
				+ dnu2p*profile->dvdy(r1,r2,r3)
				+ dnu3p*profile->dwdy(r1,r2,r3);
			sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
				* (
					dr1p*profile->ddc0dxdy(r1,r2,r3) 
					+ dr2p*profile->ddc0dydy(r1,r2,r3) 
					+ dr3p*profile->ddc0dydz(r1,r2,r3)
				  ) 
				+ nu1 * (
					dr1p*profile->ddudxdy(r1,r2,r3) 
					+ dr2p*profile->ddudydy(r1,r2,r3) 
					+ dr3p*profile->ddudydz(r1,r2,r3)
					) 
				+ nu2 * (
					dr1p*profile->ddvdxdy(r1,r2,r3) 
					+ dr2p*profile->ddvdydy(r1,r2,r3) 
					+ dr3p*profile->ddvdydz(r1,r2,r3)
					)
				+ nu3 * (
					dr1p*profile->ddwdxdy(r1,r2,r3) 
					+ dr2p*profile->ddwdydy(r1,r2,r3) 
					+ dr3p*profile->ddwdydz(r1,r2,r3)
					);
			return  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);

		case 17:
			// define source function d/ds[ dnu3/dphi ]
			sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dz(r1,r2,r3) 
				+ nu1*profile->dudz(r1,r2,r3) 
				+ nu2*profile->dvdz(r1,r2,r3)
				+ nu3*profile->dwdz(r1,r2,r3);
			sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)
				* profile->dc0dz(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3) 
				+ dnu1p*profile->dudz(r1,r2,r3) 
				+ dnu2p*profile->dvdz(r1,r2,r3)
				+ dnu3p*profile->dwdz(r1,r2,r3);
			sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
				* (
					dr1p*profile->ddc0dxdz(r1,r2,r3) 
					+ dr2p*profile->ddc0dydz(r1,r2,r3) 
					+ dr3p*profile->ddc0dzdz(r1,r2,r3)
				  ) 
				+ nu1 * (
					dr1p*profile->ddudxdz(r1,r2,r3) 
					+ dr2p*profile->ddudydz(r1,r2,r3) 
					+ dr3p*profile->ddudzdz(r1,r2,r3)
					) 
				+ nu2 * (
					dr1p*profile->ddvdxdz(r1,r2,r3) 
					+ dr2p*profile->ddvdydz(r1,r2,r3) 
					+ dr3p*profile->ddvdzdz(r1,r2,r3)
					)
				+ nu3 * (
					dr1p*profile->ddwdxdz(r1,r2,r3) 
					+ dr2p*profile->ddwdydz(r1,r2,r3) 
					+ dr3p*profile->ddwdzdz(r1,r2,r3)
					);
			return  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);





		default:
			throw new NCPA::RequestOutOfBoundsException( "Only 18 equations defined!" );
	}
}
*/

double NCPA::Acoustic3DEquationSet::calculateAmplitude( double **soln, int i ) const {
	
	double rho_i = profile->rho(soln[i][0], soln[i][1], soln[i][2]);
	double mag_nu_i = mag_nu(soln[i][0], soln[i][1], soln[i][2], soln[i][3], soln[i][4], soln[i][5]);
	double c0_i = profile->c0(soln[i][0], soln[i][1], soln[i][2]);
	double cg_mag_0 = cg_mag(soln[0][0], soln[0][1], soln[0][2], soln[0][3], soln[0][4], soln[0][5]);
	double rho_0 = profile->rho(soln[0][0], soln[0][1], soln[0][2]);
	double mag_nu_0 = mag_nu(soln[0][0], soln[0][1], soln[0][2], soln[0][3], soln[0][4], soln[0][5]);
	double c0_0 = profile->c0(soln[0][0], soln[0][1], soln[0][2]);
	double cg_mag_i = cg_mag(soln[i][0], soln[i][1], soln[i][2], soln[i][3], soln[i][4], soln[i][5]);
	double Jac = Jacobian( soln, i );

/*
	double A_num = profile->rho(soln[i][0], soln[i][1], soln[i][2])
		* mag_nu(soln[i][0], soln[i][1], soln[i][2], soln[i][3], soln[i][4], soln[i][5])
		* pow(profile->c0(soln[i][0], soln[i][1], soln[i][2]),3.0)
		* cg_mag(soln[0][0], soln[0][1], soln[0][2], soln[0][3], soln[0][4], soln[0][5])
		* cos(takeoff);

	double A_den = profile->rho(soln[0][0], soln[0][1], soln[0][2])
		* mag_nu(soln[0][0], soln[0][1], soln[0][2], soln[0][3], soln[0][4], soln[0][5])
		* pow(profile->c0(soln[0][0], soln[0][1], soln[0][2]),3.0)
		* cg_mag(soln[i][0], soln[i][1], soln[i][2], soln[i][3], soln[i][4], soln[i][5])
		* Jacobian( soln, i );
*/
	double A_num = rho_i * mag_nu_i * std::pow( c0_i, 3 ) * cg_mag_0 * cos( takeoff );
	double A_den = rho_0 * mag_nu_0 * std::pow( c0_0, 3 ) * cg_mag_i * Jac;

	/*
	std::cout << "rho_i = " << rho_i << std::endl 
		  << "mag_nu_i = " << mag_nu_i << std::endl 
		  << "c0_i = " << c0_i << std::endl
		  << "cg_mag_0 = " << cg_mag_0 << std::endl
		  << "rho_0 = " << rho_0 << std::endl
		  << "mag_nu_0 = " << mag_nu_0 << std::endl
		  << "c0_0 = " << c0_0 << std::endl
		  << "cg_mag_i = " << cg_mag_i << std::endl
		  << "Jacobian = " << Jac << std::endl;
	*/


	return 1.0 / (4.0 * Pi) * sqrt( fabs( A_num / A_den ) );
}

double NCPA::Acoustic3DEquationSet::transmissionLoss( double **soln, int k, double &refDist ) const {
	
	// iterate through ranges until you get past refDist km.  Since we expect refDist to occur very early
	// in the solution, this is probably comparably efficient to a bisection search
	int ceiling = 0;
	double range = std::sqrt(soln[ceiling][0]*soln[ceiling][0] + soln[ceiling][1]*soln[ceiling][1]);
	double bottomrange;
	while (range < refDist) {
		ceiling++;
		bottomrange = range;
		range = std::sqrt(soln[ceiling][0]*soln[ceiling][0] + soln[ceiling][1]*soln[ceiling][1]);
	}
	int index;
	
	if (std::fabs( bottomrange - refDist ) < std::fabs( range - refDist )) {
		index = ceiling-1;
		refDist = bottomrange;
	} else {
		index = ceiling;
		refDist = range;
	}
	double A0 = this->calculateAmplitude(soln, index);
	double A = this->calculateAmplitude(soln, k);
	return 20*std::log10(A/A0);
}



double NCPA::Acoustic3DEquationSet::Jacobian( double **soln, int k ) const {

	double D1 = cg1(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5])
		/ cg_mag(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5])
		* (soln[k][7]*soln[k][14] - soln[k][13]*soln[k][8]);

        double D2 = soln[k][6]
		/ cg_mag(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5])
		* (soln[k][13] * cg3(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5]) 
		- cg2(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5])
		* soln[k][14]);

	double D3 = soln[k][12]
		/ cg_mag(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5])
		* (soln[k][8] * cg2(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5]) 
		- cg3(soln[k][0],soln[k][1],soln[k][2],soln[k][3],soln[k][4],soln[k][5]) 
			* soln[k][7]);

	return D1 + D2 + D3;
}

double NCPA::Acoustic3DEquationSet::calculateTravelTime( double **soln, int steps, 
		double step_size, double azimuth ) const {

	double tau = 0.0;
	for (int i = 0; i <= steps; i++) {
		tau += step_size / cg_mag(soln[i][0],soln[i][1],soln[i][2], soln[i][3], soln[i][4], soln[i][5]);
	}
	return tau;
}

double NCPA::Acoustic3DEquationSet::calculateArrivalAzimuth( double **solution, int steps,
		int averages ) const {

	// We take the average azimuth of the last five steps
	int laststeps = (steps-1) < averages ? (steps-1) : averages;
	double mean_az = 0.0;
	for (int index = steps - laststeps + 1; index <= steps; index++) {
		double nu1 = solution[index][3];
		double nu2 = solution[index][4];
		double step_az = 90 - rad2deg( atan2( -nu2, -nu1 ) );
		step_az -= 180.0;
		while (step_az < 0) 
			step_az += 360;

		mean_az += step_az / laststeps;
	}
	return mean_az;
}

void NCPA::Acoustic3DEquationSet::results( double s, double *current_values, double *new_values ) {

	// Identify the current values by name, for ease of reading/debugging
	double r1 = current_values[ 0 ];
	double r2 = current_values[ 1 ];
	double r3 = current_values[ 2 ];
	double nu1 = current_values[ 3 ];
	double nu2 = current_values[ 4 ];
	double nu3 = current_values[ 5 ];
	double dr1t = current_values[ 6 ];
	double dr2t = current_values[ 7 ];
	double dr3t = current_values[ 8 ];
	double dnu1t = current_values[ 9 ];
	double dnu2t = current_values[ 10 ];
	double dnu3t = current_values[ 11 ];
	double dr1p = current_values[ 12 ];
	double dr2p = current_values[ 13 ];
	double dr3p = current_values[ 14 ];
	double dnu1p = current_values[ 15 ];
	double dnu2p = current_values[ 16 ];
	double dnu3p = current_values[ 17 ];

	// temp variables
	double sub_term1, sub_term2, sub_term3;
	
	// Get all the quantities from the atmospheric profile that we're going to need
	double dc0dx = profile->dc0dx(r1,r2,r3);
	double dudx = profile->dudx(r1,r2,r3);
	double dvdx = profile->dvdx(r1,r2,r3);
	double dwdx = profile->dwdx(r1,r2,r3);
	
	double ddc0dxdx = profile->ddc0dxdx(r1,r2,r3);
	double ddudxdx = profile->ddudxdx(r1,r2,r3);
	double ddvdxdx = profile->ddvdxdx(r1,r2,r3);
	double ddwdxdx = profile->ddwdxdx(r1,r2,r3);
	double ddc0dxdy = profile->ddc0dxdy(r1,r2,r3);
	double ddudxdy = profile->ddudxdy(r1,r2,r3);
	double ddvdxdy = profile->ddvdxdy(r1,r2,r3);
	double ddwdxdy = profile->ddwdxdy(r1,r2,r3);
	double ddc0dxdz = profile->ddc0dxdz(r1,r2,r3);
	double ddudxdz = profile->ddudxdz(r1,r2,r3);
	double ddvdxdz = profile->ddvdxdz(r1,r2,r3);
	double ddwdxdz = profile->ddwdxdz(r1,r2,r3);
	
	double dc0dy = profile->dc0dy(r1,r2,r3);
	double dudy = profile->dudy(r1,r2,r3);
	double dvdy = profile->dvdy(r1,r2,r3);
	double dwdy = profile->dwdy(r1,r2,r3);
	
	double ddc0dydy = profile->ddc0dydy(r1,r2,r3);
	double ddudydy = profile->ddudydy(r1,r2,r3);
	double ddvdydy = profile->ddvdydy(r1,r2,r3);
	double ddwdydy = profile->ddwdydy(r1,r2,r3);
	double ddc0dydz = profile->ddc0dydz(r1,r2,r3);
	double ddudydz = profile->ddudydz(r1,r2,r3);
	double ddvdydz = profile->ddvdydz(r1,r2,r3);
	double ddwdydz = profile->ddwdydz(r1,r2,r3);
	
	double dc0dz = profile->dc0dz(r1,r2,r3);
	double dudz = profile->dudz(r1,r2,r3);
	double dvdz = profile->dvdz(r1,r2,r3);
	double dwdz = profile->dwdz(r1,r2,r3);
	
	double ddc0dzdz = profile->ddc0dzdz(r1,r2,r3);
	double ddudzdz = profile->ddudzdz(r1,r2,r3);
	double ddvdzdz = profile->ddvdzdz(r1,r2,r3);
	double ddwdzdz = profile->ddwdzdz(r1,r2,r3);
	
	double c0 = profile->c0(r1,r2,r3);
	double c0_0 = profile->c0(0,0,profile->z0(0.0,0.0));
	double u = profile->u(r1,r2,r3);
	double v = profile->v(r1,r2,r3);
	double w = profile->w(r1,r2,r3);
	
	// precalculate some functions that call the profile, to save time
	double cg1_pc = cg1(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w);
	double cg2_pc = cg2(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w);
	double cg3_pc = cg3(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w);
	double cg_mag_pc = cg_mag(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w);
	double mag_nu_pc = mag_nu(r1,r2,r3,nu1,nu2,nu3,c0,c0_0,u,v,w);
	double du_angle_pc = du_angle(r1,r2,r3,dr1t,dr2t,dr3t,dudx,dudy,dudz);
	double dv_angle_pc = dv_angle(r1,r2,r3,dr1t,dr2t,dr3t,dvdx,dvdy,dvdz);
	double dw_angle_pc = dw_angle(r1,r2,r3,dr1t,dr2t,dr3t,dwdx,dwdy,dwdz);
	double dc_angle_pc = dc_angle(r1,r2,r3,dr1t,dr2t,dr3t,dc0dx,dc0dy,dc0dz);
	double dCg_angle_pc = dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t,
					c0, c0_0, u, v, w, dc0dx, dc0dy, dc0dz,
					dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz );
	

	// define source function dr1/ds
	//new_values[ 0 ] = cg1(r1,r2,r3,nu1,nu2,nu3) / cg_mag(r1,r2,r3,nu1,nu2,nu3);
	new_values[ 0 ] = cg1_pc / cg_mag_pc;
	
	//new_values[ 1 ] = cg2(r1,r2,r3,nu1,nu2,nu3) / cg_mag(r1,r2,r3,nu1,nu2,nu3);
	new_values[ 1 ] = cg2_pc / cg_mag_pc;
	
	//new_values[ 2 ] = cg3(r1,r2,r3,nu1,nu2,nu3) / cg_mag(r1,r2,r3,nu1,nu2,nu3);
	new_values[ 2 ] = cg3_pc / cg_mag_pc;
	
// 	new_values[ 3 ] = -1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3) 
// 		* (mag_nu(r1,r2,r3,nu1,nu2,nu3) * profile->dc0dx(r1,r2,r3) + nu1*profile->dudx(r1,r2,r3) + nu2*profile->dvdx(r1,r2,r3) + nu3*profile->dwdx(r1,r2,r3));
// 	new_values[ 4 ] = -1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3) 
// 		* (mag_nu(r1,r2,r3,nu1,nu2,nu3) * profile->dc0dy(r1,r2,r3) + nu1*profile->dudy(r1,r2,r3) + nu2*profile->dvdy(r1,r2,r3) + nu3*profile->dwdy(r1,r2,r3) );
// 	new_values[ 5 ] = -1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3) 
// 		* (mag_nu(r1,r2,r3,nu1,nu2,nu3) * profile->dc0dz(r1,r2,r3) + nu1*profile->dudz(r1,r2,r3) + nu2*profile->dvdz(r1,r2,r3) + nu3*profile->dwdz(r1,r2,r3) );
	new_values[ 3 ] = -1.0/cg_mag_pc * (mag_nu_pc * dc0dx + nu1*dudx + nu2*dvdx + nu3*dwdx );
	new_values[ 4 ] = -1.0/cg_mag_pc * (mag_nu_pc * dc0dy + nu1*dudy + nu2*dvdy + nu3*dwdy );
	new_values[ 5 ] = -1.0/cg_mag_pc * (mag_nu_pc * dc0dz + nu1*dudz + nu2*dvdz + nu3*dwdz );
	

	// define source function d/ds[ dr1/dtheta ]
// 	sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg1(r1,r2,r3,nu1,nu2,nu3);
// 	sub_term2 = du_angle(r1,r2,r3,dr1t,dr2t,dr3t) + nu1/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1t,dr2t,dr3t);
// 	sub_term3 = profile->c0(r1,r2,r3)*(dnu1t/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu1/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
// 	new_values[ 6 ] =  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = -1.0*dCg_angle_pc/pow(cg_mag_pc,2.0)*cg1_pc;
	sub_term2 = du_angle_pc + nu1/mag_nu_pc*dc_angle_pc;
	sub_term3 = c0*(dnu1t/mag_nu_pc - nu1/pow(mag_nu_pc,3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
	new_values[ 6 ] =  sub_term1 + 1.0/cg_mag_pc*(sub_term2 + sub_term3);


	// define source function d/ds[ dr2/dtheta ]
// 	sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg2(r1,r2,r3,nu1,nu2,nu3);
// 	sub_term2 = dv_angle(r1,r2,r3,dr1t,dr2t,dr3t) + nu2/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1t,dr2t,dr3t);
// 	sub_term3 = profile->c0(r1,r2,r3)*(dnu2t/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu2/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
// 	new_values[ 7 ] =  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = -1.0*dCg_angle_pc/pow(cg_mag_pc,2.0)*cg2_pc;
	sub_term2 = dv_angle_pc + nu2/mag_nu_pc*dc_angle_pc;
	sub_term3 = c0*(dnu2t/mag_nu_pc - nu2/pow(mag_nu_pc,3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
	new_values[ 7 ] =  sub_term1 + 1.0/cg_mag_pc*(sub_term2 + sub_term3);

	// define source function d/ds[ dr3/dtheta ]
// 	sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg3(r1,r2,r3,nu1,nu2,nu3);
// 	sub_term2 = dw_angle(r1,r2,r3,dr1t,dr2t,dr3t) + nu3/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1t,dr2t,dr3t);
// 	sub_term3 = profile->c0(r1,r2,r3)*(dnu3t/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu3/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
// 	new_values[ 8 ] = sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = -1.0*dCg_angle_pc/pow(cg_mag_pc,2.0)*cg3_pc;
	sub_term2 = dw_angle_pc + nu3/mag_nu_pc*dc_angle_pc;
	sub_term3 = c0*(dnu3t/mag_nu_pc - nu3/pow(mag_nu_pc,3.0)*(nu1*dnu1t + nu2*dnu2t + nu3*dnu3t));
	new_values[ 8 ] = sub_term1 + 1.0/cg_mag_pc*(sub_term2 + sub_term3);

	// define source function d/ds[ dnu1/dtheta ]
// 	sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dx(r1,r2,r3) 
// 		+ nu1*profile->dudx(r1,r2,r3) + nu2*profile->dvdx(r1,r2,r3)
// 		+ nu3*profile->dwdx(r1,r2,r3);
// 	sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*profile->dc0dx(r1,r2,r3)
// 		/ mag_nu(r1,r2,r3,nu1,nu2,nu3) 
// 		+ dnu1t*profile->dudx(r1,r2,r3) + dnu2t*profile->dvdx(r1,r2,r3)
// 		+ dnu3t*profile->dwdx(r1,r2,r3);
// 	sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
// 		* (dr1t*profile->ddc0dxdx(r1,r2,r3) + dr2t*profile->ddc0dxdy(r1,r2,r3) + dr3t*profile->ddc0dxdz(r1,r2,r3)) 
// 		+ nu1*(dr1t*profile->ddudxdx(r1,r2,r3) + dr2t*profile->ddudxdy(r1,r2,r3) + dr3t*profile->ddudxdz(r1,r2,r3)) 
// 		+ nu2*(dr1t*profile->ddvdxdx(r1,r2,r3) + dr2t*profile->ddvdxdy(r1,r2,r3) + dr3t*profile->ddvdxdz(r1,r2,r3))
// 		+ nu3*(dr1t*profile->ddwdxdx(r1,r2,r3) + dr2t*profile->ddwdxdy(r1,r2,r3) + dr3t*profile->ddwdxdz(r1,r2,r3));
// 	new_values[ 9 ] =  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = mag_nu_pc*dc0dx + nu1*dudx + nu2*dvdx + nu3*dwdx;
	sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t) * dc0dx / mag_nu_pc + dnu1t*dudx + dnu2t*dvdx + dnu3t*dwdx;
	sub_term3 = mag_nu_pc*(dr1t*ddc0dxdx + dr2t*ddc0dxdy + dr3t*ddc0dxdz) 
		+ nu1*(dr1t*ddudxdx + dr2t*ddudxdy + dr3t*ddudxdz) 
		+ nu2*(dr1t*ddvdxdx + dr2t*ddvdxdy + dr3t*ddvdxdz)
		+ nu3*(dr1t*ddwdxdx + dr2t*ddwdxdy + dr3t*ddwdxdz);
	new_values[ 9 ] =  dCg_angle_pc/pow(cg_mag_pc,2.0)*sub_term1 - 1/cg_mag_pc*(sub_term2 + sub_term3);


	// define source function d/ds[ dnu2/dtheta ]
// 	sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dy(r1,r2,r3) 
// 		+ nu1*profile->dudy(r1,r2,r3) + nu2*profile->dvdy(r1,r2,r3)
// 		+ nu3*profile->dwdy(r1,r2,r3);
// 	sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*profile->dc0dy(r1,r2,r3)
// 		/ mag_nu(r1,r2,r3,nu1,nu2,nu3) 
// 		+ dnu1t*profile->dudy(r1,r2,r3) + dnu2t*profile->dvdy(r1,r2,r3)
// 		+ dnu3t*profile->dwdy(r1,r2,r3);
// 	sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
// 		* (dr1t*profile->ddc0dxdy(r1,r2,r3) + dr2t*profile->ddc0dydy(r1,r2,r3) + dr3t*profile->ddc0dydz(r1,r2,r3)) 
// 		+ nu1*(dr1t*profile->ddudxdy(r1,r2,r3) + dr2t*profile->ddudydy(r1,r2,r3) + dr3t*profile->ddudydz(r1,r2,r3)) 
// 		+ nu2*(dr1t*profile->ddvdxdy(r1,r2,r3) + dr2t*profile->ddvdydy(r1,r2,r3) + dr3t*profile->ddvdydz(r1,r2,r3))
// 		+ nu3*(dr1t*profile->ddwdxdy(r1,r2,r3) + dr2t*profile->ddwdydy(r1,r2,r3) + dr3t*profile->ddwdydz(r1,r2,r3));
// 	new_values[ 10 ] =  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = mag_nu_pc*dc0dy + nu1*dudy + nu2*dvdy + nu3*dwdy;
	sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*dc0dy / mag_nu_pc + dnu1t*dudy + dnu2t*dvdy + dnu3t*dwdy;
	sub_term3 = mag_nu_pc*(dr1t*ddc0dxdy + dr2t*ddc0dydy + dr3t*ddc0dydz) 
		+ nu1*(dr1t*ddudxdy + dr2t*ddudydy + dr3t*ddudydz) 
		+ nu2*(dr1t*ddvdxdy + dr2t*ddvdydy + dr3t*ddvdydz)
		+ nu3*(dr1t*ddwdxdy + dr2t*ddwdydy + dr3t*ddwdydz);
	new_values[ 10 ] =  dCg_angle_pc/pow(cg_mag_pc,2.0)*sub_term1 - 1/cg_mag_pc*(sub_term2 + sub_term3);


	// define source function d/ds[ dnu3/dtheta ]
// 	sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dz(r1,r2,r3) 
// 		+ nu1*profile->dudz(r1,r2,r3) + nu2*profile->dvdz(r1,r2,r3)
// 		+ nu3*profile->dwdz(r1,r2,r3);
// 	sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*profile->dc0dz(r1,r2,r3)
// 		/ mag_nu(r1,r2,r3,nu1,nu2,nu3) 
// 		+ dnu1t*profile->dudz(r1,r2,r3) + dnu2t*profile->dvdz(r1,r2,r3)
// 		+ dnu3t*profile->dwdz(r1,r2,r3);
// 	sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
// 		* (dr1t*profile->ddc0dxdz(r1,r2,r3) + dr2t*profile->ddc0dydz(r1,r2,r3) + dr3t*profile->ddc0dzdz(r1,r2,r3)) 
// 		+ nu1*(dr1t*profile->ddudxdz(r1,r2,r3) + dr2t*profile->ddudydz(r1,r2,r3) + dr3t*profile->ddudzdz(r1,r2,r3)) 
// 		+ nu2*(dr1t*profile->ddvdxdz(r1,r2,r3) + dr2t*profile->ddvdydz(r1,r2,r3) + dr3t*profile->ddvdzdz(r1,r2,r3))
// 		+ nu3*(dr1t*profile->ddwdxdz(r1,r2,r3) + dr2t*profile->ddwdydz(r1,r2,r3) + dr3t*profile->ddwdzdz(r1,r2,r3));
// 	new_values[ 11 ] =  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1t,dr2t,dr3t,dnu1t,dnu2t,dnu3t)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = mag_nu_pc*dc0dz + nu1*dudz + nu2*dvdz + nu3*dwdz;
	sub_term2 = (nu1*dnu1t + nu2*dnu2t + nu3*dnu3t)*dc0dz / mag_nu_pc + dnu1t*dudz + dnu2t*dvdz + dnu3t*dwdz;
	sub_term3 = mag_nu_pc*(dr1t*ddc0dxdz + dr2t*ddc0dydz + dr3t*ddc0dzdz) 
		+ nu1*(dr1t*ddudxdz + dr2t*ddudydz + dr3t*ddudzdz) 
		+ nu2*(dr1t*ddvdxdz + dr2t*ddvdydz + dr3t*ddvdzdz)
		+ nu3*(dr1t*ddwdxdz + dr2t*ddwdydz + dr3t*ddwdzdz);
	new_values[ 11 ] =  dCg_angle_pc/pow(cg_mag_pc,2.0)*sub_term1 - 1/cg_mag_pc*(sub_term2 + sub_term3);

	
	// define source function d/ds[ dr1/dphi ]
// 	sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg1(r1,r2,r3,nu1,nu2,nu3);
// 	sub_term2 = du_angle(r1,r2,r3,dr1p,dr2p,dr3p) + nu1/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1p,dr2p,dr3p);
// 	sub_term3 = profile->c0(r1,r2,r3)*(dnu1p/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu1/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
// 	new_values[ 12 ] =  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = -1.0*dCg_angle_pc/pow(cg_mag_pc,2.0)*cg1_pc;
	sub_term2 = du_angle_pc + nu1/mag_nu_pc*dc_angle_pc;
	sub_term3 = c0*(dnu1p/mag_nu_pc - nu1/pow(mag_nu_pc,3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
	new_values[ 12 ] =  sub_term1 + 1.0/cg_mag_pc*(sub_term2 + sub_term3);


	// define source function d/ds[ dr2/dphi ]
// 	sub_term1 = -1.0*dCg_angle_pc(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg2(r1,r2,r3,nu1,nu2,nu3);
// 	sub_term2 = dv_angle(r1,r2,r3,dr1p,dr2p,dr3p) + nu2/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1p,dr2p,dr3p);
// 	sub_term3 = profile->c0(r1,r2,r3)*(dnu2p/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu2/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
// 	new_values[ 13 ] =  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = -1.0*dCg_angle_pc/pow(cg_mag_pc,2.0)*cg2_pc;
	sub_term2 = dv_angle_pc + nu2/mag_nu_pc*dc_angle_pc;
	sub_term3 = c0*(dnu2p/mag_nu_pc - nu2/pow(mag_nu_pc,3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
	new_values[ 13 ] =  sub_term1 + 1.0/cg_mag_pc*(sub_term2 + sub_term3);

	// define source function d/ds[ dr3/dphi ]
// 	sub_term1 = -1.0*dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*cg3(r1,r2,r3,nu1,nu2,nu3);
// 	sub_term2 = dw_angle(r1,r2,r3,dr1p,dr2p,dr3p) + nu3/mag_nu(r1,r2,r3,nu1,nu2,nu3)*dc_angle(r1,r2,r3,dr1p,dr2p,dr3p);
// 	sub_term3 = profile->c0(r1,r2,r3)*(dnu3p/mag_nu(r1,r2,r3,nu1,nu2,nu3) - nu3/pow(mag_nu(r1,r2,r3,nu1,nu2,nu3),3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
// 	new_values[ 14 ] =  sub_term1 + 1.0/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = -1.0*dCg_angle_pc/pow(cg_mag_pc,2.0)*cg3_pc;
	sub_term2 = dw_angle_pc + nu3/mag_nu_pc*dc_angle_pc;
	sub_term3 = c0*(dnu3p/mag_nu_pc - nu3/pow(mag_nu_pc,3.0)*(nu1*dnu1p + nu2*dnu2p + nu3*dnu3p));
	new_values[ 14 ] =  sub_term1 + 1.0/cg_mag_pc*(sub_term2 + sub_term3);


	// define source function d/ds[ dnu1/dphi ]
// 	sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dx(r1,r2,r3) 
// 		+ nu1*profile->dudx(r1,r2,r3) 
// 		+ nu2*profile->dvdx(r1,r2,r3)
// 		+ nu3*profile->dwdx(r1,r2,r3);
// 	sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)
// 		* profile->dc0dx(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3) 
// 		+ dnu1p*profile->dudx(r1,r2,r3) 
// 		+ dnu2p*profile->dvdx(r1,r2,r3)
// 		+ dnu3p*profile->dwdx(r1,r2,r3);
// 	sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
// 		* (
// 			dr1p*profile->ddc0dxdx(r1,r2,r3) 
// 			+ dr2p*profile->ddc0dxdy(r1,r2,r3) 
// 			+ dr3p*profile->ddc0dxdz(r1,r2,r3)
// 			) 
// 		+ nu1 * (
// 			dr1p*profile->ddudxdx(r1,r2,r3) 
// 			+ dr2p*profile->ddudxdy(r1,r2,r3) 
// 			+ dr3p*profile->ddudxdz(r1,r2,r3)
// 			) 
// 		+ nu2 * (
// 			dr1p*profile->ddvdxdx(r1,r2,r3) 
// 			+ dr2p*profile->ddvdxdy(r1,r2,r3) 
// 			+ dr3p*profile->ddvdxdz(r1,r2,r3)
// 			)
// 		+ nu3 * (
// 			dr1p*profile->ddwdxdx(r1,r2,r3) 
// 			+ dr2p*profile->ddwdxdy(r1,r2,r3) 
// 			+ dr3p*profile->ddwdxdz(r1,r2,r3)
// 			);
// 	new_values[ 15 ] =  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = mag_nu_pc*dc0dx + nu1*dudx + nu2*dvdx + nu3*dwdx;
	sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)	* dc0dx/mag_nu_pc 
		+ dnu1p*dudx + dnu2p*dvdx + dnu3p*dwdx;
	sub_term3 = mag_nu_pc * (dr1p*ddc0dxdx + dr2p*ddc0dxdy + dr3p*ddc0dxdz) 
		+ nu1 * (dr1p*ddudxdx + dr2p*ddudxdy + dr3p*ddudxdz) 
		+ nu2 * (dr1p*ddvdxdx + dr2p*ddvdxdy + dr3p*ddvdxdz)
		+ nu3 * (dr1p*ddwdxdx + dr2p*ddwdxdy + dr3p*ddwdxdz);
	new_values[ 15 ] =  dCg_angle_pc/pow(cg_mag_pc,2.0)*sub_term1 - 1/cg_mag_pc*(sub_term2 + sub_term3);



	// define source function d/ds[ dnu2/dphi ]
// 	sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dy(r1,r2,r3) 
// 		+ nu1*profile->dudy(r1,r2,r3) 
// 		+ nu2*profile->dvdy(r1,r2,r3)
// 		+ nu3*profile->dwdy(r1,r2,r3);
// 	sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)
// 		* profile->dc0dy(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3) 
// 		+ dnu1p*profile->dudy(r1,r2,r3) 
// 		+ dnu2p*profile->dvdy(r1,r2,r3)
// 		+ dnu3p*profile->dwdy(r1,r2,r3);
// 	sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
// 		* (
// 			dr1p*profile->ddc0dxdy(r1,r2,r3) 
// 			+ dr2p*profile->ddc0dydy(r1,r2,r3) 
// 			+ dr3p*profile->ddc0dydz(r1,r2,r3)
// 			) 
// 		+ nu1 * (
// 			dr1p*profile->ddudxdy(r1,r2,r3) 
// 			+ dr2p*profile->ddudydy(r1,r2,r3) 
// 			+ dr3p*profile->ddudydz(r1,r2,r3)
// 			) 
// 		+ nu2 * (
// 			dr1p*profile->ddvdxdy(r1,r2,r3) 
// 			+ dr2p*profile->ddvdydy(r1,r2,r3) 
// 			+ dr3p*profile->ddvdydz(r1,r2,r3)
// 			)
// 		+ nu3 * (
// 			dr1p*profile->ddwdxdy(r1,r2,r3) 
// 			+ dr2p*profile->ddwdydy(r1,r2,r3) 
// 			+ dr3p*profile->ddwdydz(r1,r2,r3)
// 			);
// 	new_values[ 16 ] =  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = mag_nu_pc*dc0dy + nu1*dudy + nu2*dvdy + nu3*dwdy;
	sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)	* dc0dy/mag_nu_pc 
		+ dnu1p*dudy + dnu2p*dvdy + dnu3p*dwdy;
	sub_term3 = mag_nu_pc * (dr1p*ddc0dxdy + dr2p*ddc0dydy + dr3p*ddc0dydz) 
		+ nu1 * (dr1p*ddudxdy + dr2p*ddudydy + dr3p*ddudydz) 
		+ nu2 * (dr1p*ddvdxdy + dr2p*ddvdydy + dr3p*ddvdydz)
		+ nu3 * (dr1p*ddwdxdy + dr2p*ddwdydy + dr3p*ddwdydz);
	new_values[ 16 ] =  dCg_angle_pc/pow(cg_mag_pc,2.0)*sub_term1 - 1/cg_mag_pc*(sub_term2 + sub_term3);



	// define source function d/ds[ dnu3/dphi ]
// 	sub_term1 = mag_nu(r1,r2,r3,nu1,nu2,nu3)*profile->dc0dz(r1,r2,r3) 
// 		+ nu1*profile->dudz(r1,r2,r3) 
// 		+ nu2*profile->dvdz(r1,r2,r3)
// 		+ nu3*profile->dwdz(r1,r2,r3);
// 	sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p)
// 		* profile->dc0dz(r1,r2,r3)/mag_nu(r1,r2,r3,nu1,nu2,nu3) 
// 		+ dnu1p*profile->dudz(r1,r2,r3) 
// 		+ dnu2p*profile->dvdz(r1,r2,r3)
// 		+ dnu3p*profile->dwdz(r1,r2,r3);
// 	sub_term3 = mag_nu(r1,r2,r3,nu1,nu2,nu3)
// 		* (
// 			dr1p*profile->ddc0dxdz(r1,r2,r3) 
// 			+ dr2p*profile->ddc0dydz(r1,r2,r3) 
// 			+ dr3p*profile->ddc0dzdz(r1,r2,r3)
// 			) 
// 		+ nu1 * (
// 			dr1p*profile->ddudxdz(r1,r2,r3) 
// 			+ dr2p*profile->ddudydz(r1,r2,r3) 
// 			+ dr3p*profile->ddudzdz(r1,r2,r3)
// 			) 
// 		+ nu2 * (
// 			dr1p*profile->ddvdxdz(r1,r2,r3) 
// 			+ dr2p*profile->ddvdydz(r1,r2,r3) 
// 			+ dr3p*profile->ddvdzdz(r1,r2,r3)
// 			)
// 		+ nu3 * (
// 			dr1p*profile->ddwdxdz(r1,r2,r3) 
// 			+ dr2p*profile->ddwdydz(r1,r2,r3) 
// 			+ dr3p*profile->ddwdzdz(r1,r2,r3)
// 			);
// 	new_values[ 17 ] =  dCg_angle(r1,r2,r3,nu1,nu2,nu3,dr1p,dr2p,dr3p,dnu1p,dnu2p,dnu3p)/pow(cg_mag(r1,r2,r3,nu1,nu2,nu3),2.0)*sub_term1 - 1/cg_mag(r1,r2,r3,nu1,nu2,nu3)*(sub_term2 + sub_term3);
	sub_term1 = mag_nu_pc*dc0dz + nu1*dudz + nu2*dvdz + nu3*dwdz;
	sub_term2 = (nu1*dnu1p + nu2*dnu2p + nu3*dnu3p) * dc0dz/mag_nu_pc 
		+ dnu1p*dudz + dnu2p*dvdz + dnu3p*dwdz;
	sub_term3 = mag_nu_pc * (dr1p*ddc0dxdz + dr2p*ddc0dydz + dr3p*ddc0dzdz) 
		+ nu1 * (dr1p*ddudxdz + dr2p*ddudydz + dr3p*ddudzdz) 
		+ nu2 * (dr1p*ddvdxdz + dr2p*ddvdydz + dr3p*ddvdzdz)
		+ nu3 * (dr1p*ddwdxdz + dr2p*ddwdydz + dr3p*ddwdzdz);
	new_values[ 17 ] =  dCg_angle_pc/pow(cg_mag_pc,2.0)*sub_term1 - 1/cg_mag_pc*(sub_term2 + sub_term3);



}
