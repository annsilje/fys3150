#include <iostream>
#include <cmath>
#include "modifiednewtoniangravity.h"

ModifiedNewtonianGravity::ModifiedNewtonianGravity(double G, double c) : m_G(G), m_c(c) {

}

void ModifiedNewtonianGravity::computeForces(std::vector<Particle*> particles){
	//Assume only two bodies!
	Particle* big;
	Particle* small;
	Particle* a = particles.at(0);
	Particle* b = particles.at(1);
	if (a->getName()=="Sun"){
		big = a;
		small = b;
	}
	else{
		big = b;
		small = a;
	}
	computeForce(*big, *small);
}

void ModifiedNewtonianGravity::computeForce(Particle& a, Particle& b) {

	vec3 r = b.getPosition() - a.getPosition();
	vec3 v = b.getVelocity() - a.getVelocity();
	double r_length = r.length();
	double GMM = -m_G*a.getMass()*b.getMass();
	double l = (r.cross(v)).length();
	double l2 = l*l;
	double c2 = m_c*m_c;
	double r2 = r_length*r_length;
	double r3 = r2*r_length;

	//m_potentialEnergy += GMM/r_length;
	double gr = 3*l2/(r2*c2);
	//std::cout << gr << std::endl;
	vec3 F = GMM*r/r3*(1 + gr);

	b.addForce(F.x(), F.y(), F.z());
	a.addForce(-F.x(), -F.y(), -F.z());
}

std::string ModifiedNewtonianGravity::getName() {
    return "Post-Newtonian gravity";
}
