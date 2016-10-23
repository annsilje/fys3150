#include "newtoniangravity.h"
#include <iostream>
#include <cmath>

NewtonianGravity::NewtonianGravity(double G) : m_G(G) {

}

void NewtonianGravity::computeForces(std::vector<Particle*> particles){
	for (int i=0; i<particles.size();i++){
		for (int j=i+1; j<particles.size();j++){
			Particle* a = particles.at(i);
			Particle* b = particles.at(j);
			computeForce(*a, *b);
		}
	}
}

void NewtonianGravity::computeForce(Particle& a, Particle& b) {

	vec3 r = a.getPosition() - b.getPosition();
	double r_length = r.length();
	double GMM = -m_G*a.getMass()*b.getMass();

	m_potentialEnergy += GMM/r_length;
	vec3 F_a = GMM*r/pow(r_length,3);

	a.addForce(F_a.x(), F_a.y(), F_a.z());
	b.addForce(-F_a.x(), -F_a.y(), -F_a.z());
}

std::string NewtonianGravity::getName() {
    return "Newtonian gravity";
}
