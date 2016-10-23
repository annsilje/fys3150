#include "velocityverlet.h"
#include "../system.h"

VelocityVerlet::VelocityVerlet(System* system, double dt, bool sunFixed)
    : Integrator(system, dt, sunFixed) {
}

std::string VelocityVerlet::getName() {
    return "Velocity verlet";
}

void VelocityVerlet::integrateOneStep(std::vector<Particle*> particles) {

	//First iteration only
	if (m_prev_a.size()==0){
		m_system->computeForces();
		for (int i=0; i < particles.size();i++){
			Particle* p = particles.at(i);
			m_prev_a.push_back(p->getForce()/p->getMass());
			m_prev_v.push_back(p->getVelocity());
		}
	}

	//Update position
    for (int i=0; i<particles.size(); i++) {
        Particle *p = particles.at(i);
        vec3& x = p->getPosition();
        if (!(p->getName()=="Sun" && m_sun_fixed)){
        	x += m_dt*m_prev_v.at(i) + m_dt*m_dt/2*m_prev_a.at(i); //6 FLOPs
        }
    }

    //Compute forces for the new position
    m_system->computeForces();
    //Update velocity and store the new acceleration
    for (int i=0; i<particles.size(); i++) {
        Particle *p = particles.at(i);

        vec3& v = p->getVelocity();
        vec3 a = p->getForce()/p->getMass();

        if (!(p->getName()=="Sun" && m_sun_fixed)){
        	v += m_dt/2*(a+m_prev_a[i]); //4 FLOPs
        }
        m_prev_a.at(i) = a;
        m_prev_v.at(i) = v;
    }
}
