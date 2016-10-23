#include "forwardeuler.h"
#include "../system.h"

ForwardEuler::ForwardEuler(System* system, double dt, bool sunFixed)
    : Integrator(system, dt, sunFixed) {
}

void ForwardEuler::integrateOneStep(std::vector<Particle*> particles) {
    m_system->computeForces();
    for (int i=0; i<particles.size(); i++) {
        Particle *p = particles.at(i);

        vec3& x = p->getPosition();
        vec3& v = p->getVelocity();
        vec3 a = p->getForce()/p->getMass();
        if (!(p->getName()=="Sun" && m_sun_fixed)){
        	x += m_dt*v; //2 FLOPs
        	v += m_dt*a; //2 FLOPs
        }
    }
}

std::string ForwardEuler::getName() {
    return "Forward-Euler";
}
