#pragma once
#include "integrator.h"
#include "../particle.h"
#include <string>

class VelocityVerlet : public Integrator {
private:
    bool m_firstStep = true;
    std::vector<vec3> m_prev_a;
    std::vector<vec3> m_prev_v;

public:
    VelocityVerlet(class System* system, double dt, bool sunFixed);
    std::string getName();
    void integrateOneStep(std::vector<Particle*> particles);
};
