#pragma once
#include "integrator.h"
#include "../particle.h"
#include <vector>

class ForwardEuler : public Integrator {
public:
    ForwardEuler(class System* system, double dt, bool sunFixed);
    void integrateOneStep(std::vector<Particle*> particles);
    std::string getName();
};
