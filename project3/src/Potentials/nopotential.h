#pragma once
#include "potential.h"
#include <string>

class NoPotential : public Potential {
public:
    NoPotential() {}
    std::string getName();
    void computeForces(std::vector<Particle*> particles);
    void computeForce(Particle&, Particle&);
};
