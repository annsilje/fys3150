#pragma once
#include "../particle.h"
#include "potential.h"
#include <string>
#include <vector>

class ModifiedNewtonianGravity : public Potential {
private:
    double m_G;
    double m_c;

public:
    ModifiedNewtonianGravity(double G, double c);
    void computeForces(std::vector<Particle*> particles);
    void computeForce(Particle& a, Particle& b);
    std::string getName();
};
