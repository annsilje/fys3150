#pragma once
#include "../InitialConditions/initialcondition.h"
#include <string>

class SolarSystem : public InitialCondition {
public:
	SolarSystem() {}
	SolarSystem(std::string filename);
    void setupParticles(class System& system);
    std::string getName();
};
