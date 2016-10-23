#pragma once
#include "../InitialConditions/initialcondition.h"
#include "../particle.h"
#include <vector>
#include <string>


class TwoBody : public InitialCondition {
private:
	std::string m_planet = "Earth";
public:
    TwoBody() {}
    TwoBody(std::string filename);
    TwoBody(std::string filename, std::string planet);
    void setupParticles(class System& system);
    std::string getName();
};

