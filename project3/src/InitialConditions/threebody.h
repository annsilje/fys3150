#pragma once
#include "../InitialConditions/initialcondition.h"
#include <string>
#include <sstream>
#include <iomanip>

class ThreeBody : public InitialCondition {
private:
	double m_factor = 1.0;
public:
    ThreeBody() {}
    ThreeBody(std::string filename);
    ThreeBody(std::string filename, double factor);
    void setupParticles(class System& system);
    std::string getName();
};
