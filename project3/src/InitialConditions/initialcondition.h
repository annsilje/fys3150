#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "../particle.h"

class InitialCondition {
protected:
	std::string m_filename;
	std::vector<Particle*> m_input_data;
public:
	InitialCondition(){}
    InitialCondition(std::string filename);
    void readEphemerides();
    std::string getEphFile(){return m_filename;}
    virtual void setupParticles(class System& system) = 0;
    virtual std::string getName();
    virtual ~InitialCondition();
};

