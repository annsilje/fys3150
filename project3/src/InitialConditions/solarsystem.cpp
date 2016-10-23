#include "solarsystem.h"
#include "../vec3.h"
#include "../system.h"
#include <cmath>


SolarSystem::SolarSystem(std::string filename){
	m_filename = filename;
}

void SolarSystem::setupParticles(System &system) {
	readEphemerides();
	for (int i = 0; i < m_input_data.size();i++){
		Particle* p = m_input_data.at(i);
		system.addParticle(p);
	}

}

std::string SolarSystem::getName() {
    return "Solar System";
}
