#include "twobody.h"
#include "../vec3.h"
#include "../system.h"
#include <cmath>


TwoBody::TwoBody(std::string filename){
	m_filename = filename;
}

TwoBody::TwoBody(std::string filename, std::string planet){
	m_filename = filename;
	m_planet = planet;
}

void TwoBody::setupParticles(System &system) {
	readEphemerides();
	for (int i = 0; i < m_input_data.size();i++){
		Particle* p = m_input_data.at(i);
		if (p->getName() == "Sun") system.addParticle(p);
		if (p->getName() == m_planet) system.addParticle(p);
	}
}

std::string TwoBody::getName() {
    return "Two-body";
}
