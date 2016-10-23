#include "threebody.h"
#include "../vec3.h"
#include "../system.h"
#include <cmath>


using namespace std;

ThreeBody::ThreeBody(string filename){
	m_filename = filename;
}

ThreeBody::ThreeBody(string filename, double factor){
	m_filename = filename;
	m_factor = factor;
}

void ThreeBody::setupParticles(System &system) {
	readEphemerides();
	for (int i = 0; i < m_input_data.size();i++){
		Particle* p = m_input_data.at(i);
		if (p->getName() == "Sun") system.addParticle(p);
		if (p->getName() == "Earth") system.addParticle(p);
		if (p->getName() == "Jupiter"){
			p->increaseMass(m_factor);
			system.addParticle(p);
		}
	}

}

string ThreeBody::getName() {
	stringstream ss;
	ss << "Three-body " << m_factor << "xJupiter Mass";
    return ss.str();
}
